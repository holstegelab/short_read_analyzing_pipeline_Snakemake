import argparse
import sys
import csv
import asyncio
import gzip
import os
import threading
import re
import time
import numpy
from Bio.Seq import Seq
from bam_utils import *
      
class PairedFastQReader(threading.Thread):
    def __init__(self, filea,fileb):
        self.filea = filea
        self.fileb = fileb
        self.buffer = {}
        self.splitchar = None
        self.buffer_size = 1000

        self.file_finished = False
        self.continue_reading = threading.Event()
        self.search_read = None
        self.read_found_or_buffer_full = threading.Event()

        threading.Thread.__init__(self)


    def run(self):
        with os.popen('pigz -dc ' + self.filea) as f1:
            with os.popen('pigz -dc ' + self.fileb) as f2:
                try:
                    f1 = iter(f1)
                    f2 = iter(f2)
                    while True:
                        header1 = f1.__next__().strip()
                        seq1 = f1.__next__().strip()
                        dummy = f1.__next__()
                        qual1 = f1.__next__().strip()
                        header2 = f2.__next__().strip()
                        seq2 = f2.__next__().strip()
                        dummy = f2.__next__()
                        qual2 = f2.__next__().strip()
                        assert header1.startswith('@'), 'Parsing error'
                        assert header2.startswith('@'), 'Parsing error'
                        if self.splitchar is None:
                            if ' ' in header1:
                                self.splitchar = ' '
                            elif '/' in header1:
                                self.splitchar = '/'
                            else:
                                raise RuntimeError('Cannot determine FastQ header split character')
                                
                        name1, r1 = header1.split(self.splitchar)
                        name2, r2 = header2.split(self.splitchar)
                        name1 = name1[1:] #strip '@'
                        name2 = name2[1:] #strip '@'

                        assert name1 == name2, 'FastQ input files not in same read order'

                        self.buffer[name1] = (name1, seq1, seq2, qual1, qual2)
                        
                        if name1 == self.search_read or self.search_read is True:
                            self.read_found_or_buffer_full.set()

                        if len(self.buffer) >= self.buffer_size:
                            self.read_found_or_buffer_full.set()
                            self.continue_reading.clear()
                            self.continue_reading.wait()

                except StopIteration:
                    self.read_found_or_buffer_full.set()
                    self.file_finished = True

    def popRead(self, name):
        while not name in self.buffer:
            self.search_read = name
            self.read_found_or_buffer_full.clear()
            while True:
                self.read_found_or_buffer_full.wait(timeout=5)
                if name in self.buffer:
                    break
                elif self.file_finished:
                    raise RuntimeError(f"Read {name} not found while file has ended. Is BAM file reordered?")
                elif len(self.buffer) >= self.buffer_size: #buffer full, read not found
                    self.buffer_size *=2
                    sys.stderr.write(f"Order of reads has changed between BAM and FASTQ. Growing look-back buffer to {self.buffer_size}\n")
                    sys.stderr.flush()
                    self.continue_reading.set()
                    break
                else:    
                    self.continue_reading.set()
                    #wait a bit more
                    continue

        read = self.buffer[name]
        del self.buffer[name]
        
        if len(self.buffer) < self.buffer_size:
            self.continue_reading.set()
        return read

    def retrieveRead(self):
        if len(self.buffer) < self.buffer_size:
            self.continue_reading.set()

        if not self.buffer:
            if self.file_finished:
                raise StopIteration()
                
            self.read_found_or_buffer_full.clear()
            self.search_read = True
            self.read_found_or_buffer_full.wait(timeout=15)
            if not self.buffer:
                assert self.file_finished
                raise StopIteration()

        return self.buffer.popitem()[1]





def derive_missing_sequence_tags(read, seq, qual, stats):

    if (read.flag & 0x900) != 0.0: # not a primary read
        if (read.flag & 0x800):
            stats['supplementary_alignments'] = stats.get('supplementary_alignments',0) + 1
            stats['supplementary_aligned_bp'] = stats.get('supplementary_bp',0) + read.get_aligned_read_length()
        elif (read.flag & 0x100):            
            stats['secondary_alignments'] = stats.get('secondary_alignments',0) + 1
            stats['secondary_aligned_bp'] = stats.get('secondary_bp',0) + read.get_aligned_read_length()
        return (read.cigar,{})

    
    stats['primary_reads'] = stats.get('primary_reads',0) + 1
    stats['primary_aligned_bp'] = stats.get('primary_aligned_bp',0) + read.get_aligned_read_length()

    clipping_length = len(seq) - len(read.seq)

    if clipping_length == 0: #nothing is missing w.r.t. fastq sequence
        return (read.cigar,{})
    
    #reverse to check/match sequence in record 
    tag_seq = seq[len(read.seq):]    #hardclip tags XB/XQ are w.r.t. to non-reverse-complemented sequence
    tag_qual = qual[len(read.seq):] 

    stats['size'][len(tag_seq)] = stats['size'].get(len(tag_seq),0) + 1

    if (read.flag & 0x40):
        stats['restored_read1s'] = stats.get('restored_read1s',0) + 1
        stats['restored_bp_read1'] = stats.get('restored_bp_read1',0) + clipping_length
    else:            
        stats['restored_read2s'] = stats.get('restored_read2s',0) + 1
        stats['restored_bp_read2'] = stats.get('restored_bp_read2',0) + clipping_length
    
    if read.cigar != '*':
        cigars = read.split_cigar()
        if read.is_reversed():
            rseq = str(Seq(seq).reverse_complement())
            rqual = qual[::-1]
            assert rseq[-len(read.seq):] == read.seq, f"Sequence mismatch BAM-FASTQ in fragment {read.qname}"
            assert rqual[-len(read.qual):] == read.qual, f"Quality score mismatch BAM-FASTQ in fragment {read.qname}"
            while cigars[0][0] == 'H': #remove all hard clipping cigar elements
                cigars = cigars[1:]
            cigars = [('H', clipping_length)] + cigars #re-add with correct length
            assert not any([x == 'H' for x,y in cigars[1:]]), 'Hard clipping found, but not at end of read'
        else:
            assert seq[:len(read.seq)] == read.seq, f"Sequence mismatch BAM-FASTQ in fragment {read.qname}"
            assert qual[:len(read.seq)] == read.qual, f"Quality score mismatch BAM-FASTQ in fragment {read.qname}"
            while cigars[-1][0] == 'H': #remove all hard clipping cigar elements
                cigars.pop()
            cigars.append(('H', clipping_length)) #re-add with correct length
            assert not any([x == 'H' for x,y in cigars[:-1]]), 'Hard clipping found, but not at end of read'

        ncigar = join_cigar(cigars)
    else:
        stats['restored_unaligned_reads'] = stats.get('restored_unaligned_reads',0) + 1
        ncigar = '*'

    return (ncigar, {'XB':tag_seq, 'XZ':tag_qual})        #XQ is also used by dragmap, need to rename before running revertsam
        

def process(querygroup, fastqrecord, stats):
    name, seq1, seq2, qual1, qual2 = fastqrecord
    
    stats['alignments'] = stats.get('alignments',0) + len(querygroup)
    stats['fragments'] = stats.get('fragments',0) + 1

    stats['total_bp'] = stats.get('total_bp',0) + len(seq1) + len(seq2)
   
    for read in querygroup:
        tags = []
        if (read.flag & 0x40): # r1
            ncigar, tags = derive_missing_sequence_tags(read, seq1, qual1, stats)
        elif (read.flag & 0x80):
            ncigar, tags = derive_missing_sequence_tags(read, seq2, qual2, stats)
        if tags:
            read.cigar = ncigar
            read.setTagValues(**tags)

    rg = querygroup[0].getTagValue('RG')
    return (querygroup,rg)


def fastqtosam(fastqrecord, rg, stats):
    name, seq1, seq2, qual1, qual2 = fastqrecord
    stats['restored_bp_read1'] = stats.get('restored_bp_read1') + len(seq1)
    stats['restored_bp_read2'] = stats.get('restored_bp_read2') + len(seq1)
    stats['readded_fragments'] = stats.get('readded_fragments',0) + 1


    stats['size'][len(seq1)] = stats['size'].get(len(seq1),0) + 1
    stats['size'][len(seq2)] = stats['size'].get(len(seq2),0) + 1

    flag1 = 0x1 |  0x4 | 0x8| 0x40 | 0x200  #(PAIRED | UNMAP | MUNMAP | READ1 | QCFAIL)
    flag2 = 0x1 | 0x4 | 0x8| 0x80 | 0x200  #(PAIRED | UNMAP | MUNMAP | READ2 | QCFAIL)
    querygroup = []
    read1 = BamRecord([name, str(flag1), '*', '0', '0', '*', '*', '0', '0', seq1, qual1])
    read2 = BamRecord([name, str(flag2), '*', '0', '0', '*', '*', '0', '0', seq1, qual1])
    if not rg is None:
        read1.setTagValues(RG=rg)
        read2.setTagValues(RG=rg)
    return [read1,read2]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", help='first fastq file')
    parser.add_argument("-b", help='second fastq file')
    parser.add_argument("-o", default='-', help='output sam file (default stdout)')
    parser.add_argument("-s", help='stats file')
    args = parser.parse_args()
    
    freader = PairedFastQReader(args.a, args.b)
    freader.start()
    rg = None

    alignment_counter = 0
    stats = {}

    stats['size'] = {}
    last_time = time.time()
    with sys.stdin as f:
        reader = csv.reader(f,delimiter='\t', quoting=csv.QUOTE_NONE)

        lastname = ''
        querygroup = []
        if args.o == '-':
            out = sys.stdout
        else:
            out = open(out_file, 'wt')
            
        #with open('test.out', 'w') as out_file:
        with out as out_file:
            for row in reader:
                alignment_counter += 1
                if (alignment_counter % 100000) == 0:
                    sys.stderr.write('%d reads done at %d reads/sec\r' % (alignment_counter, int(100000.0 / (time.time() - last_time))))
                    sys.stderr.flush()
                    last_time = time.time()
                if row[0][0] == '@':
                    out_file.write('\t'.join(row) + '\n')
                    continue
                row = BamRecord(row)
                if lastname == row.qname:
                    querygroup.append(row)
                else:
                    
                    if querygroup:
                        assert len(querygroup) >= 2, 'Reads are not paired'
                        fastq_record = freader.popRead(lastname)
                        res_kept, rg = process(querygroup, fastq_record, stats)
                        res_kept = '\n'.join([xrow.toSamRecord() for xrow in res_kept]) + '\n'
                        out_file.write(res_kept)

                    lastname = row.qname
                    querygroup = [row]
            if querygroup:
                assert len(querygroup) >= 2, 'Reads are not paired'
                fastq_record = freader.popRead(lastname)
                #Process
                res_kept, rg = process(querygroup, fastq_record, stats)

                res_kept = '\n'.join([xrow.toSamRecord() for xrow in res_kept]) + '\n'
                out_file.write(res_kept)

            extra_fragment_counter = 0
            sys.stderr.write('\nRe-adding removed fragments\n')
            sys.stderr.flush()
            try:
                while True:
                    x = freader.retrieveRead()
                    extra_fragment_counter += 1
                    if (extra_fragment_counter % 1000) == 0:
                        sys.stderr.write('%d fragments done\r' % extra_fragment_counter);
                        sys.stderr.flush()

                    res = fastqtosam(x, rg, stats)

                    #write out
                    res = '\n'.join([xrow.toSamRecord() for xrow in res]) + '\n'
                    out_file.write(res_kept)
            except StopIteration:
                pass
        
        #process size stats
        size_stats = numpy.zeros(max(stats['size'].keys()) + 1,dtype=int)
        for key,value in stats['size'].items():
            size_stats[key] = value
        del stats['size']

        with open(args.s,'wt') as f:            
            for key, value in stats.items():
                f.write('%s\t %d\n' % (key,value))
            f.write('\nadapter_length\tcount\n')
            for pos, s in enumerate(size_stats):
                f.write('%d\t%d\n' % (pos,s))
                

            f.flush()
            
        

