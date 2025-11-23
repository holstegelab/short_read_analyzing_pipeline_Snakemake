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
from fastq_utils import *

def check_qual(read_name, q1, q2):
    if q1 == q2:
        return

    #dragmap sometimes converts ! quality to # quality.  
    #in that case, results are not entirely lossless, but in practice it has probably no consequences. 
    if q1.replace('#', '!') == q2.replace('#','!'):
        return
    
    raise RuntimeError(f"Quality score mismatch BAM-FASTQ in fragment {read_name}: \n{q1} != \n{q2}")


def derive_missing_sequence_tags(read, seq, qual, stats):
    
    stats['primary_reads'] = stats.get('primary_reads',0) + 1
    stats['primary_aligned_bp'] = stats.get('primary_aligned_bp',0) + read.get_aligned_read_length()

    if 'S' in read.cigar:
        c = sum([e[1] for e in read.split_cigar() if e[0] == 'S'])
        stats['primary_soft_clipped_bp']  = stats.get('primary_soft_clipped_bp',0) + c

    clipping_length = len(seq) - len(read.seq)

    if clipping_length == 0: #nothing is missing w.r.t. fastq sequence
        return (read.cigar,{}, 0)
    
    #get missing sequence/qualities
    tag_seq = seq[len(read.seq):]   
    tag_qual = qual[len(read.seq):] 


    if (read.flag & 0x40):
        stats['restored_read1s'] = stats.get('restored_read1s',0) + 1
        stats['restored_bp_read1'] = stats.get('restored_bp_read1',0) + clipping_length
        stats['size_r1'][len(tag_seq)] = stats['size_r1'].get(len(tag_seq),0) + 1
    else:            
        stats['restored_read2s'] = stats.get('restored_read2s',0) + 1
        stats['restored_bp_read2'] = stats.get('restored_bp_read2',0) + clipping_length
        stats['size_r2'][len(tag_seq)] = stats['size_r2'].get(len(tag_seq),0) + 1
    
    if read.cigar != '*':
        cigars = read.split_cigar()
        if read.is_reversed():
            rseq = str(Seq(seq).reverse_complement())
            rqual = qual[::-1]
            assert rseq[-len(read.seq):] == read.seq, f"Sequence mismatch BAM-FASTQ in fragment {read.qname}: \n{read.seq} vs. \n{rseq[-len(read.seq):]}"
            check_qual(read.qname, rqual[-len(read.qual):], read.qual)

            tag_seq = str(Seq(tag_seq).reverse_complement())
            tag_qual = tag_qual[::-1]
            assert tag_seq == rseq[:-len(read.seq)]
            new_tags = {'YB': tag_seq, 'YQ': tag_qual} #YB/YQ tags for hard clips at beginnning of record

            while cigars[0][0] == 'H': #remove all hard clipping cigar elements
                cigars = cigars[1:]
            cigars = [('H', clipping_length)] + cigars #re-add with correct length
            assert not any([x == 'H' for x,y in cigars[1:]]), 'Hard clipping found, but not at end of read'
        else:
            assert seq[:len(read.seq)] == read.seq, f"Sequence mismatch BAM-FASTQ in fragment {read.qname}: \n{seq[:len(read.seq)]} != \n{read.seq}"
            check_qual(read.qname, qual[:len(read.seq)], read.qual)
            while cigars[-1][0] == 'H': #remove all hard clipping cigar elements
                cigars.pop()
            cigars.append(('H', clipping_length)) #re-add with correct length
            assert not any([x == 'H' for x,y in cigars[:-1]]), 'Hard clipping found, but not at end of read'
            new_tags = {'ZB': tag_seq, 'ZQ': tag_qual}  #ZB/ZQ tags for hard clips at end of record

        ncigar = join_cigar(cigars)
    else:
        stats['restored_unaligned_reads'] = stats.get('restored_unaligned_reads',0) + 1
        assert seq[:len(read.seq)] == read.seq, f"Sequence mismatch BAM-FASTQ in fragment {read.qname}"
        check_qual(read.qname, qual[:len(read.seq)],read.qual)
        new_tags = {'ZB': tag_seq, 'ZQ': tag_qual}
        ncigar = '*'

    return (ncigar, new_tags, len(tag_seq))                

def adapt_supplementary(read, length, stats):
    stats['supplementary_alignments'] = stats.get('supplementary_alignments',0) + 1
    stats['supplementary_aligned_bp'] = stats.get('supplementary_aligned_bp',0) + read.get_aligned_read_length()
    if length == 0:
        return

    stats['supplementary_alignments_cigar_adapted'] = stats.get('supplementary_alignments_cigar_adapted',0) + 1
    cigar = read.split_cigar(orig_orientation=True)
    
    if cigar[-1][0] == 'H':
        clip_length = cigar[-1][1]
        cigar = cigar[:-1]
    else:
        clip_length = 0
    cigar.append(('H', clip_length + length))


    if read.is_reversed():
        cigar = cigar[::-1]
    read.cigar = join_cigar(cigar)
    read.tags['XT'] = length

def process(querygroup, fastqrecord, stats):
    name, seq1, seq2, qual1, qual2 = fastqrecord
    
    stats['alignments'] = stats.get('alignments',0) + len(querygroup)
    stats['fragments'] = stats.get('fragments',0) + 1

    stats['total_bp'] = stats.get('total_bp',0) + len(seq1) + len(seq2)
    
    reads1 = process_readgroup([row for row in querygroup if row.flag & 0x40])
    reads2 = process_readgroup([row for row in querygroup if row.flag & 0x80])


    r1_ncigar, r1_tags, r1_len = derive_missing_sequence_tags(reads1['primary'], seq1, qual1, stats)
    reads1['primary'].cigar = r1_ncigar
    reads1['primary'].setTagValues(**r1_tags)
    if r1_len:
        reads1['primary'].tags['XT'] = r1_len

    r2_ncigar, r2_tags, r2_len = derive_missing_sequence_tags(reads2['primary'], seq2, qual2, stats)
    reads2['primary'].cigar = r2_ncigar
    reads2['primary'].setTagValues(**r2_tags)
    if r2_len:
        reads2['primary'].tags['XT'] = r2_len
   
    for r in reads1.get('supplementary',[]):
        adapt_supplementary(r, r1_len, stats)
    for r in reads2.get('supplementary',[]):
        adapt_supplementary(r, r2_len, stats)

    rg = querygroup[0].getTagValue('RG')

    badly_mapped = False
    #determine if record is badly mapped
    #1) both reads are unmapped
    #2) or one read is unmapped, and th other has a bad alignment score (AS tag) < 50 and < 0.5 * max(AS)
    #3) both reads have a bad alignment score (AS tag) < 50 and < 0.5 * max(AS)
    ascore1, max_ascore1 = calculate_alignment_score(reads1)
    ascore2, max_ascore2 = calculate_alignment_score(reads2)
    if(ascore1 < 50 and ascore2 < 50 and ascore1 < 0.5 * max_ascore1 and ascore2 < 0.5 * max_ascore2):
        badly_mapped = True
        if((ascore1 + ascore2) == 0): #fully unmapped
            stats['fully_unmapped_fragments'] = stats.get('fully_unmapped_fragments',0) + 1
        else:
            stats['badly_mapped_fragments'] = stats.get('badly_mapped_fragments',0) + 1


    if aligned_to_chrom(reads1, 'chrEBV') or aligned_to_chrom(reads2,'chrEBV'):
        stats['chrEBV_mapped_fragments'] = stats.get('chrEBV_mapped_fragments',0) + 1
        badly_mapped = True

    if badly_mapped:
        rgx = rg
        if rgx is not None:
            rgx = rgx.replace(' ','_').replace('@','_').replace(':','_')

        badly_mapped_fastq = (reads1['primary'].toFastqRecord(restore_seq=False, readgroup_name=rgx), reads2['primary'].toFastqRecord(restore_seq=False, readgroup_name=rgx))
    else:
        badly_mapped_fastq = None

    return (querygroup,rg, badly_mapped_fastq)



def fastqtosam(fastqrecord, rg, stats):
    name, seq1, seq2, qual1, qual2 = fastqrecord
    stats['restored_bp_read1'] = stats.get('restored_bp_read1', 0) + len(seq1)
    stats['restored_bp_read2'] = stats.get('restored_bp_read2', 0) + len(seq2)
    stats['readded_fragments'] = stats.get('readded_fragments',0) + 1


    stats['size_r1'][len(seq1)] = stats['size_r1'].get(len(seq1),0) + 1
    stats['size_r2'][len(seq2)] = stats['size_r2'].get(len(seq2),0) + 1

    flag1 = 0x1 |  0x4 | 0x8| 0x40 | 0x200  #(PAIRED | UNMAP | MUNMAP | READ1 | QCFAIL)
    flag2 = 0x1 | 0x4 | 0x8| 0x80 | 0x200  #(PAIRED | UNMAP | MUNMAP | READ2 | QCFAIL)
    querygroup = []
    read1 = BamRecord([name, str(flag1), '*', '0', '0', '*', '*', '0', '0', seq1, qual1])
    read2 = BamRecord([name, str(flag2), '*', '0', '0', '*', '*', '0', '0', seq2, qual2])
    if not rg is None:
        read1.setTagValues(RG=rg)
        read2.setTagValues(RG=rg)
    return [read1,read2]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", help='first fastq file')
    parser.add_argument("-b", help='second fastq file')
    parser.add_argument("-o", default='-', help='output sam file (default stdout)')
    parser.add_argument("-ua", help='output gz fastq file 1 for badly aligned/unmapped fragments')
    parser.add_argument("-ub", help='output gz fastq file 1 for badly aligned/unmapped fragments')
    parser.add_argument("-s", help='stats file')
    args = parser.parse_args()
    
    freader = PairedFastQReader(args.a, args.b)
    freader.start()
    rg = None

    alignment_counter = 0
    stats = {}

    stats['size_r1'] = {}
    stats['size_r2'] = {}
    last_time = time.time()
    with sys.stdin as f:
        reader = csv.reader(f,delimiter='\t', quoting=csv.QUOTE_NONE)

        lastname = ''
        querygroup = []
        if args.o == '-':
            out = sys.stdout
        else:
            out = open(args.o, 'wt')
        
        out_unmapped_1 = gzip.open(args.ua,'wt')
        out_unmapped_2 = gzip.open(args.ub,'wt')
            
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
                        assert len(querygroup) >= 2, f'Reads are not paired for {querygroup[0]}'
                        fastq_record = freader.popRead(lastname)
                        res_kept, rg, badly_mapped_fastq = process(querygroup, fastq_record, stats)
                        res_kept = '\n'.join([xrow.toSamRecord() for xrow in res_kept]) + '\n'
                        out_file.write(res_kept)
                        if badly_mapped_fastq is not None:
                            out_unmapped_1.write(badly_mapped_fastq[0])
                            out_unmapped_2.write(badly_mapped_fastq[1])

                    lastname = row.qname
                    querygroup = [row]
            if querygroup:
                assert len(querygroup) >= 2, 'Reads are not paired'
                fastq_record = freader.popRead(lastname)
                #Process
                res_kept, rg, badly_mapped_fastq = process(querygroup, fastq_record, stats)

                res_kept = '\n'.join([xrow.toSamRecord() for xrow in res_kept]) + '\n'
                out_file.write(res_kept)
                if badly_mapped_fastq is not None:
                    out_unmapped_1.write(badly_mapped_fastq[0])
                    out_unmapped_2.write(badly_mapped_fastq[1])

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
                    out_file.write(res)
            except StopIteration:
                pass
       
        out_unmapped_1.close()
        out_unmapped_2.close()

    #process size stats
    size_stats = numpy.zeros((max(max(stats['size_r1'].keys()), max(stats['size_r2'].keys()))+ 1,2),dtype=int)
    for key,value in stats['size_r1'].items():
        size_stats[key,0] = value
    for key,value in stats['size_r2'].items():
        size_stats[key,1] = value
    del stats['size_r1']
    del stats['size_r2']

    stats['primary_soft_clipped_bp_ratio'] = float(stats.get('primary_soft_clipped_bp',0.0)) / float(stats['primary_aligned_bp'] + stats.get('primary_soft_clipped_bp',0.0))
    stats['supplementary_bp_ratio'] = float(stats.get('supplementary_aligned_bp',0.0)) / float(stats['primary_aligned_bp'] + stats.get('supplementary_aligned_bp',0.0))
    sys.stderr.write(f'Storing statistics in file {args.s}\n')
    with open(args.s,'wt') as fx:            
        for key, value in stats.items():
            sys.stderr.write('%s\t%s\n' % (key,str(value)))
            sys.stderr.flush()
            if isinstance(value, float):
                fx.write('%s\t %.5f\n' % (key,value))
            else:
                fx.write('%s\t %d\n' % (key,value))
        fx.write('\nadapter_length\tcount_read1\tcount_read2\n')
        for pos in range(1,size_stats.shape[0]):
            fx.write('%d\t%d\t%d\n' % (pos,size_stats[pos,0], size_stats[pos,1]))
        fx.flush()
    sys.stderr.write('done\n')
    sys.stderr.flush()
    time.sleep(3)    
   
    import os


    if os.stat(args.s).st_size == 0:
        sys.stderr.write('Stat file is empty\n')
        sys.stderr.flush()
        sys.exit(1)

    if os.stat(args.ua).st_size == 0:
        sys.stderr.write('Badmap FQ1 file is empty\n')
        sys.stderr.flush()
        sys.exit(1)

    if os.stat(args.ub).st_size == 0:
        sys.stderr.write('Badmap FQ2 file is empty\n')
        sys.stderr.flush()
        sys.exit(1)       
