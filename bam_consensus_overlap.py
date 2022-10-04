import argparse
import sys
import csv
import asyncio
import gzip
import os
import threading
import re
import numpy
import math
import time
from Bio.Seq import Seq

from bam_utils import *


def process_readgroup(readgroup):
    result = {}
    for row in readgroup:
        if row.flag & 0x100:
            if row.flag & 0x800:
                key = 'supplementary_secondary'
            else:
                key = 'secondary'
        elif row.flag & 0x800:
            key = 'supplementary'
        else:
            key = 'primary'

        if key == 'primary':
            assert not key in result
            result['primary'] = row
        else:
            w = result.get(key,[])
            w.append(row)
            result[key] = w


    #analysis flags


    return result



def get_position_order(readgroup):
    #get order of alignments in terms of (unrevcomped) read sequence. So last alignment relates to last bases read by machine for that read

    #if remaining part of cigar ends in an insert/delete, it is pruned completely
    #'cigars': cigars in unrevcomped read order: i.e. last alignment action is cigars[-1][-1], and first alignment action is cigars[0][0]. 'H'ardclip elements are removed.

    scigars = [readgroup['primary'].split_cigar()] + [r.split_cigar() for r in readgroup.get('supplementary',[])]
    reverse = [readgroup['primary'].flag & 0x10] + [r.flag & 0x10 for r in readgroup.get('supplementary',[])]
    seq_length = readgroup['primary'].get_orig_read_length()

    read_startpos = []  
    read_endpos = []
    align_length = []
    cigars = []
    for scigar,rev in zip(scigars,reverse):
        if len(scigar) == 0:
            continue
        pos = 0
        alength = 0
        minpos = []
        maxpos = []
        #numpy.zeros(len(sequse),dtype=bool)
        for ctype, length in scigar:
            if ctype == 'S' or ctype == 'H':
                pos += length
            elif ctype == 'M' or ctype == 'I' or ctype == '=' or ctype == 'X':
                if rev:
                    minpos.append(seq_length - (pos + length))
                    maxpos.append(seq_length - pos)
                else:
                    minpos.append(pos)
                    maxpos.append(pos + length - 1)
                pos += length
                
                if ctype != 'I':
                    alength += length
                    
            elif ctype == 'D' or ctype == 'N' or ctype == 'P':
                if ctype == 'D' or ctype == 'N':
                    alength += length
                pass
        read_startpos.append(min(minpos))
        read_endpos.append(max(maxpos))
        align_length.append(alength)
        if rev:
            cigars.append(scigar[::-1])
        else:
            cigars.append(scigar)
    chrom = [readgroup['primary'].rname] + [e.rname for e in readgroup.get('supplementary',[])]
    pos = [int(readgroup['primary'].pos)] + [int(e.pos) for e in readgroup.get('supplementary',[])]
    seq = [readgroup['primary'].seq] + [e.seq for e in readgroup.get('supplementary',[])]

    #better handle seq: either merge or use hard-clip offset
    sortidx = sorted(range(len(read_startpos)), key=read_startpos.__getitem__)
    chrom_sorted = [chrom[e] for e in sortidx]
    pos_sorted = [pos[e] for e in sortidx]
    cigars_sorted = [cigars[e] for e in sortidx]
    cigars_sorted_without_hardclip = [[e for e in cigar if not e[0] == 'H'] for cigar in cigars_sorted]
    readstartpos_sorted = [read_startpos[e] for e in sortidx]
    readendpos_sorted = [read_endpos[e] for e in sortidx]
    alignlength_sorted = [align_length[e] for e in sortidx]
    endpos_sorted = [pos[e] + align_length[e] for e in sortidx]
    reverse_sorted = [reverse[e] for e in sortidx]
    seq_sorted = [seq[e] for e in sortidx]
    return {'idx':sortidx, 'chrom':chrom_sorted, 'pos':pos_sorted, 'endpos_sorted':endpos_sorted, 'seq': seq_sorted, 'readpos_start':readstartpos_sorted, 'readpos_end':readendpos_sorted, 'cigars': cigars_sorted_without_hardclip, 'align_length':alignlength_sorted, 'reverse':reverse_sorted}


def check_left_over_adapters(querygroup, stats):
    #read1 = process_readgroup([row for row in querygroup if (row.flag & 0x40)])
    #read2 = process_readgroup([row for row in querygroup if (row.flag & 0x80)])
    
    #is_unmapped1 = read1['primary'].is_unmapped()
    #is_unmapped2 = read2['primary'].is_unmapped()
    
    #pos_order1 = get_position_order(read1)
    #pos_order2 = get_position_order(read2)
    
    
    nquerygroup = querygroup
    if len(querygroup) == 2: #just primary alignments
        read1, read2 = querygroup
        if read2.flag & 0x40:
            read2,read1 = read1,read2

        #check if reads are overlapping, and map relatively uniquely, if so merge.
        if read1.rname == read2.rname and read1.rname != '*' and int(read1.mapq) >= 40 and int(read2.mapq) >= 40:
            rp1 = read1.get_read_position(orig_orientation=True)
            rp2 = read2.get_read_position(orig_orientation=True)

            overlap, ostart, ostop, inward = read_overlap(rp1,rp2)

            if overlap:
                ostart1, ostop1, start1, stop1, fstart1, fstop1, pcigar1, mcigar1, acigar1 = read1.ref_to_read_pos(ostart,ostop)
                ostart2, ostop2, start2, stop2, fstart2, fstop2, pcigar2, mcigar2, acigar2 = read2.ref_to_read_pos(ostart,ostop)

                if inward:
                    stats['inward_overlap'][ostop - ostart] += 1
                    if not (rp1['seq_length'] != stop1 or rp2['seq_length'] != stop2 or mcigar1 != mcigar2): #easy direct match
                        stats['inward_merged'] = stats.get('inward_merged',0) + 1
                        sq, mq, mismatch, total = generate_simple_consensus(read1.seq[ostart1:ostop1], read2.seq[ostart2:ostop2], read1.qual[ostart1:ostop1], read2.qual[ostart2:ostop2])
                        if read1.is_reversed():
                            newseq = read2.seq[:ostart2] + sq + read1.seq[ostop1:] 
                            newqual = read2.qual[:ostart2] + mq + read1.qual[ostop1:] 
                            assert acigar2 == ''
                            assert pcigar1 == ''
                            ncigar = pcigar2 + mcigar2 + acigar1

                            nflag = 0x2 | 0x40 | 0x10
                            pos = read2.pos
                        else:
                            newseq = read1.seq[:ostart1] + sq + read2.seq[ostop2:] 
                            newqual = read1.qual[:ostart2] + mq + read2.qual[ostop2:] 
                            assert acigar1 == ''
                            assert pcigar2 == ''
                            ncigar = pcigar1 + mcigar2 + acigar2

                            nflag = 0x2 | 0x40
                            pos = read1.pos
                        ntags = {'RG': read2.tags['RG']}
                        ntags['AS'] = read1.tags['AS'] + read2.tags['AS'] - 2 * len(sq)
                        if 'XZ' in read1.tags or 'XZ' in read2.tags:
                            ntags['XZ'] = max(read1.tags.get('XZ', int(read1.mapq)), read2.tags.get('XZ', int(read2.mapq)))

                        record = BamRecord([read1.qname, nflag, read2.rname, pos, max(int(read1.mapq), int(read2.mapq)), ncigar, '*', 0, 0, newseq,newqual])
                        record.tags = ntags
                        nquerygroup = [record]

                            #if total > 10 and mismatch / float(total) > 0.25:
                            #    print('INWARD', ostop-ostart, ostop, ostart, rp1, rp2, start1, stop1, start2, stop2, read1.cigar, read2.cigar, 'Y',pcigar1, mcigar1, acigar1, 'X',pcigar2,mcigar2, acigar2)
                            #    print(read1.mapq, read2.mapq)
                            #    print(read1.tags, read2.tags)
                            #    print((start2 - (len(read1.seq) - stop1)) * " " + read1.qual)
                            #    print((start2 - (len(read1.seq) - stop1)) * " " + read1.seq)
                            #    print(' ' * start2 + '|' * (stop2 - start2))
                            #    print(read2.seq)
                            #    print(read2.qual)




                            #print(read1.mapq, read2.mapq)
                            #print(read1.qual)
                            #print(read1.seq)
                            #print(' ' * start1 + '|' * (stop1 - start1))
                            #print((start1 - (len(read2.seq) - stop2)) * " " + read2.seq)
                            #print((start1 - (len(read2.seq) - stop2)) * " " + read2.qual)
                            
                    elif mcigar1 != mcigar2:
                        stats['inward_skipped'] = stats.get('inward_skipped',0) + 1
                        pass
                    else:
                        stats['inward_skipped'] = stats.get('inward_skipped',0) + 1
                        pass
                else:   
                    stats['outward_overlap'][ostop - ostart] += 1
                    stats['outward_skipped'] = stats.get('outward_skipped',0) + 1
                    if rp1['ref_spos'] < ostart or rp2['ref_spos'] < ostart or rp1['ref_epos'] > ostop or rp2['ref_epos'] > ostop:
                        #
                        pass
                        #print('OUTWARD', ostop-ostart, ostop, ostart, rp1, rp2, start1, stop1, start2, stop2)

    return nquerygroup
   

def process(querygroup, stats):
    
    stats['alignment_counter'] = stats.get('alignment_counter',0) + len(querygroup)
    stats['fragment_counter'] = stats.get('fragment_counter',0) + 1
    
    nquerygroup = check_left_over_adapters(querygroup, stats)

    return querygroup


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", default='-', help='output sam file (default stdout)')
    parser.add_argument("-s", help='stats file')
    args = parser.parse_args()
    
    alignment_counter = 0
    stats = {}
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
                        res_kept, rg = process(querygroup, stats)
                        res_kept = '\n'.join([xrow.toSamRecord() for xrow in res_kept]) + '\n'
                        out_file.write(res_kept)

                    lastname = row.qname
                    querygroup = [row]
            if querygroup:
                assert len(querygroup) >= 2, 'Reads are not paired'
                #Process
                res_kept, rg = process(querygroup, stats)

                res_kept = '\n'.join([xrow.toSamRecord() for xrow in res_kept]) + '\n'
                out_file.write(res_kept)

        with open(args.s,'wt') as f:            
            for key, value in stats.items():
                f.write('%s\t %d\n' % (key,value))
            f.flush()

