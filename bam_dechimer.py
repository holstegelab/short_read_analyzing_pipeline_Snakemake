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
import subprocess

from Bio.Seq import Seq
from ibidas.utils import util
from bam_utils import *



def rg_prune(readgroup, stats, stats_prefix, args):

    #prunes split alignments such that overlapping parts between the alignments are pruned from both sides
    #if remaining part of cigar ends in an insert/delete, the insert/delete is pruned completely
    assert 'supplementary' in readgroup

    stats[f"{stats_prefix}_has_supplementary_alignments"] = stats.get(f"{stats_prefix}_has_supplementary_alignments",0) + 1
    readgroup = readgroup.copy()
    reads = [readgroup['primary']] + readgroup['supplementary']
    rpos = [r.get_read_position(orig_orientation=True) for r in reads]
  

    #determine prune positions
    pruned = False
    nreadgroup = readgroup.copy()
    nreadgroup['supplementary'] = []
    for pos in range(len(rpos)):
        start = rpos[pos]['read_spos']
        stop = rpos[pos]['read_epos']
        after = False
        before = False
        for pos2 in range(len(rpos)):
            if pos == pos2:
                continue
            s = rpos[pos2]['read_spos']
            e = rpos[pos2]['read_epos']

            if s >= stop:
                after = True
            if e <= start: #no overlap
                before = True
                continue

            if s <= start: #before
                before = True
                start = e
            if e >= stop: #after
                after = True
                stop = s
            if(s > start and e < stop): #within
                if (stop - e) > (s - start):
                    start = e
                else:
                    stop = s
        start = start if start <= stop else stop

        #prune
        read = reads[pos]
        if (stop - start) < args.min_align_length: 
            pruned = True
            if read.is_primary():
                nreadgroup['primary'] = read.unmap()
                stats[f"{stats_prefix}_primary_unmapped_in_pruning"] = stats.get(f"{stats_prefix}_primary_unmapped_in_pruning",0) + 1
            else:
                #drop supplementary alignments that need to be unmapped
                stats[f"{stats_prefix}_sup_discarded_in_pruning"] = stats.get(f"{stats_prefix}_sup_discarded_in_pruning",0) + 1
            continue 
            
        if rpos[pos]['read_spos'] != start or (start > 0 and before):
            pruned = True
            read = read.clip_start(start, orig_orientation=True)
            stats[f"{stats_prefix}_start_pruned"] = stats.get(f"{stats_prefix}_start_pruned",0) + 1

        if rpos[pos]['read_epos'] != stop or (stop > 0 and after):
            pruned = True
            read = read.clip_end(stop, orig_orientation=True)
            stats[f"{stats_prefix}_end_pruned"] = stats.get(f"{stats_prefix}_end_pruned",0) + 1


        if read.is_primary():
            nreadgroup['primary'] = read
        else:
            nreadgroup['supplementary'].append(read)
            
    if not nreadgroup['supplementary']:
        del nreadgroup['supplementary']
        stats[f"{stats_prefix}_all_sup_discarded_in_pruning"] = stats.get(f"{stats_prefix}_all_sup_discarded_in_pruning",0) + 1


    return (nreadgroup,pruned)


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


def dechimer(reads1, reads2, stats, args):
    r1 = reads1['primary']
    r2 = reads2['primary']
    modified1 = modified2 = False    
    if not 'supplementary' in reads1 and not 'supplementary' in reads2: #only two alignment records
        if not 'S' in r1.cigar and not 'S' in r2.cigar:
            return (modified1,modified2) # fast path

        cigar1 = r1.split_cigar(orig_orientation=True, merge_clips=True)
        cigar2 = r2.split_cigar(orig_orientation=True, merge_clips=True)
        if not r1.is_unmapped() and not r2.is_unmapped():
            chrom1 = r1.rname
            chrom2 = r2.rname
            same_chrom = chrom1 == chrom2
            if same_chrom:
                pos1 = int(r1.pos)
                pos2 = int(r2.pos)
                diff = abs(pos1 - pos2)
            else:
                diff = 100000000000000
           
            #only primary reads --> distance > MAX_READ_DIST --> remove clipping
            if diff >= args.max_read_dist:
                if cigar1[-1][0] == 'SH':
                    r1 = r1.clip_end(cigar1[-1][1], orig_orientation=True)
                    stats[f'read1_dechimer_clip'] = stats.get('read1_dechimer_clip',0) + 1
                    modified1=True

                if cigar2[-1][0] == 'SH':
                    r2 = r2.clip_end(cigar2[-1][1], orig_orientation=True)
                    stats[f'read2_dechimer_clip'] = stats.get('read2_dechimer_clip',0) + 1
                    modified2 = True


            if args.loose_ends:
                if cigar1[0][0] == 'SH':
                    r1 = r1.clip_start(cigar1[0][1], orig_orientation=True)
                    stats[f'read1_loose_end_clip'] = stats.get('read1_loose_end_clip',0) + 1
                    modified1=True
                if cigar2[0][0] == 'SH':
                    r2 = r2.clip_start(cigar2[0][1], orig_orientation=True)
                    stats[f'read2_loose_end_clip'] = stats.get('read2_loose_end_clip',0) + 1
                    modified2 = True

        elif r1.is_unmapped() and r2.is_unmapped():
            pass
        elif r1.is_unmapped(): 
            if cigar2[-1][0] == 'SH':
                r2 = r2.clip_end(cigar2[-1][1], orig_orientation=True)
                stats[f'read2_dechimer_clip'] = stats.get('read2_dechimer_clip',0) + 1
                modified2 = True
            if args.loose_ends and cigar2[0][0] == 'SH':
                r2 = r2.clip_start(cigar2[0][1], orig_orientation=True)
                stats[f'read2_loose_end_clip'] = stats.get('read2_loose_end_clip',0) + 1
                modified2 = True
        else: #is_unmapped2 
            if cigar1[-1][0] == 'SH':
                r1 = r1.clip_end(cigar1[-1][1], orig_orientation=True)
                stats[f'read1_dechimer_clip'] = stats.get('read1_dechimer_clip',0) + 1
                modified1 = True
            if args.loose_ends and cigar1[0][0] == 'SH':
                r1 = r1.clip_start(cigar1[0][1], orig_orientation=True)
                stats[f'read1_loose_end_clip'] = stats.get('read1_loose_end_clip',0) + 1
                modified1 = True
        reads1['primary'] = r1
        reads2['primary'] = r2
    else:
        if False and  not r1.is_unmapped() and not r2.is_unmapped():
            cigar1 = r1.split_cigar(orig_orientation=True, merge_clips=True)
            cigar2 = r2.split_cigar(orig_orientation=True, merge_clips=True)

            if 'supplementary' in reads1 and not ('supplementary' in reads2) and len(reads1['supplementary']) == 1:

                
                if cigar2[-1][0] == 'SH':
                    r2 = r2.clip_end(cigar2[-1][1], orig_orientation=True)
                    stats[f'read2_dechimer_clip'] = stats.get('read2_dechimer_clip',0) + 1
                    modified2 = True

                if args.loose_ends and cigar2[0][0] == 'SH':
                    r2 = r2.clip_start(cigar2[0][1], orig_orientation=True)
                    stats[f'read2_loose_end_clip'] = stats.get('read2_loose_end_clip',0) + 1
                    modified2 = True

                #FIXME: more detailed analysis based on overlap with read1 possible
                #positions = [e.get_read_position(orig_orientation=True) for e in [reads1['primary']] + reads['supplementary']]

            if 'supplementary' in reads2 and not ('supplementary' in reads1) and len(reads2['supplementary']) == 1:
                if cigar1[-1][0] == 'SH':
                    r1 = r1.clip_end(cigar1[-1][1], orig_orientation=True)
                    stats[f'read1_dechimer_clip'] = stats.get('read1_dechimer_clip',0) + 1
                    modified1 = True

                if args.loose_end and cigar1[0][0] == 'SH':
                    r1 = r1.clip_start(cigar1[0][0], orig_orientation=True)
                    stats[f'read1_loose_end_clip'] = stats.get('read1_loose_end_clip',0) + 1
                    modified1 = True
                #FIXME: more detailed analysis based on overlap with read1 possible
                    

    return (modified1, modified2)

def record_filter(read, stats, stats_prefix, args):
    if len(read.seq) <= args.min_align_length or read.get_align_length() < args.min_align_length:
        stats[f'{stats_prefix}_min_align_length_unmap'] = stats.get(f'{stats_prefix}_min_align_length_unmap', 0) + 1
        return read.unmap()
    else:
        return read

def rg_filter(readgroup, stats, stats_prefix, args):
    readgroup['primary'] = record_filter(readgroup['primary'], stats, stats_prefix, args)

    if 'supplementary' in readgroup:
        new = [record_filter(e, stats, stats_prefix, args) for e in readgroup['supplementary']]
        new = [e for e in new if not e.is_unmapped()]

        if readgroup['primary'].is_unmapped() and len(new) > 0:
            new_primary = new.pop(0)
            stats[f'{stats_prefix}_promote_supplementary'] = stats.get(f'{stats_prefix}_promote_supplementary', 0) + 1
            new_primary = new_primary.annotate_orig_sequence(readgroup['primary'].orig_seq(orig_orientation=True),
                                                             readgroup['primary'].orig_qual(orig_orientation=True))
            new_primary.flag = new_primary.flag & (~0x800) #update flags
            readgroup['primary'] = new_primary

        if len(new) > 0:
            readgroup['supplementary'] = new
        else:
            del readgroup['supplementary']
    return readgroup           



def update_sa_tag(readgroup):
    if 'supplementary' in readgroup:
        sa_primary = readgroup['primary'].gen_sa_tag()
        sa_sups = [read.gen_sa_tag() for read in readgroup['supplementary']]
        
        if sa_sups:
            readgroup['primary'].tags['SA'] = ''.join(sa_sups)

        for pos, read in enumerate(readgroup['supplementary']):
            sa_sups_tmp = list(sa_sups)
            del sa_sups_tmp[pos]
            read.tags['SA'] = sa_primary + ''.join(sa_sups_tmp)

def update_mate_flags(reads1,reads2):
    p1 = reads1['primary']
    p2 = reads2['primary']
    

    r1_flags = 0
    r2_flags = 0

    if p1.is_unmapped():
        r2_flags |= 0x8
    if p2.is_unmapped():
        r1_flags |= 0x8
    if p1.is_reversed():
        r2_flags |= 0x20
    if p2.is_reversed():
        r1_flags |= 0x20

    if not p1.is_unmapped() and not p2.is_unmapped():
        r1_flags |= 0x2 #PROPERLY_PAIRED
        r2_flags |= 0x2 #PROPERLY_PAIRED

    cancel_flag = ~(0x8 | 0x20 | 0x2)

    reads1['primary'].flag = (reads1['primary'].flag & cancel_flag) | r1_flags
    reads2['primary'].flag = (reads2['primary'].flag & cancel_flag) | r2_flags
    for r in reads1.get('supplementary',[]):
        r.flag = (r.flag & cancel_flag) | r1_flags
    for r in reads2.get('supplementary',[]):
        r.flag = (r.flag & cancel_flag) | r2_flags

def update_pair_flag(reads1, reads2, stats, args):
    p1 = reads1['primary']
    p2 = reads2['primary']

    proper = True
    if p1.rname != p2.rname or p1.rname == "*":
        proper = False
    else:
        diff = abs(int(p1.pos) - int(p2.pos))
        if diff > args.max_read_dist:
            proper = False
    cancel_flag = ~0x2
    new_flag = 0x2 if proper else 0

    if p1.is_proper() or p2.is_proper():
        stats['proper_pair_downgrade'] = stats.get('proper_pair_downgrade',0) + 1

    p1.flag = (p1.flag & cancel_flag) | new_flag
    p2.flag = (p2.flag & cancel_flag) | new_flag



def update_mate_tags(read1, read2):
    p1 = read1['primary']
    p2 = read2['primary']

    if p2.is_unmapped():
        p1.tags.pop("MC",None)
    else:
        p1.tags['MC'] = p2.cigar
    if p1.is_unmapped():
        p2.tags.pop('MC',None)
    else:
        p2.tags['MC'] = p1.cigar

    
    p1.rnext = p2.rname
    p2.rnext = p1.rname
    p1.pnext = p2.pos
    p2.pnext = p1.pos
    

    for r in read1.get('supplementary',[]):
        r.rnext = p2.rname
        r.pnext = p2.pos
        if p2.is_unmapped():
            r.tags.pop('MC',None)
        else:
            r.tags['MC'] = p2.cigar
    for r in read2.get('supplementary',[]):
        r.rnext = p1.rname
        r.pnext = p1.pos
        if p1.is_unmapped():
            r.tags.pop("MC",None)
        else:
            r.tags['MC'] = p1.cigar


def process(querygroup, stats, args):
    
    stats['alignment_counter'] = stats.get('alignment_counter',0) + len(querygroup)
    stats['fragment_counter'] = stats.get('fragment_counter',0) + 1
    
    reads1 = process_readgroup([row for row in querygroup if row.flag & 0x40])
    reads2 = process_readgroup([row for row in querygroup if row.flag & 0x80])

    #if reads1['primary'].qname == "A01641:36:HGJYVDSX2:1:1105:10800:13745":
    #    util.debug_here()

    modified1 = False
    modified2 = False

    if 'supplementary' in reads1:
        reads1, is_pruned = rg_prune(reads1, stats, 'read1', args)
        if is_pruned:
            modified1 = True

    if 'supplementary' in reads2:
        reads2, is_pruned = rg_prune(reads2, stats, 'read2', args)
        if is_pruned:
            modified2 = True

    rmodified1,rmodified2 = dechimer(reads1,reads2, stats, args)
    modified1 = modified1 or rmodified1
    modified2 = modified2 or rmodified2

    #remove reads that have <= 1bp of sequence
    if modified1:
        reads1 = rg_filter(reads1, stats, 'read1', args)
        update_sa_tag(reads1)
    if modified2:
        reads2 = rg_filter(reads2, stats, 'read2', args)
        update_sa_tag(reads2)

    if modified1 or modified2:
        stats['fragment_modified'] = stats.get('fragment_modified',0) + 1
        update_mate_flags(reads1,reads2)
        update_mate_tags(reads1,reads2)

    if modified1 or modified2:
        update_pair_flag(reads1,reads2, stats,args)        

    if modified1 and not args.d:
        oldreads1 = process_readgroup([row for row in querygroup if (row.flag & 0x40)])
        validate(reads1, oldreads1)
    if modified2 and not args.d:
        oldreads2 = process_readgroup([row for row in querygroup if (row.flag & 0x80)])
        validate(reads2, oldreads2)

    querygroup1 = [reads1.get('primary', None)] + reads1.get('supplementary',[]) + reads1.get('secondary',[]) + reads1.get('supplementary_secondary',[])
    querygroup2 = [reads2.get('primary', None)] + reads2.get('supplementary',[]) + reads2.get('secondary',[]) + reads2.get('supplementary_secondary',[])
    querygroup = [e for e in  (querygroup1 + querygroup2) if not e is None]
    querygroup = [e for e in querygroup if not e.seq == '']


    return querygroup

def get_pos_pairs(splitted_cigar, refpos, sequence, quality):
    pairs = []

    curpos = refpos
    for stype, length in splitted_cigar:
        if stype == 'S':
            sequence = sequence[length:]
            quality = quality[length:]
        elif stype == 'M' or stype == '=' or stype == 'X':
            for i in range(length):
                pairs.append((curpos + i, sequence[i], quality[i]))
            curpos += length
            sequence = sequence[length:]
            quality = quality[length:]
        elif stype == 'I':
            for i in range(length):
                pairs.append((curpos + 0.5, sequence[i], quality[i]))
            sequence = sequence[length:]
            quality = quality[length:]
        elif stype == 'D' or stype == 'N':
            curpos += length
        elif stype == 'H' or stype == 'P':
            pass
    return pairs

def validate(readgroup, old_readgroup):
    assert 'primary' in readgroup

    reads = [readgroup['primary']] + readgroup.get('supplementary',[]) + readgroup.get('secondary',[]) + readgroup.get('supplementary_secondary',[])

    oldreads = [old_readgroup['primary']] + old_readgroup.get('supplementary',[]) + old_readgroup.get('secondary',[]) + old_readgroup.get('supplementary_secondary',[])
    oldpairs = []
    for oread in oldreads:
        cigar = oread.split_cigar()
        if cigar:
           oldpairs.extend(get_pos_pairs(cigar, int(oread.pos), oread.seq, oread.qual))
    oldpairs = set(oldpairs)
    
    had_errors = False

    #get position seq/postion ref pairs old reads
    errors = []
    for read in reads:
        errors = []
        scigar = read.split_cigar()
        if not scigar:
            if int(read.pos) != 0:
                errors.append('- Unmapped reads should  have no position')
            if not read.is_unmapped():                
                errors.append('- Mapped read without cigar')
        else: 
            if int(read.pos) == 0:
                errors.append('- Mapped reads should  have no 0 position')
            if read.is_unmapped():
                errors.append('- Unmapped read with cigar')
            #length check
            read_length = 0
            for stype, length in scigar:
                if stype == 'S' or stype == 'M' or stype == 'X' or stype == '=' or stype == 'I':
                    read_length += length
                if length == 0:
                    errors.append('- Length 0 in cigar')    

            if read_length == 0:
                errors.append('- Unmapped read without cigar = *')
                
            if read_length != len(read.seq):
                errors.append('- Length of cigar (%d) does not match length of sequence (%d)' % (read_length, len(read.seq)))

            if read_length != len(read.qual):
                errors.append('- Length of cigar (%d) does not match length of qualityscores (%d)' % (read_length, len(read.qual)))



            #check clipping
            tcigar = list(scigar)
            if tcigar and tcigar[0][0] == 'H':
                tcigar = tcigar[1:]
            if tcigar and tcigar[-1][0] == 'H':
                tcigar = tcigar[:-1]
            if tcigar and tcigar[0][0] == 'S':
                tcigar = tcigar[1:]
            if tcigar and tcigar[-1][0] == 'S':
                tcigar = tcigar[:-1]
            if not tcigar:
                errors.append(' - Cigar only has clipping: %s' % read.cigar)
            internal_clip = len([e for e in tcigar if e[0] in ['H', 'S']])
            if internal_clip:
                errors.append('- Cigar has internal clipping: %s' % read.cigar)

            
            tcigar = [e for e in tcigar if e[0] in ['M','=','X']]
            if len(tcigar) == 0:
                errors.append('- Cigar has no mapping: %s' % read.cigar)

            last_stype = ''
            for stype, length in scigar:
                if stype == last_stype:
                    errors.append('- Cigar has repeated record: %s' % read.cigar)
                last_stype = stype


            pairs = get_pos_pairs(scigar, int(read.pos), read.seq, read.qual)

            for pair in pairs:
                if pair not in oldpairs:
                    opairs = list(oldpairs)
                    opairs.sort()
                    errors.append('- New mapping relations found: %s \n%s\n%s\n' % (str(pair), str(pairs), str(opairs)))
                    break


        #checking mate/other tags?

        if errors:
            sys.stderr.write('Error found in %s:%s\n' % (read.qname, read.pos))
            sys.stderr.write('\n'.join(errors) + '\n')
            sys.stderr.flush()
            had_errors = True

    if had_errors:
        sys.stderr.write('Old read records:\n')
        for read in oldreads:
            sys.stderr.write('- ' + str(read) + '\n')
        sys.stderr.write('New read records:\n')
        for read in reads:
            sys.stderr.write('- ' + str(read) + '\n')
        sys.stderr.flush()

        sys.exit(1)
    return True
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", default='-', help='input sam file (default stdin)')
    parser.add_argument("-o", default='-', help='output sam file (default stdout)')
    parser.add_argument("-s", help='stats file')
    parser.add_argument("-d", action='store_true', default=False, help='disable validation')
    parser.add_argument("--loose_ends", action='store_true', default=False, help='Remove loose soft clipping ends at the outside of correctly paired alignments')
    parser.add_argument("--min_align_length", default=30, help='Minimum length below which alignments are discarded (default: 30)')
    parser.add_argument("--max_read_dist", default=100000, help='Maximum read distance to consider reads still regularly paired (default:100000bp)')
    args = parser.parse_args()
   

    if args.d:
        sys.stderr.write('Validation of modified records disabled\n')
        sys.stderr.flush()

    alignment_counter = 0
    stats = {}
    
    last_time = time.time()

    if args.i == '-':
        pipe_in = sys.stdin
    else:
        in_process = subprocess.Popen(['samtools','view','-h', args.i],stdout=subprocess.PIPE, universal_newlines=True)
        pipe_in = in_process.stdout

    with pipe_in as f:
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
                    sys.stderr.write(str(stats) + '\n')
                    sys.stderr.write('%d reads done at %d reads/sec\n' % (alignment_counter, int(100000.0 / (time.time() - last_time))))

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
                        assert len(querygroup) >= 2, f'Reads {[e.qname for e in querygroup]} are not paired'
                        res_kept = process(querygroup, stats, args)
                        res_kept = '\n'.join([xrow.toSamRecord() for xrow in res_kept]) + '\n'
                        out_file.write(res_kept)

                    lastname = row.qname
                    querygroup = [row]
            if querygroup:
                assert len(querygroup) >= 2, 'Reads are not paired'
                #Process
                res_kept = process(querygroup, stats, args)

                res_kept = '\n'.join([xrow.toSamRecord() for xrow in res_kept]) + '\n'
                out_file.write(res_kept)
        
        sys.stderr.write('Writing statistics\n')
        sys.stderr.flush()
        
        with open(args.s,'wt') as f:            
            for key, value in stats.items():
                f.write('%s\t %d\n' % (key,value))
            f.flush()

        sys.stderr.write('Done\n')
        sys.stderr.flush()

