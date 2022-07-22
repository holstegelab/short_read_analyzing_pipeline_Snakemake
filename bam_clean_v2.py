#!/bin/env python
import re
import csv
import sys
import gzip
import swalign

aligner = swalign.LocalAlignment(swalign.NucleotideScoringMatrix(1,-4), gap_penalty=-6)
adapter = 'AGATCGGAAGAGC'

MAX_READ_DIST=100000

looseend_counter = 0
adapter_counter = 0
pruned_counter = 0
clipped_counter = 0
modified_counter = 0
read_counter = 0
check_counter = 0
alignment_counter = 0
supclip_counter = 0

def cigar_split(cigar):
    cigar = re.split('([MIDNSHP=X]{1,1})',cigar)
    return [(ctype, int(length)) for length, ctype in zip(cigar[::2], cigar[1::2])]


def process_readgroup(readgroup):
    result = {}
    for row in readgroup:
        flag = int(row[1])
        if flag & 0x100:
            if flag & 0x800:
                key = 'supplementary_secondary'
            else:
                key = 'secondary'
        elif flag & 0x800:
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

def record_hardclip(row, only_last=False, only_first=False):
#change soft clips to hard clips
    cigar = row[5]
    reverse = int(row[1]) & 0x10
    seq = row[9]
    qual = row[10]
    tags = row[11:]
    splitted_cigar = cigar_split(cigar)

    pos = 0
    newqual = []
    newseq = []
    newcigar = []


    if not splitted_cigar:
        return (list(row), False)


    has_hardclip_begin = splitted_cigar[0][0] == 'H'
    has_hardclip_end = splitted_cigar[-1][0] == 'H'

    clipped = False
    for counter, (ctype, length) in enumerate(splitted_cigar):
        if ctype == 'S':
            if (not only_last and not only_first):
                r = ('H', length)
                clipped = True
            elif only_last and ((not reverse and counter == (len(splitted_cigar) - 1 - int(has_hardclip_end))) or (reverse and counter == int(has_hardclip_begin))):
                r = ('H', length)
                clipped = True
            elif only_first and ((reverse and counter == (len(splitted_cigar) - 1 - int(has_hardclip_end))) or (not reverse and counter == int(has_hardclip_begin))):
                r = ('H', length)
                clipped = True
            else:
                r = ('S', length)
                newseq.append(seq[pos:pos + length])
                newqual.append(qual[pos:pos + length])
            pos += length
        elif ctype == 'M' or ctype == 'I' or ctype == '=' or ctype == 'X':
            r = (ctype, length)
            newseq.append(seq[pos:pos + length])
            newqual.append(qual[pos:pos + length])
            pos += length
        elif ctype == 'H':
            
            r = (ctype, length)
            pass
        elif ctype == 'D' or ctype == 'N' or ctype == 'P':
            r = (ctype, length)
            pass
        if len(newcigar) > 0 and r[0] == newcigar[-1][0]:
            newcigar[-1] = (r[0], newcigar[-1][1] + r[1])
        else:
            newcigar.append(r)
    if pos != len(seq):
        sys.stderr.write(str((splitted_cigar, len(seq), len(qual), pos, row)) + '\n')
        sys.stderr.flush()
    assert pos == len(seq)

    newqual = ''.join(newqual)
    newseq = ''.join(newseq)


    row = list(row)
    row[9] = newseq
    row[10] = newqual 
    row[5] = ''.join(['%d%s' % (length, ctype) for ctype, length in newcigar])

    return (row,clipped)


def record_lastclip(row, bp, aligned=False):
#remove [bp] base pairs from end or record:
#- [aligned = false] clip read sequence 
#- [aligned = true] clip aligned region

    if bp <= 0:
        return row

    oldrow = list(row)

    unmapped = False
    cigar = row[5]
    reverse = int(row[1]) & 0x10
    seq = row[9]
    qual = row[10]
    tags = row[11:]
    splitted_cigar = cigar_split(cigar)

    if reverse == 0:
        splitted_cigar = splitted_cigar[::-1]
        qual = qual[::-1]
        seq = seq[::-1] # no rev comp needed as we will undo reverse

    pos_shift = 0 #shift in alignment position
    actual_bp = 0 #actual base pairs pruned from record
    newcigar = []
    hardclip_removed_bp = 0
    while bp > 0:
        ctype, length = splitted_cigar.pop(0)
        if ctype == 'S' or ctype == 'M' or ctype == 'I' or ctype == '=' or ctype == 'X':
            
            if aligned and (ctype == 'S' or ctype == 'I'):
                prune = length
            else:
                prune = min(length, bp)
                bp -= prune

            seq = seq[prune:]
            qual = qual[prune:]
            actual_bp += prune
            if length > prune:
                newcigar.append('%d%s' % (length - prune, ctype))

            if ctype == 'M' or ctype == '=' or ctype == 'X':
                pos_shift += prune
        
        elif ctype == 'H' or ctype == 'D' or ctype == 'N' or ctype == 'P':
            if ctype == 'D' or ctype == 'N':
                if aligned:
                    bp -= length
                pos_shift += length
            elif ctype == 'H':
                hardclip_removed_bp += length


    if len([e for e in newcigar if not e[0] in ['S','H','D','P','I']]) == 0 and len([e for e in splitted_cigar if not e[0] in ['S','H', 'D','I','P']]) == 0:
        unmapped = True

    newcigar = ['%dH' % (actual_bp + hardclip_removed_bp)] + newcigar
    newcigar.extend(['%d%s' % (length, ctype) for ctype, length in splitted_cigar])

    if reverse == 0:
        newcigar = newcigar[::-1]
        seq = seq[::-1]
        qual = qual[::-1]

    newcigar = ''.join(newcigar)
    
    row = list(row)
    row[9] = seq
    row[10] = qual 
    if reverse:
        row[3] = str(int(row[3]) + pos_shift)

    if unmapped:
        newcigar = '*'
        row[3] = '0'
        row[9] = oldrow[9] #if unmapped, set back old sequence
        row[10] = oldrow[10] #if unmapped, set back old sequence
        

    row[5] = newcigar
    return row


def get_orig_read_length(cigar):
    pos = 0
    for (ctype, length) in cigar:
        if ctype == 'S' or ctype == 'M' or ctype == 'I' or ctype == '=' or ctype == 'X' or ctype == 'H':
            pos += length
    return pos 


def rg_hardclip(readgroup):
    readgroup = readgroup.copy()
    readgroup['primary'],clipped = record_hardclip(readgroup['primary'])
    if 'supplementary' in readgroup:
        sclipped = []
        for r in readgroup['supplementary']:
            s, is_clipped = record_hardclip(r)
            clipped = clipped or is_clipped
            sclipped.append(s)

        readgroup['supplementary'] = sclipped
    return (readgroup,clipped)


def rg_prune(readgroup):
    #prunes split alignments such that overlapping parts between the alignments are pruned from both sides
    #if remaining part of cigar ends in an insert/delete, it is pruned completely
    if 'is_pruned' in readgroup:
        return readgroup
    readgroup = readgroup.copy()
    
    cigars = [cigar_split(readgroup['primary'][5])] + [cigar_split(r[5]) for r in readgroup['supplementary']]
    reverse = [int(readgroup['primary'][1]) & 0x10] + [int(r[1]) & 0x10 for r in readgroup['supplementary']]
     
    #part of read that is in-use
    inuse = []  #filter per alignment
    sequse = [0] * get_orig_read_length(cigars[0])
    for scigar,rev in zip(cigars,reverse):
        pos = 0
        w = []
        #numpy.zeros(len(sequse),dtype=bool)
        for ctype, length in scigar:
            if ctype == 'S':
                pos += length
            elif ctype == 'M' or ctype == 'I' or ctype == '=' or ctype == 'X':
                if rev:
                    for p in range(pos, pos+length):
                        sequse[-(p + 1)] += 1
                else:
                    for p in range(pos, pos+length):
                        sequse[p]+=1
                w.extend(range(pos,  pos+length))
                pos += length
            elif ctype == 'H' or ctype == 'D' or ctype == 'N' or ctype == 'P':
                pass
        if rev:
            w = [pos - e - 1 for e in w[::-1]]
        inuse.append(w)


    
    #determine first and last part of alignment that is aligned more than once
    shift = []
    pruning_required=False
    for rev,iu in zip(reverse,inuse):
        usecount = [sequse[a] for a in iu]
        begin_shift = 0
        end_shift = 0
        drop = False


        while begin_shift < len(usecount) and usecount[begin_shift] > 1: 
            begin_shift += 1
        reverse_usecount = usecount[::-1]
        while end_shift < len(usecount) and reverse_usecount[end_shift] > 1: 
            end_shift += 1

        #drop if whole read sequence is double used
        drop = (usecount[begin_shift] > 1) and (reverse_usecount[end_shift] > 1)


        if begin_shift > 0 or end_shift > 0:
            pruning_required=True

        if rev:
            begin_shift,end_shift = end_shift,begin_shift
        shift.append((begin_shift, end_shift, drop))


    #calculate new cigars and new position per alignment
    newcigars = []
    newpos = []
   
    for (begin_shift,end_shift, drop),rev,scigar in zip(shift, reverse, cigars):
        if drop:
            newcigars.append('*')
            newpos.append(0)
            continue

        begin_pos = 0
        begin_shift_full = begin_shift

        while scigar[0][0] == 'S' or scigar[0][0] == 'H':
            if scigar[0][0] == 'S':
                begin_shift_full += scigar[0][1]
            scigar = scigar[1:]
 
        while begin_shift and scigar:
            ctype, length = scigar.pop(0)
            assert ctype != 'S' and ctype != 'H'  #otherwise it should be a drop
            if ctype == 'M' or ctype == '=' or ctype == 'X':
                if length <= begin_shift:
                    begin_pos += length
                    begin_shift -= length
                else:
                    remaining_length = length - begin_shift
                    scigar.insert(0,(ctype, remaining_length))
                    begin_pos += begin_shift
                    begin_shift = 0
            elif ctype == 'I':
                if length <= begin_shift:
                    begin_shift -= length
                else:
                    remaining_length = length - begin_shift
                    begin_shift_full += remaining_length #remove rest of insert also
                    begin_shift = 0 #remove insert
                    if remaining_length > 0:
                        sys.stderr.write('Removing begin insert\n')
                        sys.stderr.flush()
            elif ctype == 'D' or ctype == 'N': #not part of sequence, so is pruned
                begin_pos += length #affects begin position
            elif ctype == 'P':
                pass #not part of sequence anymore

          

        end_shift_full = end_shift
        while scigar[-1][0] == 'S' or scigar[-1][0] == 'H':
            if scigar[-1][0] == 'S':
                end_shift_full += scigar[-1][1]
            scigar = scigar[:-1]


        while end_shift and scigar:
            ctype, length = scigar.pop()
            assert ctype != 'S' and ctype != 'H'  #otherwise it should be a drop
            if ctype == 'M' or ctype == '=' or ctype == 'X':
                if length <= end_shift:
                    end_shift -= length
                else:
                    remaining_length = length - end_shift
                    scigar.append((ctype, remaining_length))
                    end_shift = 0
            elif ctype == 'I':
                if length <= end_shift:
                    end_shift -= length
                else:
                    remaining_length = length - end_shift
                    end_shift_full += remaining_length
                    end_shift = 0
                    if remaining_length > 0:
                        sys.stderr.write('Removing end insert\n')
                        sys.stderr.flush()
            elif ctype == 'D' or ctype == 'N' or ctype == 'P':
                pass #not part of sequence anymore

        
        newcigar = ""

        #remove begin deletes/inserts
        while scigar:
            if scigar[0][0] == 'D' or scigar[0][0] == 'N':
                begin_pos += scigar[0][1]
                sys.stderr.write('Removing begin delete\n')
                scigar.pop(0)
            elif scigar[0][0] == 'I':
                begin_shift_full += scigar[0][1]
                sys.stderr.write('Removing begin insert\n')
                scigar.pop(0)
            else:
                break

        #remove end deletes/inserts
        while scigar:
            if scigar[-1][0] == 'D' or scigar[-1][0] == 'N':            
                sys.stderr.write('Removing end delete\n')
                scigar.pop()
            elif scigar[-1][0] == 'I':
                end_shift_full += scigar[-1][1]
                sys.stderr.write('Removing end insert\n')
                scigar.pop()
            else:
                break
    
        if begin_shift_full > 0:
            beginpart= '%dS' % begin_shift_full
        else:
            beginpart = ''
        if end_shift_full > 0:
            lastpart = '%dS' % end_shift_full
        else:
            lastpart = ''

        if len([e for e in scigar if not e[0] in ['S','H','D','P','I']]) == 0:
            newcigars.append('*') #unmapped
            newpos.append(0)
        else:
            for ctype, length in scigar:
                newcigar += '%d%s' % (length, ctype)
            newcigar = beginpart + newcigar + lastpart
            newcigars.append(newcigar)
            newpos.append(begin_pos)

        

    #update records
    new_supplementary = []
    for rg,nc,np in zip(readgroup['supplementary'], newcigars[1:], newpos[1:]):
        if nc == '*':
            continue #unmapped supplementary alignments are removed
        rg = list(rg)
        rg[3] = str(int(rg[3]) + np)
        rg[5] = nc
        new_supplementary.append(rg)
    if new_supplementary:
        readgroup['supplementary'] = new_supplementary
    else:
        del readgroup['supplementary']

    readgroup['primary'] = list(readgroup['primary'])
    readgroup['primary'][5] = newcigars[0]
    if newcigars[0] == '*':
        readgroup['primary'][3] = '0'
    else:
        readgroup['primary'][3] = str(int(readgroup['primary'][3]) + newpos[0])

    readgroup['is_pruned'] = pruning_required
    return (readgroup,pruning_required)

def get_sa_tag(row):
    nm = '0'
    for elem in row[11:]:
        if elem.startswith('NM'):
            nm = elem[5:]
            break
    strand = '-' if (int(row[1]) & 0x10) else '+'
    return '%s,%s,%s,%s,%s,%s;' % (row[2], row[3], strand, row[5], row[4], nm)


def set_sa_tag(row, sa_tag):
    row = list(row)
    sa_tag_set = False
    for pos, elem in enumerate(row[11:]):
        if elem.startswith('SA'):
            if sa_tag:
                row[11 + pos] = 'SA:Z:' + sa_tag
                sa_tag_set = True
            else:
                del row[11 + pos]
            break
    if sa_tag and not sa_tag_set:
        row.append('SA:Z:' + sa_tag)
    return row

def update_sa_tag(readgroup):
    readgroup = readgroup.copy()

    if not 'supplementary' in readgroup:
        return readgroup

    sa_primary = get_sa_tag(readgroup['primary'])
    sa_sups = [get_sa_tag(rg) for rg in readgroup['supplementary']]

    readgroup['primary'] = set_sa_tag(readgroup['primary'], ''.join(sa_sups))

    new_supplementary = []
    for pos, rg in enumerate(readgroup['supplementary']):
        sa_sups_tmp = list(sa_sups)
        del sa_sups_tmp[pos]
        rg = set_sa_tag(rg, sa_primary + ''.join(sa_sups_tmp))
        new_supplementary.append(rg)

    if not new_supplementary:
        del readgroup['supplementary']
    else:
        readgroup['supplementary'] = new_supplementary
    return readgroup


def update_mate_tag(read1, read2):
    read1 = read1.copy()
    read2 = read2.copy()
    read1['primary'] = list(read1['primary'])
    read2['primary'] = list(read2['primary'])

    pos1 = read1['primary'][3]
    pos2 = read2['primary'][3]

    mtag1 = 'MC:Z:' + read1['primary'][5]
    mtag2 = 'MC:Z:' + read2['primary'][5]

    read1['primary'][7] = pos2
    read2['primary'][7] = pos1
    
    for pos, tag in enumerate(read1['primary'][11:]): 
        if tag.startswith('MC'):
            read1['primary'][pos + 11] = mtag2

    for pos, tag in enumerate(read2['primary'][11:]): 
        if tag.startswith('MC'):
            read2['primary'][pos + 11] = mtag1


    if 'supplementary' in read1:
        new_sup = []
        for row in read1['supplementary']:
            row = list(row)
            if mtag2 != 'MC:Z:*':
                row[7] = pos2
                for pos, tag in enumerate(row[11:]): 
                    if tag.startswith('MC'):
                        row[pos + 11] = mtag2
            else:
                row[7] = row[3]
            new_sup.append(row)
        read1['supplementary'] = new_sup
    
    if 'supplementary' in read2:
        new_sup = []
        for row in read2['supplementary']:
            row = list(row)
            if mtag1 != 'MC:Z:*':
                row[7] = pos1
                for pos, tag in enumerate(row[11:]): 
                    if tag.startswith('MC'):
                        row[pos + 11] = mtag1
            else:
                row[7] = row[3]
            new_sup.append(row)
        read2['supplementary'] = new_sup
    return (read1, read2)


def get_position_order(readgroup):
    #get order of alignments in terms of (unrevcomped) read sequence. So last alignment relates to last bases read by machine for that read

    #if remaining part of cigar ends in an insert/delete, it is pruned completely
    #'cigars': cigars in unrevcomped read order: i.e. last alignment action is cigars[-1][-1], and first alignment action is cigars[0][0]. 'H'ardclip elements are removed.

    scigars = [cigar_split(readgroup['primary'][5])] + [cigar_split(r[5]) for r in readgroup.get('supplementary',[])]
    reverse = [int(readgroup['primary'][1]) & 0x10] + [int(r[1]) & 0x10 for r in readgroup.get('supplementary',[])]
    seq_length = len(readgroup['primary'][10]) 
    
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
            if ctype == 'S':
                pos += length
            elif ctype == 'M' or ctype == 'I' or ctype == '=' or ctype == 'X':
                if rev:
                    minpos.append(seq_length - (pos + length))
                    maxpos.append(seq_length - pos)
                else:
                    minpos.append(pos)
                    maxpos.append(pos + length - 1)
                pos += length
                
                if type != 'I':
                    alength += length
                    
            elif ctype == 'H' or ctype == 'D' or ctype == 'N' or ctype == 'P':
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
    chrom = [readgroup['primary'][2]] + [e[2] for e in readgroup.get('supplementary',[])]
    pos = [int(readgroup['primary'][3])] + [int(e[3]) for e in readgroup.get('supplementary',[])]
    seq = [readgroup['primary'][9]] + [e[9] for e in readgroup.get('supplementary',[])]


    sortidx = sorted(range(len(read_startpos)), key=read_startpos.__getitem__)
    chrom_sorted = [chrom[e] for e in sortidx]
    pos_sorted = [pos[e] for e in sortidx]
    cigars_sorted = [cigars[e] for e in sortidx]
    cigars_sorted_without_hardclip = [[e for e in cigar if not e[0] == 'H'] for cigar in cigars_sorted]
    readstartpos_sorted = [read_startpos[e] for e in sortidx]
    readendpos_sorted = [read_endpos[e] for e in sortidx]
    alignlength_sorted = [align_length[e] for e in sortidx]
    reverse_sorted = [reverse[e] for e in sortidx]
    seq_sorted = [seq[e] for e in sortidx]
    return {'idx':sortidx, 'chrom':chrom_sorted, 'pos':pos_sorted, 'seq': seq_sorted, 'readpos_start':readstartpos_sorted, 'readpos_end':readendpos_sorted, 'cigars': cigars_sorted_without_hardclip, 'align_length':alignlength_sorted, 'reverse':reverse_sorted}


def clip_first_elem(readgroup, pos_order):
    pos = pos_order['idx'][0]
    readgroup = readgroup.copy()
    if pos == 0:
        readgroup['primary'],clipped = record_hardclip(readgroup['primary'], only_first=True)
    else:
        readgroup['supplementary'] = readgroup['supplementary']
        readgroup['supplementary'][pos - 1],clipped = record_hardclip(readgroup['supplementary'][pos - 1], only_first=True)
    return readgroup

def clip_last_elem(readgroup, pos_order):
    pos = pos_order['idx'][-1]
    readgroup = readgroup.copy()
    if pos == 0:
        readgroup['primary'],clipped = record_hardclip(readgroup['primary'], only_last=True)
    else:
        readgroup['supplementary'] = readgroup['supplementary']
        readgroup['supplementary'][pos - 1],clipped = record_hardclip(readgroup['supplementary'][pos - 1], only_last=True)
    return readgroup


def record_filter(read):
    if read[5] == '*':
        return read

    
    if len(read[9]) <= 1:
        read = list(read)
        read[5] = '*'
        read[3] = '0'
    return read
        

def rg_filter(readgroup):
    readgroup['primary'] = record_filter(readgroup['primary'])

    if 'supplementary' in readgroup:
        new = [record_filter(e) for e in readgroup['supplementary']]
        new = [e for e in new if not e[5] == '*']

        if readgroup['primary'][5] == '*' and len(new) > 0:
            new_primary = new.pop(0)
            new_primary[1] = str(int(new_primary[1]) & (~0x800)) #update flags
            readgroup['primary'] = new_primary

        if len(new) > 0:
            readgroup['supplementary'] = new
        else:
            del readgroup['supplementary']

    return readgroup           


def is_unmapped(read):
    return (read[5] == '*')

def process_querygroup(querygroup, loose_ends, adapters):
    global looseend_counter, adapter_counter, pruned_counter, clipped_counter, modified_counter, read_counter, check_counter, supclip_counter, alignment_counter

    read_counter += 2
    alignment_counter += len(querygroup)

    #fast path for normal pairing
    if len(querygroup) == 2:
        flag1 = int(querygroup[0][1])
        flag2 = int(querygroup[1][1])
        read_length = len(querygroup[0][10])

        if not ((flag1 & 0x900) or (flag2 & 0x900)): #not supplementary nor secondary
            chrom1 = querygroup[0][2]
            chrom2 = querygroup[1][2]
            same_chrom = chrom1 == chrom2

            cigar1 = querygroup[0][5]
            cigar2 = querygroup[1][5]

            

            if cigar1 == '*' or cigar2 == '*':
                #one is unmapped, no fast path for now
                pass
            else:
                #remove hardclip and digits
                cigar1 = [e for e in cigar1 if not (e.isdigit() or e == 'H')] 
                cigar2 = [e for e in cigar2 if not (e.isdigit() or e == 'H')]

                if flag1 & 0x10:
                    cigar1 = cigar1[::-1]

                if flag2 & 0x10:
                    cigar2 = cigar2[::-1]


                if same_chrom:
                    pos1 = int(querygroup[0][3])
                    pos2 = int(querygroup[1][3])
                    diff = abs(pos1 - pos2)
                else:
                    diff = 1000000000000000

                check = False
                #reads are overlapping, check for possible  adapter read through
                if adapters and diff <= 20 and (cigar1[-1] == 'S' and cigar2[-1] == 'S'): 
                    check = True # possible adapter read through
                
                #check for possible chimeric switch at the beginning or the end of the fragment:    (S?R1) --  (R2S?)
                if loose_ends and (cigar1[0] == 'S' or cigar2[0] == 'S'):  #check if the end of the reads (near adapter) have no soft-clipping
                    check = True
                
                #reads are far apart, chimeric switch is located either in reads or in 'insert' fragment bewteen reads
                #check if the reads have a soft-clip at the non-adapter end   (R1S?) -- (S?R2)
                if diff >= MAX_READ_DIST and (cigar1[-1] == 'S' or cigar2[-1] == 'S'):
                    check = True

                if not check:
                    return (querygroup,[])

    check_counter += 2

    read1 = process_readgroup([row for row in querygroup if (int(row[1]) & 0x40)])
    read2 = process_readgroup([row for row in querygroup if (int(row[1]) & 0x80)])
   
    
    is_unmapped1 = is_unmapped(read1['primary'])
    is_unmapped2 = is_unmapped(read2['primary'])
    
    pos_order1 = get_position_order(read1)
    pos_order2 = get_position_order(read2)

    modified1 = False
    modified2 = False

    #STEP1: REMOVING READTHROUGH ADAPTERS

    if adapters and not is_unmapped1 and not is_unmapped2 and len(pos_order1['cigars']) == len(pos_order2['cigars']) and  pos_order1['cigars'][-1][-1][0] == 'S' and pos_order2['cigars'][-1][-1][0] == 'S':
        #ADAPTER:  ****-----------------****
        #R1:           -----------------*>
        #R2:         <*----------------- 
        #Startpos:     ^
        #Features:
        #1) Alignment start position (^) is very close together (should be same, except for soft clip threshold errors)
        #2) Alignment length (-) is very similar
        #3) Read alignment orientation is dissimilar
        #4) soft clip at read ends

        #only primary reads --> distance > 100kb --> remove clipping
        chrom1 = pos_order1['chrom'][-1]
        chrom2 = pos_order2['chrom'][-1]
        pos1 = pos_order1['pos'][-1]
        pos2 = pos_order2['pos'][-1]

        #ADAPTER FIXME: only primary reads --> distance < read length --> remove clipping
        reverse1 = pos_order1['reverse'][-1]
        reverse2 = pos_order2['reverse'][-1]


        if (reverse1 != reverse2) and chrom1 == chrom2 and abs(pos1 - pos2) < 20 and abs(pos_order1['align_length'][-1] - pos_order2['align_length'][-1]) < 20:
            #align adapter and find start point
            prune1 = pos_order1['cigars'][-1][-1][1]
            prune2 = pos_order2['cigars'][-1][-1][1]
            seq1 = pos_order1['seq'][-1]
            seq2 = pos_order2['seq'][-1]
            max_prune1 = len(seq1) - pos_order1['readpos_start'][-1]
            max_prune2 = len(seq2) - pos_order2['readpos_start'][-1]
            modified_prune = True

            while modified_prune:
                modified_prune = False
                if reverse1:
                    testseq1 = swalign.revcomp(seq1[:prune1])
                else:
                    testseq1 = seq1[-prune1:]

                if reverse2:
                    testseq2 = swalign.revcomp(seq2[:prune2])
                else:
                    testseq2 = seq2[-prune2:]
                
                align1 = aligner.align(testseq1, adapter)
                align2 = aligner.align(testseq2, adapter)

                if (align1.q_pos - align1.r_pos) > 0 and prune1 < max_prune1:
                    prune1 = min(prune1 + (align1.q_pos - align1.r_pos),max_prune1)
                    modified_prune = True
                if (align2.q_pos - align2.r_pos) > 0 and prune2 < max_prune2:
                    prune2 = min(prune2 + (align2.q_pos - align2.r_pos),max_prune2)
                    modified_prune = True
            
            prune1 = min(prune1 + (align1.q_pos - align1.r_pos), max_prune1)
            prune2 = min(prune2 + (align2.q_pos - align2.r_pos), max_prune2)

            score1 = align1.score
            score2 = align2.score

            if align1.r_end == len(testseq1) and align1.q_end < len(adapter):
                score1 += len(adapter) - align1.q_end
            if align2.r_end == len(testseq2) and align2.q_end < len(adapter):
                score2 += len(adapter) - align2.q_end
            

            if score1 >= 9 and score2 >= 9:
                p1 = pos_order1['idx'][-1]
                p2 = pos_order2['idx'][-1]
                if p1 == 0:
                    read1['primary'] = record_lastclip(read1['primary'],prune1)
                else:
                    read1['supplementary'][p1-1] = record_lastclip(read1['supplementary'][p1-1],prune1)
                
                if p2 == 0:
                    read2['primary'] = record_lastclip(read2['primary'],prune2)
                else:
                    read2['supplementary'][p2-1] = record_lastclip(read2['supplementary'][p2-1],prune2)
                
                modified1 = True
                modified2 = True
                adapter_counter += 2

                #update read stats
                is_unmapped1 = is_unmapped(read1['primary'])
                is_unmapped2 = is_unmapped(read2['primary'])
        
                pos_order1 = get_position_order(read1)
                pos_order2 = get_position_order(read2)

    
    ####STEP2:  PRUNING multi-aligned read regions
    if 'supplementary' in read1:
        read1,is_pruned = rg_prune(read1)
        if is_pruned:
            is_unmapped1 = is_unmapped(read1['primary'])
            pos_order1 = get_position_order(read1)
            modified1 = True
            pruned_counter += 1

    if 'supplementary' in read2:
        read2,is_pruned = rg_prune(read2)
        if is_pruned:
            is_unmapped2 = is_unmapped(read2['primary'])
            pos_order2 = get_position_order(read2)
            modified2 = True
            pruned_counter += 1




    #STEP 3: REMOVING soft-clips due to chimeric reads and optionally loose-end soft-clips
    if not 'supplementary' in read1 and not 'supplementary' in read2: #only two alignment records
        if not is_unmapped1 and not is_unmapped2:
            cigar1 = pos_order1['cigars'][0]
            cigar2 = pos_order2['cigars'][0]

            chrom1 = read1['primary'][2]
            chrom2 = read2['primary'][2]
            same_chrom = chrom1 == chrom2
            if same_chrom:
                pos1 = int(read1['primary'][3])
                pos2 = int(read2['primary'][3])
                diff = abs(pos1 - pos2)
            else:
                diff = 100000000000000
           
            #only primary reads --> distance > MAX_READ_DIST --> remove clipping
            if diff >= MAX_READ_DIST and (cigar1[-1][0] == 'S' or cigar2[-1][0] == 'S'):
                if cigar1[-1][0] == 'S':
                    read1 = clip_last_elem(read1, pos_order1)
                    modified1=True
                    clipped_counter += 1

                if cigar2[-1][0] == 'S':
                    read2 = clip_last_elem(read2, pos_order2)
                    modified2 = True
                    clipped_counter += 1


            if loose_ends and ((cigar1[0][0] == 'S') or (cigar2[0][0] == 'S')):
                if cigar1[0][0] == 'S':
                    read1 = clip_first_elem(read1, pos_order1)
                    modified1=True
                    looseend_counter +=1
                if cigar2[0][0] == 'S':
                    read2 = clip_first_elem(read2, pos_order2)
                    modified2 = True
                    looseend_counter +=1

        elif is_unmapped1 and is_unmapped2:
            pass
        elif is_unmapped1: 
            if pos_order2['cigars'][-1][-1][0] == 'S':
                read2 = clip_last_elem(read2, pos_order2)
                modified2 = True
                clipped_counter += 1
            if loose_ends and pos_order2['cigars'][0][0][0] == 'S':
                read2 = clip_first_elem(read2, pos_order2)
                looseend_counter +=1
                modified2 = True
        else: #is_unmapped2 
            if pos_order1['cigars'][-1][-1][0] == 'S':
                read1 = clip_last_elem(read1, pos_order1)
                modified1 = True
                clipped_counter += 1
            if loose_ends and pos_order1['cigars'][0][0][0] == 'S':
                read1 = clip_first_elem(read1, pos_order1)
                looseend_counter +=1
                modified1 = True
    else:
        if not is_unmapped1 and not is_unmapped2:
            reverse1 = pos_order1['reverse'][-1]
            reverse2 = pos_order2['reverse'][-1]
            chrom1 = pos_order1['chrom'][-1]
            chrom2 = pos_order2['chrom'][-1]
            alignment_start1 = int(pos_order1['pos'][-1])
            alignment_start2 = int(pos_order2['pos'][-1])

            alignment_stop1 = alignment_start1 + pos_order1['align_length'][-1]
            alignment_stop2 = alignment_start2 + pos_order2['align_length'][-1]

            if 'supplementary' in read1 and not ('supplementary' in read2) and len(read1['supplementary']) == 1:
                if pos_order2['cigars'][0][-1][0] == 'S' and chrom1 == chrom2 and reverse1 != reverse2:  
                    if reverse1 and alignment_start1 >= alignment_start2 and alignment_stop1 < (alignment_stop2 + pos_order2['cigars'][0][-1][1]):
                        prune = alignment_stop2  - alignment_stop1 
                        if prune > 0:
                            read2['primary'] = record_lastclip(read2['primary'],prune,aligned=True)
                        else:
                            read2['primary'],clipped = record_hardclip(read2['primary'],only_last=True)
                        modified2 =True
                        clipped_counter += 1

                    elif reverse2 and alignment_stop2 >= alignment_stop1 and alignment_start1 > (alignment_start2 - pos_order2['cigars'][0][-1][1]):
                        prune = alignment_start1 - alignment_start2 
                        if prune > 0:
                            read2['primary'] = record_lastclip(read2['primary'],prune,aligned=True)
                        else:
                            read2['primary'],clipped = record_hardclip(read2['primary'], only_last=True)
                        modified2 =True
                        clipped_counter += 1

            if 'supplementary' in read2 and not ('supplementary' in read1) and len(read2['supplementary']) == 1:
                if pos_order1['cigars'][0][-1][0] == 'S' and chrom1 == chrom2 and reverse1 != reverse2:

                    if reverse2 and alignment_start2 >= alignment_start1 and alignment_stop2 < (alignment_stop1 + pos_order1['cigars'][0][-1][1]):
                        prune = alignment_stop1  - alignment_stop2
                        if prune > 0:
                            read1['primary'] = record_lastclip(read1['primary'],prune, aligned=True)
                        else:
                            read1['primary'],clipped = record_hardclip(read1['primary'], only_last=True)
                        modified1 =True
                        clipped_counter += 1
                        
                    elif reverse1 and alignment_stop1 >= alignment_stop2 and alignment_start2 > (alignment_start1 - pos_order1['cigars'][0][-1][1]):
                        prune = alignment_start2  - alignment_start1
                        if prune > 0:
                            read1['primary'] = record_lastclip(read1['primary'],prune, aligned=True)
                        else:
                            read1['primary'],clipped = record_hardclip(read1['primary'], only_last=True)
                        modified1 =True
                        clipped_counter += 1

        if 'supplementary' in read1:
            read1,clipped = rg_hardclip(read1)
            if clipped:
                supclip_counter +=1
                modified1 =True
        if 'supplementary' in read2:
            read2,clipped = rg_hardclip(read2)
            if clipped:
                supclip_counter +=1
                modified2 =True

    #remove reads that have <= 1bp of sequence
    if modified1:
        read1 = rg_filter(read1)
    if modified2:
        read2 = rg_filter(read2)

    
    
    if 'supplementary' in read1:
        read1 = update_sa_tag(read1)

    if 'supplementary' in read2:
        read2 = update_sa_tag(read2)


    if modified1 or modified2:
        read1, read2 = update_mate_tag(read1, read2)
    if modified1:
        oldread1 = process_readgroup([row for row in querygroup if (int(row[1]) & 0x40)])
        validate(read1, oldread1)
    if modified2:
        oldread2 = process_readgroup([row for row in querygroup if (int(row[1]) & 0x80)])
        validate(read2, oldread2)

    querygroup1 = [read1.get('primary', None)] + read1.get('supplementary',[]) + read1.get('secondary',[]) + read1.get('supplementary_secondary',[])
    querygroup2 = [read2.get('primary', None)] + read2.get('supplementary',[]) + read2.get('secondary',[]) + read2.get('supplementary_secondary',[])
    querygroup = [e for e in  (querygroup1 + querygroup2) if not e is None]
    querygroup = [e for e in querygroup if not e[9] == '']


    if modified1:
        modified_counter += 1
    if modified2:
        modified_counter += 2


    return (querygroup,[])



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
        cigar = oread[5]
        if cigar == '*':
            continue
        scigar = cigar_split(cigar)

        oldpairs.extend(get_pos_pairs(scigar, int(oread[3]), oread[9], oread[10]))
    oldpairs = set(oldpairs)
    
    had_errors = False

    #get position seq/postion ref pairs old reads
    errors = []
    for read in reads:
        errors = []
        cigar = read[5]
        if cigar == '*':
            if int(read[3]) != 0:
                errors.append('- Unmapped reads should  have no position')
            if len(read[9]) > 0:
                errors.append('- Unmapped reads should  have no sequence')
            if len(read[10]) > 0:
                errors.append('- Unmapped reads should  have no qualityscores')
            continue
        
        scigar = cigar_split(read[5])
        
        #length check
        read_length = 0
        for stype, length in scigar:
            if stype == 'S' or stype == 'M' or stype == 'X' or stype == '=' or stype == 'I':
                read_length += length
            if length == 0:
                errors.append('- Length 0 in cigar')    

        if read_length == 0:
            errors.append('- Unmapped read without cigar = *')
            
        if read_length != len(read[9]):
            errors.append('- Length of cigar (%d) does not match length of sequence (%d)' % (read_length, len(read[9])))

        if read_length != len(read[10]):
            errors.append('- Length of cigar (%d) does not match length of qualityscores (%d)' % (read_length, len(read[10])))



        #check clipping
        tcigar = list(scigar)
        if tcigar[0][0] == 'H':
            tcigar = tcigar[1:]
        if tcigar[-1][0] == 'H':
            tcigar = tcigar[:-1]
        if tcigar[0][0] == 'S':
            tcigar = tcigar[1:]
        if tcigar[-1][0] == 'S':
            tcigar = tcigar[:-1]
        internal_clip = len([e for e in tcigar if e[0] in ['H', 'S']])
        if internal_clip:
            errors.append('- Cigar has internal clipping: %s' % read[5])

        
        if len(tcigar) == 0:
            tcigar = [e for e in tcigar if e[0] in ['M','=','X']]
            errors.append('- Cigar has no mapping: %s' % read[5])

        last_stype = ''
        for stype, length in scigar:
            if stype == last_stype:
                errors.append('- Cigar has repeated record: %s' % read[5])
            last_stype = stype


        pairs = get_pos_pairs(scigar, int(read[3]), read[9], read[10])

        for pair in pairs:
            if pair not in oldpairs:
                opairs = list(oldpairs)
                opairs.sort()
                errors.append('- New mapping relations found: %s \n%s\n%s\n' % (str(pair), str(pairs), str(opairs)))
                break


        #checking mate/other tags?

        if errors:
            sys.stderr.write('Error found in %s:%s\n' % (read[0], read[3]))
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

        raise
    return True
    




if __name__ == '__main__':
    #with open('test.sam','r') as f:
    with sys.stdin as f:
        reader = csv.reader(f,delimiter='\t', quoting=csv.QUOTE_NONE)
        if len(sys.argv) > 1:
            loose_ends = bool(sys.argv[1])
            adapters = bool(sys.argv[2])
        else:
            loose_ends = True
            adapters = True
        sys.stderr.write('Loose ends: %s\n' % str(loose_ends))
        sys.stderr.write('Adapters: %s\n' % str(adapters))
        lastname = ''
        querygroup = []
        #with open('test.out', 'w') as out_file:
        with sys.stdout as out_file:
            counter = 0
            for row in reader:
                counter += 1
                if (counter % 10000) == 0:
                    sys.stderr.write(str(counter) + '\n')
                    sys.stderr.flush()

                if row[0][0] == '@':
                    out_file.write('\t'.join(row) + '\n')
                    #filter_file.write('\t'.join(row) + '\n')
                    continue
                name = row[0]
                if lastname == name:
                    querygroup.append(row)
                else:
                    if querygroup:
                        assert len(querygroup) >= 2
                        res_kept, res_filtered = process_querygroup(querygroup, loose_ends, adapters)
                        #ry:
                        #   res_kept, res_filtered = process_querygroup(querygroup, loose_ends, adapters)
                        #xcept Exception as e:
                        #   sys.stdout.write(str(e))
                        #   sys.stdout.flush()
                        #   raise e

                        if len(res_kept) < 2 and len(res_kept) < len(querygroup):
                            sys.stderr.write(str(querygroup) + '\n\n')
                            sys.stderr.write(str(res_kept) + '\n\n')
                            sys.stderr.flush()
                            raise
                        res_kept = '\n'.join(['\t'.join(xrow) for xrow in res_kept]) + '\n'
                        out_file.write(res_kept)
                        if res_filtered:
                            res_filtered = '\n'.join(['\t'.join(xrow) for xrow in res_filtered]) + '\n'
                            #filter_file.write(res_filtered)
                    lastname = name
                    querygroup = [row]
            if querygroup and len(querygroup) >= 2: 
                res_kept, res_filtered = process_querygroup(querygroup, loose_ends, adapters)
                res_kept = '\n'.join(['\t'.join(xrow) for xrow in res_kept]) + '\n'
                out_file.write(res_kept)
                if res_filtered:
                    res_filtered = '\n'.join(['\t'.join(xrow) for xrow in res_filtered]) + '\n'
                    #filter_file.write(res_filtered)
        sys.stderr.write('Read counts: %d\n' %read_counter)
        sys.stderr.write('Alignment counts: %d\n' % alignment_counter)
        sys.stderr.write('Reads detail checked: %d\n' % check_counter)
        sys.stderr.write('Reads modified: %d\n' % modified_counter)
        sys.stderr.write('- Loose ends: %d\n' %looseend_counter)
        sys.stderr.write('- Adapters: %d\n' % adapter_counter)
        sys.stderr.write('- Hardclipped: %d\n' %clipped_counter)
        sys.stderr.write('- Supplementary clips: %d\n' % supclip_counter)
        sys.stderr.write('- Pruned: %d\n' %pruned_counter)
        sys.stderr.flush()

