#!/bin/env python
import re
import numpy
import csv
import sys
import gzip
import re

#Poly-A, Illumina, PCR primer illumina/nextera 1, PCR primer illumina/nextera 2, Nextera transposon 
adapter_finder = re.compile('AAAAAAAA|AGATCGGAAGAGC|AATGATACGGCGACCACC|CAAGCAGAAGACGGCA|AGATGTGTATAAGAGACAG')


def detect_adapter(seq):
    poly_a = 0
    illumina = 0
    pcr_illumina_1 = 0
    pcr_illumina_2 = 0
    nextera = 0

    res = adapter_finder.findall(seq)
    for e in res:
        if e == 'AAAAAAAA':
            poly_a += 1
        elif e == 'AGATCGGAAGAGC':
            illumina += 1
        elif e == 'AATGATACGGCGACCACC':
            pcr_illumina_1 += 1
        elif e == 'CAAGCAGAAGACGGCA':
            pcr_illumina_2 += 1
        else:
            nextera += 1

    return (poly_a, illumina, pcr_illumina_1, pcr_illumina_2, nextera)

def soft_clipped_region(splitted_cigar):
    clipped_region = 0
    aligned_region = 0
    inserted_region = 0
    deleted_region = 0
    for ctype, length in splitted_cigar:
        if ctype == 'S':
            clipped_region += length
        elif ctype == 'M' or ctype == '=':
            aligned_region += length
        elif ctype == 'D':
            deleted_region += length
        elif ctype == 'I':
            inserted_region += length
    return (clipped_region, aligned_region, inserted_region, deleted_region)


def cigar_filter(splitted_cigar, qual):
    res = numpy.zeros(len(qual),dtype=bool)
    pos = 0
    for ctype, length in splitted_cigar:
        if ctype == 'S':
            res[pos:pos+length] = True
            pos += length
        elif ctype == 'M' or ctype == 'I' or ctype == '=' or ctype == 'X':
            pos += length
        elif ctype == 'H' or ctype == 'D' or ctype == 'N' or ctype == 'P':
            pass
    assert len(qual) == pos
    return res

   

def toqual(qual):
    return numpy.array([ord(e) + 33 for e in qual],dtype='uint8')

def split_seq(splitted_cigar, seq):
    res_clipped = []
    res_nonclipped = []
    pos = 0
    for ctype, length in splitted_cigar:
        if ctype == 'S':
            res_clipped.append(seq[pos:pos+length])
            pos += length
        elif ctype == 'M' or ctype == 'I' or ctype == '=' or ctype == 'X':
            res_nonclipped.append(seq[pos:pos+length])
            pos += length
        elif ctype == 'H' or ctype == 'D' or ctype == 'N' or ctype == 'P':
            pass
    assert pos == len(seq)
    return (numpy.concatenate(res_clipped), numpy.concatenate(res_nonclipped))

def cigar_split(cigar):
    cigar = re.split('([MIDNSHP=X]{1,1})',cigar)
    return [(ctype, int(length)) for length, ctype in zip(cigar[::2], cigar[1::2])]


def split_flag(flag):
    f = "{0:012b}".format(int(flag))
    f = int(flag)

    res =[]

    return tuple([f  & 1<<i for i in range(12)])


def average_quality_by_cycle(data):
    reverse, qual = data.Get(_.reverse, _.qual)()
    maxlen = max([len(e) for e in qual])

    qual = [q if not r else q[::-1] for (r,q) in zip(reverse,qual)]

    res = []
    for i in range(0, maxlen):
        res.append(numpy.mean([e[i] for e in qual if len(e) > i]))
    return numpy.array(res)
        


def clipped_qual_by_cycle(data, threshold, threshold_softclip):
    read1, reverse, qual, cigar = data.Get(_.read1, _.reverse, _.qual, _.cigar)()
    reverse = reverse.copy()
    reverse[~read1] = ~reverse[~read1]

    cigar_filters = [cigar_filter(c,q) for c,q in zip(cigar, qual)]
    maxlen = max([len(e) for e in qual])

    qual = [q if not r else q[::-1] for (r,q) in zip(reverse,qual)]
    cigar_filters = [c if not r else c[::-1] for r,c  in zip(reverse,cigar_filters)]

    res_count = numpy.zeros(maxlen, dtype=int)
    res_value = numpy.zeros(maxlen,dtype=float)

    res_count_notclipped = numpy.zeros(maxlen,dtype=int)
    res_value_notclipped = numpy.zeros(maxlen,dtype=float)


    res_filter = numpy.zeros(maxlen,dtype=int)
    
    filter_read_count = 0
    nofilter_read_count = 0

    total_map_area = 0
    total_trimmed_area = 0
    total_clipped_area = 0 
    total_filtered_area = 0
    for q, cgf,ci in zip(qual, cigar_filters, cigar):

        f = numpy.zeros(len(q),dtype=bool)
        startpos = 0
        endpos = 0
        while (startpos < len(q) and q[startpos] < threshold):
            startpos += 1
        f[:startpos] = True
        while (endpos < len(q) and (q[-(endpos + 1)] < threshold)):
            endpos += 1
        if endpos > 0 or q[-1] < threshold:
            f[-(endpos + 1):] = True
        tcgf = cgf.copy()
        tcgf[f] = False
        ntcgf = ~cgf
        ntcgf[f] = False
        
        
        if (ntcgf).sum() < threshold_softclip:
            filter_read_count += 1
            total_filtered_area += ntcgf.sum()
            continue
        total_clipped_area += tcgf.sum()
        total_trimmed_area += (~cgf).sum() - ntcgf.sum()

        nofilter_read_count += 1
        
        res_filter += f
        res_value[tcgf] += q[tcgf]
        res_count[tcgf] += 1

        res_count_notclipped[ntcgf] += 1
        res_value_notclipped[ntcgf] += q[ntcgf]
        total_map_area += ntcgf.sum()
    
    return {'average_qual_clipped':(res_value / res_count),
            'average_qual_notclipped':(res_value_notclipped / res_count_notclipped), 
            'read_fraction_clipped':res_count / float(nofilter_read_count), 
            'read_fraction_trimmed':(res_filter / float(nofilter_read_count)), 
            'read_fraction_filtered':(filter_read_count / (filter_read_count + float(nofilter_read_count))),
            'total_mapped_area':total_map_area,
            'total_filtered_area_mapped': total_filtered_area,
            'total_trimmed_area_mapped':total_trimmed_area,
            'total_clipped_area':total_clipped_area}



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

def record_hardclip(row, update_mate_tag):
    cigar = row[5]
    seq = row[9]
    qual = row[10]
    tags = row[11:]
    splitted_cigar = cigar_split(cigar)

    pos = 0
    newqual = []
    newseq = []
    newcigar = ''
    for ctype, length in splitted_cigar:
        if ctype == 'S':
            newcigar += '%dH' % length
            pos += length
        elif ctype == 'M' or ctype == 'I' or ctype == '=' or ctype == 'X':
            newcigar += '%d%s' % (length, ctype)
            newseq.append(seq[pos:pos + length])
            newqual.append(qual[pos:pos + length])
            pos += length
        elif ctype == 'H' or ctype == 'D' or ctype == 'N' or ctype == 'P':
            newcigar += '%d%s' % (length, ctype)
            pass
    if pos != len(seq):
        sys.stderr.write(str((splitted_cigar, len(seq), len(qual), pos, row)) + '\n')
        sys.stderr.flush()
    assert pos == len(seq)

    newqual = ''.join(newqual)
    newseq = ''.join(newseq)

    row = list(row)
    row[9] = newseq
    row[10] = newqual 
    row[5] = newcigar

    newtags = []
    for tag in tags:
        if tag.startswith('SA'):
            newrecs = []
            for srec in tag[5:].split(';')[:-1]:
                w = list(srec.split(','))
                sc = cigar_split(w[3]) #cigar tag
                w[3] = ''.join(['%d%s' % (length, 'H' if ctype == 'S' else ctype) for ctype,length in sc])
                newrecs.append(','.join(w))
            newtags.append(tag[:5] + ';'.join(newrecs) + ';')
        elif update_mate_tag and tag.startswith('MC'):
            sc = cigar_split(update_mate_tag)
            newtags.append(tag[:5] + ''.join(['%d%s' % (length, 'H' if ctype == 'S' else ctype) for ctype,length in sc]))
        else:
            newtags.append(tag)
    
    #FIXME: adapt SA tag, MC tag, 
    row = row[:11] + newtags
    return row


def rg_hardclip(readgroup, update_mate_tag):
    readgroup = readgroup.copy()
    readgroup['primary'] = record_hardclip(readgroup['primary'], update_mate_tag)
    readgroup['supplementary'] = [record_hardclip(r, update_mate_tag) for r in readgroup['supplementary']]
    return readgroup


def rg_prune(readgroup):
    readgroup = readgroup.copy()
    
    cigars = [cigar_split(readgroup['primary'][5])] + [cigar_split(r[5]) for r in readgroup['supplementary']]
    reverse = [int(readgroup['primary'][1]) & 0x10] + [int(r[1]) & 0x10 for r in readgroup['supplementary']]
    
    #part of read that is in-use
    inuse = []  #filter per alignment
    sequse = numpy.zeros(len(readgroup['primary'][9]),dtype=int)   #total aligncount per base
    for scigar,rev in zip(cigars,reverse):
        pos = 0
        w = numpy.zeros(len(sequse),dtype=bool)
        for ctype, length in scigar:
            if ctype == 'S':
                pos += length
            elif ctype == 'M' or ctype == 'I' or ctype == '=' or ctype == 'X':
                if rev:
                    sequse[::-1][pos:pos + length] += 1
                    w[::-1][pos:pos + length] = True
                else:
                    sequse[pos:pos + length] += 1
                    w[pos:pos + length] = True
                pos += length
            elif ctype == 'H' or ctype == 'D' or ctype == 'N' or ctype == 'P':
                pass
        inuse.append(w)


    
    #determine first and last part of alignment that is aligned more than once
    shift = []
    for rev,iu in zip(reverse,inuse):
        usecount = sequse[iu]
        begin_shift = 0
        end_shift = 0
        drop = False

        while usecount[begin_shift] > 1 and begin_shift < len(usecount): 
            begin_shift += 1
        while usecount[::-1][end_shift] > 1 and end_shift < len(usecount): 
            end_shift += 1
        drop = (begin_shift == len(usecount))

        if rev:
            begin_shift,end_shift = end_shift,begin_shift
        shift.append((begin_shift, end_shift, drop))


    #calculate niew cigars and new position per alignment
    newcigars = []
    newpos = []
    
    for (begin_shift,end_shift, drop),rev,scigar in zip(shift, reverse, cigars):
        if drop:
            newcigars.append('')
            newpos.append(None)
            continue

        begin_pos = 0
        begin_shift_full = begin_shift
        if scigar[0][0] == 'S':
            begin_shift_full += scigar[0][1]
            scigar = scigar[1:]
 
        while begin_shift:
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
                    begin_shift_full += remaining_length 
                    begin_shift = 0 #remove insert
                    if remaining_length > 0:
                        sys.stderr.write('Removing begin insert\n')
                        sys.stderr.flush()
            elif ctype == 'D' or ctype == 'N' or ctype == 'P':
                pass #not part of sequence anymore

        if begin_shift_full > 0:
            beginpart= '%dS' % begin_shift_full
        else:
            beginpart = ''

          

        end_shift_full = end_shift
        if scigar[-1][0] == 'S':
            end_shift_full += scigar[-1][1]
            scigar = scigar[:-1]

        while end_shift:
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

        if end_shift_full > 0:
            lastpart = '%dS' % end_shift_full
        else:
            lastpart = ''
        
        newcigar = ""
        for ctype, length in scigar:
            newcigar += '%d%s' % (length, ctype)

        newcigar = beginpart + newcigar + lastpart

        newcigars.append(newcigar)
        newpos.append(begin_pos)

    readgroup['primary'] = list(readgroup['primary'])
    readgroup['primary'][5] = newcigars[0]
    readgroup['primary'][3] = str(int(readgroup['primary'][3]) + newpos[0])

    #update records
    new_supplementary = []
    for rg,nc,np in zip(readgroup['supplementary'], newcigars[1:], newpos[1:]):
        if not nc:
            continue #FIXME --> FIXME to filter file? and handling in update_sa_tag
        rg = list(rg)
        rg[3] = str(int(rg[3]) + np)
        rg[5] = nc
        new_supplementary.append(rg)
    readgroup['supplementary'] = new_supplementary


    return readgroup

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
    for pos, elem in enumerate(row[11:]):
        if elem.startswith('SA'):
            if sa_tag:
                row[11 + pos] = 'SA:Z:' + sa_tag
            else:
                del row[11 + pos]
            break
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

def process_querygroup(querygroup):
    #fast path for normal pairing
    if len(querygroup) == 2:
        flag1 = int(querygroup[0][1])
        flag2 = int(querygroup[1][1])

        if not ((flag1 & 0x900) or (flag2 & 0x900)):
            return (querygroup,[])

    read1 = process_readgroup([row for row in querygroup if (int(row[1]) & 0x40)])
    read2 = process_readgroup([row for row in querygroup if (int(row[1]) & 0x80)])

    if not ('supplementary' in read1 or 'supplementary' in read2):
        return (querygroup,[])


    update_mate_tag_r1 = ''
    update_mate_tag_r2 = ''

    if 'supplementary' in read1:
        read1_pruned = rg_prune(read1)
        update_mate_tag_r1 = read1_pruned['primary'][5]
    if 'supplementary' in read2:
        read2_pruned = rg_prune(read2)
        update_mate_tag_r2 = read2_pruned['primary'][5]
    
    if 'supplementary' in read1:
        read1 = rg_hardclip(read1_pruned, update_mate_tag_r2)
        read1 = update_sa_tag(read1)
    if 'supplementary' in read2:
        read2 = rg_hardclip(read2_pruned, update_mate_tag_r1)
        read2 = update_sa_tag(read2)


    querygroup1 = [read1.get('primary', None)] + read1.get('supplementary',[]) + read1.get('secondary',[]) + read1.get('supplementary_secondary',[])
    querygroup2 = [read2.get('primary', None)] + read2.get('supplementary',[]) + read2.get('secondary',[]) + read2.get('supplementary_secondary',[])
    querygroup = [e for e in  (querygroup1 + querygroup2) if not e is None]

    return (querygroup,[])




if __name__ == '__main__':
    with sys.stdin as f:
        reader = csv.reader(f,delimiter='\t', quoting=csv.QUOTE_NONE)
        mode = sys.argv[1]

        if mode == 'clips':
            if len(sys.argv) > 2:
                xlen= int(sys.argv[2])
            else:
                xlen = 25
            with sys.stdout as f2:
                for row in reader:
                    if row[0][0] == '@':
                        continue
                    name = row[0]
                    flag = row[1]
                    cigar = row[5]
                    seq = row[9]

                    paired, properpair, unmapped, mate_unmapped, reverse, mate_reverse, read1, read2, secondary, qcfail, duplicate, supplementary = split_flag(flag)
                    if duplicate:
                        continue
                    scigar = cigar_split(cigar)
                    softclipped, aligned, inserted, deleted = soft_clipped_region(scigar)
                    if softclipped < numpy.abs(xlen):
                        continue
        
                    res = ''
                    pos = 0
                    for ctype, length in scigar:
                        if ctype == 'S':
                            res = res + seq[pos:pos+length]
                        elif ctype == 'M' or ctype == '=' or ctype == 'X':
                            if xlen < 0:
                                res = res + ' [' + seq[pos:pos+length] + '] '
                            pos += length
                        elif ctype == 'I':
                            if xlen < 0:
                               res = res + '<' + seq[pos:pos+length] + '>'
                            pos += length
                        elif ctype == 'H' or ctype == 'D' or ctype == 'N' or ctype == 'P':
                            pass
                    f2.write('%s\t%d\t%d\t%d\t%d\t%d\t%s\n' % (name, softclipped, aligned, inserted, deleted,int(supplementary),res))
        elif mode == 'stats':
            total_count = 0
            unmapped_count = 0
            mqual20_count = 0
            secondary_count=0
            supplementary_count = 0
            duplicate_count = 0
            duplicate_supplement_count = 0

            softclipped_base_count = 0
            aligned_base_count = 0
            inserted_base_count = 0
            deleted_base_count = 0
            base_count = 0

            soft_clipped_50 = 0
            soft_clipped_60 = 0
            soft_clipped_70 = 0

            aligned_50 = 0
            aligned_60 = 0
            aligned_70 = 0

            secondary_count = 0
            supplementary_count = 0
            sup_diffchrom_count = 0

            poly_a = 0
            illumina_adapter = 0
            pcr_illumina_adapter_1 = 0
            pcr_illumina_adapter_2 = 0
            nextera = 0


            for row in reader:
                if row[0][0] == '@':
                    continue
                flag = row[1]
                my_contig = row[2]
                mqual = row[4]
                cigar = row[5]
                mate_contig = row[6]
                seq = row[9]
                tags = row[11:]
                total_count += 1
                if cigar == '*':
                    unmapped_count += 1
                    continue

                mqual = int(mqual)
                scigar = cigar_split(cigar)

                soft_clipped, aligned, inserted, deleted = soft_clipped_region(scigar)
                softclipped_base_count += soft_clipped
                aligned_base_count += aligned
                inserted_base_count += inserted
                deleted_base_count += deleted
                base_count += len(seq)
                
                paired, properpair, unmapped, mate_unmapped, reverse, mate_reverse, read1, read2, secondary, qcfail, duplicate, supplementary = split_flag(flag)

                
                if (aligned + inserted) <= 50:
                    soft_clipped_50 += soft_clipped
                    aligned_50 += aligned

                if (aligned + inserted) <= 60:
                    soft_clipped_60 += soft_clipped
                    aligned_60 += aligned

                if (aligned + inserted) <= 70:
                    soft_clipped_70 += soft_clipped
                    aligned_70 += aligned

                if mqual >= 20:
                    mqual20_count += 1

                secondary_count += int(secondary > 0)
                supplementary_count += int(supplementary > 0)
                if duplicate:
                    duplicate_count += 1
                    duplicate_supplement_count += int(supplementary > 0)
               
               
                pa, il, pcril1, pcril2, n = detect_adapter(seq)
                poly_a += (pa > 0)
                illumina_adapter += (il > 0)
                pcr_illumina_adapter_1 += (pcril1 > 0)
                pcr_illumina_adapter_2 += (pcril2 > 0)
                nextera += (n > 0)

                if supplementary > 0:
                    suptag = None
                    for tag in tags:
                        if tag.startswith('SA:Z'):
                            suptag = tag
                            break
                    assert not suptag is None
                    prim_contig = suptag.split(':')[2].split(',')[0]
                    if prim_contig != my_contig:
                        sup_diffchrom_count += 1

            with sys.stdout as f:
                w = csv.writer(f,delimiter='\t')
                w.writerow(('#unmapped_ratio','mqual20_ratio','secondary_ratio', 'supplementary_ratio','sup_diffchrom_ratio','duplicate_ratio','duplicate_supplement_ratio','soft_clipped_bp_ratio', 'aligned_bp_ratio', 'inserted_bp_ratio', 'deleted_bp_ratio', 'total_bp',\
                                'soft_clipped_bp_ratio_filter_50', 'soft_clipped_bp_ratio_filter_60', 'soft_clipped_bp_ratio_filter_70',
                                'aligned_bp_ratio_filter_50', 'aligned_bp_ratio_filter_60', 'aligned_bp_ratio_filter_70', 'poly_a', 'illumina_adapter', 'pcr_adapter_1', 'pcr_adapter_2', 'nextera'
                                ))
                
                w.writerow((float(unmapped_count)/ float(total_count), mqual20_count / float(total_count), secondary_count / float(total_count), supplementary_count / float(total_count),float(sup_diffchrom_count) / max(float(supplementary_count),1.0), duplicate_count / float(total_count), duplicate_supplement_count / max(float(supplementary_count),1.0),\
                        softclipped_base_count/float(base_count), aligned_base_count/float(base_count), inserted_base_count/float(base_count), deleted_base_count/float(base_count), base_count,
                        (softclipped_base_count - soft_clipped_50)/float(base_count), (softclipped_base_count - soft_clipped_60) / float(base_count), (softclipped_base_count - soft_clipped_70)/float(base_count) ,
                        (aligned_base_count - aligned_50)/float(base_count), (aligned_base_count - aligned_60) / float(base_count), (aligned_base_count - aligned_70)/float(base_count) ,
                        poly_a / float(total_count), illumina_adapter / float(total_count), pcr_illumina_adapter_1 / float(total_count), pcr_illumina_adapter_2 / float(total_count), nextera / float(total_count)
                        ))
                    
                    
        elif mode == 'clean':
            lastname = ''
            querygroup = []
            with sys.stdout as out_file:
                with gzip.open(sys.argv[2],'wt') as filter_file:
                    counter = 0
                    for row in reader:
                        counter += 1
                        if (counter % 10000):
                            sys.stderr.write(str(counter) + '\n')
                            sys.stderr.flush()

                        if row[0][0] == '@':
                            out_file.write('\t'.join(row) + '\n')
                            filter_file.write('\t'.join(row) + '\n')
                            continue
                        name = row[0]
                        if lastname == name:
                            querygroup.append(row)
                        else:
                            if querygroup:
                                res_kept, res_filtered = process_querygroup(querygroup)
                                res_kept = '\n'.join(['\t'.join(xrow) for xrow in res_kept]) + '\n'
                                out_file.write(res_kept)
                                res_filtered = '\n'.join(['\t'.join(xrow) for xrow in res_filtered]) + '\n'
                                filter_file.write(res_filtered)
                            lastname = name
                            querygroup = [row]

