import sys
import re
import os
import copy

from Bio.Seq import Seq

convert_tag_type = {'i':int, 'f': float, 'A':str, 'Z':str}
extract_tag_type = {int:'i', float:'f', str:'Z'}

MAX_READ_LENGTH = 1000
MIN_ALIGN_LENGTH = 30
class BamRecord:
    def __init__(self, row):
        self.row = row
        self.qname = row[0]
        self.flag = int(row[1])
        self.rname = row[2]
        self.pos = row[3]
        self.mapq = row[4]
        self.cigar = row[5]
        self.rnext = row[6]
        self.pnext = row[7]
        self.tlen = row[8]
        self.seq = row[9]
        self.qual = row[10]
        tags = [e.split(':') for e in row[11:]]
        self.tags = {e[0]: convert_tag_type[e[1]](':'.join(e[2:])) for e in tags}

    def __repr__(self):
        return self.toSamRecord()

    def toSamRecord(self):
        tags = ['%s:%s:%s' % (tag, extract_tag_type[tag_value.__class__], tag_value)  for tag, tag_value in self.tags.items()]
        return '\t'.join([self.qname, str(self.flag), self.rname, self.pos, self.mapq, self.cigar, self.rnext, self.pnext, self.tlen, self.seq, self.qual] + tags)

    def getTagValue(self, name, default=''):
        return self.tags.get(name, default)

    def setTagValue(self,name, value):
        self.tags[name] = value

    def setTagValues(self, **kwargs):
        self.tags.update(kwargs)
   
    def is_primary(self):
        return not bool(self.flag & 0x900)

    def is_supplementary(self):
        return bool(self.flag & 0x800)

    def is_unmapped(self):
        if self.flag & 0x4:
            return True
        else:
            return False
    
    def is_reversed(self):
        return bool(self.flag & 0x10)

    def unmap(self):
        if not self.is_unmapped():
            read = self.copy()
            read.cigar = '*'
            read.pos = '0'
            read.flag = read.flag | 0x4
            read.mapq = '0'
            read.rname = '*'
            if read.is_reversed():
                read.seq = str(Seq(read.orig_seq()).reverse_complement())
                read.qual = read.orig_qual()[::-1]
                read.flag = read.flag & (~0x10)
            else:
                read.seq = read.orig_seq()
                read.qual = read.orig_qual()

            for k in ['YB','YQ','ZB','ZQ']:
                if k in read.tags:
                    del read.tags[k]

            return read
        return self            

            
    
    def gen_sa_tag(self):
        nm = self.tags.get('NM',0)
        strand = '-' if self.is_reversed() else '+'
        return '%s,%s,%s,%s,%s,%d;' % (self.rname, self.pos, strand, self.cigar, self.mapq, nm)
   

    def get_aligned_read_length(self):
        if self.cigar == '*':
            return 0
        if not 'S' in self.cigar and not 'H' in self.cigar:
            return len(self.seq)
        pos = 0
        for (ctype, length) in self.split_cigar():
            if ctype == 'M' or ctype == 'I' or ctype == '=' or ctype == 'X':
                pos += length
        return pos

    def get_orig_read_length(self):
        if self.cigar == '*' or not 'H' in self.cigar:
            return len(self.seq) 
        else:
            pos = 0
            for (ctype, length) in self.split_cigar():
                if ctype == 'S' or ctype == 'M' or ctype == 'I' or ctype == '=' or ctype == 'X' or ctype == 'H':
                    pos += length
            return pos 

    def split_cigar(self, orig_orientation=False, merge_clips=False):
        #orig_orientation=True:  return cigar that is un-reversed
        #merge_clips: merge S and H clips into a unified SH clip.
        if self.cigar == '*':
            return []
        cigar = re.split('([MIDNSHP=X]{1,1})',self.cigar)
        res = [(ctype, int(length)) for length, ctype in zip(cigar[::2], cigar[1::2])]
        if orig_orientation and self.is_reversed():
            res = res[::-1]
        if merge_clips:
            pos = 0
            while res[-1][0] in ['S','H']:
                pos += res.pop()[1]
            res.append(('SH',pos))
            pos = 0
            while res[0][0] in ['S','H']:
                pos += res[0][1]
                res = res[1:]
            res = [('SH',pos)] + res
        return res
            
    
    def get_read_position(self, orig_orientation=False):
        #get alignment start and end position on read and alignment start and end position on reference.
        #if orig_orientation is True, read positions are w.r.t. to original read orientation (reverse-complement undone)
        scigar = self.split_cigar()
        reverse = self.is_reversed()
        seq_length = self.get_orig_read_length()

        if not scigar or self.pos == '*':
            return {'read_spos':0, 'read_epos':seq_length, 'ref_spos':0, 'ref_epos':0, 'rname':self.rname, 'reverse':reverse, 'seq_length':seq_length}

        pos = 0
        alength = 0
        minpos = []
        maxpos = []
        #numpy.zeros(len(sequse),dtype=bool)
        for ctype, length in scigar:
            if ctype == 'S' or ctype == 'H':
                pos += length
            elif ctype == 'M' or ctype == 'I' or ctype == '=' or ctype == 'X':
                if reverse and orig_orientation:
                    minpos.append(seq_length - (pos + length))
                    maxpos.append(seq_length - pos)
                else:
                    minpos.append(pos)
                    maxpos.append(pos + length)
                pos += length
                
                if ctype != 'I':
                    alength += length
                    
            elif ctype == 'D' or ctype == 'N' or ctype == 'P':
                if ctype == 'D' or ctype == 'N':
                    alength += length
                pass
        read_startpos = min(minpos)
        read_endpos = max(maxpos)
        ref_startpos = int(self.pos)
        ref_endpos = int(self.pos) + alength

        return {'read_spos':read_startpos, 'read_epos':read_endpos, 'ref_spos':ref_startpos, 'ref_epos':ref_endpos, 'rname':self.rname, 'reverse':reverse, 'seq_length':seq_length}


    def clip_start(self, read_start_pos, orig_orientation=True):
        if not orig_orientation and self.is_reversed():
            return self.clip_end(self.get_orig_read_length() - read_start_pos, orig_orientation=True)

        scigar = self.split_cigar()
        scigar = scigar[::-1] if self.is_reversed() else scigar

        scigar, ref_pos_shift, read_pos = prune_cigar_end(scigar[::-1], read_start_pos)

        if all([e[0] in ['S','H'] for e in scigar]):
            return self.unmap()
        else:
            scigar = [('H', read_pos)] + scigar[::-1]
            read = self.copy()
            if not self.is_reversed():
                read.pos = str(int(self.pos) + ref_pos_shift)
                read.cigar = join_cigar(scigar)
                read.tags['YB'] = read.orig_seq()[:read_pos]
                read.tags['YQ'] = read.orig_qual()[:read_pos]
            else:
                read.cigar = join_cigar(scigar[::-1])
                read.tags['ZB'] = read.orig_seq()[-read_pos:]
                read.tags['ZQ'] = read.orig_qual()[-read_pos:]
        return read
    
    def clip_end(self, read_end_pos, orig_orientation=True):
        if not orig_orientation and self.is_reversed():
            return self.clip_start(self.get_orig_read_length() - read_end_pos, orig_orientation=True)
        
        scigar = self.split_cigar()
        scigar = scigar[::-1] if self.is_reversed() else scigar

        scigar, ref_pos_shift, read_pos = prune_cigar_end(scigar, self.get_orig_read_length() - read_end_pos)

        if all([e[0] in ['S','H'] for e in scigar]):
            return self.unmap()
        else:
            scigar = scigar + [('H', read_pos)]
            read = self.copy()
            if not self.is_reversed():
                read.cigar = join_cigar(scigar)
                read.tags['ZB'] = read.orig_seq()[-read_pos:]
                read.tags['ZQ'] = read.orig_qual()[-read_pos:]
            else:
                read.pos = str(int(self.pos) + ref_pos_shift)
                read.cigar = join_cigar(scigar[::-1])
                read.tags['YB'] = read.orig_seq()[:read_pos]
                read.tags['YQ'] = read.orig_qual()[:read_pos]
        return read

    def copy(self):
        read = copy.copy(self)
        read.tags = dict(self.tags)
        return read


    def annotate_orig_sequence(self, seq, qual):
        read = self.copy()

        scigar = read.split_cigar()
        reverse = read.is_reversed()

        if reverse:
            seq = str(Seq(seq).reverse_complement())
            qual = qual[::-1]
        
        startpos = 0
        endpos = 0
        if scigar[0][0] == 'H':
            startpos = scigar[0][1]
        if scigar[-1][0] == 'H':
            endpos = scigar[-1][1]
        assert seq[startpos:self.get_orig_read_length() -endpos]  == read.seq
        assert qual[startpos:self.get_orig_read_length() - endpos] == read.qual
        
        if startpos:
            read.tags['YB'] = seq[:startpos]
            read.tags['YQ'] = qual[:startpos]
        if endpos:
            read.tags['ZB'] = seq[-endpos:]
            read.tags['ZQ'] = qual[-endpos:]
        return read


    def orig_seq(self, orig_orientation=False):
        res = ''.join([self.tags.get('YB',''), self.seq, self.tags.get('ZB','')])
        if orig_orientation and self.is_reversed():
            res = str(Seq(res).reverse_complement())
        return res            
    
    def orig_qual(self, orig_orientation=False):
        res = ''.join([self.tags.get('YQ',''), self.qual, self.tags.get('ZQ','')])
        if orig_orientation and self.is_reversed():
            res = res[::-1]
        return res


        

    def ref_to_read_pos(self, start, stop):
        # get read positions and cigar corresponding to a reference range
        #FIXME: check handling of hard clipping
        scigar = self.split_cigar()
        reverse = self.is_reversed()

        if not scigar or self.pos == '*':
            return (0, 0, True, True)

        readpos = 0
        refpos = int(self.pos)
        
        read_start = None
        read_stop = None

        fuzzy_start = False #ref start position is within deletion event
        fuzzy_stop = False #ref stop position is within deletion evven


        pcigar = []
        mcigar = []
        acigar = []
        scigar = scigar[::-1] #reverse, so we can use it as a a stack to pop from and append to
        while scigar:
            ctype, length = scigar.pop()
            if ctype == 'S' or ctype == 'I':
                readpos += length
            elif ctype == 'M' or ctype == '=' or ctype == 'X':
                if refpos+length >= start and read_start is None:
                    read_start = readpos + (start - refpos)
                    rem_length = (readpos + length - read_start)
                    if start - refpos:
                        pcigar.append((ctype, (start - refpos)))
                    scigar.append((ctype, rem_length))
                    readpos = read_start
                    refpos = start
                    break
                else:
                    readpos += length
                    refpos += length
            elif ctype == 'D' or ctype == 'N':
                if refpos+length >= start and read_start is None:
                    read_start = readpos
                    fuzzy_start = True
                    rem_length = (refpos + length) - start
                    if start - refpos:
                        pcigar.append((ctype, (start - refpos)))
                    scigar.append((ctype, rem_length))
                    refpos = start
                    break
                else:
                    refpos += length
            elif ctype == 'P' or ctype == 'H':
                pass
            pcigar.append((ctype, length))

        while scigar:
            ctype, length = scigar.pop()
            if ctype == 'S' or ctype == 'I':
                readpos += length
            elif ctype == 'M' or ctype == '=' or ctype == 'X':
                if refpos + length >= stop and read_stop is None:
                    read_stop = readpos + (stop - refpos)
                    rem_length = read_stop - readpos
                    if rem_length:
                       mcigar.append((ctype, rem_length))
                    scigar.append((ctype, (length - rem_length)))
                    refpos = stop
                    break
                else: 
                    readpos += length
                    refpos += length
            elif ctype == 'D' or ctype == 'N':
                if refpos +  length >= stop and read_stop is None:
                    read_stop = readpos
                    fuzzy_stop = True
                    rem_length = stop - refpos
                    if rem_length:
                        mcigar.append((ctype, rem_length))
                    scigar.append((ctype, (refpos + length - stop)))
                    refpos = stop
                    break
                else:
                    refpos += length
            elif ctype == 'P' or ctype =='H':
                pass
            if length:                
                mcigar.append((ctype, length))

        acigar = [e for e in scigar[::-1] if e[1] != 0]
        pcigar = join_cigar(pcigar)
        mcigar = join_cigar(mcigar)                
        acigar = join_cigar(acigar)
        assert not read_stop is None,str((scigar, refpos, start, stop))

        assert not read_start is None, str((scigar, self.pos, start, stop))
        if reverse:
            seq_length = self.get_orig_read_length()
            xread_start, xread_stop = seq_length - read_stop, seq_length - read_start
        else:
            xread_start = read_start
            xread_stop = read_stop
        return (read_start, read_stop, xread_start, xread_stop, fuzzy_start, fuzzy_stop, pcigar, mcigar, acigar)                
 
def read_overlap(rpos1, rpos2):
    if rpos1['rname'] != rpos2['rname'] or rpos1['reverse'] == rpos2['reverse']: #they need to read into one another
        return (False, 0,0, True)
    
    start_overlap = max(rpos1['ref_spos'], rpos2['ref_spos'])
    end_overlap = min(rpos1['ref_epos'], rpos2['ref_epos'])
    if not rpos1['reverse']:
        inward = rpos1['ref_spos'] < rpos2['ref_spos']
    else:
        inward = rpos2['ref_spos'] < rpos1['ref_spos']

    return (start_overlap <= end_overlap, start_overlap, end_overlap, inward)

def join_cigar(split_cigar):
    return ''.join(['%d%s' % (length,ctype) for ctype,length in split_cigar])

def generate_simple_consensus(seq1, seq2, qual1, qual2):
    mq = []
    sq = []
    mismatch = 0
    match = 0
    for s1,s2,q1,q2 in zip(seq1, seq2, qual1, qual2):
        e1 = 10.0**(-(ord(q1) - 33) / 10.0)
        e2 = 10.0**(-(ord(q2) - 33) / 10.0)
        if s1 == s2:
            match += 1
            ef = e1 * e2
            sf = s1
        else:
            mismatch += 1
            if e1 < e2:
                ef = (1 - e1) * e2
                sf = s1
            elif e1 == e2:
                sf = 'N'
                ef = 1.0
            else:
                ef = (1 - e2) * e1
                sf = s2
        q = chr(max(min(int(round(-(math.log10(ef) * 10.0) + 33)),73),2))
        sq.append(sf)
        mq.append(q)
    return (''.join(sq), ''.join(mq), mismatch, mismatch + match)


def prune_cigar_end(scigar, pos):
    read_pos = 0
    ref_pos = 0
    #pop till we are at read pos, then readd clip
    while scigar and read_pos < pos:
        ctype, length = scigar.pop()
        if ctype == 'S' or ctype == 'H':
            read_pos += length
        elif ctype == 'M' or ctype == "I" or ctype == '=' or ctype == 'X':
            if read_pos + length <= pos: #complete removal
                read_pos += length
                if not ctype == 'I':
                    ref_pos += length

            else:    #partial removal
                remove_length = pos - read_pos
                if not ctype == 'I':
                    ref_pos += remove_length
                read_pos += remove_length

                scigar.append((ctype, length - remove_length))
        elif ctype == 'D' or ctype == 'N':
            ref_pos += length
        elif ctype == 'P':
            pass

    while scigar:
        ctype,length = scigar.pop()
        if ctype == 'D' or ctype == 'N':
            ref_pos += length
        elif ctype == 'P':
            pass
        elif ctype == 'I':
            ref_pos += length
            read_pos += length
        else:
            scigar.append((ctype,length))
            break

    return (scigar, ref_pos, read_pos)
   
