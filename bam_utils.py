import sys
import re
import os


convert_tag_type = {'i':int, 'f': float, 'A':str, 'Z':str}
extract_tag_type = {int:'i', float:'f', str:'Z'}

MAX_READ_LENGTH = 1000

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

    def toSamRecord(self):
        tags = ['%s:%s:%s' % (tag, extract_tag_type[tag_value.__class__], tag_value)  for tag, tag_value in self.tags.items()]
        return '\t'.join([self.qname, str(self.flag), self.rname, self.pos, self.mapq, self.cigar, self.rnext, self.pnext, self.tlen, self.seq, self.qual] + tags)

    def getTagValue(self, name, default=''):
        return self.tags.get(name, default)

    def setTagValue(self,name, value):
        self.tags[name] = value

    def setTagValues(self, **kwargs):
        self.tags.update(kwargs)
    
    def is_unmapped(self):
        return (self.cigar == '*')
    
    def is_reversed(self):
        return bool(self.flag & 0x10)
    

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

    def split_cigar(self):
        if self.cigar == '*':
            return []
        cigar = re.split('([MIDNSHP=X]{1,1})',self.cigar)
        return [(ctype, int(length)) for length, ctype in zip(cigar[::2], cigar[1::2])]
    
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
            if ctype == 'S':
                pos += length
            elif ctype == 'M' or ctype == 'I' or ctype == '=' or ctype == 'X':
                if reverse and orig_orientation:
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
        read_startpos = min(minpos)
        read_endpos = max(maxpos)
        ref_startpos = int(self.pos)
        ref_endpos = int(self.pos) + alength

        return {'read_spos':read_startpos, 'read_epos':read_endpos, 'ref_spos':ref_startpos, 'ref_epos':ref_endpos, 'rname':self.rname, 'reverse':reverse, 'seq_length':seq_length}

    def ref_to_read_pos(self, start, stop):
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
        scigar = scigar[::-1]
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
        q = chr(min(int(round(-(math.log10(ef) * 10.0) + 33)),73))
        sq.append(sf)
        mq.append(q)
    return (''.join(sq), ''.join(mq), mismatch, mismatch + match)

   
