# bam_dechimer Python ↔ C line-by-line algorithmic comparison

This report walks through the Python implementation (`scripts/bam_dechimer.py` + helpers in `scripts/bam_utils.py`) and the C port (`scripts/bam_dechimer.c`) and, for each algorithmic step, shows the relevant Python and C snippets followed by a conclusion on logical equivalence.

The comparison focuses on the core pipeline:
- Read grouping and per-fragment processing
- Supplementary pruning (`rg_prune`)
- Primary pair dechimering (`dechimer`) for the simple 2-record case
- Record filtering and promotion (`record_filter`, `rg_filter`)
- SA tag recomputation (`update_sa_tag`)
- Mate flag and tag updates (`update_mate_flags`, `update_mate_tags`)
- Proper-pair flag recomputation (`update_pair_flag`)
- Clipping helpers (`clip_start`, `clip_end`), read span computation (`get_read_position`)

Where helpful, small behavioral differences are highlighted explicitly.

---

## Top-level streaming, grouping, and stats writing

- **Python** (excerpts)
```python
# bam_dechimer.py
if args.i == '-':
    pipe_in = sys.stdin
else:
    in_process = subprocess.Popen(['samtools','view','-h', args.i], stdout=subprocess.PIPE, universal_newlines=True)
    pipe_in = in_process.stdout
...
reader = csv.reader(pipe_in, delimiter='\t', quoting=csv.QUOTE_NONE)
...
# strip /1,/2 by BamRecord ctor
row = BamRecord(row)
if lastname == row.qname:
    querygroup.append(row)
else:
    if querygroup:
        res_kept = process(querygroup, stats, args)
        out_file.write('\n'.join([x.toSamRecord() for x in res_kept]) + '\n')
    lastname = row.qname
    querygroup = [row]
...
# final group + stats
with open(args.s,'wt') as sf:
    for k,v in stats.items():
        sf.write('%s\t %d\n' % (k, v))
```

- **C**
```c
// bam_dechimer.c
if(strcmp(in,"-")){
    size_t n=strlen(in)+32; char *cmd=malloc(n); snprintf(cmd,n,"samtools view -h %s", in);
    pipe_proc = popen(cmd, "r"); free(cmd); inF = pipe_proc;
} else { inF = stdin; }
...
// grouping by qname, stripping /1 or /2
split_tab(tmp,f,&nf);
size_t L=strlen(f[0]); if(L>2 && f[0][L-2]=='/' && (f[0][L-1]=='1'||f[0][L-1]=='2')) f[0][L-2]='\0';
...
parse_bam_line_full(pline, &group_rows[ng]) // builds BamRow
...
process_fragment(outF, &sm, group_rows, ng, ...);
...
FILE *sf=fopen(stats, "w");
for(int i=0;i<sm.n;i++) fprintf(sf, "%s\t %lld\n", sm.keys[i], sm.vals[i]);
```

- **Conclusion**
- Equivalent behavior for streaming and grouping by qname.
- Headers are passed through in both.
- Stats file format matches.

---

## Core data model mapping

- **Python**
```python
class BamRecord:
    # fields: qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags
    def is_unmapped(self): return bool(self.flag & 0x4)
    def is_reversed(self): return bool(self.flag & 0x10)
    def is_supplementary(self): return bool(self.flag & 0x800)
    def get_align_length(self):  # M,=,X only
        pos=0
        for t,l in self.split_cigar():
            if t in ['M','=','X']: pos += l
        return pos
```

- **C**
```c
typedef struct { ... int flag; char *cigar; char *seq; TagList tags; } BamRow;
static int is_unmapped(const BamRow *r){ return (r->flag & 0x4)!=0; }
static int is_reversed(const BamRow *r){ return (r->flag & 0x10)!=0; }
static int is_supplementary(const BamRow *r){ return (r->flag & 0x800)!=0; }
static int get_aligned_length_MeqX(const char *c){
    // counts M,=,X only
}
```

- **Conclusion**
- Field mapping and helpers are equivalent.
- Alignment length for filtering uses M/= /X only in both.

---

## Helper: get_read_position (orig read orientation)

- **Python** (`bam_utils.py:get_read_position`)
```python
def get_read_position(self, orig_orientation=False):
    scigar = self.split_cigar()
    reverse = self.is_reversed()
    seq_length = self.get_orig_read_length()
    if not scigar or self.pos == '*':
        return { 'read_spos':0, 'read_epos':seq_length, 'ref_spos':0, 'ref_epos':0, ...}
    pos=0; alength=0; minpos=[]; maxpos=[]
    for ctype, length in scigar:
        if ctype in ['S','H']: pos += length
        elif ctype in ['M','I','=','X']:
            if reverse and orig_orientation:
                minpos.append(seq_length - (pos + length))
                maxpos.append(seq_length - pos)
            else:
                minpos.append(pos)
                maxpos.append(pos + length)
            pos += length
            if ctype != 'I': alength += length
        elif ctype in ['D','N','P']:
            if ctype in ['D','N']: alength += length
    return {'read_spos':min(minpos), 'read_epos':max(maxpos), 'ref_spos':int(self.pos), 'ref_epos':int(self.pos)+alength, ...}
```

- **C** (`get_read_position`)
```c
static void get_read_position(const BamRow *r, int orig_orientation,
    int *read_spos, int *read_epos, int *ref_spos, int *ref_epos){
  int reverse=is_reversed(r);
  int seq_length=get_orig_read_length_from_cigar(r->cigar, strlen(r->seq));
  if(n==0 || r->pos==0){ *read_spos=0; *read_epos=seq_length; *ref_spos=0; *ref_epos=0; return; }
  int pos=0, alength=0, minpos=1e9, maxpos=0;
  for(each cigar tok){
    if t in S/H: pos += l;
    else if t in M/I/=/X: // compute s,e possibly reversed
      if(reverse && orig_orientation){ s = seq_length - (pos + l); e = seq_length - pos; }
      else { s = pos; e = pos + l; }
      update minpos/maxpos; pos += l; if(t!='I') alength += l;
    else if t in D/N/P: if t in D/N: alength += l;
  }
  *read_spos=minpos==1e9?0:minpos; *read_epos=maxpos; *ref_spos=r->pos; *ref_epos=r->pos + alength;
}
```

- **Conclusion**
- Equivalent logic. Python uses half-open endpoint internally (pos+length); C mirrors that.

---

## Helper: prune_cigar_end and clipping

- **Python** (`bam_utils.py:prune_cigar_end`)
```python
def prune_cigar_end(scigar, pos):
    # pops from end until 'pos' is reached in read space
    # returns (new_scigar, ref_pos_shift, read_pos, hard_clip_start)
```

- **C** (`prune_cigar_end` + `clip_start`/`clip_end`)
```c
static void prune_cigar_end(CigarVec *scigar, int pos, int *ref_pos_out, int *read_pos_out, int *hard_clip_start){ ... }
static BamRow clip_start(const BamRow *in, int read_start_pos, int orig_orientation){ ... }
static BamRow clip_end(const BamRow *in, int read_end_pos, int orig_orientation){ ... }
```

- Key behaviors compared per path:
  - Handle reversed reads by operating in original orientation for clipping.
  - If result would be only S/H, unmap the read (restoring original SEQ/QUAL and removing YB/YQ/ZB/ZQ).
  - Prepend/append hard clip element to CIGAR and adjust tags:
    - Forward `clip_start`: set YB/YQ from prefix, trim seq/qual; advance POS by reference shift.
    - Forward `clip_end`: set ZB/ZQ from suffix, trim seq/qual; POS unchanged.
    - Reverse `clip_start`: set ZB/ZQ, trim from end; CIGAR reversed, POS unchanged.
    - Reverse `clip_end`: set YB/YQ, trim from start; advance POS by reference shift.
  - If a leading hard-clip existed and no stored prefix/postfix tags, subtract `hard_clip_start` from counted read_pos (both ports implement this adjustment).

- **Conclusion**
- C mirrors Python’s clipping semantics and tag handling closely for all four orientation/direction cases.

---

## Function: rg_prune (prune overlapping supplementary alignments)

- **Python** (`bam_dechimer.py:rg_prune`)
```python
stats[f"{prefix}_has_supplementary_alignments"] += 1
reads = [primary] + supplementary
rpos = [r.get_read_position(orig_orientation=True) for r in reads]
...
# For each read i, shrink [start, stop] by other reads' spans
if (stop - start) < args.min_align_length:
    if read.is_primary(): primary = read.unmap(); stats[f"{prefix}_primary_unmapped_in_pruning"] += 1
    else: stats[f"{prefix}_sup_discarded_in_pruning"] += 1
    continue
if rpos[i]['read_spos'] != start or (start > 0 and before):
    read = read.clip_start(start, orig_orientation=True); stats[f"{prefix}_start_pruned"] += 1
if rpos[i]['read_epos'] != stop or (stop > 0 and after):
    read = read.clip_end(stop, orig_orientation=True); stats[f"{prefix}_end_pruned"] += 1
...
if not nreadgroup['supplementary']:
    del nreadgroup['supplementary']; stats[f"{prefix}_all_sup_discarded_in_pruning"] += 1
```

- **C** (`rg_prune`)
```c
smap_inc(sm, "{prefix}_has_supplementary_alignments", 1);
// compute start/stop from get_read_position(..., orig=1)
...
if(keep_len < min_align_length){
  if(r==g->primary){ bam_replace_inplace(r, unmap_record(r,1,1)); smap_inc(sm, "{prefix}_primary_unmapped_in_pruning", 1); }
  else { smap_inc(sm, "{prefix}_sup_discarded_in_pruning", 1); }
  continue;
}
if(start[i] != s || (s > 0 && before[i])) { r=clip_start(...); smap_inc(sm, "{prefix}_start_pruned", 1); }
if(stop[i] != e  || (e > 0 && after[i])) { r=clip_end(...);   smap_inc(sm, "{prefix}_end_pruned", 1); }
...
if(nns==0){ smap_inc(sm, "{prefix}_all_sup_discarded_in_pruning", 1); g->supplementary=NULL; g->nsup=0; }
```

- **Conclusion**
- Equivalent pruning logic, including the “before/after” triggers to force terminal clipping even if boundaries did not numerically change.

---

## Function: dechimer (simple two-primary case)

- **Python** (`bam_dechimer.py:dechimer`)
```python
if not 'supplementary' in reads1 and not 'supplementary' in reads2:
    if 'S' not in r1.cigar and 'S' not in r2.cigar: return False, False
    cigar1 = r1.split_cigar(orig_orientation=True, merge_clips=True)
    cigar2 = r2.split_cigar(orig_orientation=True, merge_clips=True)
    if not r1.is_unmapped() and not r2.is_unmapped():
        same_chrom = (r1.rname == r2.rname)
        diff = abs(int(r1.pos) - int(r2.pos)) if same_chrom else huge
        if diff >= args.max_read_dist:
            if cigar1[-1][0] == 'SH': r1 = r1.clip_end(cigar1[-1][1], orig_orientation=True); stats['read1_dechimer_clip'] += 1; modified1=True
            if cigar2[-1][0] == 'SH': r2 = r2.clip_end(cigar2[-1][1], orig_orientation=True); stats['read2_dechimer_clip'] += 1; modified2=True
        if args.loose_ends:
            if cigar1[0][0] == 'SH': r1 = r1.clip_start(cigar1[0][1], orig_orientation=True); stats['read1_loose_end_clip'] += 1; modified1=True
            if cigar2[0][0] == 'SH': r2 = r2.clip_start(cigar2[0][1], orig_orientation=True); stats['read2_loose_end_clip'] += 1; modified2=True
    elif r1.is_unmapped() and r2.is_unmapped(): pass
    elif r1.is_unmapped():
        if cigar2[-1][0] == 'SH': ...
        if args.loose_ends and cigar2[0][0] == 'SH': ...
    else:
        if cigar1[-1][0] == 'SH': ...
        if args.loose_ends and cigar1[0][0] == 'SH': ...
    reads1['primary']=r1; reads2['primary']=r2
```

- **C** (`dechimer`)
```c
if(r1->nsup==0 && r2->nsup==0){
    if(!strchr(p1->cigar,'S') && !strchr(p2->cigar,'S')) return;
    int f1=0,b1=0,f2=0,b2=0; cigar_front_back_clip(p1,&f1,&b1); cigar_front_back_clip(p2,&f2,&b2);
    if(!is_unmapped(p1) && !is_unmapped(p2)){
        long long diff = (strcmp(p1->rname,p2->rname)==0) ? llabs((long long)p1->pos - (long long)p2->pos) : huge;
        if(diff >= max_read_dist){ if(b1>0){ p1=clip_end(..., orl1-b1, 1); smap_inc(..."read1_dechimer_clip"); *mod1=1; }
                                   if(b2>0){ p2=clip_end(..., orl2-b2, 1); smap_inc(..."read2_dechimer_clip"); *mod2=1; } }
        if(loose_ends){ if(f1>0){ p1=clip_start(..., f1,1); smap_inc(..."read1_loose_end_clip"); *mod1=1; }
                         if(f2>0){ p2=clip_start(..., f2,1); smap_inc(..."read2_loose_end_clip"); *mod2=1; } }
    } else if(is_unmapped(p1) && is_unmapped(p2)){ /* nothing */ }
    else if(is_unmapped(p1)){ if(b2>0){ ... } if(loose_ends && f2>0){ ... } }
    else { if(b1>0){ ... } if(loose_ends && f1>0){ ... } }
}
```

- **Conclusion**
- Equivalent behavior. Python’s merged `SH` tokens correspond to C’s summed leading/trailing S/H counts.

---

## Function: record_filter (min alignment length / very short reads)

- **Python**
```python
def record_filter(read, stats, prefix, args):
    if len(read.seq) <= args.min_align_length or read.get_align_length() < args.min_align_length:
        stats[f'{prefix}_min_align_length_unmap'] += 1
        return read.unmap()
    else:
        return read
```

- **C**
```c
static BamRow record_filter(const BamRow *in, StatsMap *sm, const char *prefix, int min_align_length){
    if(strlen(in->seq) <= min_align_length || get_aligned_length_MeqX(in->cigar) < min_align_length){
        smap_inc(sm, "{prefix}_min_align_length_unmap", 1);
        return unmap_record(in,1,1);
    }
    return bam_copy(in);
}
```

- **Conclusion**
- Equivalent logic and thresholds. Both use only M/= /X for aligned length.

---

## Function: rg_filter (apply record_filter, possibly promote a supplementary)

- **Python**
```python
def rg_filter(readgroup, stats, prefix, args):
    readgroup['primary'] = record_filter(...)
    if 'supplementary' in readgroup:
        new = [record_filter(e, ...) for e in readgroup['supplementary']]
        new = [e for e in new if not e.is_unmapped()]
        if readgroup['primary'].is_unmapped() and len(new) > 0:
            new_primary = new.pop(0)
            stats[f'{prefix}_promote_supplementary'] += 1
            new_primary = new_primary.annotate_orig_sequence(
                readgroup['primary'].orig_seq(orig_orientation=True),
                readgroup['primary'].orig_qual(orig_orientation=True))
            new_primary.flag = new_primary.flag & (~0x800)
            readgroup['primary'] = new_primary
        readgroup['supplementary'] = new if len(new)>0 else del
```

- **C**
```c
static void rg_filter(RG *g, StatsMap *sm, const char *prefix, int min_align_length){
    BamRow np = record_filter(g->primary, sm, prefix, min_align_length); bam_replace_inplace(g->primary, np);
    if(g->nsup>0){
        // apply record_filter to sups and drop unmapped
        if(is_unmapped(g->primary) && nk>0){
            BamRow *np1=keep[0];
            char *os=orig_seq(g->primary,1,1,1), *oq=orig_qual(g->primary,1,1,1);
            BamRow ann=annotate_orig_sequence(np1, os, oq); ann.flag &= (~0x800);
            bam_replace_inplace(np1, ann); g->primary = np1; smap_inc(sm, "{prefix}_promote_supplementary", 1);
            // remaining kept sups shifted
        }
        // update g->supplementary to kept list or clear
    }
}
```

- **Conclusion**
- Equivalent behavior, including promotion with annotation of YB/YQ/ZB/ZQ based on original primary’s full sequence and quality.

---

## Function: update_sa_tag

- **Python**
```python
def update_sa_tag(readgroup):
    if 'supplementary' in readgroup:
        sa_primary = readgroup['primary'].gen_sa_tag()
        sa_sups = [read.gen_sa_tag() for read in readgroup['supplementary']]
        if sa_sups:
            readgroup['primary'].tags['SA'] = ''.join(sa_sups)
        for pos, read in enumerate(readgroup['supplementary']):
            sa_sups_tmp = list(sa_sups); del sa_sups_tmp[pos]
            read.tags['SA'] = sa_primary + ''.join(sa_sups_tmp)
```

- **C**
```c
static void update_sa_tag(RG *rg){
    // primary gets concat of all sups; each sup gets primary + other sups
    if(rg->nsup<=0 || !rg->primary) return;
    char *concat=...; tlist_set(&rg->primary->tags, "SA", 'Z', concat or remove);
    for(int pos=0; pos<rg->nsup; pos++){ build primary + other sups; tlist_set(&rg->supplementary[pos]->tags, "SA", 'Z', buf); }
}
```

- **Conclusion**
- Equivalent SA tag construction.

---

## Function: update_mate_flags

- **Python**
```python
def update_mate_flags(reads1,reads2):
    r1_flags=r2_flags=0
    if p1.is_unmapped(): r2_flags |= 0x8
    if p2.is_unmapped(): r1_flags |= 0x8
    if p1.is_reversed(): r2_flags |= 0x20
    if p2.is_reversed(): r1_flags |= 0x20
    if not p1.is_unmapped() and not p2.is_unmapped(): r1_flags|=r2_flags|=0x2
    cancel = ~(0x8|0x20|0x2)
    primary and sups := (flag & cancel) | r*_flags
```

- **C**
```c
static void update_mate_flags(RG *r1, RG *r2){ ... same bit logic ... apply to primary and all sups }
```

- **Conclusion**
- Equivalent bit-twiddling and propagation to supplementary records.

---

## Function: update_mate_tags (MC, RNEXT/PNEXT)

- **Python**
```python
def update_mate_tags(read1, read2):
    if p2.is_unmapped(): p1.tags.pop('MC',None) else: p1.tags['MC']=p2.cigar
    if p1.is_unmapped(): p2.tags.pop('MC',None) else: p2.tags['MC']=p1.cigar
    p1.rnext=p2.rname; p2.rnext=p1.rname; p1.pnext=p2.pos; p2.pnext=p1.pos
    # propagate MC/rnext/pnext to sups accordingly
```

- **C**
```c
static void update_mate_tags(RG *r1, RG *r2){
    if(is_unmapped(p2)) tlist_del(&p1->tags, "MC"); else tlist_set(&p1->tags, "MC", 'Z', p2->cigar);
    if(is_unmapped(p1)) tlist_del(&p2->tags, "MC"); else tlist_set(&p2->tags, "MC", 'Z', p1->cigar);
    p1->rnext=sdup(p2->rname); p2->rnext=sdup(p1->rname);
    p1->pnext=sdup(fmt(p2->pos)); p2->pnext=sdup(fmt(p1->pos));
    // propagate to sups
}
```

- **Conclusion**
- Equivalent mate tag and mate position updates.
- Both ports only update mates when the fragment was modified (in their respective `process` blocks).

---

## Function: update_pair_flag (properly paired)

- **Python**
```python
def update_pair_flag(reads1, reads2, stats, args):
    proper=True
    if p1.rname != p2.rname or p1.rname == "*": proper=False
    else:
        diff = abs(int(p1.pos) - int(p2.pos))
        if diff > args.max_read_dist: proper=False
    cancel_flag = ~0x2; new_flag = 0x2 if proper else 0
    if p1.is_proper() or p2.is_proper():
        stats['proper_pair_downgrade'] = stats.get('proper_pair_downgrade',0) + 1
    p1.flag = (p1.flag & cancel_flag) | new_flag
    p2.flag = (p2.flag & cancel_flag) | new_flag
```

- **C**
```c
static void update_pair_flag(RG *r1, RG *r2, StatsMap *sm, int max_read_dist){
    int proper=1; int old_proper = ((p1->flag & 0x2) || (p2->flag & 0x2));
    if(strcmp(p1->rname,p2->rname)!=0 || strcmp(p1->rname,"*")==0) proper=0;
    else { int diff = abs(p1->pos - p2->pos); if(diff > max_read_dist) proper=0; }
    int nf = proper? 0x2: 0;
    if(old_proper && !proper) smap_inc(sm, "proper_pair_downgrade", 1);
    p1->flag = (p1->flag & ~0x2) | nf; p2->flag = (p2->flag & ~0x2) | nf;
}
```

- **Conclusion**
- Both recompute proper pairing with the same condition and reset bit 0x2 on both reads.
- The Python implementation increments `proper_pair_downgrade` whenever either read was flagged proper at entry, regardless of whether the recomputed pair remains proper. The C version increments only on actual downgrade (old proper -> new not proper). This is a deliberate, small divergence that more precisely counts true downgrades.

---

## Function: process/process_fragment (per fragment pipeline)

- **Python**
```python
def process(querygroup, stats, args):
    stats['alignment_counter'] += len(querygroup)
    stats['fragment_counter'] += 1
    reads1 = process_readgroup([row for row in querygroup if row.flag & 0x40])
    reads2 = process_readgroup([row for row in querygroup if row.flag & 0x80])
    modified1=modified2=False
    if 'supplementary' in reads1: reads1, p = rg_prune(...); modified1 |= p
    if 'supplementary' in reads2: reads2, p = rg_prune(...); modified2 |= p
    rmodified1, rmodified2 = dechimer(reads1, reads2, stats, args); modified1|=rmodified1; modified2|=rmodified2
    if modified1: reads1 = rg_filter(...); update_sa_tag(reads1)
    if modified2: reads2 = rg_filter(...); update_sa_tag(reads2)
    if modified1 or modified2:
        stats['fragment_modified'] += 1
        update_mate_flags(reads1,reads2)
        update_mate_tags(reads1,reads2)
        update_pair_flag(reads1,reads2, stats,args)
    # (optional) validate
    # rebuild and return ordered list; filter out empty-seq
```

- **C**
```c
static void process_fragment(FILE *out, StatsMap *sm, BamRow *group, int ng, ...){
    smap_inc(sm, "alignment_counter", ng); smap_inc(sm, "fragment_counter", 1);
    split rows into read1/read2 arrays; process_readgroup_assign(...)
    int mod1=0, mod2=0;
    if(r1.nsup>0) mod1 |= rg_prune(&r1, ...);
    if(r2.nsup>0) mod2 |= rg_prune(&r2, ...);
    int d1=0,d2=0; dechimer(&r1,&r2, ...,&d1,&d2); mod1|=d1; mod2|=d2;
    if(mod1){ rg_filter(&r1, ...); update_sa_tag(&r1); }
    if(mod2){ rg_filter(&r2, ...); update_sa_tag(&r2); }
    if(mod1 || mod2){ smap_inc(sm, "fragment_modified", 1); update_mate_flags(&r1,&r2); update_mate_tags(&r1,&r2); update_pair_flag(&r1,&r2, ...); }
    // validation intentionally omitted
    // build output list in same order; print only if seq not empty
}
```

- **Conclusion**
- Equivalent pipeline ordering and gating on modification. C intentionally skips the optional validation step.

---

## Summary of differences

- **Proper-pair downgrade counting**
  - Python increments `proper_pair_downgrade` if either read had bit 0x2 set on entry, even if the recomputed pair remains proper.
  - C increments only when old proper transitions to new not-proper, which more precisely tracks true downgrades.

- **Validation**
  - Python performs invariant checks when `-d` is not provided.
  - C omits runtime validation (commented as intentionally skipped).

- **All other core behaviors**
  - Supplementary pruning, terminal clipping triggers, min alignment length unmapping, promotion of a supplementary to primary, SA/MC/tag updates, mate flag propagation, and the simple two-primary dechimer logic are functionally equivalent.

---

## Conclusion

Except for the intentionally different `proper_pair_downgrade` statistic and the omission of validation in the C port, the C implementation reproduces the Python algorithm’s logic step-by-step for the present codebase:
- Clipping and tag maintenance match across orientations and ends.
- Pruning and unmapping thresholds match (aligned length M/= /X only).
- Mate flags/tags are only recalculated when a fragment was modified.
- Output ordering and filtering of empty-SEQ records match.
