#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

/* Minimal skeleton for bam_dechimer: argument parsing, SAM streaming, grouping by qname, and stats writing.
Full dechimer logic (prune/clip/update tags) will be added incrementally. */

static ssize_t safe_getline(char **lineptr, size_t *n, FILE *stream){
    ssize_t r = getline(lineptr, n, stream);
    if(r<=0) return r;
    if((*lineptr)[r-1]=='\n'){ (*lineptr)[r-1]='\0'; --r; }
    if(r>0 && (*lineptr)[r-1]=='\r'){ (*lineptr)[r-1]='\0'; --r; }
    return r;
}

static void split_tab(char *s, char **f, int *nf){
    int i=0; char *p=s; f[i++]=p;
    while(*p){ if(*p=='\t'){ *p='\0'; f[i++]=p+1; } p++; if(i>=256) break; }
    *nf=i;
}

static char* sdup(const char *s){ if(!s) return NULL; size_t n=strlen(s); char *p=(char*)malloc(n+1); memcpy(p,s,n+1); return p; }

/* Simple string->int map preserving insertion order */
typedef struct { char **keys; long long *vals; int n, cap; } StatsMap;
static void smap_init(StatsMap *m){ m->keys=NULL; m->vals=NULL; m->n=0; m->cap=0; }
static void smap_free(StatsMap *m){ if(!m) return; for(int i=0;i<m->n;i++) free(m->keys[i]); free(m->keys); free(m->vals); m->keys=NULL; m->vals=NULL; m->n=m->cap=0; }
static void smap_inc(StatsMap *m, const char *key, long long by){
    for(int i=0;i<m->n;i++){ if(strcmp(m->keys[i],key)==0){ m->vals[i]+=by; return; } }
    if(m->n==m->cap){ m->cap = m->cap? m->cap*2:64; m->keys=(char**)realloc(m->keys,sizeof(char*)*m->cap); m->vals=(long long*)realloc(m->vals,sizeof(long long)*m->cap); }
    m->keys[m->n] = sdup(key); m->vals[m->n] = by; m->n++;
}

/* Tag utilities */
typedef struct { char key[3]; char type; char *val; } Tag;
typedef struct { Tag *v; int n; int cap; } TagList;
static void tlist_init(TagList *tl){ tl->v=NULL; tl->n=0; tl->cap=0; }
static void tlist_free(TagList *tl){ if(!tl) return; for(int i=0;i<tl->n;i++){ free(tl->v[i].val); } free(tl->v); tl->v=NULL; tl->n=tl->cap=0; }
static int tlist_find(const TagList *tl, const char *key){ for(int i=0;i<tl->n;i++){ if(tl->v[i].key[0]==key[0] && tl->v[i].key[1]==key[1]) return i; } return -1; }
static const char* tlist_get(const TagList *tl, const char *key){ int i=tlist_find(tl,key); return i>=0? tl->v[i].val : NULL; }
static void tlist_del(TagList *tl, const char *key){ int i=tlist_find(tl,key); if(i>=0){ free(tl->v[i].val); for(int j=i+1;j<tl->n;j++) tl->v[j-1]=tl->v[j]; tl->n--; } }
static void tlist_set(TagList *tl, const char *key, char type, const char *val){ int i=tlist_find(tl,key); if(i<0){ if(tl->n==tl->cap){ tl->cap = tl->cap? tl->cap*2:8; tl->v=(Tag*)realloc(tl->v,sizeof(Tag)*tl->cap); } tl->v[tl->n].key[0]=key[0]; tl->v[tl->n].key[1]=key[1]; tl->v[tl->n].key[2]='\0'; tl->v[tl->n].type=type; tl->v[tl->n].val=sdup(val?val:""); tl->n++; } else { free(tl->v[i].val); tl->v[i].type=type; tl->v[i].val=sdup(val?val:""); } }
static int parse_tag_field(const char *s, char key[3], char *type, char **val){ size_t L=strlen(s); if(L<5) return 0; if(s[2]!=':'||s[4]!=':') return 0; key[0]=s[0]; key[1]=s[1]; key[2]='\0'; *type=s[3]; *val=sdup(s+5); return 1; }

/* CIGAR helpers */
typedef struct { char t; int l; } CigTok;
static int has_char(const char *s, char c){ if(!s) return 0; for(const char *p=s; *p; ++p){ if(*p==c) return 1; } return 0; }
static int cigar_parse(const char *c, CigTok *v, int maxv){ int n=0, val=0; if(!c||c[0]=='*') return 0; for(const char *p=c; *p; ++p){ if(isdigit((unsigned char)*p)){ val=val*10+(*p-'0'); } else { if(n<maxv){ v[n].t=*p; v[n].l=val; n++; } val=0; } } return n; }
static char* cigar_join(const CigTok *v, int n){ size_t sz=1; for(int i=0;i<n;i++){ int tmp=v[i].l; int digits=1; while(tmp>=10){digits++; tmp/=10;} sz += digits+1; } char *s=(char*)malloc(sz); size_t k=0; for(int i=0;i<n;i++){ k += sprintf(s+k, "%d%c", v[i].l, v[i].t); } s[k]='\0'; return s; }
static void cigar_reverse(CigTok *a, int n){ for(int i=0,j=n-1;i<j;i++,j--){ CigTok t=a[i]; a[i]=a[j]; a[j]=t; } }
static int get_orig_read_length_from_cigar(const char *c, int seqlen){ if(!c||c[0]=='*' || !has_char(c,'H')) return seqlen; CigTok a[256]; int n=cigar_parse(c,a,256); int pos=0; for(int i=0;i<n;i++){ char t=a[i].t; int l=a[i].l; if(t=='S'||t=='M'||t=='I'||t=='='||t=='X'||t=='H') pos += l; } return pos; }
static int get_aligned_length_MeqX(const char *c){ if(!c||c[0]=='*') return 0; CigTok a[256]; int n=cigar_parse(c,a,256); int pos=0; for(int i=0;i<n;i++){ char t=a[i].t; int l=a[i].l; if(t=='M'||t=='='||t=='X') pos += l; } return pos; }

/* prune_cigar_end */
typedef struct { CigTok *v; int n; int cap; } CigarVec;
static void cvec_from(const CigTok *a, int n, CigarVec *cv){ cv->v=(CigTok*)malloc(sizeof(CigTok)*n); memcpy(cv->v,a,sizeof(CigTok)*n); cv->n=n; cv->cap=n; }
static void cvec_free(CigarVec *cv){ if(cv && cv->v){ free(cv->v); cv->v=NULL; } }
static void cvec_push(CigarVec *cv, char t, int l){ if(cv->n==cv->cap){ cv->cap=cv->cap? cv->cap*2:16; cv->v=(CigTok*)realloc(cv->v,sizeof(CigTok)*cv->cap); } cv->v[cv->n].t=t; cv->v[cv->n].l=l; cv->n++; }
static void prune_cigar_end(CigarVec *scigar, int pos, int *ref_pos_out, int *read_pos_out, int *hard_clip_start){
    int read_pos=0, ref_pos=0; *hard_clip_start=0;
    if(scigar->n>0 && scigar->v[scigar->n-1].t=='H') *hard_clip_start=scigar->v[scigar->n-1].l;
    while(scigar->n>0 && read_pos < pos){
        CigTok c = scigar->v[--scigar->n];
        if(c.t=='S' || c.t=='H') read_pos += c.l;
        else if(c.t=='M' || c.t=='I' || c.t=='=' || c.t=='X'){
            if(read_pos + c.l <= pos){ read_pos += c.l; if(c.t!='I') ref_pos += c.l; }
            else { int rem = pos - read_pos; if(c.t!='I') ref_pos += rem; read_pos += rem; scigar->v[scigar->n++] = (CigTok){c.t, c.l - rem}; }
        } else if(c.t=='D' || c.t=='N') ref_pos += c.l;
        else if(c.t=='P'){}
    }
    while(scigar->n>0){
        CigTok c = scigar->v[--scigar->n];
        if(c.t=='D' || c.t=='N') ref_pos += c.l;
        else if(c.t=='P'){}
        else if(c.t=='I') read_pos += c.l;
        else { scigar->v[scigar->n++] = c; break; }
    }
    *ref_pos_out = ref_pos; *read_pos_out = read_pos;
}

/* BamRow parsing/printing */
typedef struct {
    char *qname; int flag; char *rname; int pos; char *mapq; char *cigar; char *rnext; char *pnext; char *tlen; char *seq; char *qual;
    int read1, read2, unmapped, reversed, secondary, supplementary;
    TagList tags;
} BamRow;

static int parse_bam_line_full(char *line, BamRow *out){ char *f[256]; int nf=0; split_tab(line,f,&nf); if(nf<11) return 0; size_t L=strlen(f[0]); if(L>2 && f[0][L-2]=='/' && (f[0][L-1]=='1'||f[0][L-1]=='2')) f[0][L-2]='\0';
    out->qname=sdup(f[0]); out->flag=(int)strtol(f[1],NULL,10); out->rname=sdup(f[2]); out->pos=(int)strtol(f[3],NULL,10); out->mapq=sdup(f[4]); out->cigar=sdup(f[5]); out->rnext=sdup(f[6]); out->pnext=sdup(f[7]); out->tlen=sdup(f[8]); out->seq=sdup(f[9]); out->qual=sdup(f[10]);
    out->read1=(out->flag&0x40)!=0; out->read2=(out->flag&0x80)!=0; out->unmapped=(out->flag&0x4)!=0; out->reversed=(out->flag&0x10)!=0; out->secondary=(out->flag&0x100)!=0; out->supplementary=(out->flag&0x800)!=0; tlist_init(&out->tags);
    for(int i=11;i<nf;i++){ char key[3]; char tp; char *val=NULL; if(parse_tag_field(f[i], key, &tp, &val)){ tlist_set(&out->tags, key, tp, val); free(val); } }
    return 1;
}
static void print_bam(FILE *out, const BamRow *r){ fprintf(out, "%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s", r->qname, r->flag, r->rname, r->pos, r->mapq, r->cigar, r->rnext, r->pnext, r->tlen, r->seq, r->qual); for(int i=0;i<r->tags.n;i++) fprintf(out, "\t%.*s:%c:%s", 2, r->tags.v[i].key, r->tags.v[i].type, r->tags.v[i].val); fputc('\n', out); }
static void bam_free(BamRow *r){ if(!r) return; free(r->qname); free(r->rname); free(r->mapq); free(r->cigar); free(r->rnext); free(r->pnext); free(r->tlen); free(r->seq); free(r->qual); tlist_free(&r->tags); memset(r,0,sizeof(*r)); }
static BamRow bam_copy(const BamRow *r){ BamRow x=*r; x.qname=sdup(r->qname); x.rname=sdup(r->rname); x.mapq=sdup(r->mapq); x.cigar=sdup(r->cigar); x.rnext=sdup(r->rnext); x.pnext=sdup(r->pnext); x.tlen=sdup(r->tlen); x.seq=sdup(r->seq); x.qual=sdup(r->qual); tlist_init(&x.tags); for(int i=0;i<r->tags.n;i++){ tlist_set(&x.tags, r->tags.v[i].key, r->tags.v[i].type, r->tags.v[i].val); } return x; }
static void bam_replace_inplace(BamRow *dst, BamRow src){ bam_free(dst); *dst = src; }
static int is_unmapped(const BamRow *r){ return (r->flag & 0x4)!=0; }
static int is_reversed(const BamRow *r){ return (r->flag & 0x10)!=0; }
static int is_supplementary(const BamRow *r){ return (r->flag & 0x800)!=0; }

/* Orig seq/qual assembly (YB/ZB + seq + ZB postfix) */
static inline void rev_inplace(char *s){ if(!s) return; size_t i=0, j=strlen(s); if(j==0) return; j--; while(i<j){ char t=s[i]; s[i]=s[j]; s[j]=t; ++i; --j; } }
static inline char rc_base(char c){ switch(c){ case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; case 'a': return 't'; case 'c': return 'g'; case 'g': return 'c'; case 't': return 'a'; default: return c; } }
static void revcomp_inplace(char *s){ if(!s) return; size_t i=0, j=strlen(s); if(j==0) return; j--; while(i<j){ char a=s[i], b=s[j]; s[i]=rc_base(b); s[j]=rc_base(a); ++i; --j; } if(i==j) s[i]=rc_base(s[i]); }
static char* orig_seq(const BamRow *r, int orig_orientation, int include_prefix, int include_postfix){ const char *yb=tlist_get(&r->tags,"YB"); const char *zb=tlist_get(&r->tags,"ZB"); size_t L1= include_prefix && yb? strlen(yb):0; size_t L2=strlen(r->seq); size_t L3= include_postfix && zb? strlen(zb):0; size_t N=L1+L2+L3; char *s=(char*)malloc(N+1); size_t k=0; if(L1){ memcpy(s+k,yb,L1); k+=L1; } memcpy(s+k,r->seq,L2); k+=L2; if(L3){ memcpy(s+k,zb,L3); k+=L3; } s[k]='\0'; if(orig_orientation && is_reversed(r)) revcomp_inplace(s); return s; }
static char* orig_qual(const BamRow *r, int orig_orientation, int include_prefix, int include_postfix){ const char *yq=tlist_get(&r->tags,"YQ"); const char *zq=tlist_get(&r->tags,"ZQ"); size_t L1= include_prefix && yq? strlen(yq):0; size_t L2=strlen(r->qual); size_t L3= include_postfix && zq? strlen(zq):0; size_t N=L1+L2+L3; char *s=(char*)malloc(N+1); size_t k=0; if(L1){ memcpy(s+k,yq,L1); k+=L1; } memcpy(s+k,r->qual,L2); k+=L2; if(L3){ memcpy(s+k,zq,L3); k+=L3; } s[k]='\0'; if(orig_orientation && is_reversed(r)) rev_inplace(s); return s; }

/* clip_start/clip_end need get_read_position and unmap_record; add simple versions */
static BamRow unmap_record(const BamRow *in, int restore_seq, int orig_orientation){ BamRow r=bam_copy(in); free(r.cigar); r.cigar=sdup("*"); r.pos=0; free(r.rname); r.rname=sdup("*"); free(r.mapq); r.mapq=sdup("0"); r.flag &= (~0x10); r.flag |= 0x4; if(restore_seq){ char *s=orig_seq(in, orig_orientation, 1,1); char *q=orig_qual(in, orig_orientation, 1,1); free(r.seq); r.seq=s; free(r.qual); r.qual=q; tlist_del(&r.tags,"YB"); tlist_del(&r.tags,"YQ"); tlist_del(&r.tags,"ZB"); tlist_del(&r.tags,"ZQ"); } return r; }

/* Compute read start/end in original read orientation (half-open) and reference range */
static void get_read_position(const BamRow *r, int orig_orientation, int *read_spos, int *read_epos, int *ref_spos, int *ref_epos){
    CigTok a[512]; int n=cigar_parse(r->cigar,a,512); int reverse=is_reversed(r); int seq_length=get_orig_read_length_from_cigar(r->cigar,(int)strlen(r->seq));
    if(n==0 || r->pos==0){ *read_spos=0; *read_epos=seq_length; *ref_spos=0; *ref_epos=0; return; }
    int pos=0; int alength=0; int minpos=(int)1e9; int maxpos=0;
    for(int i=0;i<n;i++){
        char t=a[i].t; int l=a[i].l;
        if(t=='S' || t=='H'){ pos += l; }
        else if(t=='M' || t=='I' || t=='=' || t=='X'){
            int s = pos; int e = pos + l;
            if(reverse && orig_orientation){ s = seq_length - (pos + l); e = seq_length - pos; }
            if(s < minpos) minpos = s;
            if(e > maxpos) maxpos = e;
            pos += l;
            if(t!='I') alength += l;
        } else if(t=='D' || t=='N' || t=='P'){ if(t=='D' || t=='N') alength += l; }
    }
    *read_spos = minpos==(int)1e9? 0: minpos; *read_epos = maxpos; *ref_spos = r->pos; *ref_epos = r->pos + alength;
}

/* Clip from start/end in original orientation, mirroring Python BamRecord.clip_start/clip_end */
static BamRow clip_start(const BamRow *in, int read_start_pos, int orig_orientation){
    if(is_unmapped(in)) return bam_copy(in);
    if(!orig_orientation && is_reversed(in)){
        int orl=get_orig_read_length_from_cigar(in->cigar,(int)strlen(in->seq));
        return clip_start(in, orl - read_start_pos, 1);
    }
    CigTok a[512]; int n=cigar_parse(in->cigar,a,512); if(is_reversed(in)) cigar_reverse(a,n);
    /* operate on reversed copy to pop from end */
    CigTok arev[512]; for(int i=0;i<n;i++) arev[i]=a[i]; cigar_reverse(arev,n);
    CigarVec cv; cvec_from(arev,n,&cv); int ref_shift=0, read_pos=0, hard_clip_start=0;
    prune_cigar_end(&cv, read_start_pos, &ref_shift, &read_pos, &hard_clip_start);
    int only_clip=1; for(int i=0;i<cv.n;i++){ if(!(cv.v[i].t=='S'||cv.v[i].t=='H')){ only_clip=0; break; } }
    if(only_clip){ BamRow u=unmap_record(in,1,1); cvec_free(&cv); return u; }
    CigarVec nx; nx.v=NULL; nx.n=0; nx.cap=0; cvec_push(&nx,'H', read_pos); for(int i=cv.n-1;i>=0;i--) cvec_push(&nx, cv.v[i].t, cv.v[i].l);
    if(hard_clip_start && !tlist_get(&in->tags,"YB") && !tlist_get(&in->tags,"ZB")) read_pos -= hard_clip_start;
    BamRow r=bam_copy(in);
    if(!is_reversed(in)){
        char *s=orig_seq(&r,0,1,0); char *q=orig_qual(&r,0,1,0);
        r.pos = in->pos + ref_shift; free(r.cigar); r.cigar=cigar_join(nx.v,nx.n);
        if(!is_supplementary(&r) && read_pos>0){ char *yb=(char*)malloc(read_pos+1); memcpy(yb,s,read_pos); yb[read_pos]='\0'; tlist_set(&r.tags,"YB",'Z',yb); free(yb); char *yq=(char*)malloc(read_pos+1); memcpy(yq,q,read_pos); yq[read_pos]='\0'; tlist_set(&r.tags,"YQ",'Z',yq); free(yq); } else { tlist_del(&r.tags,"YB"); tlist_del(&r.tags,"YQ"); }
        int Ls=(int)strlen(s); if(read_pos>Ls) read_pos=Ls; memmove(s, s+read_pos, Ls-read_pos+1); memmove(q, q+read_pos, (int)strlen(q)-read_pos+1); free(r.seq); r.seq=s; free(r.qual); r.qual=q;
    } else {
        char *s=orig_seq(&r,0,0,1); char *q=orig_qual(&r,0,0,1);
        cigar_reverse(nx.v,nx.n); free(r.cigar); r.cigar=cigar_join(nx.v,nx.n);
        int Ls=(int)strlen(s); if(read_pos>Ls) read_pos=Ls; if(!is_supplementary(&r) && read_pos>0){ char *zb=(char*)malloc(read_pos+1); memcpy(zb,s+Ls-read_pos,read_pos); zb[read_pos]='\0'; tlist_set(&r.tags,"ZB",'Z',zb); free(zb); char *zq=(char*)malloc(read_pos+1); memcpy(zq,q+Ls-read_pos,read_pos); zq[read_pos]='\0'; tlist_set(&r.tags,"ZQ",'Z',zq); free(zq); } else { tlist_del(&r.tags,"ZB"); tlist_del(&r.tags,"ZQ"); }
        s[Ls-read_pos]='\0'; q[Ls-read_pos]='\0'; free(r.seq); r.seq=s; free(r.qual); r.qual=q;
    }
    cvec_free(&cv); cvec_free(&nx); return r;
}

static BamRow clip_end(const BamRow *in, int read_end_pos, int orig_orientation){
    if(is_unmapped(in)) return bam_copy(in);
    if(!orig_orientation && is_reversed(in)){
        int orl=get_orig_read_length_from_cigar(in->cigar,(int)strlen(in->seq));
        return clip_end(in, orl - read_end_pos, 1);
    }
    CigTok a[512]; int n=cigar_parse(in->cigar,a,512); if(is_reversed(in)) cigar_reverse(a,n);
    int orl=get_orig_read_length_from_cigar(in->cigar,(int)strlen(in->seq)); int from_end = orl - read_end_pos;
    CigarVec cv; cvec_from(a,n,&cv); int ref_shift=0, read_pos=0, hard_clip_start=0;
    prune_cigar_end(&cv, from_end, &ref_shift, &read_pos, &hard_clip_start);
    int only_clip=1; for(int i=0;i<cv.n;i++){ if(!(cv.v[i].t=='S'||cv.v[i].t=='H')){ only_clip=0; break; } }
    if(only_clip){ BamRow u=unmap_record(in,1,1); cvec_free(&cv); return u; }
    CigarVec nx; nx.v=NULL; nx.n=0; nx.cap=0; for(int i=0;i<cv.n;i++) cvec_push(&nx, cv.v[i].t, cv.v[i].l); cvec_push(&nx,'H', read_pos);
    if(hard_clip_start && !tlist_get(&in->tags,"YB") && !tlist_get(&in->tags,"ZB")) read_pos -= hard_clip_start;
    BamRow r=bam_copy(in);
    if(!is_reversed(in)){
        char *s=orig_seq(&r,0,0,1); char *q=orig_qual(&r,0,0,1);
        free(r.cigar); r.cigar=cigar_join(nx.v,nx.n);
        int Ls=(int)strlen(s); if(read_pos>Ls) read_pos=Ls; if(!is_supplementary(&r) && read_pos>0){ char *zb=(char*)malloc(read_pos+1); memcpy(zb,s+Ls-read_pos,read_pos); zb[read_pos]='\0'; tlist_set(&r.tags,"ZB",'Z',zb); free(zb); char *zq=(char*)malloc(read_pos+1); memcpy(zq,q+Ls-read_pos,read_pos); zq[read_pos]='\0'; tlist_set(&r.tags,"ZQ",'Z',zq); free(zq); } else { tlist_del(&r.tags,"ZB"); tlist_del(&r.tags,"ZQ"); }
        s[Ls-read_pos]='\0'; q[Ls-read_pos]='\0'; free(r.seq); r.seq=s; free(r.qual); r.qual=q;
    } else {
        char *s=orig_seq(&r,0,1,0); char *q=orig_qual(&r,0,1,0);
        r.pos = in->pos + ref_shift; cigar_reverse(nx.v,nx.n); free(r.cigar); r.cigar=cigar_join(nx.v,nx.n);
        if(!is_supplementary(&r) && read_pos>0){ char *yb=(char*)malloc(read_pos+1); memcpy(yb,s,read_pos); yb[read_pos]='\0'; tlist_set(&r.tags,"YB",'Z',yb); free(yb); char *yq=(char*)malloc(read_pos+1); memcpy(yq,q,read_pos); yq[read_pos]='\0'; tlist_set(&r.tags,"YQ",'Z',yq); free(yq); } else { tlist_del(&r.tags,"YB"); tlist_del(&r.tags,"YQ"); }
        int Ls=(int)strlen(s); if(read_pos>Ls) read_pos=Ls; memmove(s, s+read_pos, Ls-read_pos+1); memmove(q, q+read_pos, (int)strlen(q)-read_pos+1); free(r.seq); r.seq=s; free(r.qual); r.qual=q;
    }
    cvec_free(&cv); cvec_free(&nx); return r;
}

/* annotate_orig_sequence: set YB/YQ/ZB/ZQ based on hardclip lengths */
static BamRow annotate_orig_sequence(const BamRow *in, const char *seq, const char *qual){
    BamRow r=bam_copy(in);
    CigTok a[512]; int n=cigar_parse(r.cigar,a,512);
    int rev=is_reversed(&r);
    char *s=sdup(seq); char *q=sdup(qual);
    if(rev){ revcomp_inplace(s); rev_inplace(q); }
    int startpos=0, endpos=0; if(n>0 && a[0].t=='H') startpos=a[0].l; if(n>0 && a[n-1].t=='H') endpos=a[n-1].l;
    if(startpos>0){ char *yb=(char*)malloc(startpos+1); memcpy(yb,s,startpos); yb[startpos]='\0'; tlist_set(&r.tags,"YB",'Z',yb); free(yb); char *yq=(char*)malloc(startpos+1); memcpy(yq,q,startpos); yq[startpos]='\0'; tlist_set(&r.tags,"YQ",'Z',yq); free(yq); }
    if(endpos>0){ int Ls=(int)strlen(s); char *zb=(char*)malloc(endpos+1); memcpy(zb,s+Ls-endpos,endpos); zb[endpos]='\0'; tlist_set(&r.tags,"ZB",'Z',zb); free(zb); int Lq=(int)strlen(q); char *zq=(char*)malloc(endpos+1); memcpy(zq,q+Lq-endpos,endpos); zq[endpos]='\0'; tlist_set(&r.tags,"ZQ",'Z',zq); free(zq); }
    free(s); free(q); return r;
}

/* Readgroup structure and helpers */
typedef struct {
    BamRow *primary;
    BamRow **supplementary; int nsup;
    BamRow **secondary; int nsec;
    BamRow **supplementary_secondary; int nss;
} RG;

static void rg_init(RG *g){ memset(g,0,sizeof(*g)); }
static void rg_free(RG *g){ if(!g) return; free(g->supplementary); free(g->secondary); free(g->supplementary_secondary); memset(g,0,sizeof(*g)); }

static void process_readgroup_assign(BamRow **rows, int n, RG *out){ rg_init(out); for(int i=0;i<n;i++){ BamRow *row=rows[i]; if(row->secondary){ if(row->supplementary){ out->supplementary_secondary=(BamRow**)realloc(out->supplementary_secondary, sizeof(BamRow*)*(out->nss+1)); out->supplementary_secondary[out->nss++]=row; } else { out->secondary=(BamRow**)realloc(out->secondary, sizeof(BamRow*)*(out->nsec+1)); out->secondary[out->nsec++]=row; } } else if(row->supplementary){ out->supplementary=(BamRow**)realloc(out->supplementary, sizeof(BamRow*)*(out->nsup+1)); out->supplementary[out->nsup++]=row; } else { out->primary=row; } } }

static char* gen_sa_tag(const BamRow *r){
    const char *nm = tlist_get(&r->tags, "NM");
    int nmv = nm? atoi(nm): 0;
    char strand = is_reversed(r)? '-': '+';
    
    // Prefer original MAPQ from XQ tag if present
    const char *xq = tlist_get(&r->tags, "XQ");
    const char *mq = xq ? xq : (r->mapq ? r->mapq : "0");
    
    size_t need = strlen(r->rname) + strlen(r->cigar) + strlen(mq) + 64;
    char *s = (char*)malloc(need);
    sprintf(s, "%s,%d,%c,%s,%s,%d;",
        r->rname, r->pos, strand, r->cigar, mq, nmv);
        return s;
    }
    
    static void update_sa_tag(RG *rg){ if(!rg || rg->nsup<=0 || !rg->primary) return; // primary SA = concat all sups
        char *concat=NULL; size_t clen=0; for(int i=0;i<rg->nsup;i++){ char *z=gen_sa_tag(rg->supplementary[i]); size_t zl=strlen(z); concat=(char*)realloc(concat, clen+zl+1); memcpy(concat+clen,z,zl); clen+=zl; concat[clen]='\0'; free(z); }
        if(clen>0) tlist_set(&rg->primary->tags, "SA", 'Z', concat); else tlist_del(&rg->primary->tags, "SA");
        for(int pos=0;pos<rg->nsup;pos++){
            char *p=gen_sa_tag(rg->primary); size_t pl=strlen(p); size_t tot=pl; for(int j=0;j<rg->nsup;j++){ if(j==pos) continue; char *z=gen_sa_tag(rg->supplementary[j]); tot += strlen(z); free(z); }
            char *buf=(char*)malloc(tot+1); size_t off=0; memcpy(buf+off,p,pl); off+=pl; for(int j=0;j<rg->nsup;j++){ if(j==pos) continue; char *z=gen_sa_tag(rg->supplementary[j]); size_t zl=strlen(z); memcpy(buf+off,z,zl); off+=zl; free(z); } buf[off]='\0'; tlist_set(&rg->supplementary[pos]->tags, "SA", 'Z', buf); free(buf); free(p);
        }
        free(concat);
    }
    
    static BamRow record_filter(const BamRow *in, StatsMap *sm, const char *prefix, int min_align_length){
        int seq_len = (int)strlen(in->seq);
        int aln_len = get_aligned_length_MeqX(in->cigar);
        if(seq_len <= min_align_length || aln_len < min_align_length){
            char key[128];
            snprintf(key,sizeof(key), "%s_min_align_length_unmap", prefix);
            smap_inc(sm, key, 1);
            BamRow u=unmap_record(in,1,1);
            return u;
        }
        BamRow cp=bam_copy(in);
        return cp;
    }
    
    static void rg_filter(RG *g, StatsMap *sm, const char *prefix, int min_align_length){ 
        
        if(!g || !g->primary) return; 
        BamRow *orig_primary = g->primary; 
        BamRow np = record_filter(orig_primary, sm, prefix, min_align_length); 
        bam_replace_inplace(orig_primary, np); 
        
        if(g->nsup>0){ // filter sups
            BamRow **keep=(BamRow**)malloc(sizeof(BamRow*)*g->nsup);
            int nk=0; 
            for(int i=0;i<g->nsup;i++)
            { 
                BamRow *s=g->supplementary[i]; 
                BamRow ns = record_filter(s, sm, prefix, min_align_length); 
                bam_replace_inplace(s, ns); 
                if(!is_unmapped(s)) keep[nk++]=s; 
            }
            if(is_unmapped(orig_primary) && nk>0)
            { 
                BamRow *np1=keep[0]; 
                char *os=orig_seq(g->primary,1,1,1); 
                char *oq=orig_qual(g->primary,1,1,1); 
                BamRow ann=annotate_orig_sequence(np1, os, oq); 
                free(os); free(oq); 
                
                ann.flag &= (~0x800); 
                bam_replace_inplace(np1, ann); 
                g->primary = np1; 
                if(orig_primary && orig_primary != np1)
                { 
                    BamRow cleared = bam_copy(orig_primary);
                    if(cleared.seq){
                        cleared.seq[0] = '\0';  /* mark as "logically removed" */
                    }
                    bam_replace_inplace(orig_primary, cleared);
                } 
                char k2[128]; 
                snprintf(k2,sizeof(k2), "%s_promote_supplementary", prefix); 
                smap_inc(sm, k2, 1); 
                for(int i=1;i<nk;i++){ 
                    keep[i-1]=keep[i]; 
                } 
                nk--; 
            }
            
            free(g->supplementary); 
            if(nk>0){ 
                g->supplementary=(BamRow**)malloc(sizeof(BamRow*)*nk); 
                memcpy(g->supplementary, keep, sizeof(BamRow*)*nk); g->nsup=nk; 
            } 
            else 
            { 
                g->supplementary=NULL; g->nsup=0; 
            }
            free(keep);
        }
    }

// update mate flags for primary and supplementary
static void update_mate_flags(RG *r1, RG *r2){
    BamRow *p1=r1->primary, *p2=r2->primary;
    /* Ignore "dead" primaries that were cleared after supplementary promotion */
    if(p1 && p1->seq && p1->seq[0]=='\0') p1 = NULL;
    if(p2 && p2->seq && p2->seq[0]=='\0') p2 = NULL;

    int r1f=0, r2f=0;
    int p1_unmapped = (p1==NULL) || is_unmapped(p1);
    int p2_unmapped = (p2==NULL) || is_unmapped(p2);
    int p1_rev = (p1!=NULL) && is_reversed(p1);
    int p2_rev = (p2!=NULL) && is_reversed(p2);

    if(p1_unmapped) r2f |= 0x8;
    if(p2_unmapped) r1f |= 0x8;
    if(p1_rev) r2f |= 0x20;
    if(p2_rev) r1f |= 0x20;
    if(!p1_unmapped && !p2_unmapped){ r1f |= 0x2; r2f |= 0x2; }

    int cancel = ~(0x8|0x20|0x2);
    if(p1) p1->flag = (p1->flag & cancel) | r1f;
    if(p2) p2->flag = (p2->flag & cancel) | r2f;
    for(int i=0;i<r1->nsup;i++){ r1->supplementary[i]->flag = (r1->supplementary[i]->flag & cancel) | r1f; }
    for(int i=0;i<r2->nsup;i++){ r2->supplementary[i]->flag = (r2->supplementary[i]->flag & cancel) | r2f; }
}

static void update_pair_flag(RG *r1, RG *r2, StatsMap *sm, int max_read_dist){ BamRow *p1=r1->primary, *p2=r2->primary; int proper=1; int old_proper = ((p1->flag & 0x2) || (p2->flag & 0x2)); if(strcmp(p1->rname,p2->rname)!=0 || strcmp(p1->rname,"*")==0) proper=0; else { int diff = abs(p1->pos - p2->pos); if(diff > max_read_dist) proper=0; } int cancel=~0x2; int nf = proper? 0x2: 0; if(old_proper) smap_inc(sm, "proper_pair_downgrade", 1); p1->flag = (p1->flag & cancel) | nf; p2->flag = (p2->flag & cancel) | nf; }

static void update_mate_tags(RG *r1, RG *r2)
{ 
    BamRow *p1=r1->primary, *p2=r2->primary; 
    /* Ignore "dead" primaries cleared after supplementary promotion */
    if(p1 && p1->seq && p1->seq[0] == '\0') p1 = NULL;
    if(p2 && p2->seq && p2->seq[0] == '\0') p2 = NULL;

    /* Primary MC tags */
    if(p1){
        if(!p2 || is_unmapped(p2)) tlist_del(&p1->tags, "MC");
        else tlist_set(&p1->tags, "MC", 'Z', p2->cigar);
    }
    if(p2){
        if(!p1 || is_unmapped(p1)) tlist_del(&p2->tags, "MC");
        else tlist_set(&p2->tags, "MC", 'Z', p1->cigar);
    }

    /* rnext/pnext only when both primaries are valid */
    if(p1 && p2){
        if(p1->rnext) free(p1->rnext); 
        p1->rnext=sdup(p2->rname); 
        if(p2->rnext) free(p2->rnext);     
        p2->rnext=sdup(p1->rname); 

        char buf[64];
        snprintf(buf,sizeof(buf), "%d", p2->pos); 
        if(p1->pnext) free(p1->pnext); 
        p1->pnext=sdup(buf); 
        snprintf(buf,sizeof(buf), "%d", p1->pos); 
        if(p2->pnext) free(p2->pnext);     
        p2->pnext=sdup(buf);
    }

    /* Supplementary records on each side: mate is the other side's primary */
    for(int i=0;i<r1->nsup;i++){ 
        BamRow *s=r1->supplementary[i]; 
        if(p2){
            if(s->rnext) free(s->rnext); 
            s->rnext=sdup(p2->rname); 
            if(s->pnext) free(s->pnext); 
            char buf[64]; snprintf(buf,sizeof(buf), "%d", p2->pos); 
            s->pnext=sdup(buf); 
            if(is_unmapped(p2)) tlist_del(&s->tags, "MC"); 
            else tlist_set(&s->tags, "MC", 'Z', p2->cigar); 
        } else {
            tlist_del(&s->tags, "MC");
        }
    }
    for(int i=0;i<r2->nsup;i++){ 
        BamRow *s=r2->supplementary[i]; 
        if(p1){
            if(s->rnext) free(s->rnext); 
            s->rnext=sdup(p1->rname); 
            if(s->pnext) free(s->pnext); 
            char buf[64]; snprintf(buf,sizeof(buf), "%d", p1->pos); 
            s->pnext=sdup(buf); 
            if(is_unmapped(p1)) tlist_del(&s->tags, "MC"); 
            else tlist_set(&s->tags, "MC", 'Z', p1->cigar); 
        } else {
            tlist_del(&s->tags, "MC");
        }
    }
}

static void cigar_front_back_clip(const BamRow *r, int *front, int *back){ *front=0; *back=0; CigTok a[512]; int n=cigar_parse(r->cigar,a,512); if(n<=0) return; if(is_reversed(r)) cigar_reverse(a,n); int i=0; while(i<n && (a[i].t=='S'||a[i].t=='H')){ *front += a[i].l; i++; } int j=n-1; while(j>=0 && (a[j].t=='S'||a[j].t=='H')){ *back += a[j].l; j--; } }

static void dechimer(RG *r1, RG *r2, StatsMap *sm, int loose_ends, int max_read_dist, int *mod1, int *mod2){ *mod1=0; *mod2=0; BamRow *p1=r1->primary, *p2=r2->primary; if(r1->nsup==0 && r2->nsup==0){ if(!strchr(p1->cigar,'S') && !strchr(p2->cigar,'S')) return; int f1=0,b1=0,f2=0,b2=0; cigar_front_back_clip(p1,&f1,&b1); cigar_front_back_clip(p2,&f2,&b2); if(!is_unmapped(p1) && !is_unmapped(p2)){
    int same = strcmp(p1->rname,p2->rname)==0; long long diff = same? llabs((long long)p1->pos - (long long)p2->pos) : 100000000000000LL;
    if(diff >= max_read_dist){
        if(b1>0){
            int orl1=get_orig_read_length_from_cigar(p1->cigar,(int)strlen(p1->seq));
            BamRow nx=clip_end(p1, orl1 - b1, 1);
            bam_replace_inplace(p1, nx);
            smap_inc(sm, "read1_dechimer_clip", 1);
            *mod1=1;
        }
        if(b2>0){
            int orl2=get_orig_read_length_from_cigar(p2->cigar,(int)strlen(p2->seq));
            BamRow nx=clip_end(p2, orl2 - b2, 1);
            bam_replace_inplace(p2, nx);
            smap_inc(sm, "read2_dechimer_clip", 1);
            *mod2=1;
        }
    }
    if(loose_ends){
        if(f1>0){
            BamRow nx=clip_start(p1,f1,1);
            bam_replace_inplace(p1, nx);
            smap_inc(sm, "read1_loose_end_clip", 1);
            *mod1=1;
        }
        if(f2>0){
            BamRow nx=clip_start(p2,f2,1);
            bam_replace_inplace(p2, nx);
            smap_inc(sm, "read2_loose_end_clip", 1);
            *mod2=1;
        }
    }
} else if(is_unmapped(p1) && is_unmapped(p2)){
    // nothing
} else if(is_unmapped(p1)){
    if(b2>0){ int orl2=get_orig_read_length_from_cigar(p2->cigar,(int)strlen(p2->seq)); BamRow nx=clip_end(p2, orl2 - b2, 1); bam_replace_inplace(p2, nx); smap_inc(sm, "read2_dechimer_clip", 1); *mod2=1; }
    if(loose_ends && f2>0){ BamRow nx=clip_start(p2,f2,1); bam_replace_inplace(p2, nx); smap_inc(sm, "read2_loose_end_clip", 1); *mod2=1; }
} else {
    if(b1>0){ int orl1=get_orig_read_length_from_cigar(p1->cigar,(int)strlen(p1->seq)); BamRow nx=clip_end(p1, orl1 - b1, 1); bam_replace_inplace(p1, nx); smap_inc(sm, "read1_dechimer_clip", 1); *mod1=1; }
    if(loose_ends && f1>0){ BamRow nx=clip_start(p1,f1,1); bam_replace_inplace(p1, nx); smap_inc(sm, "read1_loose_end_clip", 1); *mod1=1; }
}
}
}

static int rg_prune(RG *g, StatsMap *sm, const char *prefix, int min_align_length){ if(!g || g->nsup<=0 || !g->primary) return 0; smap_inc(sm, (strcmp(prefix,"read1")==0? "read1_has_supplementary_alignments" : "read2_has_supplementary_alignments"), 1); int pruned=0; int n=g->nsup+1; BamRow **reads=(BamRow**)malloc(sizeof(BamRow*)*n); reads[0]=g->primary; for(int i=0;i<g->nsup;i++) reads[1+i]=g->supplementary[i]; int *start=(int*)malloc(sizeof(int)*n); int *stop=(int*)malloc(sizeof(int)*n); int *before=(int*)calloc(n,sizeof(int)); int *after=(int*)calloc(n,sizeof(int)); for(int i=0;i<n;i++){ int rs,re, rfs,rfe; get_read_position(reads[i], 1, &rs, &re, &rfs, &rfe); start[i]=rs; stop[i]=re; }
BamRow **nsupv=(BamRow**)malloc(sizeof(BamRow*)*g->nsup); int nns=0;
for(int i=0;i<n;i++){
    int s=start[i];
    int e=stop[i];
    before[i]=0;
    after[i]=0;
    for(int j=0;j<n;j++){
        if(i == j) continue;
        int sj = start[j];
        int ej = stop[j];
        if(sj >= e){ after[i] = 1; }
        if(ej <= s){ before[i] = 1; continue; }
        if(sj <= s){ before[i] = 1; s = ej; }
        if(ej >= e){ after[i] = 1; e = sj; }
        if(sj > s && ej < e){ if((e - ej) > (sj - s)) s = ej; else e = sj; }
    }
    if(s > e) s = e;
    int keep_len = e - s;
    BamRow *r=reads[i];
    if(keep_len < min_align_length){
        pruned=1;
        if(r==g->primary){
            BamRow u=unmap_record(r,1,1);
            bam_replace_inplace(r, u);
            char key[128]; snprintf(key,sizeof(key), "%s_primary_unmapped_in_pruning", prefix);
            smap_inc(sm, key, 1);
        } else {
            char key2[128]; snprintf(key2,sizeof(key2), "%s_sup_discarded_in_pruning", prefix);
            smap_inc(sm, key2, 1);
        }
        continue;
    }
    if(start[i] != s || (s > 0 && before[i])){ pruned=1; BamRow nx=clip_start(r, s, 1); bam_replace_inplace(r, nx); char kk[128]; snprintf(kk,sizeof(kk), "%s_start_pruned", prefix); smap_inc(sm, kk, 1); }
    if(stop[i] != e  || (e > 0 && after[i])){ pruned=1; BamRow nx2=clip_end(r, e, 1); bam_replace_inplace(r, nx2); char kk2[128]; snprintf(kk2,sizeof(kk2), "%s_end_pruned", prefix); smap_inc(sm, kk2, 1); }
    if(r!=g->primary){ nsupv[nns++]=r; }
}
if(nns==0){ char k3[128]; snprintf(k3,sizeof(k3), "%s_all_sup_discarded_in_pruning", prefix); smap_inc(sm, k3, 1); free(g->supplementary); g->supplementary=NULL; g->nsup=0; } else { free(g->supplementary); g->supplementary=(BamRow**)malloc(sizeof(BamRow*)*nns); memcpy(g->supplementary, nsupv, sizeof(BamRow*)*nns); g->nsup=nns; }
free(nsupv); free(reads); free(start); free(stop); free(before); free(after); return pruned; }

static void process_fragment(FILE *out, StatsMap *sm, BamRow *group, int ng, int loose_ends, int min_align_length, int max_read_dist, int disable_validation){ smap_inc(sm, "alignment_counter", ng); smap_inc(sm, "fragment_counter", 1);
    // split by read1/read2
    int n1=0,n2=0; for(int i=0;i<ng;i++){ if(group[i].flag & 0x40) n1++; else if(group[i].flag & 0x80) n2++; }
    BamRow **r1arr=(BamRow**)malloc(sizeof(BamRow*)*n1); BamRow **r2arr=(BamRow**)malloc(sizeof(BamRow*)*n2); int i1=0,i2=0; for(int i=0;i<ng;i++){ if(group[i].flag & 0x40) r1arr[i1++]=&group[i]; else if(group[i].flag & 0x80) r2arr[i2++]=&group[i]; }
    if(n1==0 || n2==0){
        fprintf(stderr, "Error: Reads are not paired for qname %s (group size=%d)\n", ng>0? group[0].qname: "", ng);
        fflush(stderr);
        exit(1);
    }
    RG r1, r2; process_readgroup_assign(r1arr,n1,&r1); process_readgroup_assign(r2arr,n2,&r2);
    int mod1=0, mod2=0;
    if(r1.nsup>0){ mod1 |= rg_prune(&r1, sm, "read1", min_align_length); }
    if(r2.nsup>0){ mod2 |= rg_prune(&r2, sm, "read2", min_align_length); }
    int d1=0,d2=0; dechimer(&r1,&r2, sm, loose_ends, max_read_dist, &d1,&d2); mod1 |= d1; mod2 |= d2;
    if(mod1){ rg_filter(&r1, sm, "read1", min_align_length); update_sa_tag(&r1); }
    if(mod2){ rg_filter(&r2, sm, "read2", min_align_length); update_sa_tag(&r2); }
    if(mod1 || mod2){ smap_inc(sm, "fragment_modified", 1); }
    if(mod1 || mod2){
        update_mate_flags(&r1,&r2);
        update_mate_tags(&r1,&r2);
        update_pair_flag(&r1,&r2, sm, max_read_dist);
    }
    if((mod1 || mod2) && !disable_validation){ /* validation skipped in C port for now */ }
    // build output: primary + sup + secondary + supp_secondary for read1 then read2
    BamRow *olist[4096]; int no=0;
    if(r1.primary) {
        olist[no++]=r1.primary;
    }
    for(int i=0;i<r1.nsup;i++) {
        olist[no++]=r1.supplementary[i];
    }
    for(int i=0;i<r1.nsec;i++) {
        olist[no++]=r1.secondary[i];
    }
    for(int i=0;i<r1.nss;i++) {
        olist[no++]=r1.supplementary_secondary[i];
    }
    if(r2.primary) {
        olist[no++]=r2.primary;
    }
    for(int i=0;i<r2.nsup;i++) {
        olist[no++]=r2.supplementary[i];
    }
    for(int i=0;i<r2.nsec;i++) {
        olist[no++]=r2.secondary[i];
    }
    for(int i=0;i<r2.nss;i++) {
        olist[no++]=r2.supplementary_secondary[i];
    }
    for(int i=0;i<no;i++){ if(olist[i] && olist[i]->seq && olist[i]->seq[0]!='\0') print_bam(out, olist[i]); }
    rg_free(&r1); rg_free(&r2); free(r1arr); free(r2arr);
}

static void usage(const char *p){
    fprintf(stderr, "Usage: %s -i <in.sam|bam|-] -o <out.sam|-> -s stats.tsv [--loose_ends] [--min_align_length N] [--max_read_dist N] [-d]\n", p);
}

int main(int argc, char **argv){
    const char *in="-"; const char *outp="-"; const char *stats=NULL; int disable_validation=0;
    int loose_ends = 0; int min_align_length = 30; int max_read_dist = 100000;
    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"-i") && i+1<argc) in=argv[++i];
        else if(!strcmp(argv[i],"-o") && i+1<argc) outp=argv[++i];
        else if(!strcmp(argv[i],"-s") && i+1<argc) stats=argv[++i];
        else if(!strcmp(argv[i],"-d")) disable_validation=1;
        else if(!strcmp(argv[i],"--loose_ends")) loose_ends=1;
        else if(!strcmp(argv[i],"--min_align_length") && i+1<argc) min_align_length=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--max_read_dist") && i+1<argc) max_read_dist=atoi(argv[++i]);
        else { usage(argv[0]); return 1; }
    }
    if(!stats){ usage(argv[0]); fprintf(stderr, "Missing -s stats file\n"); return 1; }
    (void)disable_validation; (void)loose_ends; (void)min_align_length; (void)max_read_dist;
    
    /* Open input */
    FILE *inF = NULL; FILE *outF = NULL; FILE *pipe_proc = NULL;
    if(strcmp(in,"-")){
        /* stream BAM via samtools view -h */
        size_t n = strlen(in) + 32; char *cmd=(char*)malloc(n); snprintf(cmd,n,"samtools view -h '%s'", in);
        pipe_proc = popen(cmd, "r"); free(cmd);
        if(!pipe_proc){ perror("popen samtools"); return 1; }
        inF = pipe_proc;
    } else {
        inF = stdin;
    }
    if(strcmp(outp,"-")){
        outF = fopen(outp, "w"); if(!outF){ perror("open out"); if(pipe_proc) pclose(pipe_proc); return 1; }
    } else {
        outF = stdout;
    }
    
    StatsMap sm; smap_init(&sm);
    
    char *line=NULL; size_t cap=0; ssize_t rl;
    char *cur_q=NULL;
    BamRow *group_rows=NULL; int ng=0, ng_cap=0;
    
    long long alignment_counter=0; // for progress only
    
    double last_time = (double)clock() / CLOCKS_PER_SEC;
    while((rl = safe_getline(&line, &cap, inF)) > 0){
        if(line[0]=='@'){
            /* header lines pass-through */
            fputs(line, outF); fputc('\n', outF);
            continue;
        }
        alignment_counter++;
        /* parse first column (qname); we need it to group */
        char *tmp=sdup(line); char *f[64]; int nf=0; split_tab(tmp,f,&nf);
        if(nf<1){ free(tmp); continue; }
        /* strip trailing /1 or /2 for grouping */
        size_t L=strlen(f[0]); if(L>2 && f[0][L-2]=='/' && (f[0][L-1]=='1'||f[0][L-1]=='2')) f[0][L-2]='\0';
        if(!cur_q){ cur_q=sdup(f[0]); }
        if(strcmp(cur_q,f[0])!=0){
            /* process previous group */
            if(ng>0){
                process_fragment(outF, &sm, group_rows, ng, loose_ends, min_align_length, max_read_dist, disable_validation);
                for(int i=0;i<ng;i++){ bam_free(&group_rows[i]); }
                free(group_rows); group_rows=NULL; ng=0; ng_cap=0;
            }
            free(cur_q); cur_q=sdup(f[0]);
        }
        free(tmp);
        if(ng==ng_cap){ ng_cap = ng_cap? ng_cap*2: 16; group_rows=(BamRow*)realloc(group_rows, sizeof(BamRow)*ng_cap); }
        /* parse current line into BamRow */
        char *pline=sdup(line);
        if(!parse_bam_line_full(pline, &group_rows[ng])){ free(pline); continue; }
        free(pline);
        ng++;
        if((alignment_counter % 100000)==0){
            double now = (double)clock() / CLOCKS_PER_SEC;
            double dt = now - last_time;
            long long rate = dt > 0 ? (long long)(100000.0 / dt) : -1;
            fprintf(stderr, "alignment_counter=%lld reads_per_sec=%lld\n", alignment_counter, rate);
            last_time = now;
        }
    }
    if(ng>0){
        process_fragment(outF, &sm, group_rows, ng, loose_ends, min_align_length, max_read_dist, disable_validation);
        for(int i=0;i<ng;i++){ bam_free(&group_rows[i]); }
        free(group_rows);
    }
    
    if(outF && outF!=stdout) fclose(outF);
    if(pipe_proc) pclose(pipe_proc);
    
    fprintf(stderr, "Writing statistics\n");
    FILE *sf=fopen(stats, "w"); if(!sf){ perror("open stats"); smap_free(&sm); return 1; }
    for(int i=0;i<sm.n;i++) fprintf(sf, "%s\t %lld\n", sm.keys[i], sm.vals[i]);
    fflush(sf); fclose(sf);
    fprintf(stderr, "Done\n");
    
    smap_free(&sm);
    free(line); free(cur_q);
    return 0;
}