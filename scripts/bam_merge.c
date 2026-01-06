#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <zlib.h>
#include <stdint.h>
#include <ctype.h>
#include <stdarg.h>
#include <sys/stat.h>

 

static ssize_t safe_getline(char **lineptr, size_t *n, FILE *stream){
    ssize_t r = getline(lineptr, n, stream);
    if(r<=0) return r;
    if(r>0 && (*lineptr)[r-1]=='\n'){ (*lineptr)[r-1]='\0'; --r; }
    if(r>0 && (*lineptr)[r-1]=='\r'){ (*lineptr)[r-1]='\0'; --r; }
    return r;
}

static inline char qual_norm(char c){ return (c=='#') ? '!' : c; }
static inline char rc_base(char c){
    switch(c){
        case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
        case 'a': return 't'; case 'c': return 'g'; case 'g': return 'c'; case 't': return 'a';
        default: return c;
    }
}
static void rev_inplace(char *s, size_t n){ size_t i=0, j = n? n-1:0; while(i<j){ char t=s[i]; s[i]=s[j]; s[j]=t; i++; j--; } }
static void revcomp_inplace(char *s, size_t n){ size_t i=0, j = n? n-1:0; while(i<j){ char ci=rc_base(s[i]); char cj=rc_base(s[j]); s[i]=cj; s[j]=ci; i++; j--; } if(n%2==1) s[i]=rc_base(s[i]); }

static void usage(const char *p){
    fprintf(stderr,"Usage: %s -a fq1.gz -b fq2.gz -ua bad1.gz -ub bad2.gz -s stats.tsv [-o out|-]\n",p);
}

typedef struct {
    char *raw;
    char *qname;
    int flag;
    char *rname;
    char *cigar;
    char *seq;
    char *qual;
    int read1, read2, unmapped, reversed, secondary, supplementary;
    int ascore;
    char *rg;
} BamRow;

static char* sdup(const char *s){ size_t n=strlen(s); char *p=(char*)malloc(n+1); memcpy(p,s,n+1); return p; }

static void split_tab(char *s, char **f, int *nf){
    int i=0; char *p=s; f[i++]=p;
    while(*p){ if(*p=='\t'){ *p='\0'; f[i++]=p+1; } p++; if(i>=64) break; }
    *nf=i;
}

static BamRow parse_bam_line(char *line){
    BamRow r; memset(&r,0,sizeof(r));
    r.raw = sdup(line);
    char *f[64]; int nf=0; split_tab(line,f,&nf);
    if(nf<11){ free(line); return r; }
    r.qname=sdup(f[0]);
    {
        size_t L = strlen(r.qname);
        if(L>2 && r.qname[L-2]=='/' && (r.qname[L-1]=='1' || r.qname[L-1]=='2')){
            r.qname[L-2] = '\0';
        }
    }
    r.flag=(int)strtol(f[1],NULL,10); r.rname=sdup(f[2]); r.cigar=sdup(f[5]);
    r.seq=sdup(f[9]); r.qual=sdup(f[10]);
    r.read1=(r.flag&0x40)!=0; r.read2=(r.flag&0x80)!=0; r.unmapped=(r.flag&0x4)!=0; r.reversed=(r.flag&0x10)!=0;
    r.secondary=(r.flag&0x100)!=0; r.supplementary=(r.flag&0x800)!=0; r.ascore=0; r.rg=NULL;
    for(int i=11;i<nf;i++){
        if(!strncmp(f[i],"AS:i:",5)) r.ascore=(int)strtol(f[i]+5,NULL,10);
        else if(!strncmp(f[i],"RG:Z:",5)) r.rg=sdup(f[i]+5);
    }
    free(line);
    return r;
}

static void free_bamrow(BamRow *r){
    if(!r) return;
    free(r->raw); free(r->qname); free(r->rname); free(r->cigar);
    free(r->seq); free(r->qual); if(r->rg) free(r->rg);
    memset(r,0,sizeof(*r));
}

// --- CIGAR helpers and RG sanitization ---
static int cigar_soft_bp(const char *c){ if(!c||c[0]=='*') return 0; int n=0,v=0; for(const char *p=c; *p; ++p){ if(isdigit((unsigned char)*p)) v=v*10+(*p-'0'); else { if(*p=='S') n+=v; v=0; } } return n; }
static int cigar_aligned_read_bp(const char *c, int seq_len){ if(!c||c[0]=='*'){ return 0; } int v=0, n=0; int hasSH=0; for(const char *p=c; *p; ++p){ if(isdigit((unsigned char)*p)) v=v*10+(*p-'0'); else { if(*p=='S' || *p=='H') { hasSH=1; } if(*p=='M'||*p=='I'||*p=='='||*p=='X') n+=v; v=0; } } if(!hasSH) return seq_len; return n; }
static void sanitize_rg(char *s){ if(!s) return; for(char *p=s; *p; ++p){ if(*p==' '||*p=='@'||*p==':') *p='_'; } }

static void gz_write_fastq(gzFile gz, const char *qname, const char *rg, int readnr, int reversed, const char *seq, const char *qual){
    if(!qname||!seq||!qual) return;
    if(rg && *rg){
        gzprintf(gz, "@%s:%s/%d\n", qname, rg, readnr);
    } else {
        gzprintf(gz, "@%s/%d\n", qname, readnr);
    }
    size_t L = strlen(seq);
    char *s = sdup(seq);
    char *q = sdup(qual);
    if(reversed){ revcomp_inplace(s, L); rev_inplace(q, L); }
    gzputs(gz, s); gzputs(gz, "\n+\n"); gzputs(gz, q); gzputs(gz, "\n");
    free(s); free(q);
}

// --- FastqPair and simple hash ---
typedef struct { char *name,*seq1,*seq2,*qual1,*qual2; } FastqPair;
typedef struct Entry { char *key; FastqPair *val; struct Entry *next; } Entry;
typedef struct { Entry **b; size_t nb; size_t n; } Hash;

static uint64_t fnv1a_k(const char *s){ uint64_t h=0xcbf29ce484222325ULL; const unsigned char *p=(const unsigned char*)s; while(*p){ h^=*p++; h*=0x100000001b3ULL; } return h; }
static void hash_init(Hash *h, size_t nb){ h->nb= nb?nb:4096; h->n=0; h->b=(Entry**)calloc(h->nb,sizeof(Entry*)); }
static void hash_rehash(Hash *h, size_t new_nb){
    if(!h||!h->b||new_nb<=h->nb) return;
    Entry **nb = (Entry**)calloc(new_nb, sizeof(Entry*));
    for(size_t i=0;i<h->nb;i++){
        Entry *e = h->b[i];
        while(e){
            Entry *nx = e->next;
            uint64_t hv = fnv1a_k(e->key);
            size_t j = hv % new_nb;
            e->next = nb[j];
            nb[j] = e;
            e = nx;
        }
    }
    free(h->b);
    h->b = nb;
    h->nb = new_nb;
}
static inline void hash_maybe_grow(Hash *h){ if(h->n > (h->nb*3)/4) hash_rehash(h, h->nb*2); }
static void fastqpair_free(FastqPair *v){ if(!v) return; free(v->name); free(v->seq1); free(v->seq2); free(v->qual1); free(v->qual2); free(v); }
static void hash_free(Hash *h){ if(!h||!h->b) return; for(size_t i=0;i<h->nb;i++){ Entry *e=h->b[i]; while(e){ Entry*nx=e->next; free(e->key); fastqpair_free(e->val); free(e); e=nx; } } free(h->b); h->b=NULL; h->nb=0; h->n=0; }
static void hash_put(Hash *h, const char *key, FastqPair *v){ uint64_t hv=fnv1a_k(key); size_t i=hv%h->nb; Entry *e=(Entry*)malloc(sizeof(Entry)); e->key=sdup(key); e->val=v; e->next=h->b[i]; h->b[i]=e; h->n++; hash_maybe_grow(h); }
 
static FastqPair* hash_pop(Hash *h, const char *key){ if(!h||!h->b) return NULL; uint64_t hv=fnv1a_k(key); size_t i=hv%h->nb; Entry *prev=NULL,*e=h->b[i]; while(e){ if(strcmp(e->key,key)==0){ if(prev) prev->next=e->next; else h->b[i]=e->next; FastqPair* v=e->val; free(e->key); free(e); h->n--; return v; } prev=e; e=e->next; } return NULL; }

// Helpers for assembling pairs by name
static Entry* hash_find_entry(Hash *h, const char *key){ if(!h||!h->b) return NULL; uint64_t hv=fnv1a_k(key); size_t i=hv%h->nb; Entry *e=h->b[i]; while(e){ if(strcmp(e->key,key)==0) return e; e=e->next; } return NULL; }
static FastqPair* hash_get_or_insert(Hash *h, const char *key){ Entry *e = hash_find_entry(h,key); if(e) return e->val; FastqPair *v=(FastqPair*)calloc(1, sizeof(FastqPair)); v->name=sdup(key); hash_put(h, key, v); return v; }
static int fqpair_is_complete(FastqPair *v){ return v && v->seq1 && v->qual1 && v->seq2 && v->qual2; }
static FastqPair* hash_pop_any_complete(Hash *h){ if(!h||!h->b) return NULL; for(size_t i=0;i<h->nb;i++){ Entry *prev=NULL,*e=h->b[i]; while(e){ if(fqpair_is_complete(e->val)){ if(prev) prev->next=e->next; else h->b[i]=e->next; FastqPair *v=e->val; free(e->key); free(e); h->n--; return v; } prev=e; e=e->next; } } return NULL; }

/* moved fq_retrieve_any below FastqCtx definition */

// Print two unmapped SAM records from a FastqPair, with optional RG
static void print_unmapped_sam_pair(FILE *out, const FastqPair *fq, const char *rg){
    int flag1 = 0x1 | 0x4 | 0x8 | 0x40 | 0x200;
    int flag2 = 0x1 | 0x4 | 0x8 | 0x80 | 0x200;
    if(rg && *rg){
        fprintf(out, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tRG:Z:%s\n", fq->name, flag1, fq->seq1, fq->qual1, rg);
        fprintf(out, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tRG:Z:%s\n", fq->name, flag2, fq->seq2, fq->qual2, rg);
    } else {
        fprintf(out, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", fq->name, flag1, fq->seq1, fq->qual1);
        fprintf(out, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", fq->name, flag2, fq->seq2, fq->qual2);
    }
}

// --- FASTQ reader via pigz -dc ---
typedef struct { FILE *fa; FILE *fb; Hash buf; size_t max_buf; int eof_a; int eof_b; int qual_shift; } FastqCtx;

static FILE* pigz_open(const char *path){
    char *cmd=NULL; size_t n = strlen(path)+20; cmd = (char*)malloc(n);
    snprintf(cmd, n, "pigz -dc '%s'", path);
    FILE *fp = popen(cmd, "r");
    free(cmd);
    return fp;
}

static void fqctx_init(FastqCtx *fc, const char *a, const char *b){
    fc->fa = pigz_open(a);
    fc->fb = pigz_open(b);
    hash_init(&fc->buf, 8192);
    fc->max_buf = 200000; fc->eof_a=0; fc->eof_b=0;
    fc->qual_shift = 0;
    if(fc->fa) setvbuf(fc->fa, NULL, _IOFBF, 1<<20);
    if(fc->fb) setvbuf(fc->fb, NULL, _IOFBF, 1<<20);
}
static void fqctx_free(FastqCtx *fc){ if(fc->fa) pclose(fc->fa); if(fc->fb) pclose(fc->fb); hash_free(&fc->buf); }

static char* trim_newline(char *s){ if(!s) return s; size_t l=strlen(s); while(l>0 && (s[l-1]=='\n'||s[l-1]=='\r')) s[--l]='\0'; return s; }
static char* extract_qname(const char *hdr){
    // hdr begins with '@'
    const char *p = hdr[0]=='@'? hdr+1: hdr;
    size_t i=0; while(p[i] && p[i] != ' ' && p[i] != '\t' && p[i] != '\n' && p[i] != '\r') i++;
    char *name = (char*)malloc(i+1); memcpy(name,p,i); name[i]='\0';
    size_t L=strlen(name);
    if(L>2 && name[L-2]=='/' && (name[L-1]=='1'||name[L-1]=='2')) name[L-2]='\0';
    return name;
}

static int read_fastq(FILE *f, char **name, char **seq, char **qual){
    static char *h=NULL,*s=NULL,*p=NULL,*q=NULL; static size_t ch=0,cs=0,cp=0,cq=0; ssize_t r;
    if(!f) return 0;
    r=safe_getline(&h,&ch,f); if(r<=0) return 0; trim_newline(h);
    r=safe_getline(&s,&cs,f); if(r<=0) return 0; trim_newline(s);
    r=safe_getline(&p,&cp,f); if(r<=0) return 0;
    r=safe_getline(&q,&cq,f); if(r<=0) return 0; trim_newline(q);
    *name = extract_qname(h);
    *seq = s;
    *qual = q;
    return 1;
}

static void qual_shift_inplace(char *q, int shift){
    if(!q || shift==0) return;
    for(char *p=q; *p; ++p){ *p = (char)((int)(unsigned char)(*p) + shift); }
}

static void fqctx_apply_shift(FastqCtx *fc, int shift){
    if(!fc) return;
    if(fc->qual_shift == shift) return;
    fc->qual_shift = shift;
    if(!fc->buf.b) return;
    for(size_t i=0;i<fc->buf.nb;i++){
        Entry *e = fc->buf.b[i];
        while(e){
            FastqPair *v = e->val;
            if(v){
                if(v->qual1) qual_shift_inplace(v->qual1, shift);
                if(v->qual2) qual_shift_inplace(v->qual2, shift);
            }
            e = e->next;
        }
    }
}

static void fq_store_pair(FastqCtx *fc, const char *name, const char *s1, const char *s2, const char *q1, const char *q2){
    FastqPair *v = hash_get_or_insert(&fc->buf, name);
    if(s1 && q1){ if(v->seq1) free(v->seq1); if(v->qual1) free(v->qual1); v->seq1=sdup(s1); v->qual1=sdup(q1); if(fc && fc->qual_shift) qual_shift_inplace(v->qual1, fc->qual_shift); }
    if(s2 && q2){ if(v->seq2) free(v->seq2); if(v->qual2) free(v->qual2); v->seq2=sdup(s2); v->qual2=sdup(q2); if(fc && fc->qual_shift) qual_shift_inplace(v->qual2, fc->qual_shift); }
}

static int fq_pop_by_name(FastqCtx *fc, const char *target, FastqPair **out){
    Entry *e = hash_find_entry(&fc->buf, target);
    if(e && fqpair_is_complete(e->val)){
        *out = hash_pop(&fc->buf, target);
        return 1;
    }
    while(1){
        char *n1=NULL,*s1=NULL,*q1=NULL; char *n2=NULL,*s2=NULL,*q2=NULL;
        int ok1 = read_fastq(fc->fa,&n1,&s1,&q1);
        if(!ok1){ free(n1); *out=NULL; return 0; }
        fq_store_pair(fc, n1, s1, NULL, q1, NULL);
        free(n1);
        int ok2 = read_fastq(fc->fb,&n2,&s2,&q2);
        if(!ok2){ free(n2); *out=NULL; return 0; }
        fq_store_pair(fc, n2, NULL, s2, NULL, q2);
        free(n2);
        e = hash_find_entry(&fc->buf, target);
        if(e && fqpair_is_complete(e->val)){
            *out = hash_pop(&fc->buf, target);
            return 1;
        }
        if(fc->buf.n > fc->max_buf){ /* optional: could warn about memory pressure */ }
    }
}

// Retrieve any available complete FASTQ pair from the buffers, reading forward if needed.
// Returns 1 and sets *out when a pair is available; returns 0 when no more pairs can be read.
static int fq_retrieve_any(FastqCtx *fc, FastqPair **out){
    if(!fc || !out) return 0;
    // First try to pop any complete pair already buffered
    FastqPair *v = hash_pop_any_complete(&fc->buf);
    if(v){ *out = v; return 1; }
    // Otherwise, keep reading until we form at least one complete pair or hit EOF
    while(1){
        char *n1=NULL,*s1=NULL,*q1=NULL; char *n2=NULL,*s2=NULL,*q2=NULL;
        int ok1 = read_fastq(fc->fa,&n1,&s1,&q1);
        if(!ok1){ fc->eof_a=1; free(n1); break; }
        fq_store_pair(fc, n1, s1, NULL, q1, NULL);
        free(n1);
        int ok2 = read_fastq(fc->fb,&n2,&s2,&q2);
        if(!ok2){ fc->eof_b=1; free(n2); break; }
        fq_store_pair(fc, n2, NULL, s2, NULL, q2);
        free(n2);
        v = hash_pop_any_complete(&fc->buf);
        if(v){ *out = v; return 1; }
    }
    // Drain any remaining complete pairs after EOF
    v = hash_pop_any_complete(&fc->buf);
    if(v){ *out = v; return 1; }
    *out = NULL;
    return 0;
}

typedef struct { char t; int l; } CigTok;
static int cigar_parse(const char *c, CigTok *v, int maxv){ int n=0, val=0; if(!c||c[0]=='*') return 0; for(const char *p=c; *p; ++p){ if(isdigit((unsigned char)*p)) val = val*10 + (*p-'0'); else { if(n<maxv){ v[n].t=*p; v[n].l=val; n++; } val=0; } } return n; }
static char* cigar_join(const CigTok *v, int n){ size_t sz=1; for(int i=0;i<n;i++){ int tmp=v[i].l; int digits=1; while(tmp>=10){digits++; tmp/=10;} sz += digits+1; } char *s=(char*)malloc(sz); size_t k=0; for(int i=0;i<n;i++){ k += sprintf(s+k, "%d%c", v[i].l, v[i].t); } s[k]='\0'; return s; }
static char* cigar_add_hardclip(const char *c, int add_len, int at_front, int reversed){
    (void)reversed; // parameter currently not used
    if(add_len<=0) return sdup(c?c:"*");
    if(!c||c[0]=='*'){ return sdup("*"); }
    CigTok a[128]; int n=cigar_parse(c,a,128); if(n<=0) return sdup(c);
    // Work in original orientation as Python does for supplementary handling; for primaries we only need front/back
    if(at_front){ // remove leading H
        int i0=0; while(i0<n && a[i0].t=='H'){ i0++; }
        int m=n-i0; for(int i=0;i<m;i++){ a[i]=a[i+i0]; }
        // prepend H
        for(int i=m;i>0;i--){ a[i]=a[i-1]; }
        a[0].t='H'; a[0].l=add_len; n=m+1;
    } else {
        // remove trailing H
        int m=n; while(m>0 && a[m-1].t=='H'){ m--; }
        a[m].t='H'; a[m].l=add_len; n=m+1;
    }
    return cigar_join(a,n);
}

// Extend existing leading/trailing H clip by add_len
static char* cigar_extend_hardclip(const char *c, int add_len, int at_front){
    if(add_len<=0) return sdup(c?c:"*");
    if(!c || c[0]=='*') return sdup("*");
    CigTok a[256]; int n=cigar_parse(c,a,256); if(n<=0) return sdup(c);
    if(at_front){
        int i0=0; int clip=0; while(i0<n && a[i0].t=='H'){ clip += a[i0].l; i0++; }
        int m=n-i0; for(int i=0;i<m;i++){ a[i]=a[i+i0]; }
        for(int i=m;i>0;i--){ a[i]=a[i-1]; }
        a[0].t='H'; a[0].l=clip + add_len; n=m+1;
    } else {
        int m=n; int clip=0; while(m>0 && a[m-1].t=='H'){ clip += a[m-1].l; m--; }
        a[m].t='H'; a[m].l=clip + add_len; n=m+1;
    }
    return cigar_join(a,n);
}

// --- Qual check allowing # vs ! ---
static int qual_equal_norm(const char *a, const char *b, int len){ for(int i=0;i<len;i++){ char ca= a[i]=='#'?'!':a[i]; char cb= b[i]=='#'?'!':b[i]; if(ca!=cb) return 0; } return 1; }

static int qual_equal_norm_shift(const char *a, const char *b, int len, int shift){
    for(int i=0;i<len;i++){
        char ca = a[i]=='#'?'!':a[i];
        char cb = b[i]=='#'?'!':b[i];
        ca = (char)((int)(unsigned char)ca + shift);
        if(ca!=cb) return 0;
    }
    return 1;
}

// --- Print SAM with optional new cigar and appended tags ---
static void sam_print_qname_only(FILE *out, const char *raw, const char *new_qname){
    const char *tab = strchr(raw, '\t');
    if(tab && new_qname){
        fputs(new_qname, out);
        fputc('\t', out);
        fputs(tab+1, out);
        fputc('\n', out);
    } else {
        fputs(raw, out);
        fputc('\n', out);
    }
}
static void sam_print_with_mods(FILE *out, const char *raw, const char *new_qname, const char *new_cigar,
                                const char **extra_tags, int ntags){
    char *tmp = sdup(raw);
    char *f[128]; int nf=0; split_tab(tmp,f,&nf);
    if(nf<11){ fprintf(out, "%s\n", raw); free(tmp); return; }
    // print field 0 (QNAME) possibly replaced, then fields 1..4
    fputs(new_qname? new_qname: f[0], out);
    for(int i=1;i<5;i++){ fputs("\t",out); fputs(f[i],out); }
    fputs("\t",out);
    fputs(new_cigar? new_cigar: f[5], out);
    for(int i=6;i<11;i++){ fputs("\t",out); fputs(f[i],out); }
    for(int i=11;i<nf;i++){ fputs("\t",out); fputs(f[i],out); }
    for(int i=0;i<ntags;i++){ fputs("\t",out); fputs(extra_tags[i],out); }
    fputc('\n',out);
    free(tmp);
}

// --- Derive missing sequence tags for primary ---
typedef struct {
    long long primary_reads;
    long long primary_aligned_bp;
    long long primary_soft_bp;
    long long restored_bp_r1;
    long long restored_bp_r2;
    long long restored_read1s;
    long long restored_read2s;
    long long supplementary_aligned_bp;
    long long supplementary_alignments;
    long long supplementary_alignments_cigar_adapted;
    long long fragments;
    long long alignments;
    long long total_bp;
    long long badly_mapped_fragments;
    long long fully_unmapped_fragments;
    long long chrEBV_mapped_fragments;
    long long readded_fragments;
    long long restored_unaligned_reads;
    int size_r1_cap, size_r2_cap; int *size_r1; int *size_r2;
} Stats;

static void stats_init(Stats *st){ memset(st,0,sizeof(*st)); st->size_r1_cap=4096; st->size_r2_cap=4096; st->size_r1=(int*)calloc(st->size_r1_cap,sizeof(int)); st->size_r2=(int*)calloc(st->size_r2_cap,sizeof(int)); }
static void stats_free(Stats *st){ free(st->size_r1); free(st->size_r2); }
static void stats_add_size(int **arr, int *cap, int len){ if(len<=0) return; if(len>=*cap){ int ncap=*cap; while(len>=ncap) ncap*=2; int *nx=(int*)realloc(*arr, ncap*sizeof(int)); memset(nx+*cap,0,(ncap-*cap)*sizeof(int)); *arr=nx; *cap=ncap; } (*arr)[len] += 1; }

static int derive_missing_sequence_tags_primary(const BamRow *br, const char *fqseq, const char *fqqual,
                                               char **new_cigar_out, const char **tags_out, int *ntags_out,
                                               Stats *st){
    // returns restored length
    *new_cigar_out=NULL; *ntags_out=0;
    if(!br || !fqseq || !fqqual) return 0;
    st->primary_reads += 1;
    st->primary_aligned_bp += cigar_aligned_read_bp(br->cigar, (int)strlen(br->seq));
    if(strchr(br->cigar,'S')) st->primary_soft_bp += cigar_soft_bp(br->cigar);

    int Lbam = (int)strlen(br->seq);
    int Lfq = (int)strlen(fqseq);
    int clip = Lfq - Lbam;
    if(clip<=0){ return 0; }
    const char *qname = br->qname;
    // Special handling for unmapped primaries: treat as forward orientation like Python implementation
    if(br->cigar && br->cigar[0]=='*'){
        // forward: check prefix equals BAM seq/qual
        if(Lbam>0){
            if(strncmp(fqseq, br->seq, Lbam)!=0){ fprintf(stderr,"Sequence mismatch BAM-FASTQ in fragment %s (unmapped)\n", qname); exit(1); }
            if(!qual_equal_norm(fqqual, br->qual, Lbam)){ fprintf(stderr,"Quality mismatch BAM-FASTQ in fragment %s (unmapped)\n", qname); exit(1); }
        }
        if(clip>0){
            // tag is suffix
            const char *tseq = fqseq + Lbam; const char *tqual = fqqual + Lbam;
            char *zbtag; char *zqtag; size_t tbl=strlen(tseq), tql=strlen(tqual);
            zbtag=(char*)malloc(5+tbl+1); sprintf(zbtag,"ZB:Z:%s", tseq);
            zqtag=(char*)malloc(5+tql+1); sprintf(zqtag,"ZQ:Z:%s", tqual);
            tags_out[(*ntags_out)++] = zbtag; tags_out[(*ntags_out)++] = zqtag;
            // keep CIGAR as '*'
            *new_cigar_out = sdup("*");
        }
    }
    else if(br->reversed){
        // reverse-complement fastq to align to BAM orientation
        char *rseq = sdup(fqseq); revcomp_inplace(rseq, (size_t)Lfq);
        char *rqual = sdup(fqqual); rev_inplace(rqual, (size_t)Lfq);
        // check suffix equals BAM seq/qual
        if(strncmp(rseq + (Lfq-Lbam), br->seq, Lbam)!=0){ fprintf(stderr,"Sequence mismatch BAM-FASTQ in fragment %s (rev)\n", qname); exit(1); }
        if(!qual_equal_norm(rqual + (Lfq-Lbam), br->qual, Lbam)){ fprintf(stderr,"Quality mismatch BAM-FASTQ in fragment %s (rev)\n", qname); exit(1); }
        // tag is reversed prefix
        rseq[Lfq-Lbam]='\0'; rqual[Lfq-Lbam]='\0';
        // build tags YB/YQ
        char *yb = sdup(rseq); yb[Lfq-Lbam]='\0';
        char *yq = sdup(rqual); yq[Lfq-Lbam]='\0';
        char *ybtag; char *yqtag; size_t ybl=strlen(yb), yql=strlen(yq);
        ybtag=(char*)malloc(5+ybl+1); sprintf(ybtag,"YB:Z:%s", yb);
        yqtag=(char*)malloc(5+yql+1); sprintf(yqtag,"YQ:Z:%s", yq);
        tags_out[(*ntags_out)++] = ybtag; tags_out[(*ntags_out)++] = yqtag;
        // adjust cigar: add H at front
        *new_cigar_out = cigar_add_hardclip(br->cigar, clip, 1, 0);
        free(yb); free(yq); free(rseq); free(rqual);
    } else {
        // forward: check prefix equals BAM seq/qual
        if(strncmp(fqseq, br->seq, Lbam)!=0){ fprintf(stderr,"Sequence mismatch BAM-FASTQ in fragment %s (fwd)\n", qname); exit(1); }
        if(!qual_equal_norm(fqqual, br->qual, Lbam)){ fprintf(stderr,"Quality mismatch BAM-FASTQ in fragment %s (fwd)\n", qname); exit(1); }
        // tag is suffix
        const char *tseq = fqseq + Lbam; const char *tqual = fqqual + Lbam;
        char *zbtag; char *zqtag; size_t tbl=strlen(tseq), tql=strlen(tqual);
        zbtag=(char*)malloc(5+tbl+1); sprintf(zbtag,"ZB:Z:%s", tseq);
        zqtag=(char*)malloc(5+tql+1); sprintf(zqtag,"ZQ:Z:%s", tqual);
        tags_out[(*ntags_out)++] = zbtag; tags_out[(*ntags_out)++] = zqtag;
        // adjust cigar: add H at end
        *new_cigar_out = cigar_add_hardclip(br->cigar, clip, 0, 0);
    }
    // XT tag only when clip>0
    if(clip>0){ char *xt=(char*)malloc(32); sprintf(xt,"XT:i:%d", clip); tags_out[(*ntags_out)++] = xt; }
    // stats
    if(br->read1){ st->restored_bp_r1 += clip; st->restored_read1s += 1; stats_add_size(&st->size_r1, &st->size_r1_cap, clip); }
    else if(br->read2){ st->restored_bp_r2 += clip; st->restored_read2s += 1; stats_add_size(&st->size_r2, &st->size_r2_cap, clip); }
    if(br->cigar[0]=='*') st->restored_unaligned_reads += 1;
    return clip;
}

// Flush a group of BAM rows that share the same qname
static void flush_group(FILE *out, gzFile g1, gzFile g2,
                        BamRow *group, int *ng_ptr, char **cur_q_ptr,
                        long long *primary_aligned_bp, long long *primary_soft_bp, long long *supplementary_aligned_bp,
                        int *wrote_badmap_ptr,
                        FastqCtx *fc, Stats *st, char **last_rg_ptr){
    int ng = *ng_ptr;
    if(ng==0) return;
    int as1=0, as2=0; int max1=0, max2=0; int ebv=0;
    BamRow *p1=NULL, *p2=NULL;
    // identify primaries and compute stats
    for(int i=0;i<ng;i++){
        BamRow *br=&group[i];
        if(!br->secondary && !br->supplementary){
            // primary (mapped or unmapped)
            if(br->read1){ p1=br; max1=(int)strlen(br->seq); }
            if(br->read2){ p2=br; max2=(int)strlen(br->seq); }
            if(!br->unmapped){
                *primary_aligned_bp += cigar_aligned_read_bp(br->cigar, (int)strlen(br->seq));
                *primary_soft_bp += cigar_soft_bp(br->cigar);
            }
        } else if(!br->unmapped && br->supplementary){
            *supplementary_aligned_bp += cigar_aligned_read_bp(br->cigar, (int)strlen(br->seq));
            if(st) st->supplementary_alignments += 1;
        }
        // sum AS for primary + supplementary only, exclude secondary
        if(!br->secondary && !br->unmapped){
            if(br->read1) as1 += br->ascore;
            if(br->read2) as2 += br->ascore;
        }
        if(!br->secondary && !br->unmapped && strcmp(br->rname,"chrEBV")==0) ebv=1;
    }
    // mirror bam_merge.py: both reads low absolute and relative AS; when max==0, relative check is false
    int badly_lowAS = 0;
    double thr1 = 0.5 * (double)max1;
    double thr2 = 0.5 * (double)max2;
    if( (as1<50 && as2<50) && ( (max1>0 ? (as1 < thr1) : 0) && (max2>0 ? (as2 < thr2) : 0) ) ) badly_lowAS=1;
    int badly = badly_lowAS || ebv;
    if(ebv){ if(st) st->chrEBV_mapped_fragments += 1; }

    // Retrieve corresponding FASTQ pair and compute fragment-level stats and restoration for primaries
    if(fc){
        FastqPair *fq=NULL; if(!fq_pop_by_name(fc, group[0].qname, &fq)){ fq=NULL; }
         if(fq){ st->fragments += 1; st->alignments += ng; st->total_bp += (long long)strlen(fq->seq1) + (long long)strlen(fq->seq2);
            if(fc->qual_shift == 0 && (p1 || p2)){
                int ok0 = 1;
                int ok64 = 1;
                int tested = 0;
                int L1 = p1 ? (int)strlen(p1->qual) : 0;
                int L2 = p2 ? (int)strlen(p2->qual) : 0;
                if(p1){
                    tested = 1;
                    if(p1->reversed){
                        for(int i=0;i<L1;i++){
                            char ca = fq->qual1[L1-1-i]=='#'?'!':fq->qual1[L1-1-i];
                            char cb = p1->qual[i]=='#'?'!':p1->qual[i];
                            if(ca!=cb){ ok0=0; break; }
                        }
                        if(ok64){
                            for(int i=0;i<L1;i++){
                                char ca = fq->qual1[L1-1-i]=='#'?'!':fq->qual1[L1-1-i];
                                char cb = p1->qual[i]=='#'?'!':p1->qual[i];
                                ca = (char)((int)(unsigned char)ca - 31);
                                if(ca!=cb){ ok64=0; break; }
                            }
                        }
                    } else {
                        if(!qual_equal_norm(fq->qual1, p1->qual, L1)) ok0=0;
                        if(ok64 && !qual_equal_norm_shift(fq->qual1, p1->qual, L1, -31)) ok64=0;
                    }
                }
                if(p2 && (ok0 || ok64)){
                    tested = 1;
                    if(p2->reversed){
                        for(int i=0;i<L2;i++){
                            char ca = fq->qual2[L2-1-i]=='#'?'!':fq->qual2[L2-1-i];
                            char cb = p2->qual[i]=='#'?'!':p2->qual[i];
                            if(ca!=cb){ ok0=0; break; }
                        }
                        if(ok64){
                            for(int i=0;i<L2;i++){
                                char ca = fq->qual2[L2-1-i]=='#'?'!':fq->qual2[L2-1-i];
                                char cb = p2->qual[i]=='#'?'!':p2->qual[i];
                                ca = (char)((int)(unsigned char)ca - 31);
                                if(ca!=cb){ ok64=0; break; }
                            }
                        }
                    } else {
                        if(!qual_equal_norm(fq->qual2, p2->qual, L2)) ok0=0;
                        if(ok64 && !qual_equal_norm_shift(fq->qual2, p2->qual, L2, -31)) ok64=0;
                    }
                }
                if(tested && !ok0 && ok64){
                    fqctx_apply_shift(fc, -31);
                    qual_shift_inplace(fq->qual1, -31);
                    qual_shift_inplace(fq->qual2, -31);
                }
            }
            // precompute primary modifications
            const char *tags1[8]; int nt1=0; char *nc1=NULL; int clip1=0;
            const char *tags2[8]; int nt2=0; char *nc2=NULL; int clip2=0;
            if(p1){ clip1 = derive_missing_sequence_tags_primary(p1, fq->seq1, fq->qual1, &nc1, tags1, &nt1, st); }
            if(p2){ clip2 = derive_missing_sequence_tags_primary(p2, fq->seq2, fq->qual2, &nc2, tags2, &nt2, st); }
            // track last RG from primaries
            if(last_rg_ptr){ const char *rgp = (p1 && p1->rg)? p1->rg : ((p2 && p2->rg)? p2->rg : NULL); if(rgp){ if(*last_rg_ptr) free(*last_rg_ptr); *last_rg_ptr = sdup(rgp); } }
            // print in original order
            for(int i=0;i<ng;i++){
                BamRow *br=&group[i];
                if(br==p1){ if(nt1>0 || nc1){ sam_print_with_mods(out, p1->raw, p1->qname, nc1?nc1:NULL, tags1, nt1); } else { sam_print_qname_only(out, p1->raw, p1->qname); } }
                else if(br==p2){ if(nt2>0 || nc2){ sam_print_with_mods(out, p2->raw, p2->qname, nc2?nc2:NULL, tags2, nt2); } else { sam_print_qname_only(out, p2->raw, p2->qname); } }
                else {
                    // adapt supplementary
                    if(br->supplementary && !br->unmapped){
                        int add = br->read1? clip1 : (br->read2? clip2 : 0);
                        if(add>0){
                            char *nc_sup = cigar_extend_hardclip(br->cigar, add, br->reversed?1:0);
                            char *xt = (char*)malloc(32); sprintf(xt,"XT:i:%d", add);
                            const char *tagsup[1]; tagsup[0]=xt;
                            sam_print_with_mods(out, br->raw, br->qname, nc_sup, tagsup, 1);
                            free(xt); free(nc_sup);
                            if(st) st->supplementary_alignments_cigar_adapted += 1;
                        } else {
                            sam_print_qname_only(out, br->raw, br->qname);
                        }
                    } else {
                        sam_print_qname_only(out, br->raw, br->qname);
                    }
                }
            }
            // write badmap FASTQ if needed using SAM sequences (unrestored), oriented to original by reversed flag
            if(badly && p1 && p2){
                char *rg = NULL; if(p1->rg) rg=sdup(p1->rg); if(rg) sanitize_rg(rg);
                gz_write_fastq(g1, p1->qname, rg, 1, p1->reversed, p1->seq, p1->qual);
                gz_write_fastq(g2, p2->qname, rg, 2, p2->reversed, p2->seq, p2->qual);
                if(rg) free(rg);
                if(wrote_badmap_ptr) *wrote_badmap_ptr += 1;
                if(st){
                    if((as1+as2)==0) st->fully_unmapped_fragments += 1;
                    else if(badly_lowAS) st->badly_mapped_fragments += 1;
                }
            }
            for(int k=0;k<nt1;k++){ free((void*)tags1[k]); }
            if(nc1) free(nc1);
            for(int k=0;k<nt2;k++){ free((void*)tags2[k]); }
            if(nc2) free(nc2);
            fastqpair_free(fq);
        } else {
            for(int i=0;i<ng;i++) sam_print_qname_only(out, group[i].raw, group[i].qname);
        }
    } else {
        for(int i=0;i<ng;i++) sam_print_qname_only(out, group[i].raw, group[i].qname);
    }


    // free group and reset
    for(int i=0;i<ng;i++) free_bamrow(&group[i]);
    *ng_ptr = 0;
    if(*cur_q_ptr){ free(*cur_q_ptr); *cur_q_ptr=NULL; }
}

int main(int argc,char**argv){
    const char *fq1=NULL,*fq2=NULL,*ua=NULL,*ub=NULL,*stats=NULL,*outp="-";
    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"-a")&&i+1<argc) fq1=argv[++i];
        else if(!strcmp(argv[i],"-b")&&i+1<argc) fq2=argv[++i];
        else if(!strcmp(argv[i],"-ua")&&i+1<argc) ua=argv[++i];
        else if(!strcmp(argv[i],"-ub")&&i+1<argc) ub=argv[++i];
        else if(!strcmp(argv[i],"-s")&&i+1<argc) stats=argv[++i];
        else if(!strcmp(argv[i],"-o")&&i+1<argc) outp=argv[++i];
        else { usage(argv[0]); return 1; }
    }
    if(!fq1||!fq2||!ua||!ub||!stats){ usage(argv[0]); return 1; }

    FILE*out=stdout; if(strcmp(outp,"-")) { out=fopen(outp,"w"); if(!out){perror("open out"); return 1;} }
    setvbuf(stdin, NULL, _IOFBF, 1<<20);
    if(out==stdout) setvbuf(stdout, NULL, _IOFBF, 1<<20); else setvbuf(out, NULL, _IOFBF, 1<<20);
    gzFile g1=gzopen(ua,"wb"); if(!g1){perror("gzopen ua"); return 1;}
    gzFile g2=gzopen(ub,"wb"); if(!g2){perror("gzopen ub"); return 1;}
    gzbuffer(g1, 1<<20);
    gzbuffer(g2, 1<<20);
    FastqCtx fq; fqctx_init(&fq, fq1, fq2);
    Stats st; stats_init(&st);

    long long primary_aligned_bp=0, primary_soft_bp=0, supplementary_aligned_bp=0;

    char *line=NULL; size_t cap=0; ssize_t rl;
    char *cur_q=NULL;
    BamRow *group=NULL; int ng=0; int gcap=0; int wrote_badmap=0; char *last_rg=NULL;

    while((rl=safe_getline(&line,&cap,stdin))>0){
        if(line[0]=='@'){ fprintf(out, "%s\n", line); continue; }
        char *cpy=sdup(line);
        BamRow br=parse_bam_line(cpy); // takes ownership of cpy
        if(!cur_q){ cur_q=sdup(br.qname); }
        if(!group || ng>=gcap){ gcap = gcap? gcap*2 : 64; group=(BamRow*)realloc(group, gcap*sizeof(BamRow)); }
        if(strcmp(cur_q, br.qname)==0){ group[ng++]=br; }
        else { flush_group(out, g1, g2, group, &ng, &cur_q, &primary_aligned_bp, &primary_soft_bp, &supplementary_aligned_bp, &wrote_badmap, &fq, &st, &last_rg); cur_q=sdup(br.qname); group[ng++]=br; }
    }
    flush_group(out, g1, g2, group, &ng, &cur_q, &primary_aligned_bp, &primary_soft_bp, &supplementary_aligned_bp, &wrote_badmap, &fq, &st, &last_rg);

    // Drain leftover FASTQs: re-add as unmapped SAM (do not count in fragments/alignments/total_bp)
    FastqPair *left=NULL; while(fq_retrieve_any(&fq, &left)){
        stats_add_size(&st.size_r1, &st.size_r1_cap, (int)strlen(left->seq1));
        stats_add_size(&st.size_r2, &st.size_r2_cap, (int)strlen(left->seq2));
        st.restored_bp_r1 += (long long)strlen(left->seq1);
        st.restored_bp_r2 += (long long)strlen(left->seq2);
        st.readded_fragments += 1;
        print_unmapped_sam_pair(out, left, last_rg);
        fastqpair_free(left); left=NULL;
    }

    if(out!=stdout) fclose(out);
    gzclose(g1); gzclose(g2);
    fqctx_free(&fq);
    if(last_rg) free(last_rg);

    double p_soft_ratio = 0.0; if((primary_aligned_bp + primary_soft_bp) > 0) p_soft_ratio = (double)primary_soft_bp / (double)(primary_aligned_bp + primary_soft_bp);
    double supp_ratio = 0.0; if((primary_aligned_bp + supplementary_aligned_bp) > 0) supp_ratio = (double)supplementary_aligned_bp / (double)(primary_aligned_bp + supplementary_aligned_bp);
    fprintf(stderr, "Storing statistics in file %s\n", stats);
    // Also print all stats to stderr
    fprintf(stderr, "alignments\t %lld\n", st.alignments);
    fprintf(stderr, "fragments\t %lld\n", st.fragments);
    fprintf(stderr, "total_bp\t %lld\n", st.total_bp);
    fprintf(stderr, "primary_reads\t %lld\n", st.primary_reads);
    fprintf(stderr, "primary_aligned_bp\t %lld\n", primary_aligned_bp);
    fprintf(stderr, "primary_soft_clipped_bp\t %lld\n", primary_soft_bp);
    fprintf(stderr, "restored_read1s\t %lld\n", st.restored_read1s);
    fprintf(stderr, "restored_bp_read1\t %lld\n", st.restored_bp_r1);
    fprintf(stderr, "restored_read2s\t %lld\n", st.restored_read2s);
    fprintf(stderr, "restored_bp_read2\t %lld\n", st.restored_bp_r2);
    fprintf(stderr, "restored_unaligned_reads\t %lld\n", st.restored_unaligned_reads);
    fprintf(stderr, "supplementary_alignments\t %lld\n", st.supplementary_alignments);
    fprintf(stderr, "supplementary_alignments_cigar_adapted\t %lld\n", st.supplementary_alignments_cigar_adapted);
    fprintf(stderr, "supplementary_aligned_bp\t %lld\n", supplementary_aligned_bp);
    fprintf(stderr, "chrEBV_mapped_fragments\t %lld\n", st.chrEBV_mapped_fragments);
    fprintf(stderr, "badly_mapped_fragments\t %lld\n", st.badly_mapped_fragments);
    fprintf(stderr, "fully_unmapped_fragments\t %lld\n", st.fully_unmapped_fragments);
    fprintf(stderr, "readded_fragments\t %lld\n", st.readded_fragments);
    fprintf(stderr, "primary_soft_clipped_bp_ratio\t %.5f\n", p_soft_ratio);
    fprintf(stderr, "supplementary_bp_ratio\t %.5f\n", supp_ratio);

    FILE*sf=fopen(stats,"w"); if(!sf){perror("stats"); return 1; }
    fprintf(sf, "alignments\t %lld\n", st.alignments);
    fprintf(sf, "fragments\t %lld\n", st.fragments);
    fprintf(sf, "total_bp\t %lld\n", st.total_bp);
    fprintf(sf, "primary_reads\t %lld\n", st.primary_reads);
    fprintf(sf, "primary_aligned_bp\t %lld\n", primary_aligned_bp);
    fprintf(sf, "primary_soft_clipped_bp\t %lld\n", primary_soft_bp);
    fprintf(sf, "restored_read1s\t %lld\n", st.restored_read1s);
    fprintf(sf, "restored_bp_read1\t %lld\n", st.restored_bp_r1);
    fprintf(sf, "restored_read2s\t %lld\n", st.restored_read2s);
    fprintf(sf, "restored_bp_read2\t %lld\n", st.restored_bp_r2);
    fprintf(sf, "restored_unaligned_reads\t %lld\n", st.restored_unaligned_reads);
    fprintf(sf, "supplementary_alignments\t %lld\n", st.supplementary_alignments);
    fprintf(sf, "supplementary_alignments_cigar_adapted\t %lld\n", st.supplementary_alignments_cigar_adapted);
    fprintf(sf, "supplementary_aligned_bp\t %lld\n", supplementary_aligned_bp);
    fprintf(sf, "chrEBV_mapped_fragments\t %lld\n", st.chrEBV_mapped_fragments);
    fprintf(sf, "badly_mapped_fragments\t %lld\n", st.badly_mapped_fragments);
    fprintf(sf, "fully_unmapped_fragments\t %lld\n", st.fully_unmapped_fragments);
    fprintf(sf, "readded_fragments\t %lld\n", st.readded_fragments);
    fprintf(sf, "primary_soft_clipped_bp_ratio\t %.5f\n", p_soft_ratio);
    fprintf(sf, "supplementary_bp_ratio\t %.5f\n", supp_ratio);
    // adapter size table
    fprintf(sf, "\nadapter_length\tcount_read1\tcount_read2\n");
    int mx = 0;
    for(int i=st.size_r1_cap-1;i>0;i--){ if(st.size_r1[i]>0){ if(i>mx) mx=i; break; } }
    for(int i=st.size_r2_cap-1;i>0;i--){ if(st.size_r2[i]>0){ if(i>mx) mx=i; break; } }
    for(int i=1;i<=mx;i++){
        int c1 = (i<st.size_r1_cap)? st.size_r1[i] : 0;
        int c2 = (i<st.size_r2_cap)? st.size_r2[i] : 0;
        fprintf(sf, "%d\t%d\t%d\n", i, c1, c2);
    }
    fflush(sf); fclose(sf);
    fprintf(stderr, "done\n");

    // Non-empty checks for s, ua, ub
    struct stat stS, stA, stB;
    if(stat(stats, &stS)!=0 || stS.st_size == 0){ fprintf(stderr, "Stat file is empty\n"); return 1; }
    if(stat(ua, &stA)!=0 || stA.st_size == 0){ fprintf(stderr, "Badmap FQ1 file is empty\n"); return 1; }
    if(stat(ub, &stB)!=0 || stB.st_size == 0){ fprintf(stderr, "Badmap FQ2 file is empty\n"); return 1; }
    stats_free(&st);
    return 0;
}
