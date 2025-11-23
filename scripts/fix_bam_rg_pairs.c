#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <errno.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

typedef struct { char *old_id; char *new_id; char *pu_norm; char *pu_orig; } RgMap;

typedef struct { RgMap *v; size_t n; size_t cap; } RgVec;

static void *xmalloc(size_t n){ void *p = malloc(n); if(!p){ fprintf(stderr, "OOM %zu\n", n); exit(1);} return p; }
static char *xstrdup(const char *s){ if(!s) return NULL; size_t n=strlen(s); char *p=xmalloc(n+1); memcpy(p,s,n+1); return p; }
static void rgvec_push(RgVec *rv, RgMap m){ if(rv->n==rv->cap){ rv->cap = rv->cap? rv->cap*2: 16; rv->v = (RgMap*)realloc(rv->v, rv->cap*sizeof(RgMap)); if(!rv->v){ fprintf(stderr,"OOM rgvec\n"); exit(1);} } rv->v[rv->n++] = m; }
static void rgmap_free(RgMap *m){ if(!m) return; free(m->old_id); free(m->new_id); free(m->pu_norm); free(m->pu_orig); memset(m,0,sizeof(*m)); }
static void rgvec_free(RgVec *rv){ for(size_t i=0;i<rv->n;i++) rgmap_free(&rv->v[i]); free(rv->v); rv->v=NULL; rv->n=rv->cap=0; }

// Strict normalization: drop the PENULTIMATE token when it equals "1" or "3".
// Example: TM014_4_1_GATCAG -> TM014_4_GATCAG (keeps index token at end).
// Requires at least two underscores.
static char *normalize_pu_last_drop(const char *pu){
    if(!pu) return NULL;
    char *tmp = xstrdup(pu);
    char *last = strrchr(tmp, '_');
    if(!last){
        free(tmp);
        return NULL; // needs at least one underscore before final token
    }
    char *final_token = last + 1;
    *last = '\0';
    char *pen = strrchr(tmp, '_');
    if(!pen){
        free(tmp);
        return NULL; // requires penultimate token
    }
    char *pen_token = pen + 1;
    if(!(strcmp(pen_token, "1")==0 || strcmp(pen_token, "3")==0)){
        free(tmp);
        return NULL;
    }
    *pen = '\0';
    const char *prefix = tmp;
    size_t prefix_len = strlen(prefix);
    size_t final_len = strlen(final_token);
    size_t need = prefix_len + (prefix_len>0 && final_len>0 ? 1 : 0) + final_len + 1;
    char *res = (char*)xmalloc(need);
    size_t k = 0;
    if(prefix_len>0){
        memcpy(res+k, prefix, prefix_len);
        k += prefix_len;
        if(final_len>0){
            res[k++] = '_';
        }
    }
    if(final_len>0){
        memcpy(res+k, final_token, final_len);
        k += final_len;
    }
    res[k] = '\0';
    free(tmp);
    return res;
}

static void parse_rg_fields(const char *line, char **ID, char **PU){ *ID=*PU=NULL; const char *p=line; if(strncmp(p,"@RG\t",4)==0) p+=4; char *tmp=xstrdup(p); char *save=NULL; for(char *tok=strtok_r(tmp,"\t\n\r",&save); tok; tok=strtok_r(NULL,"\t\n\r",&save)){ if(strncmp(tok,"ID:",3)==0) *ID=xstrdup(tok+3); else if(strncmp(tok,"PU:",3)==0) *PU=xstrdup(tok+3); } free(tmp); }

static char* rebuild_hdr_without_rg(const char *text, size_t l_text){ if(!text||l_text==0) return xstrdup(""); char *buf=(char*)xmalloc(l_text+2); size_t k=0; const char *p=text; const char *end=text+l_text; while(p<end){ const char *nl=memchr(p,'\n',(size_t)(end-p)); size_t len=nl?(size_t)(nl-p):(size_t)(end-p); if(!(len>=3 && p[0]=='@' && p[1]=='R' && p[2]=='G')){ memcpy(buf+k,p,len); k+=len; buf[k++]='\n'; } if(!nl) break; p=nl+1; } buf[k]='\0'; return buf; }

static void build_mappings_from_hdr(const char *text, size_t l_text, RgVec *maps){
    const char *p=text; const char *end=text+l_text;
    while(p<end){
        const char *nl=memchr(p,'\n',(size_t)(end-p)); size_t len=nl?(size_t)(nl-p):(size_t)(end-p);
        if(len>=3 && p[0]=='@' && p[1]=='R' && p[2]=='G'){
            char *line=(char*)xmalloc(len+1); memcpy(line,p,len); line[len]='\0';
            char *id=NULL,*pu=NULL; parse_rg_fields(line,&id,&pu);
            if(id){
                const char *src = pu? pu: id;
                char *norm = normalize_pu_last_drop(src);
                if(!norm){
                    norm = xstrdup(src);
                }
                RgMap m; memset(&m,0,sizeof(m)); m.old_id=id; m.new_id=xstrdup(norm); m.pu_norm=norm; m.pu_orig=pu; rgvec_push(maps,m);
            } else { free(id); free(pu); }
            free(line);
        }
        if(!nl){
            break;
        }
        p = nl + 1;
    }
}

static int uniq_seen(const char **arr, int n, const char *s){ for(int i=0;i<n;i++){ if(strcmp(arr[i],s)==0) return 1; } return 0; }

static char* build_new_rg_block(const RgVec *maps){ // emit one RG per unique new_id
    size_t cap=maps->n*256+16; char *buf=(char*)xmalloc(cap); size_t k=0; const char *emitted[1024]; int ne=0; for(size_t i=0;i<maps->n;i++){ const char *nid=maps->v[i].new_id; if(uniq_seen(emitted,ne,nid)) continue; emitted[ne++]=nid; const char *pu = maps->v[i].pu_norm?maps->v[i].pu_norm:""; int w = snprintf(buf+k,cap-k,"@RG\tID:%s\tPU:%s\n", nid, pu); if(w<0) w=0; k += (size_t)w; if(k+256>cap){ cap*=2; buf=(char*)realloc(buf,cap); }
 }
 buf[k]='\0'; return buf; }

static const char* map_rg_id(const RgVec *maps, const char *old){ if(!old) return NULL; for(size_t i=0;i<maps->n;i++){ if(maps->v[i].old_id && strcmp(maps->v[i].old_id,old)==0) return maps->v[i].new_id; } return NULL; }

static int infer_readnr_from_qname(const char *q){ if(!q) return 0; const char *slash = strrchr(q,'/'); if(slash && slash[1]){ if(slash[1]=='1') return 1; if(slash[1]=='2' || slash[1]=='3') return 2; }
 // also handle '#x/1' within
 const char *hash=strrchr(q,'#'); if(hash){ const char *sl=strchr(hash,'/'); if(sl && sl[1]){ if(sl[1]=='1') return 1; if(sl[1]=='2' || sl[1]=='3') return 2; } }
 return 0; }

static void qname_fix_inplace(char *q){ if(!q) return; size_t n=strlen(q); for(size_t i=0;i+1<n;i++){ if(q[i]=='/' && q[i+1]=='3'){ q[i+1]='2'; } } }

static int main_impl(const char *inpath, const char *outpath, const char *ref, int threads, int no_rename){ htsFile *in=hts_open(inpath, "r"); if(!in){ fprintf(stderr,"Cannot open %s\n", inpath); return 1; } if(ref && *ref) hts_set_fai_filename(in, ref); if(threads>1) hts_set_threads(in, threads);
 bam_hdr_t *hdr = sam_hdr_read(in); if(!hdr){ fprintf(stderr,"Failed to read header\n"); hts_close(in); return 1; }
 RgVec maps={0}; build_mappings_from_hdr(hdr->text, (size_t)hdr->l_text, &maps);
 char *base = rebuild_hdr_without_rg(hdr->text, (size_t)hdr->l_text); char *rgblk = build_new_rg_block(&maps);
 size_t newlen = strlen(base) + strlen(rgblk) + 1; char *newtext=(char*)xmalloc(newlen); newtext[0]='\0'; strcat(newtext, base); strcat(newtext, rgblk);
 free(hdr->text); hdr->text=newtext; hdr->l_text=(int)strlen(newtext);
 free(base); free(rgblk);
    // choose mode by extension
    const char *mode = "wb";
    const char *ext = strrchr(outpath, '.');
    if(ext){
        if(strcmp(ext, ".cram")==0){
            mode = "wc";
        } else if(strcmp(ext, ".sam")==0){
            mode = "w";
        }
    }
 htsFile *out=hts_open(outpath, mode); if(!out){ fprintf(stderr,"Cannot open %s\n", outpath); rgvec_free(&maps); bam_hdr_destroy(hdr); hts_close(in); return 1; }
 if(mode[1]=='c' && ref && *ref){
     hts_set_fai_filename(out, ref);
 }
 if(threads>1){
     hts_set_threads(out, threads);
 }
 if(sam_hdr_write(out, hdr)<0){ fprintf(stderr,"Failed to write header\n"); hts_close(out); rgvec_free(&maps); bam_hdr_destroy(hdr); hts_close(in); return 1; }
 bam1_t *b = bam_init1(); if(!b){ fprintf(stderr,"OOM bam1\n"); hts_close(out); rgvec_free(&maps); bam_hdr_destroy(hdr); hts_close(in); return 1; }
 while(sam_read1(in, hdr, b) >= 0){
 uint8_t *rg = bam_aux_get(b, "RG"); const char *oldrg = rg? bam_aux2Z(rg): NULL; const char *newrg = map_rg_id(&maps, oldrg); if(newrg){ bam_aux_update_str(b, "RG", (uint32_t)strlen(newrg)+1, newrg); }
 char *qname = bam_get_qname(b); int rn = infer_readnr_from_qname(qname);
 if(rn==0){
     fprintf(stderr, "ERROR: cannot infer read number from QNAME: %s\n", qname?qname:"(null)");
     bam_destroy1(b); hts_close(out); rgvec_free(&maps); bam_hdr_destroy(hdr); hts_close(in);
     return 2;
 }
 uint16_t f = b->core.flag;
 if(!no_rename){ qname_fix_inplace(qname); }
 f |= 0x1;
 f &= ~(0x40|0x80);
 if(rn==1){ f |= 0x40; }
 else /* rn==2 */ { f |= 0x80; }
 b->core.flag = f;
 if(sam_write1(out, hdr, b) < 0){ fprintf(stderr,"Write failed\n"); bam_destroy1(b); hts_close(out); rgvec_free(&maps); bam_hdr_destroy(hdr); hts_close(in); return 1; }
 }
 bam_destroy1(b); hts_close(out); rgvec_free(&maps); bam_hdr_destroy(hdr); hts_close(in); return 0; }

static void usage(const char *p){ fprintf(stderr, "Usage: %s -i in.bam -o out.bam [-r ref.fa] [--threads N] [--no-rename]\n", p); }

int main(int argc, char **argv){ const char *inpath=NULL,*outpath=NULL,*ref=NULL; int threads=1; int no_rename=0; for(int i=1;i<argc;i++){ if(!strcmp(argv[i],"-i")&&i+1<argc) inpath=argv[++i]; else if(!strcmp(argv[i],"-o")&&i+1<argc) outpath=argv[++i]; else if(!strcmp(argv[i],"-r")&&i+1<argc) ref=argv[++i]; else if(!strcmp(argv[i],"--threads")&&i+1<argc) threads=atoi(argv[++i]); else if(!strcmp(argv[i],"--no-rename")) no_rename=1; else { usage(argv[0]); return 1; } }
 if(!inpath||!outpath){ usage(argv[0]); return 1; }
 return main_impl(inpath,outpath,ref,threads,no_rename);
}
