#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

static inline uint64_t fnv1a64_init(void) {
    return 0xcbf29ce484222325ULL;
}
static inline uint64_t fnv1a64_update(uint64_t h, const char *data, size_t len) {
    const uint64_t prime = 0x100000001b3ULL;
    for (size_t i = 0; i < len; ++i) { h ^= (unsigned char)data[i]; h *= prime; }
    return h;
}
static inline uint64_t fnv1a64_update_1(uint64_t h, unsigned char c) {
    const uint64_t prime = 0x100000001b3ULL; h ^= c; h *= prime; return h;
}
static inline char qual_norm(char c) { return (c == '#') ? '!' : c; }
static inline char comp_base(char c) {
    switch (c) {
        case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
        case 'a': return 't'; case 'c': return 'g'; case 'g': return 'c'; case 't': return 'a';
        default: return c;
    }
}

static uint64_t hash_seq_with_orientation(const char *yb, size_t yb_len,
                                          const char *seq, size_t seq_len,
                                          const char *zb, size_t zb_len,
                                          int reversed) {
    uint64_t h = fnv1a64_init();
    if (!reversed) {
        if (yb_len) h = fnv1a64_update(h, yb, yb_len);
        if (seq_len) h = fnv1a64_update(h, seq, seq_len);
        if (zb_len) h = fnv1a64_update(h, zb, zb_len);
    } else {
        for (ssize_t i = (ssize_t)zb_len - 1; i >= 0; --i) h = fnv1a64_update_1(h, (unsigned char)comp_base(zb[i]));
        for (ssize_t i = (ssize_t)seq_len - 1; i >= 0; --i) h = fnv1a64_update_1(h, (unsigned char)comp_base(seq[i]));
        for (ssize_t i = (ssize_t)yb_len - 1; i >= 0; --i) h = fnv1a64_update_1(h, (unsigned char)comp_base(yb[i]));
    }
    return h;
}

static uint64_t hash_qual_with_orientation(const char *yq, size_t yq_len,
                                           const char *qual, size_t qual_len,
                                           const char *zq, size_t zq_len,
                                           int reversed) {
    uint64_t h = fnv1a64_init();
    if (!reversed) {
        for (size_t i = 0; i < yq_len; ++i) h = fnv1a64_update_1(h, (unsigned char)qual_norm(yq[i]));
        for (size_t i = 0; i < qual_len; ++i) h = fnv1a64_update_1(h, (unsigned char)qual_norm(qual[i]));
        for (size_t i = 0; i < zq_len; ++i) h = fnv1a64_update_1(h, (unsigned char)qual_norm(zq[i]));
    } else {
        for (ssize_t i = (ssize_t)zq_len - 1; i >= 0; --i) h = fnv1a64_update_1(h, (unsigned char)qual_norm(zq[i]));
        for (ssize_t i = (ssize_t)qual_len - 1; i >= 0; --i) h = fnv1a64_update_1(h, (unsigned char)qual_norm(qual[i]));
        for (ssize_t i = (ssize_t)yq_len - 1; i >= 0; --i) h = fnv1a64_update_1(h, (unsigned char)qual_norm(yq[i]));
    }
    return h;
}

static PyObject *py_bam_stats(PyObject *self, PyObject *args, PyObject *kwargs) {
    (void)self;
    static char *kwlist[] = {"path", "threads", "reference", NULL};
    const char *path = NULL;
    int threads = 1;
    const char *reference = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s|is", kwlist, &path, &threads, &reference)) {
        return NULL;
    }

    htsFile *fp = hts_open(path, "r");
    if (!fp) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to open BAM/CRAM");
        return NULL;
    }
    if (threads > 1) {
        hts_set_threads(fp, threads);
    }
    if (reference && reference[0] != '\0') {
        hts_set_fai_filename(fp, reference);
    }

    bam_hdr_t *hdr = sam_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        PyErr_SetString(PyExc_RuntimeError, "Failed to read header");
        return NULL;
    }

    bam1_t *b = bam_init1();
    if (!b) {
        bam_hdr_destroy(hdr);
        hts_close(fp);
        PyErr_SetString(PyExc_RuntimeError, "Failed to allocate BAM record");
        return NULL;
    }

    unsigned long long c1_seq = 0, c1_qual = 0, c2_seq = 0, c2_qual = 0;
    long long nrow1 = 0, nrow2 = 0, nbases1 = 0, nbases2 = 0;

    char *seq = NULL;
    char *qual = NULL;
    int seq_cap = 0;
    int qual_cap = 0;

    const char *nt_table = "=ACMGRSVTWYHKDBN"; // htslib seq_nt16_str, but local to avoid linkage

    while (sam_read1(fp, hdr, b) >= 0) {
        uint16_t flag = b->core.flag;
        if ((flag & 0x100) || (flag & 0x800)) continue; // skip secondary/supplementary

        int is_read1 = (flag & 0x40) ? 1 : 0;
        int reversed = (flag & 0x10) ? 1 : 0;
        int unmapped = (flag & 0x4) ? 1 : 0;

        // sequence
        uint8_t *seqptr = bam_get_seq(b);
        int lseq = b->core.l_qseq;
        // convert to bases
        if ((int)((size_t)lseq + 1) > seq_cap) {
            char *tmp = (char*)realloc(seq, (size_t)lseq + 1);
            if (!tmp) { free(seq); free(qual); bam_destroy1(b); bam_hdr_destroy(hdr); hts_close(fp); PyErr_NoMemory(); return NULL; }
            seq = tmp;
            seq_cap = (int)((size_t)lseq + 1);
        }
        for (int i = 0; i < lseq; ++i) {
            seq[i] = nt_table[bam_seqi(seqptr, i)];
        }
        seq[lseq] = '\0';

        // qualities
        uint8_t *qualptr = bam_get_qual(b);
        if ((int)((size_t)lseq + 1) > qual_cap) {
            char *tmpq = (char*)realloc(qual, (size_t)lseq + 1);
            if (!tmpq) { free(seq); free(qual); bam_destroy1(b); bam_hdr_destroy(hdr); hts_close(fp); PyErr_NoMemory(); return NULL; }
            qual = tmpq;
            qual_cap = (int)((size_t)lseq + 1);
        }
        if (qualptr && qualptr[0] != 0xFF) {
            for (int i = 0; i < lseq; ++i) { qual[i] = (char)(qualptr[i] + 33); }
        } else {
            // missing quals -> fill with '!' like SAM convention
            for (int i = 0; i < lseq; ++i) { qual[i] = '!'; }
        }
        qual[lseq] = '\0';

        // Optional tags YB/ZB/YQ/ZQ
        uint8_t *yb_tag = bam_aux_get(b, "YB");
        const char *yb = NULL; size_t yb_len = 0;
        if (yb_tag) { yb = bam_aux2Z(yb_tag); if (yb) yb_len = strlen(yb); }
        uint8_t *zb_tag = bam_aux_get(b, "ZB");
        const char *zb = NULL; size_t zb_len = 0;
        if (zb_tag) { zb = bam_aux2Z(zb_tag); if (zb) zb_len = strlen(zb); }
        uint8_t *yq_tag = bam_aux_get(b, "YQ");
        const char *yq = NULL; size_t yq_len = 0;
        if (yq_tag) { yq = bam_aux2Z(yq_tag); if (yq) yq_len = strlen(yq); }
        uint8_t *zq_tag = bam_aux_get(b, "ZQ");
        const char *zq = NULL; size_t zq_len = 0;
        if (zq_tag) { zq = bam_aux2Z(zq_tag); if (zq) zq_len = strlen(zq); }

        int do_restore = (!unmapped) || yb_len > 0 || zb_len > 0 || yq_len > 0 || zq_len > 0;

        uint64_t hs = fnv1a64_init();
        uint64_t hq = fnv1a64_init();
        long long total_bases = 0;
        if (do_restore) {
            hs = hash_seq_with_orientation(yb ? yb : "", yb_len, seq, (size_t)lseq, zb ? zb : "", zb_len, reversed);
            hq = hash_qual_with_orientation(yq ? yq : "", yq_len, qual, (size_t)lseq, zq ? zq : "", zq_len, reversed);
            total_bases = (long long)(yb_len + (size_t)lseq + zb_len);
        } else {
            hs = fnv1a64_update(hs, seq, (size_t)lseq);
            for (int i = 0; i < lseq; ++i) hq = fnv1a64_update_1(hq, (unsigned char)qual_norm(qual[i]));
            total_bases = (long long)lseq;
        }

        if (is_read1) {
            c1_seq ^= hs; c1_qual ^= hq; nrow1 += 1; nbases1 += total_bases;
        } else {
            c2_seq ^= hs; c2_qual ^= hq; nrow2 += 1; nbases2 += total_bases;
        }

        // buffers reused across iterations
    }

    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    hts_close(fp);
    free(seq);
    free(qual);

    return Py_BuildValue("KKLLKKLL",
                         (unsigned long long)c1_seq, (unsigned long long)c1_qual, (long long)nrow1, (long long)nbases1,
                         (unsigned long long)c2_seq, (unsigned long long)c2_qual, (long long)nrow2, (long long)nbases2);
}

static PyMethodDef FastcheckHTSMethods[] = {
    {"bam_stats", (PyCFunction)py_bam_stats, METH_VARARGS | METH_KEYWORDS, "Compute read checksums/stats directly from BAM/CRAM using htslib (bam_stats(path, threads=1, reference=None))"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef fastcheck_hts_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "fastcheck_hts",
    .m_doc = NULL,
    .m_size = -1,
    .m_methods = FastcheckHTSMethods
};

PyMODINIT_FUNC PyInit_fastcheck_hts(void) {
    return PyModule_Create(&fastcheck_hts_module);
}
