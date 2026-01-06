#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

static inline uint64_t fnv1a64_init(void) {
    return 0xcbf29ce484222325ULL;
}

static ssize_t safe_getline(char **lineptr, size_t *n, FILE *stream);
static inline uint64_t fnv1a64_update(uint64_t h, const char *data, size_t len);
static inline uint64_t fnv1a64_update_1(uint64_t h, unsigned char c);
static inline char qual_norm(char c);

typedef struct {
    char *q1;
    char *q2;
    ssize_t l1;
    ssize_t l2;
} QualPair;

static inline void qual_minmax_update(const char *q, ssize_t len, int *minq, int *maxq) {
    for (ssize_t i = 0; i < len; ++i) {
        unsigned char c = (unsigned char)qual_norm(q[i]);
        if ((int)c < *minq) *minq = (int)c;
        if ((int)c > *maxq) *maxq = (int)c;
    }
}

static inline int determine_qual_shift(int minq, int maxq) {
    if (minq == INT_MAX) return 0;
    if (minq < 64) return 0;
    if (maxq > 74) return -31;
    return 0;
}

static inline uint64_t hash_qual_line_with_shift(const char *q, ssize_t len, int shift) {
    uint64_t h = fnv1a64_init();
    for (ssize_t i = 0; i < len; ++i) {
        unsigned char c = (unsigned char)qual_norm(q[i]);
        c = (unsigned char)((int)c + shift);
        h = fnv1a64_update_1(h, c);
    }
    return h;
}

static PyObject *py_fastq_stats_interleaved(PyObject *self, PyObject *args) {
    (void)self;
    const char *path;
    if (!PyArg_ParseTuple(args, "s", &path)) {
        return NULL;
    }

    FILE *f = NULL;
    int need_pclose = 0;
    if (strcmp(path, "-") == 0) {
        f = stdin;
        need_pclose = 0;
    } else {
        char cmd[4096];
        if (snprintf(cmd, sizeof(cmd), "pigz -dc '%s'", path) >= (int)sizeof(cmd)) {
            PyErr_SetString(PyExc_RuntimeError, "Command too long");
            return NULL;
        }
        f = popen(cmd, "r");
        if (!f) {
            PyErr_SetFromErrno(PyExc_OSError);
            return NULL;
        }
        need_pclose = 1;
    }

    char *h1=NULL,*s1=NULL,*p1=NULL,*q1=NULL;
    char *h2=NULL,*s2=NULL,*p2=NULL,*q2=NULL;
    size_t n1=0,n2=0,n3=0,n4=0;
    size_t m1=0,m2=0,m3=0,m4=0;

    uint64_t c1_seq = 0, c1_qual = 0, c2_seq = 0, c2_qual = 0;
    long long nrow = 0;
    long long nbases1 = 0, nbases2 = 0;

    int have_splitchar = 0;
    char splitchar = ' ';
    int layout = 0;

    int qual_shift = 0;
    int qual_shift_determined = 0;
    int minq = INT_MAX;
    int maxq = 0;
    QualPair *qbuf = NULL;
    size_t qbuf_n = 0;
    size_t qbuf_cap = 0;

    while (1) {
        ssize_t rh1 = safe_getline(&h1, &n1, f);
        if (rh1 <= 0) break;
        if (h1[0] != '@') {
            if (need_pclose) pclose(f);
            free(h1); free(s1); free(p1); free(q1);
            free(h2); free(s2); free(p2); free(q2);
            PyErr_SetString(PyExc_RuntimeError, "Parsing error: FASTQ headers must start with '@'");
            return NULL;
        }
        ssize_t rs1 = safe_getline(&s1, &n2, f);
        if (rs1 < 0) break;
        ssize_t r3 = safe_getline(&p1, &n3, f);
        if (r3 < 0) break;
        if (layout == 0) {
            if (p1[0] == '+') layout = 4;
            else if (p1[0] == '@') layout = 2;
            else {
                if (need_pclose) pclose(f);
                free(h1); free(s1); free(p1); free(q1);
                free(h2); free(s2); free(p2); free(q2);
                PyErr_SetString(PyExc_RuntimeError, "Cannot determine FASTQ layout (expected '+' or '@')");
                return NULL;
            }
        }
        if (!have_splitchar) {
            if (strchr(h1, ' ')) { splitchar = ' '; have_splitchar = 1; }
            else if (strchr(h1, '/')) { splitchar = '/'; have_splitchar = 1; }
            else { splitchar = ' '; have_splitchar = 1; }
        }
        ssize_t rs2 = -1;
        if (layout == 4) {
            ssize_t rq1_ = safe_getline(&q1, &n4, f);
            if (rq1_ < 0) break;
            ssize_t rh2 = safe_getline(&h2, &m1, f);
            if (rh2 <= 0) break;
            rs2 = safe_getline(&s2, &m2, f);
            if (rs2 < 0) break;
            ssize_t rp2_ = safe_getline(&p2, &m3, f);
            if (rp2_ < 0) break;
            if (p2[0] != '+') {
                if (need_pclose) pclose(f);
                free(h1); free(s1); free(p1); free(q1);
                free(h2); free(s2); free(p2); free(q2);
                PyErr_SetString(PyExc_RuntimeError, "Malformed FASTQ: expected '+' line for read2");
                return NULL;
            }
            ssize_t rq2_ = safe_getline(&q2, &m4, f);
            if (rq2_ < 0) break;
        } else {
            h2 = p1; p1 = NULL; m1 = n3; n3 = 0;
            rs2 = safe_getline(&s2, &m2, f);
            if (rs2 < 0) break;
        }
        char *split1 = strchr(h1, splitchar);
        char *split2 = strchr(h2, splitchar);
        if (split1) *split1 = '\0';
        if (split2) *split2 = '\0';
        const char *name1 = h1 + 1;
        const char *name2 = h2 + 1;
        if (strcmp(name1, name2) != 0) {
            if (need_pclose) pclose(f);
            free(h1); free(s1); free(p1); free(q1);
            free(h2); free(s2); free(p2); free(q2);
            PyErr_SetString(PyExc_RuntimeError, "FastQ input records not in same read order (interleaved)");
            return NULL;
        }
        uint64_t hseq1 = fnv1a64_init();
        hseq1 = fnv1a64_update(hseq1, s1, (size_t)rs1);
        c1_seq ^= hseq1;
        uint64_t hseq2 = fnv1a64_init();
        hseq2 = fnv1a64_update(hseq2, s2, (size_t)rs2);
        c2_seq ^= hseq2;
        if (layout == 4) {
            if (!qual_shift_determined) {
                qual_minmax_update(q1, (ssize_t)strlen(q1), &minq, &maxq);
                qual_minmax_update(q2, (ssize_t)strlen(q2), &minq, &maxq);
                if (qbuf_n == qbuf_cap) {
                    size_t ncap = qbuf_cap ? (qbuf_cap * 2) : 256;
                    QualPair *nx = (QualPair *)realloc(qbuf, ncap * sizeof(QualPair));
                    if (!nx) {
                        if (need_pclose) pclose(f);
                        free(h1); free(s1); free(p1); free(q1);
                        free(h2); free(s2); free(p2); free(q2);
                        for (size_t j = 0; j < qbuf_n; ++j) { free(qbuf[j].q1); free(qbuf[j].q2); }
                        free(qbuf);
                        PyErr_NoMemory();
                        return NULL;
                    }
                    qbuf = nx;
                    qbuf_cap = ncap;
                }
                ssize_t lq1 = (ssize_t)strlen(q1);
                ssize_t lq2 = (ssize_t)strlen(q2);
                qbuf[qbuf_n].q1 = (char *)malloc((size_t)lq1 + 1);
                qbuf[qbuf_n].q2 = (char *)malloc((size_t)lq2 + 1);
                if (!qbuf[qbuf_n].q1 || !qbuf[qbuf_n].q2) {
                    if (need_pclose) pclose(f);
                    free(h1); free(s1); free(p1); free(q1);
                    free(h2); free(s2); free(p2); free(q2);
                    if (qbuf[qbuf_n].q1) free(qbuf[qbuf_n].q1);
                    if (qbuf[qbuf_n].q2) free(qbuf[qbuf_n].q2);
                    for (size_t j = 0; j < qbuf_n; ++j) { free(qbuf[j].q1); free(qbuf[j].q2); }
                    free(qbuf);
                    PyErr_NoMemory();
                    return NULL;
                }
                memcpy(qbuf[qbuf_n].q1, q1, (size_t)lq1 + 1);
                memcpy(qbuf[qbuf_n].q2, q2, (size_t)lq2 + 1);
                qbuf[qbuf_n].l1 = lq1;
                qbuf[qbuf_n].l2 = lq2;
                qbuf_n += 1;

                if (qbuf_n >= 2000) {
                    qual_shift = determine_qual_shift(minq, maxq);
                    qual_shift_determined = 1;
                    for (size_t j = 0; j < qbuf_n; ++j) {
                        c1_qual ^= hash_qual_line_with_shift(qbuf[j].q1, qbuf[j].l1, qual_shift);
                        c2_qual ^= hash_qual_line_with_shift(qbuf[j].q2, qbuf[j].l2, qual_shift);
                        free(qbuf[j].q1);
                        free(qbuf[j].q2);
                    }
                    qbuf_n = 0;
                }
            } else {
                c1_qual ^= hash_qual_line_with_shift(q1, (ssize_t)strlen(q1), qual_shift);
                c2_qual ^= hash_qual_line_with_shift(q2, (ssize_t)strlen(q2), qual_shift);
            }
        }
        nbases1 += (long long)rs1;
        nbases2 += (long long)rs2;
        nrow += 1;
    }

    if (need_pclose) pclose(f);
    if (layout == 4 && !qual_shift_determined && qbuf_n > 0) {
        qual_shift = determine_qual_shift(minq, maxq);
        for (size_t j = 0; j < qbuf_n; ++j) {
            c1_qual ^= hash_qual_line_with_shift(qbuf[j].q1, qbuf[j].l1, qual_shift);
            c2_qual ^= hash_qual_line_with_shift(qbuf[j].q2, qbuf[j].l2, qual_shift);
            free(qbuf[j].q1);
            free(qbuf[j].q2);
        }
        qbuf_n = 0;
    }
    free(qbuf);
    free(h1); free(s1); free(p1); free(q1);
    free(h2); free(s2); free(p2); free(q2);

    return Py_BuildValue("KKKKLLL",
                         (unsigned long long)c1_seq, (unsigned long long)c1_qual,
                         (unsigned long long)c2_seq, (unsigned long long)c2_qual,
                         (long long)nrow, (long long)nbases1, (long long)nbases2);
}

static PyObject *py_fastq_stats(PyObject *self, PyObject *args) {
    (void)self;
    const char *file1, *file2;
    if (!PyArg_ParseTuple(args, "ss", &file1, &file2)) {
        return NULL;
    }

    char cmd1[4096];
    char cmd2[4096];
    if (snprintf(cmd1, sizeof(cmd1), "pigz -dc '%s'", file1) >= (int)sizeof(cmd1) ||
        snprintf(cmd2, sizeof(cmd2), "pigz -dc '%s'", file2) >= (int)sizeof(cmd2)) {
        PyErr_SetString(PyExc_RuntimeError, "Command too long");
        return NULL;
    }

    FILE *f1 = popen(cmd1, "r");
    if (!f1) {
        PyErr_SetFromErrno(PyExc_OSError);
        return NULL;
    }
    FILE *f2 = popen(cmd2, "r");
    if (!f2) {
        pclose(f1);
        PyErr_SetFromErrno(PyExc_OSError);
        return NULL;
    }

    char *h1 = NULL, *s1 = NULL, *p1 = NULL, *q1 = NULL;
    char *h2 = NULL, *s2 = NULL, *p2 = NULL, *q2 = NULL;
    size_t n1 = 0, n2 = 0, n3 = 0, n4 = 0;
    size_t m1 = 0, m2 = 0, m3 = 0, m4 = 0;

    uint64_t c1_seq = 0, c1_qual = 0, c2_seq = 0, c2_qual = 0;
    long long nrow = 0;
    long long nbases1 = 0, nbases2 = 0;

    int have_splitchar = 0;
    char splitchar = ' ';

    int qual_shift = 0;
    int qual_shift_determined = 0;
    int minq = INT_MAX;
    int maxq = 0;
    QualPair *qbuf = NULL;
    size_t qbuf_n = 0;
    size_t qbuf_cap = 0;

    while (1) {
        ssize_t rh1 = safe_getline(&h1, &n1, f1);
        ssize_t rh2 = safe_getline(&h2, &m1, f2);
        if (rh1 <= 0 || rh2 <= 0) break; // EOF
        if (h1[0] != '@' || h2[0] != '@') {
            pclose(f1); pclose(f2);
            free(h1); free(h2); free(s1); free(s2); free(p1); free(p2); free(q1); free(q2);
            PyErr_SetString(PyExc_RuntimeError, "Parsing error: FASTQ headers must start with '@'");
            return NULL;
        }
        if (!have_splitchar) {
            if (strchr(h1, ' ')) { splitchar = ' '; have_splitchar = 1; }
            else if (strchr(h1, '/')) { splitchar = '/'; have_splitchar = 1; }
            else {
                pclose(f1); pclose(f2);
                free(h1); free(h2); free(s1); free(s2); free(p1); free(p2); free(q1); free(q2);
                PyErr_SetString(PyExc_RuntimeError, "Cannot determine FastQ header split character");
                return NULL;
            }
        }
        ssize_t rs1 = safe_getline(&s1, &n2, f1);
        ssize_t rs2 = safe_getline(&s2, &m2, f2);
        if (rs1 < 0 || rs2 < 0) break;
        ssize_t rp1 = safe_getline(&p1, &n3, f1);
        ssize_t rp2 = safe_getline(&p2, &m3, f2);
        if (rp1 < 0 || rp2 < 0) break;
        ssize_t rq1 = safe_getline(&q1, &n4, f1);
        ssize_t rq2 = safe_getline(&q2, &m4, f2);
        if (rq1 < 0 || rq2 < 0) break;

        char *split1 = strchr(h1, splitchar);
        char *split2 = strchr(h2, splitchar);
        if (!split1 || !split2) {
            pclose(f1); pclose(f2);
            free(h1); free(h2); free(s1); free(s2); free(p1); free(p2); free(q1); free(q2);
            PyErr_SetString(PyExc_RuntimeError, "Malformed FASTQ header");
            return NULL;
        }
        *split1 = '\0';
        *split2 = '\0';
        const char *name1 = h1 + 1; // skip '@'
        const char *name2 = h2 + 1;
        if (strcmp(name1, name2) != 0) {
            pclose(f1); pclose(f2);
            free(h1); free(h2); free(s1); free(s2); free(p1); free(p2); free(q1); free(q2);
            PyErr_SetString(PyExc_RuntimeError, "FastQ input files not in same read order");
            return NULL;
        }
        // compute checksums
        uint64_t hseq1 = fnv1a64_init();
        hseq1 = fnv1a64_update(hseq1, s1, (size_t)rs1);
        c1_seq ^= hseq1;

        if (!qual_shift_determined) {
            qual_minmax_update(q1, rq1, &minq, &maxq);
            qual_minmax_update(q2, rq2, &minq, &maxq);
            if (qbuf_n == qbuf_cap) {
                size_t ncap = qbuf_cap ? (qbuf_cap * 2) : 256;
                QualPair *nx = (QualPair *)realloc(qbuf, ncap * sizeof(QualPair));
                if (!nx) {
                    pclose(f1); pclose(f2);
                    free(h1); free(h2); free(s1); free(s2); free(p1); free(p2); free(q1); free(q2);
                    for (size_t j = 0; j < qbuf_n; ++j) { free(qbuf[j].q1); free(qbuf[j].q2); }
                    free(qbuf);
                    PyErr_NoMemory();
                    return NULL;
                }
                qbuf = nx;
                qbuf_cap = ncap;
            }
            qbuf[qbuf_n].q1 = (char *)malloc((size_t)rq1 + 1);
            qbuf[qbuf_n].q2 = (char *)malloc((size_t)rq2 + 1);
            if (!qbuf[qbuf_n].q1 || !qbuf[qbuf_n].q2) {
                pclose(f1); pclose(f2);
                free(h1); free(h2); free(s1); free(s2); free(p1); free(p2); free(q1); free(q2);
                if (qbuf[qbuf_n].q1) free(qbuf[qbuf_n].q1);
                if (qbuf[qbuf_n].q2) free(qbuf[qbuf_n].q2);
                for (size_t j = 0; j < qbuf_n; ++j) { free(qbuf[j].q1); free(qbuf[j].q2); }
                free(qbuf);
                PyErr_NoMemory();
                return NULL;
            }
            memcpy(qbuf[qbuf_n].q1, q1, (size_t)rq1 + 1);
            memcpy(qbuf[qbuf_n].q2, q2, (size_t)rq2 + 1);
            qbuf[qbuf_n].l1 = rq1;
            qbuf[qbuf_n].l2 = rq2;
            qbuf_n += 1;

            if (qbuf_n >= 2000) {
                qual_shift = determine_qual_shift(minq, maxq);
                qual_shift_determined = 1;
                for (size_t j = 0; j < qbuf_n; ++j) {
                    c1_qual ^= hash_qual_line_with_shift(qbuf[j].q1, qbuf[j].l1, qual_shift);
                    c2_qual ^= hash_qual_line_with_shift(qbuf[j].q2, qbuf[j].l2, qual_shift);
                    free(qbuf[j].q1);
                    free(qbuf[j].q2);
                }
                qbuf_n = 0;
            }
        } else {
            c1_qual ^= hash_qual_line_with_shift(q1, rq1, qual_shift);
        }

        uint64_t hseq2 = fnv1a64_init();
        hseq2 = fnv1a64_update(hseq2, s2, (size_t)rs2);
        c2_seq ^= hseq2;

        if (qual_shift_determined) {
            c2_qual ^= hash_qual_line_with_shift(q2, rq2, qual_shift);
        }

        nbases1 += (long long)rs1;
        nbases2 += (long long)rs2;
        nrow += 1;
    }

    pclose(f1);
    pclose(f2);

    if (!qual_shift_determined && qbuf_n > 0) {
        qual_shift = determine_qual_shift(minq, maxq);
        for (size_t j = 0; j < qbuf_n; ++j) {
            c1_qual ^= hash_qual_line_with_shift(qbuf[j].q1, qbuf[j].l1, qual_shift);
            c2_qual ^= hash_qual_line_with_shift(qbuf[j].q2, qbuf[j].l2, qual_shift);
            free(qbuf[j].q1);
            free(qbuf[j].q2);
        }
        qbuf_n = 0;
    }
    free(qbuf);

    free(h1); free(h2); free(s1); free(s2); free(p1); free(p2); free(q1); free(q2);

    return Py_BuildValue("KKKKLLL", (unsigned long long)c1_seq, (unsigned long long)c1_qual,
                         (unsigned long long)c2_seq, (unsigned long long)c2_qual,
                         (long long)nrow, (long long)nbases1, (long long)nbases2);
}

static inline uint64_t fnv1a64_update(uint64_t h, const char *data, size_t len) {
    const uint64_t prime = 0x100000001b3ULL;
    for (size_t i = 0; i < len; ++i) {
        h ^= (unsigned char)data[i];
        h *= prime;
    }
    return h;
}

static inline uint64_t fnv1a64_update_1(uint64_t h, unsigned char c) {
    const uint64_t prime = 0x100000001b3ULL;
    h ^= c;
    h *= prime;
    return h;
}

static inline char qual_norm(char c) {
    return (c == '#') ? '!' : c;
}

static inline char comp_base(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'a': return 't';
        case 'c': return 'g';
        case 'g': return 'c';
        case 't': return 'a';
        default: return c;
    }
}

static ssize_t safe_getline(char **lineptr, size_t *n, FILE *stream) {
    ssize_t r = getline(lineptr, n, stream);
    if (r <= 0) return r;
    if (r > 0 && (*lineptr)[r-1] == '\n') {
        (*lineptr)[r-1] = '\0';
        --r;
        if (r > 0 && (*lineptr)[r-1] == '\r') {
            (*lineptr)[r-1] = '\0';
            --r;
        }
    }
    return r;
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

static PyObject *py_sam_stats(PyObject *self, PyObject *args) {
    (void)self;
    const char *path;
    if (!PyArg_ParseTuple(args, "s", &path)) return NULL;

    FILE *fp = NULL;
    int need_pclose = 0;
    if (strcmp(path, "-") == 0) {
        fp = stdin;
        need_pclose = 0;
    } else {
        char cmd[4096];
        if (snprintf(cmd, sizeof(cmd), "samtools view -h '%s'", path) >= (int)sizeof(cmd)) {
            PyErr_SetString(PyExc_RuntimeError, "Command too long");
            return NULL;
        }
        fp = popen(cmd, "r");
        if (!fp) {
            PyErr_SetFromErrno(PyExc_OSError);
            return NULL;
        }
        need_pclose = 1;
    }

    char *line = NULL;
    size_t cap = 0;

    unsigned long long c1_seq = 0, c1_qual = 0, c2_seq = 0, c2_qual = 0;
    long long nrow1 = 0, nrow2 = 0, nbases1 = 0, nbases2 = 0;

    while (1) {
        ssize_t r = safe_getline(&line, &cap, fp);
        if (r <= 0) break;
        if (line[0] == '@') continue;

        // Split by tab
        char *fields[64];
        int nf = 0;
        char *p = line;
        while (*p && nf < 64) {
            fields[nf++] = p;
            char *t = strchr(p, '\t');
            if (!t) break;
            *t = '\0';
            p = t + 1;
        }
        if (nf < 11) continue; // malformed

        long flag = strtol(fields[1], NULL, 10);
        if ((flag & 0x100) || (flag & 0x800)) continue;

        int is_read1 = (flag & 0x40) ? 1 : 0;
        int reversed = (flag & 0x10) ? 1 : 0;
        int unmapped = (flag & 0x4) ? 1 : 0;

        const char *seq = fields[9];
        const char *qual = fields[10];
        size_t seq_len = strlen(seq);
        size_t qual_len = strlen(qual);

        const char *yb = ""; size_t yb_len = 0;
        const char *yq = ""; size_t yq_len = 0;
        const char *zb = ""; size_t zb_len = 0;
        const char *zq = ""; size_t zq_len = 0;

        for (int i = 11; i < nf; ++i) {
            if (strncmp(fields[i], "YB:Z:", 5) == 0) { yb = fields[i] + 5; yb_len = strlen(yb); }
            else if (strncmp(fields[i], "YQ:Z:", 5) == 0) { yq = fields[i] + 5; yq_len = strlen(yq); }
            else if (strncmp(fields[i], "ZB:Z:", 5) == 0) { zb = fields[i] + 5; zb_len = strlen(zb); }
            else if (strncmp(fields[i], "ZQ:Z:", 5) == 0) { zq = fields[i] + 5; zq_len = strlen(zq); }
        }

        int do_restore = (!unmapped) || yb_len > 0 || zb_len > 0 || yq_len > 0 || zq_len > 0;
        uint64_t hs = fnv1a64_init();
        uint64_t hq = fnv1a64_init();
        long long total_bases = 0;
        if (do_restore) {
            hs = hash_seq_with_orientation(yb, yb_len, seq, seq_len, zb, zb_len, reversed);
            hq = hash_qual_with_orientation(yq, yq_len, qual, qual_len, zq, zq_len, reversed);
            total_bases = (long long)(yb_len + seq_len + zb_len);
        } else {
            hs = fnv1a64_update(hs, seq, seq_len);
            for (size_t i = 0; i < qual_len; ++i) hq = fnv1a64_update_1(hq, (unsigned char)qual_norm(qual[i]));
            total_bases = (long long)seq_len;
        }

        if (is_read1) {
            c1_seq ^= hs;
            c1_qual ^= hq;
            nrow1 += 1;
            nbases1 += total_bases;
        } else {
            c2_seq ^= hs;
            c2_qual ^= hq;
            nrow2 += 1;
            nbases2 += total_bases;
        }
    }

    if (need_pclose) pclose(fp);
    free(line);

    return Py_BuildValue("KKLLKKLL",
                         (unsigned long long)c1_seq, (unsigned long long)c1_qual, (long long)nrow1, (long long)nbases1,
                         (unsigned long long)c2_seq, (unsigned long long)c2_qual, (long long)nrow2, (long long)nbases2);
}

static PyMethodDef FastcheckMethods[] = {
    {"fastq_stats", py_fastq_stats, METH_VARARGS, "Compute paired FASTQ checksums and stats"},
    {"fastq_stats_interleaved", py_fastq_stats_interleaved, METH_VARARGS, "Compute paired FASTQ checksums and stats from a single interleaved FASTQ (file or '-')"},
    {"sam_stats",   py_sam_stats,   METH_VARARGS, "Compute checksums from SAM/BAM (via samtools view -h)"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef fastcheckmodule = {
    PyModuleDef_HEAD_INIT,
    .m_name = "fastcheck",
    .m_doc = NULL,
    .m_size = -1,
    .m_methods = FastcheckMethods,
    .m_slots = NULL,
    .m_traverse = NULL,
    .m_clear = NULL,
    .m_free = NULL
};

PyMODINIT_FUNC PyInit_fastcheck(void) {
    return PyModule_Create(&fastcheckmodule);
}
