import argparse
import sys
import csv
import subprocess
from typing import Tuple

from fastcheck_hts_loader import ensure_fastcheck_hts
from fastcheck_loader import ensure_fastcheck as ensure_fastcheck_sam


def _log(msg: str):
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


def load_stats(path: str) -> dict:
    res = {}
    with open(path, 'rt') as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if not parts:
                continue
            key = parts[0]
            try:
                val_str = parts[1].strip()
                val = int(val_str)
            except Exception:
                val = 0
            res[key] = val
    return res


def sam_stats_fallback(path: str) -> Tuple[int, int, int, int, int, int, int, int]:
    has_fc, fastcheck = ensure_fastcheck_sam(logger=_log)
    if has_fc and fastcheck is not None:
        return fastcheck.sam_stats(path)
    # ultimate fallback: stream SAM text via samtools and compute in Python (slower)
    # import locally to avoid cost if not needed
    from bam_utils import BamRecord

    if path == '-':
        pipe_in = sys.stdin
    else:
        in_process = subprocess.Popen(['samtools', 'view', '-h', path], stdout=subprocess.PIPE, universal_newlines=True)
        pipe_in = in_process.stdout

    def fnv1a64(data: str) -> int:
        FNV_OFFSET = 0xcbf29ce484222325
        FNV_PRIME = 0x100000001b3
        h = FNV_OFFSET
        for ch in data:
            h ^= ord(ch)
            h = (h * FNV_PRIME) & 0xFFFFFFFFFFFFFFFF
        return h

    def qual_seqaware(seq: str, qual: str) -> str:
        # Canonicalize quality at N-bases so BAM/CRAM and FASTQ remain comparable
        # when aligners normalize N quality to '#'.
        out = []
        for b, q in zip(seq, qual):
            if b in ('N', 'n'):
                out.append('!')
            else:
                out.append('!' if q == '#' else q)
        return ''.join(out)

    fastq1_checksum_seq = 0
    fastq1_checksum_qual = 0
    fastq1_nrow = 0
    fastq1_nbases = 0
    fastq2_checksum_seq = 0
    fastq2_checksum_qual = 0
    fastq2_nrow = 0
    fastq2_nbases = 0

    with pipe_in as f:
        reader = csv.reader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        for row in reader:
            if not row:
                continue
            if row[0][0] == '@':
                continue
            rec = BamRecord(row)
            if rec.flag & 0x100 or rec.flag & 0x800:
                continue
            rec = rec.unmap(True, orig_orientation=True)
            cseq = fnv1a64(rec.seq)
            cqual = fnv1a64(qual_seqaware(rec.seq, rec.qual))
            if rec.flag & 0x40:
                fastq1_checksum_seq ^= cseq
                fastq1_checksum_qual ^= cqual
                fastq1_nrow += 1
                fastq1_nbases += len(rec.seq)
            else:
                fastq2_checksum_seq ^= cseq
                fastq2_checksum_qual ^= cqual
                fastq2_nrow += 1
                fastq2_nbases += len(rec.seq)
    return (fastq1_checksum_seq, fastq1_checksum_qual, fastq1_nrow, fastq1_nbases,
            fastq2_checksum_seq, fastq2_checksum_qual, fastq2_nrow, fastq2_nbases)


def _debug_verbose_scan(path: str, expect_stats: dict) -> None:
    """Per-record debug scan on a saved BAM file.

    Computes both naive (non-N-aware, old behaviour) and N-aware quality checksums
    per record and compares them against the stored FASTQ-side 'expect' values.
    Logs the first 20 reads that have N bases whose quality required canonicalization.
    The two aggregate values at the end tell you immediately:
      - naive matches expect  → the FASTQ stats TSV is stale (regenerate it)
      - N-aware matches expect → bug is on FASTQ stats side only
      - neither matches       → something unrelated to N-quality normalization is wrong
    """
    import csv as _csv
    from bam_utils import BamRecord

    _FNV_OFFSET = 0xcbf29ce484222325
    _FNV_PRIME  = 0x100000001b3

    def _fnv1a64(data: str) -> int:
        h = _FNV_OFFSET
        for ch in data:
            h ^= ord(ch)
            h = (h * _FNV_PRIME) & 0xFFFFFFFFFFFFFFFF
        return h

    def _qual_naive(seq: str, qual: str) -> str:
        return qual.replace('#', '!')

    def _qual_seqaware(seq: str, qual: str) -> str:
        out = []
        for b, q in zip(seq, qual):
            if b in ('N', 'n'):
                out.append('!')
            else:
                out.append('!' if q == '#' else q)
        return ''.join(out)

    c1_naive = c1_aware = 0
    c2_naive = c2_aware = 0
    n1_canon = n2_canon = 0
    logged = 0
    MAX_LOG = 20

    _log('[DEBUG] starting per-record verbose scan ...')
    proc = subprocess.Popen(['samtools', 'view', '-h', path],
                            stdout=subprocess.PIPE, universal_newlines=True)
    with proc.stdout as f:
        reader = _csv.reader(f, delimiter='\t', quoting=_csv.QUOTE_NONE)
        for row in reader:
            if not row or row[0][0] == '@':
                continue
            rec = BamRecord(row)
            if rec.flag & 0x100 or rec.flag & 0x800:
                continue
            rec = rec.unmap(True, orig_orientation=True)

            q_naive = _qual_naive(rec.seq, rec.qual)
            q_aware = _qual_seqaware(rec.seq, rec.qual)
            is_r1 = bool(rec.flag & 0x40)

            if q_naive != q_aware:
                if is_r1:
                    n1_canon += 1
                else:
                    n2_canon += 1
                if logged < MAX_LOG:
                    n_pos = [
                        (i, rec.seq[i], rec.qual[i])
                        for i in range(len(rec.seq))
                        if rec.seq[i] in ('N', 'n') and rec.qual[i] not in ('!', '#')
                    ]
                    pos_str = ', '.join(
                        f"pos={i} base={b} qual={q}(phred={ord(q)-33})"
                        for i, b, q in n_pos[:5]
                    )
                    _log(f"  [DEBUG] {'R1' if is_r1 else 'R2'} {rec.qname}: "
                         f"{len(n_pos)} N-pos with non-trivial qual: {pos_str}")
                    logged += 1

            h_naive = _fnv1a64(q_naive)
            h_aware = _fnv1a64(q_aware)
            if is_r1:
                c1_naive ^= h_naive
                c1_aware ^= h_aware
            else:
                c2_naive ^= h_naive
                c2_aware ^= h_aware

    exp1 = expect_stats.get('compare_fastq1_checksum_qual', 0)
    exp2 = expect_stats.get('compare_fastq2_checksum_qual', 0)
    _log(f'[DEBUG] reads with N-canonicalization: R1={n1_canon}  R2={n2_canon}')
    _log(f'[DEBUG] R1 qual naive   = {hex(c1_naive)}  expect={hex(exp1)}  match={c1_naive == exp1}')
    _log(f'[DEBUG] R1 qual N-aware = {hex(c1_aware)}  expect={hex(exp1)}  match={c1_aware == exp1}')
    _log(f'[DEBUG] R2 qual naive   = {hex(c2_naive)}  expect={hex(exp2)}  match={c2_naive == exp2}')
    _log(f'[DEBUG] R2 qual N-aware = {hex(c2_aware)}  expect={hex(exp2)}  match={c2_aware == exp2}')
    if c1_naive == exp1 and c1_aware != exp1:
        _log('[DEBUG] HINT R1: naive matches expect → FASTQ stats TSV was built with old code; regenerate it')
    elif c2_naive == exp2 and c2_aware != exp2:
        _log('[DEBUG] HINT R2: naive matches expect → FASTQ stats TSV was built with old code; regenerate it')
    elif c1_aware != exp1 or c2_aware != exp2:
        _log('[DEBUG] HINT: neither naive nor N-aware matches expect for at least one readgroup; '
             'mismatch is NOT purely N-quality normalization — check for other quality differences')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', required=True, help='input BAM/CRAM path or - for stdin')
    ap.add_argument('--fastq-stats', required=True, help='precomputed FASTQ stats TSV from fastq_stats.py')
    ap.add_argument('-s', required=True, help='output stats TSV')
    ap.add_argument('-c', required=True, help='output check file (Match) when all comparisons hold')
    ap.add_argument('--threads', type=int, default=2, help='htslib BGZF threads (default: 2)')
    ap.add_argument('--reference', default='', help='FAI reference for CRAM decoding (optional)')
    ap.add_argument('--ignore-qual-checksum-diff', action='store_true',
                    help='Ignore read1/read2 quality checksum mismatches when deciding pass/fail')
    ap.add_argument('--debug', action='store_true',
                    help='When a qual checksum fails and input is a file path, run a verbose '
                         'per-record scan: logs reads with N-canonicalization applied, and computes '
                         'both naive and N-aware checksums to diagnose stale TSV vs algorithm bugs')
    args = ap.parse_args()

    stats = load_stats(args.fastq_stats)

    has_hts, fastcheck_hts = ensure_fastcheck_hts(logger=_log)
    if has_hts and fastcheck_hts is not None:
        try:
            (b1_seq, b1_qual, b1_nrow, b1_nbases,
             b2_seq, b2_qual, b2_nrow, b2_nbases) = fastcheck_hts.bam_stats(args.i, threads=args.threads, reference=args.reference)
        except Exception as e:
            _log(f"fastcheck_hts failed: {e}; falling back to SAM path")
            (b1_seq, b1_qual, b1_nrow, b1_nbases,
             b2_seq, b2_qual, b2_nrow, b2_nbases) = sam_stats_fallback(args.i)
    else:
        _log('fastcheck_hts unavailable; using SAM path fallback')
        (b1_seq, b1_qual, b1_nrow, b1_nbases,
         b2_seq, b2_qual, b2_nrow, b2_nbases) = sam_stats_fallback(args.i)

    stats.update({
        'fastq1_checksum_seq': b1_seq,
        'fastq1_checksum_qual': b1_qual,
        'fastq1_nrow': b1_nrow,
        'fastq1_nbases': b1_nbases,
        'fastq2_checksum_seq': b2_seq,
        'fastq2_checksum_qual': b2_qual,
        'fastq2_nrow': b2_nrow,
        'fastq2_nbases': b2_nbases,
    })

    stats['seq1_compare'] = stats['fastq1_checksum_seq'] == stats.get('compare_fastq1_checksum_seq', 0)
    stats['seq2_compare'] = stats['fastq2_checksum_seq'] == stats.get('compare_fastq2_checksum_seq', 0)
    stats['qual1_compare_raw'] = stats['fastq1_checksum_qual'] == stats.get('compare_fastq1_checksum_qual', 0)
    stats['qual2_compare_raw'] = stats['fastq2_checksum_qual'] == stats.get('compare_fastq2_checksum_qual', 0)
    stats['qual_checks_ignored'] = bool(args.ignore_qual_checksum_diff)
    stats['qual1_compare'] = True if args.ignore_qual_checksum_diff else stats['qual1_compare_raw']
    stats['qual2_compare'] = True if args.ignore_qual_checksum_diff else stats['qual2_compare_raw']
    stats['nrow1_compare'] = stats['fastq1_nrow'] == stats.get('compare_fastq_nrow', 0)
    stats['nrow2_compare'] = stats['fastq2_nrow'] == stats.get('compare_fastq_nrow', 0)
    stats['compare_all'] = stats['seq1_compare'] and stats['qual1_compare'] and stats['seq2_compare'] and stats['qual2_compare'] and stats['nrow1_compare'] and stats['nrow2_compare']

    with open(args.s, 'wt') as f:
        for k, v in stats.items():
            f.write(f"{k}\t {int(v)}\n")
        f.flush()

    if stats['compare_all']:
        with open(args.c, 'wt') as f:
            f.write('Match\n')

    def _hex(x: int) -> str:
        try:
            return hex(int(x))
        except Exception:
            return str(x)

    _log('BAM/CRAM (htslib) vs FASTQ comparison:')
    _log(f"  read1_nrow: bam={b1_nrow} expect={stats.get('compare_fastq_nrow', 0)} match={stats['nrow1_compare']}")
    _log(f"  read2_nrow: bam={b2_nrow} expect={stats.get('compare_fastq_nrow', 0)} match={stats['nrow2_compare']}")
    _log(f"  read1_nbases: bam={b1_nbases} expect={stats.get('compare_fastq_nbases1', 0)}")
    _log(f"  read2_nbases: bam={b2_nbases} expect={stats.get('compare_fastq_nbases2', 0)}")
    _log(f"  read1_seq_checksum: bam={_hex(b1_seq)} expect={_hex(stats.get('compare_fastq1_checksum_seq', 0))} match={stats['seq1_compare']}")
    if args.ignore_qual_checksum_diff:
        _log("  quality checksum policy: IGNORE mismatches (ERF mode)")
    _log(f"  read1_qual_checksum: bam={_hex(b1_qual)} expect={_hex(stats.get('compare_fastq1_checksum_qual', 0))} raw_match={stats['qual1_compare_raw']} effective_match={stats['qual1_compare']}")
    _log(f"  read2_seq_checksum: bam={_hex(b2_seq)} expect={_hex(stats.get('compare_fastq2_checksum_seq', 0))} match={stats['seq2_compare']}")
    _log(f"  read2_qual_checksum: bam={_hex(b2_qual)} expect={_hex(stats.get('compare_fastq2_checksum_qual', 0))} raw_match={stats['qual2_compare_raw']} effective_match={stats['qual2_compare']}")

    if not stats['compare_all']:
        sys.stdout.write('Checksums do not match\n')
        if args.debug:
            if args.i == '-':
                _log('[DEBUG] --debug requires a file path for -i (not stdin); '
                     'save the BAM to disk and re-run with the file path to get per-record diagnostics')
            else:
                _debug_verbose_scan(args.i, stats)
    else:
        sys.stdout.write('Checksums match\n')
    sys.stdout.flush()


if __name__ == '__main__':
    main()
