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
            cqual = fnv1a64(rec.qual.replace('#', '!'))
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', required=True, help='input BAM/CRAM path or - for stdin')
    ap.add_argument('--fastq-stats', required=True, help='precomputed FASTQ stats TSV from fastq_stats.py')
    ap.add_argument('-s', required=True, help='output stats TSV')
    ap.add_argument('-c', required=True, help='output check file (Match) when all comparisons hold')
    ap.add_argument('--threads', type=int, default=2, help='htslib BGZF threads (default: 2)')
    ap.add_argument('--reference', default='', help='FAI reference for CRAM decoding (optional)')
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
    stats['qual1_compare'] = stats['fastq1_checksum_qual'] == stats.get('compare_fastq1_checksum_qual', 0)
    stats['seq2_compare'] = stats['fastq2_checksum_seq'] == stats.get('compare_fastq2_checksum_seq', 0)
    stats['qual2_compare'] = stats['fastq2_checksum_qual'] == stats.get('compare_fastq2_checksum_qual', 0)
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
    _log(f"  read1_qual_checksum: bam={_hex(b1_qual)} expect={_hex(stats.get('compare_fastq1_checksum_qual', 0))} match={stats['qual1_compare']}")
    _log(f"  read2_seq_checksum: bam={_hex(b2_seq)} expect={_hex(stats.get('compare_fastq2_checksum_seq', 0))} match={stats['seq2_compare']}")
    _log(f"  read2_qual_checksum: bam={_hex(b2_qual)} expect={_hex(stats.get('compare_fastq2_checksum_qual', 0))} match={stats['qual2_compare']}")

    if not stats['compare_all']:
        sys.stdout.write('Checksums do not match\n')
    else:
        sys.stdout.write('Checksums match\n')
    sys.stdout.flush()


if __name__ == '__main__':
    main()
