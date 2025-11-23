import argparse
import sys
import time
from fastcheck_loader import ensure_fastcheck

def _log(msg: str):
    sys.stderr.write(msg + '\n')
    sys.stderr.flush()

_HAS_FASTCHECK, fastcheck = ensure_fastcheck(logger=_log)

from fastq_utils import PairedFastQReaderSimple
from hash_utils import fnv1a64, normalize_qual

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--f1')
    ap.add_argument('--f2')
    ap.add_argument('--interleaved', action='store_true', help='read interleaved FASTQ from stdin or --input')
    ap.add_argument('--input', default='-', help='interleaved FASTQ path or - for stdin (only with --interleaved)')
    ap.add_argument('-s', required=True, help='output stats TSV')
    args = ap.parse_args()

    if args.interleaved:
        if _HAS_FASTCHECK:
            c1_seq, c1_qual, c2_seq, c2_qual, c_nrow, c_nbases1, c_nbases2 = fastcheck.fastq_stats_interleaved(args.input)
            stats = {
                'compare_fastq1_checksum_seq': c1_seq,
                'compare_fastq1_checksum_qual': c1_qual,
                'compare_fastq2_checksum_seq': c2_seq,
                'compare_fastq2_checksum_qual': c2_qual,
                'compare_fastq_nrow': c_nrow,
                'compare_fastq_nbases1': c_nbases1,
                'compare_fastq_nbases2': c_nbases2,
            }
        else:
            compare_fastq1_checksum_seq = 0
            compare_fastq2_checksum_seq = 0
            compare_fastq1_checksum_qual = 0
            compare_fastq2_checksum_qual = 0
            compare_fastq_nrow = 0
            compare_fastq_nbases1 = 0
            compare_fastq_nbases2 = 0
            f = sys.stdin if args.input == '-' else open(args.input, 'rt')
            try:
                while True:
                    h1 = f.readline()
                    if not h1:
                        break
                    s1 = f.readline().rstrip('\n')
                    p1 = f.readline()
                    q1 = f.readline().rstrip('\n')
                    h2 = f.readline()
                    s2 = f.readline().rstrip('\n')
                    p2 = f.readline()
                    q2 = f.readline().rstrip('\n')
                    if not s2:
                        break
                    q1 = normalize_qual(q1)
                    q2 = normalize_qual(q2)
                    compare_fastq1_checksum_seq ^= fnv1a64(s1)
                    compare_fastq2_checksum_seq ^= fnv1a64(s2)
                    compare_fastq1_checksum_qual ^= fnv1a64(q1)
                    compare_fastq2_checksum_qual ^= fnv1a64(q2)
                    compare_fastq_nrow += 1
                    compare_fastq_nbases1 += len(s1)
                    compare_fastq_nbases2 += len(s2)
            finally:
                if f is not sys.stdin:
                    f.close()
            stats = {
                'compare_fastq1_checksum_seq': compare_fastq1_checksum_seq,
                'compare_fastq1_checksum_qual': compare_fastq1_checksum_qual,
                'compare_fastq2_checksum_seq': compare_fastq2_checksum_seq,
                'compare_fastq2_checksum_qual': compare_fastq2_checksum_qual,
                'compare_fastq_nrow': compare_fastq_nrow,
                'compare_fastq_nbases1': compare_fastq_nbases1,
                'compare_fastq_nbases2': compare_fastq_nbases2,
            }
    else:
        assert args.f1 and args.f2, '--f1 and --f2 are required when not using --interleaved'
        if _HAS_FASTCHECK:
            c1_seq, c1_qual, c2_seq, c2_qual, c_nrow, c_nbases1, c_nbases2 = fastcheck.fastq_stats(args.f1, args.f2)
            stats = {
                'compare_fastq1_checksum_seq': c1_seq,
                'compare_fastq1_checksum_qual': c1_qual,
                'compare_fastq2_checksum_seq': c2_seq,
                'compare_fastq2_checksum_qual': c2_qual,
                'compare_fastq_nrow': c_nrow,
                'compare_fastq_nbases1': c_nbases1,
                'compare_fastq_nbases2': c_nbases2,
            }
        else:
            freader = PairedFastQReaderSimple(args.f1, args.f2)
            compare_fastq1_checksum_seq = 0
            compare_fastq2_checksum_seq = 0
            compare_fastq1_checksum_qual = 0
            compare_fastq2_checksum_qual = 0
            compare_fastq_nrow = 0
            compare_fastq_nbases1 = 0
            compare_fastq_nbases2 = 0
            for _, seq1, seq2, qual1, qual2 in freader.retrieveRead():
                qual1 = normalize_qual(qual1)
                qual2 = normalize_qual(qual2)
                compare_fastq1_checksum_seq ^= fnv1a64(seq1)
                compare_fastq2_checksum_seq ^= fnv1a64(seq2)
                compare_fastq1_checksum_qual ^= fnv1a64(qual1)
                compare_fastq2_checksum_qual ^= fnv1a64(qual2)
                compare_fastq_nrow += 1
                compare_fastq_nbases1 += len(seq1)
                compare_fastq_nbases2 += len(seq2)
            stats = {
                'compare_fastq1_checksum_seq': compare_fastq1_checksum_seq,
                'compare_fastq1_checksum_qual': compare_fastq1_checksum_qual,
                'compare_fastq2_checksum_seq': compare_fastq2_checksum_seq,
                'compare_fastq2_checksum_qual': compare_fastq2_checksum_qual,
                'compare_fastq_nrow': compare_fastq_nrow,
                'compare_fastq_nbases1': compare_fastq_nbases1,
                'compare_fastq_nbases2': compare_fastq_nbases2,
            }

    with open(args.s, 'wt') as f:
        for k, v in stats.items():
            f.write(f"{k}\t {int(v)}\n")
        f.flush()

    # minimal logging
    def _hex(x: int) -> str:
        try:
            return hex(int(x))
        except Exception:
            return str(x)
    if args.interleaved:
        _log(f"FASTQ interleaved stats: nrow={int(stats.get('compare_fastq_nrow',0))} nbases1={int(stats.get('compare_fastq_nbases1',0))} nbases2={int(stats.get('compare_fastq_nbases2',0))}")
        _log(f"  read1_seq_checksum={_hex(stats.get('compare_fastq1_checksum_seq',0))} read1_qual_checksum={_hex(stats.get('compare_fastq1_checksum_qual',0))}")
        _log(f"  read2_seq_checksum={_hex(stats.get('compare_fastq2_checksum_seq',0))} read2_qual_checksum={_hex(stats.get('compare_fastq2_checksum_qual',0))}")
    else:
        _log(f"FASTQ paired stats: nrow={int(stats.get('compare_fastq_nrow',0))} nbases1={int(stats.get('compare_fastq_nbases1',0))} nbases2={int(stats.get('compare_fastq_nbases2',0))}")
        _log(f"  read1_seq_checksum={_hex(stats.get('compare_fastq1_checksum_seq',0))} read1_qual_checksum={_hex(stats.get('compare_fastq1_checksum_qual',0))}")
        _log(f"  read2_seq_checksum={_hex(stats.get('compare_fastq2_checksum_seq',0))} read2_qual_checksum={_hex(stats.get('compare_fastq2_checksum_qual',0))}")

if __name__ == '__main__':
    main()
