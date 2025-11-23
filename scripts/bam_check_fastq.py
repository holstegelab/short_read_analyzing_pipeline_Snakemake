import argparse
import sys
import csv
import gzip

try:
    import fastcheck
    _HAS_FASTCHECK = True
except ModuleNotFoundError:
    sys.stderr.write('fastcheck module not found; falling back to Python implementation\n')
    sys.stderr.flush()
    _HAS_FASTCHECK = False
except Exception:
    import traceback
    traceback.print_exc(file=sys.stderr)
    sys.stderr.write('fastcheck import failed; falling back to Python implementation\n')
    sys.stderr.flush()
    _HAS_FASTCHECK = False

import time
import subprocess
from bam_utils import *
from fastq_utils import *


def process(querygroup, stats, args):
    
    stats['alignment_counter'] = stats.get('alignment_counter',0) + len(querygroup)
    stats['fragment_counter'] = stats.get('fragment_counter',0) + 1

    reads1 = process_readgroup([row for row in querygroup if row.flag & 0x40])
    reads2 = process_readgroup([row for row in querygroup if row.flag & 0x80])

    read1 = reads1['primary']
    read2 = reads2['primary']

    return [read1, read2]

   

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", default='-', help='input sam file (default stdin)')
    parser.add_argument("--f1", help='first fastq file')
    parser.add_argument("--f2", help='second fastq file')
    parser.add_argument("-s", help='stats file')
    parser.add_argument("-c", help='check file')
    args = parser.parse_args()
    

    stats = {}
    alignment_counter = 0
    last_time = time.time()

    if _HAS_FASTCHECK:
        c1_seq, c1_qual, c2_seq, c2_qual, c_nrow, c_nbases1, c_nbases2 = fastcheck.fastq_stats(args.f1, args.f2)
        fastq1_checksum_seq = 0
        fastq1_checksum_qual = 0
        fastq1_nrow = 0
        fastq1_nbases = 0
        fastq2_checksum_seq = 0
        fastq2_checksum_qual = 0
        fastq2_nrow = 0
        fastq2_nbases = 0

        compare_fastq1_checksum_seq = c1_seq
        compare_fastq2_checksum_seq = c2_seq
        compare_fastq1_checksum_qual = c1_qual
        compare_fastq2_checksum_qual = c2_qual
        compare_fastq_nrow = c_nrow
        compare_fastq_nbases1 = c_nbases1
        compare_fastq_nbases2 = c_nbases2

        b1_seq, b1_qual, b1_nrow, b1_nbases, b2_seq, b2_qual, b2_nrow, b2_nbases = fastcheck.sam_stats(args.i)
        fastq1_checksum_seq = b1_seq
        fastq1_checksum_qual = b1_qual
        fastq1_nrow = b1_nrow
        fastq1_nbases = b1_nbases
        fastq2_checksum_seq = b2_seq
        fastq2_checksum_qual = b2_qual
        fastq2_nrow = b2_nrow
        fastq2_nbases = b2_nbases

        stats = {'fastq1_checksum_seq':fastq1_checksum_seq, 'fastq1_checksum_qual':fastq1_checksum_qual, 'fastq1_nrow':fastq1_nrow, 'fastq2_checksum_seq':fastq2_checksum_seq, 'fastq2_checksum_qual':fastq2_checksum_qual, 'fastq2_nrow':fastq2_nrow}
        stats.update({'fastq1_nbases':fastq1_nbases, 'fastq2_nbases':fastq2_nbases})
        stats.update({'compare_fastq1_checksum_seq':compare_fastq1_checksum_seq, 'compare_fastq1_checksum_qual':compare_fastq1_checksum_qual, 'compare_fastq2_checksum_seq':compare_fastq2_checksum_seq, 'compare_fastq2_checksum_qual':compare_fastq2_checksum_qual})
        stats.update({'compare_fastq_nrow':compare_fastq_nrow, 'compare_fastq_nbases1':compare_fastq_nbases1, 'compare_fastq_nbases2':compare_fastq_nbases2})
    else:
        freader = PairedFastQReaderSimple(args.f1, args.f2)
        if args.i == '-':
            pipe_in = sys.stdin
        else:
            in_process = subprocess.Popen(['samtools','view','-h', args.i],stdout=subprocess.PIPE, universal_newlines=True)
            pipe_in = in_process.stdout

        with pipe_in as f:
            reader = csv.reader(f,delimiter='\t', quoting=csv.QUOTE_NONE)

            lastname = ''
            fastq1_checksum_seq = 0
            fastq1_checksum_qual = 0
            fastq1_nrow = 0
            fastq1_nbases = 0
            fastq2_checksum_seq = 0
            fastq2_checksum_qual = 0
            fastq2_nrow = 0
            fastq2_nbases = 0

            compare_fastq1_checksum_seq = 0
            compare_fastq2_checksum_seq = 0
            compare_fastq1_checksum_qual = 0
            compare_fastq2_checksum_qual = 0
            compare_fastq_nrow = 0
            compare_fastq_nbases1 = 0
            compare_fastq_nbases2 = 0

            fastq_counter = 0
            for fread_name, fread_seq1, fread_seq2, fread_qual1, fread_qual2 in freader.retrieveRead():
                if (fastq_counter % 100000) == 0:
                    sys.stderr.write('%d fastq reads done at %d reads/sec\n' % (alignment_counter, int(100000.0 / (time.time() - last_time))))
                fastq_counter += 1
                fread_qual1 = fread_qual1.replace('#', '!')
                fread_qual2 = fread_qual2.replace('#', '!')
                compare_fastq1_checksum_seq ^= hash(fread_seq1)
                compare_fastq1_checksum_qual ^= hash(fread_qual1)
                compare_fastq2_checksum_seq ^= hash(fread_seq2)
                compare_fastq2_checksum_qual ^= hash(fread_qual2)
                compare_fastq_nrow += 1
                compare_fastq_nbases1 += len(fread_seq1)
                compare_fastq_nbases2 += len(fread_seq2)

            sys.stderr.write('End of fastq loop\n')
            sys.stderr.flush()
            for row in reader:
                alignment_counter += 1
                if (alignment_counter % 100000) == 0:
                    sys.stderr.write('%d reads done at %d reads/sec\n' % (alignment_counter, int(100000.0 / (time.time() - last_time))))
                    sys.stderr.flush()
                    last_time = time.time()
                if row[0][0] == '@':
                    continue
                row = BamRecord(row)
                if row.flag & 0x100 or row.flag & 0x800:
                    continue
                row, cseq, cqual = row.toFastqChecksum()
                if row.flag & 0x40:
                    fastq1_checksum_seq ^= cseq
                    fastq1_checksum_qual ^= cqual
                    fastq1_nrow += 1
                    fastq1_nbases += len(row.seq)
                else:
                    fastq2_checksum_seq ^= cseq
                    fastq2_checksum_qual ^= cqual
                    fastq2_nrow += 1
                    fastq2_nbases += len(row.seq)
            sys.stderr.write('End of BAM loop\n')
            sys.stderr.flush()

            stats = {'fastq1_checksum_seq':fastq1_checksum_seq, 'fastq1_checksum_qual':fastq1_checksum_qual, 'fastq1_nrow':fastq1_nrow, 'fastq2_checksum_seq':fastq2_checksum_seq, 'fastq2_checksum_qual':fastq2_checksum_qual, 'fastq2_nrow':fastq2_nrow}
            stats.update({'fastq1_nbases':fastq1_nbases, 'fastq2_nbases':fastq2_nbases})
            stats.update({'compare_fastq1_checksum_seq':compare_fastq1_checksum_seq, 'compare_fastq1_checksum_qual':compare_fastq1_checksum_qual, 'compare_fastq2_checksum_seq':compare_fastq2_checksum_seq, 'compare_fastq2_checksum_qual':compare_fastq2_checksum_qual})
            stats.update({'compare_fastq_nrow':compare_fastq_nrow, 'compare_fastq_nbases1':compare_fastq_nbases1, 'compare_fastq_nbases2':compare_fastq_nbases2})

    for k,v in stats.items():
        sys.stderr.write('%s\t%d\n' % (k,v))
    sys.stderr.flush()

    stats['seq1_compare'] = stats['fastq1_checksum_seq'] == stats['compare_fastq1_checksum_seq']
    stats['qual1_compare'] = stats['fastq1_checksum_qual'] == stats['compare_fastq1_checksum_qual']
    stats['seq2_compare'] = stats['fastq2_checksum_seq'] == stats['compare_fastq2_checksum_seq']
    stats['qual2_compare'] = stats['fastq2_checksum_qual'] == stats['compare_fastq2_checksum_qual']
    stats['nrow1_compare'] = stats['fastq1_nrow'] == stats['compare_fastq_nrow']
    stats['nrow2_compare'] = stats['fastq2_nrow'] == stats['compare_fastq_nrow']
    stats['compare_all'] = stats['seq1_compare'] and stats['qual1_compare'] and stats['seq2_compare'] and stats['qual2_compare'] and stats['nrow1_compare'] and stats['nrow2_compare']

    sys.stderr.write('Writing statistics\n')
    sys.stderr.flush()
    
    with open(args.s,'wt') as f:
        for key, value in stats.items():
            f.write('%s\t %d\n' % (key,value))
        f.flush()

    sys.stderr.write('Done\n')
    sys.stderr.flush()

    if stats['compare_all']:
        with open(args.c,'wt') as f:
            f.write('Match\n')

    if(not stats['compare_all']):
        sys.stdout.write('Checksums do not match\n')
    else:
        sys.stdout.write('Checksums match\n')
    sys.stdout.flush()
