import argparse
import sys
import csv
import gzip

import time
import subprocess
from bam_utils import *
      
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
    args = parser.parse_args()
    

    stats = {}
    alignment_counter = 0
    last_time = time.time()

    if args.i == '-':
        pipe_in = sys.stdin
    else:
        in_process = subprocess.Popen(['samtools','view','-h', args.i],stdout=subprocess.PIPE, universal_newlines=True)
        pipe_in = in_process.stdout

    with pipe_in as f:
        reader = csv.reader(f,delimiter='\t', quoting=csv.QUOTE_NONE)

        lastname = ''
        querygroup = []

        with gzip.open(args.f1,'w') as out_file1, gzip.open(args.f2,'w') as out_file2:
            for row in reader:
                alignment_counter += 1
                if (alignment_counter % 100000) == 0:
                    sys.stderr.write(str(stats) + '\n')
                    sys.stderr.write('%d reads done at %d reads/sec\n' % (alignment_counter, int(100000.0 / (time.time() - last_time))))

                    sys.stderr.flush()
                    last_time = time.time()
                if row[0][0] == '@':
                    continue
                row = BamRecord(row)
                if lastname == row.qname:
                    querygroup.append(row)
                else:
                    
                    if querygroup:
                        assert len(querygroup) >= 2, f'Reads {[e.qname for e in querygroup]} are not paired'
                        res_kept = process(querygroup, stats, args)
                        for xrow in res_kept:
                            if xrow.flag & 0x40:
                                out_file1.write(xrow.toFastqRecord().encode('utf-8'))
                            else:
                                out_file2.write(xrow.toFastqRecord().encode("utf-8"))
                        
                    lastname = row.qname
                    querygroup = [row]

            if querygroup:
                assert len(querygroup) >= 2, f'Reads {[e.qname for e in querygroup]} are not paired'
                res_kept = process(querygroup, stats, args)
                for xrow in res_kept:
                    if xrow.flag & 0x40:
                        out_file1.write(xrow.toFastqRecord().encode('utf-8'))
                    else:
                        out_file2.write(xrow.toFastqRecord().encode('utf-8'))
        
        sys.stderr.write('Writing statistics\n')
        sys.stderr.flush()
        
        with open(args.s,'wt') as f:            
            for key, value in stats.items():
                f.write('%s\t %d\n' % (key,value))
            f.flush()

        sys.stderr.write('Done\n')
        sys.stderr.flush()