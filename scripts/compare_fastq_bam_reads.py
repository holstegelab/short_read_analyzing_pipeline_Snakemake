#!/usr/bin/env python3
"""
Compare FASTQ and BAM files to identify which reads are duplicated or extra in BAM.

This script helps debug rare mismatches between FASTQ input and BAM output
by identifying specific read names that differ.
"""

import argparse
import sys
import os
import subprocess
from collections import Counter


def parse_fastq_readnames(fastq_file):
    """Extract read names from a FASTQ file, counting occurrences."""
    readnames = Counter()
    
    cmd = ['pigz', '-dc', fastq_file] if fastq_file.endswith('.gz') else ['cat', fastq_file]
    
    with subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True
    ) as proc:
        line_num = 0
        for line in proc.stdout:
            if line_num % 4 == 0:  # Header line
                header = line.rstrip('\n')
                if not header.startswith('@'):
                    snippet = header[:200]
                    raise RuntimeError(
                        "FASTQ parse error: expected header line starting with '@' "
                        f"but got: {snippet}\n"
                        f"File: {fastq_file}\n"
                        f"At line: {line_num + 1}\n"
                        "This often indicates a truncated/corrupted .gz file or a decompression error. "
                        "Run: pigz -t <file.gz> (or gzip -t) to verify integrity."
                    )
                # Parse read name (strip @ and everything after space or /)
                name = header[1:]
                if ' ' in name:
                    name = name.split(' ')[0]
                if name.endswith('/1') or name.endswith('/2'):
                    name = name[:-2]
                readnames[name] += 1
            line_num += 1

        if line_num % 4 != 0:
            raise RuntimeError(
                f"FASTQ parse error: file ended mid-record (line count {line_num} not divisible by 4).\n"
                f"File: {fastq_file}\n"
                "This often indicates truncation/corruption."
            )

        rc = proc.wait()
        stderr_txt = (proc.stderr.read() or '').strip()
        if rc != 0:
            msg = (
                f"Decompression failed for {fastq_file} (exit code {rc}).\n"
                "Stderr (first 2000 chars):\n"
                f"{stderr_txt[:2000]}\n"
                "\nLikely cause: corrupted .gz stream (CRC mismatch) or truncated file. "
                "Verify with: pigz -t <file.gz> or gzip -t <file.gz>."
            )
            raise RuntimeError(msg)
    
    return readnames


def parse_bam_readnames(bam_file, read_number=None):
    """
    Extract read names from a BAM file, counting occurrences.
    Only counts primary alignments (excludes secondary 0x100 and supplementary 0x800).
    
    read_number: 1 for read1 (0x40), 2 for read2 (0x80), None for both
    """
    readnames = Counter()
    
    cmd = ['samtools', 'view', bam_file]
    
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True) as proc:
        for line in proc.stdout:
            if line.startswith('@'):
                continue
            fields = line.split('\t')
            qname = fields[0]
            flag = int(fields[1])
            
            # Skip secondary (0x100) and supplementary (0x800) alignments
            if flag & 0x100 or flag & 0x800:
                continue
            
            # Filter by read number if specified
            if read_number == 1 and not (flag & 0x40):
                continue
            if read_number == 2 and not (flag & 0x80):
                continue
            
            # Strip /1 or /2 suffix if present
            if qname.endswith('/1') or qname.endswith('/2'):
                qname = qname[:-2]
            
            readnames[qname] += 1
    
    return readnames


def compare_readnames(fastq_counts, bam_counts, label=""):
    """Compare read name counts between FASTQ and BAM."""
    
    all_names = set(fastq_counts.keys()) | set(bam_counts.keys())
    
    extra_in_bam = []      # In BAM but not in FASTQ
    missing_in_bam = []    # In FASTQ but not in BAM
    duplicated_in_bam = [] # More copies in BAM than FASTQ
    fewer_in_bam = []      # Fewer copies in BAM than FASTQ
    
    for name in all_names:
        fq_count = fastq_counts.get(name, 0)
        bam_count = bam_counts.get(name, 0)
        
        if fq_count == 0 and bam_count > 0:
            extra_in_bam.append((name, bam_count))
        elif fq_count > 0 and bam_count == 0:
            missing_in_bam.append((name, fq_count))
        elif bam_count > fq_count:
            duplicated_in_bam.append((name, fq_count, bam_count))
        elif bam_count < fq_count:
            fewer_in_bam.append((name, fq_count, bam_count))
    
    return {
        'extra_in_bam': extra_in_bam,
        'missing_in_bam': missing_in_bam,
        'duplicated_in_bam': duplicated_in_bam,
        'fewer_in_bam': fewer_in_bam,
        'label': label
    }


def print_results(results, output_file=None):
    """Print comparison results."""
    
    out = open(output_file, 'w') if output_file else sys.stdout
    
    label = results['label']
    prefix = f"[{label}] " if label else ""
    
    out.write(f"\n{'='*60}\n")
    out.write(f"{prefix}COMPARISON RESULTS\n")
    out.write(f"{'='*60}\n\n")
    
    # Extra in BAM (not in FASTQ)
    out.write(f"{prefix}EXTRA IN BAM (not in FASTQ): {len(results['extra_in_bam'])}\n")
    for name, count in results['extra_in_bam'][:50]:  # Limit output
        out.write(f"  {name}\t(appears {count}x in BAM)\n")
    if len(results['extra_in_bam']) > 50:
        out.write(f"  ... and {len(results['extra_in_bam']) - 50} more\n")
    out.write("\n")
    
    # Missing in BAM
    out.write(f"{prefix}MISSING IN BAM (in FASTQ but not BAM): {len(results['missing_in_bam'])}\n")
    for name, count in results['missing_in_bam'][:50]:
        out.write(f"  {name}\t(appears {count}x in FASTQ)\n")
    if len(results['missing_in_bam']) > 50:
        out.write(f"  ... and {len(results['missing_in_bam']) - 50} more\n")
    out.write("\n")
    
    # Duplicated in BAM
    out.write(f"{prefix}DUPLICATED IN BAM (more copies than FASTQ): {len(results['duplicated_in_bam'])}\n")
    for name, fq_count, bam_count in results['duplicated_in_bam'][:50]:
        out.write(f"  {name}\t(FASTQ: {fq_count}, BAM: {bam_count}, extra: {bam_count - fq_count})\n")
    if len(results['duplicated_in_bam']) > 50:
        out.write(f"  ... and {len(results['duplicated_in_bam']) - 50} more\n")
    out.write("\n")
    
    # Fewer in BAM
    out.write(f"{prefix}FEWER IN BAM (fewer copies than FASTQ): {len(results['fewer_in_bam'])}\n")
    for name, fq_count, bam_count in results['fewer_in_bam'][:50]:
        out.write(f"  {name}\t(FASTQ: {fq_count}, BAM: {bam_count}, missing: {fq_count - bam_count})\n")
    if len(results['fewer_in_bam']) > 50:
        out.write(f"  ... and {len(results['fewer_in_bam']) - 50} more\n")
    out.write("\n")
    
    # Summary
    total_extra = sum(c for _, c in results['extra_in_bam'])
    total_missing = sum(c for _, c in results['missing_in_bam'])
    total_dup_extra = sum(bam - fq for _, fq, bam in results['duplicated_in_bam'])
    total_fewer = sum(fq - bam for _, fq, bam in results['fewer_in_bam'])
    
    out.write(f"{prefix}SUMMARY:\n")
    out.write(f"  Total extra reads in BAM: {total_extra}\n")
    out.write(f"  Total missing reads in BAM: {total_missing}\n")
    out.write(f"  Total duplicated read instances in BAM: {total_dup_extra}\n")
    out.write(f"  Total fewer read instances in BAM: {total_fewer}\n")
    out.write(f"  Net difference (BAM - FASTQ): {total_extra - total_missing + total_dup_extra - total_fewer}\n")
    
    if output_file:
        out.close()


def main():
    parser = argparse.ArgumentParser(
        description='Compare FASTQ and BAM files to find extra/duplicated reads in BAM'
    )
    parser.add_argument('--f1', required=True, help='First FASTQ file (R1)')
    parser.add_argument('--f2', help='Second FASTQ file (R2), optional')
    parser.add_argument('-i', '--bam', required=True, help='BAM file to compare')
    parser.add_argument('-o', '--output', help='Output file (default: stdout)')
    parser.add_argument('--read1-only', action='store_true', 
                        help='Only compare read1 (useful for single-end or debugging)')
    parser.add_argument('--by-fragment', action='store_true',
                        help='Compare by fragment (read pair) instead of individual reads')
    
    args = parser.parse_args()
    
    sys.stderr.write("Parsing FASTQ file 1...\n")
    fastq1_counts = parse_fastq_readnames(args.f1)
    sys.stderr.write(f"  Found {sum(fastq1_counts.values())} reads, {len(fastq1_counts)} unique names\n")
    
    if args.f2 and not args.read1_only:
        sys.stderr.write("Parsing FASTQ file 2...\n")
        fastq2_counts = parse_fastq_readnames(args.f2)
        sys.stderr.write(f"  Found {sum(fastq2_counts.values())} reads, {len(fastq2_counts)} unique names\n")
    else:
        fastq2_counts = None
    
    if args.by_fragment:
        # Compare by fragment - use read names as fragments
        sys.stderr.write("Parsing BAM file (fragment mode)...\n")
        bam_counts = Counter()
        
        cmd = ['samtools', 'view', args.bam]
        seen_fragments = {}  # Track which reads we've seen for each fragment
        
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True) as proc:
            for line in proc.stdout:
                if line.startswith('@'):
                    continue
                fields = line.split('\t')
                qname = fields[0]
                flag = int(fields[1])
                
                if flag & 0x100 or flag & 0x800:
                    continue
                
                if qname.endswith('/1') or qname.endswith('/2'):
                    qname = qname[:-2]
                
                is_read1 = bool(flag & 0x40)
                key = (qname, is_read1)
                
                if key not in seen_fragments:
                    seen_fragments[key] = True
                    # Count as fragment when we see either read
                    if is_read1:
                        bam_counts[qname] += 1
        
        sys.stderr.write(f"  Found {sum(bam_counts.values())} fragments in BAM\n")
        
        results = compare_readnames(fastq1_counts, bam_counts, "fragments")
        print_results(results, args.output)
        
    else:
        # Compare by individual reads
        sys.stderr.write("Parsing BAM file (read1)...\n")
        bam1_counts = parse_bam_readnames(args.bam, read_number=1)
        sys.stderr.write(f"  Found {sum(bam1_counts.values())} read1 entries\n")
        
        results1 = compare_readnames(fastq1_counts, bam1_counts, "read1")
        print_results(results1, args.output)
        
        if fastq2_counts:
            sys.stderr.write("Parsing BAM file (read2)...\n")
            bam2_counts = parse_bam_readnames(args.bam, read_number=2)
            sys.stderr.write(f"  Found {sum(bam2_counts.values())} read2 entries\n")
            
            results2 = compare_readnames(fastq2_counts, bam2_counts, "read2")
            print_results(results2, args.output)


if __name__ == '__main__':
    main()
