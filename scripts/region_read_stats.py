#!/usr/bin/env python3
import argparse
import csv
import re
import subprocess
import sys
from typing import Dict, Optional


REWRITE = {
    '1st_fragments': 'first_fragments',
}


def _sanitize(name: str) -> str:
    name = name.lower()
    name = re.sub(r"[^a-z0-9]+", "_", name)
    name = re.sub(r"_+", "_", name)
    name = name.strip("_")
    return REWRITE.get(name, name)


def _run_samtools_stats(bam: str, region: Optional[str], bed: Optional[str], threads: int, exclude_flags: Optional[int]) -> Dict[str, str]:
    if region and bed:
        raise ValueError("Provide at most one of --region or --bed.")

    cmd = ["samtools", "stats", "-@", str(threads)]
    if exclude_flags is not None:
        cmd.extend(["-F", str(exclude_flags)])
    if bed:
        cmd.extend(["-t", bed])
    cmd.append(bam)
    if region:
        cmd.append(region)

    completed = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        sys.stderr.write(completed.stderr)
        sys.stderr.flush()
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")

    stats: Dict[str, str] = {}
    for line in completed.stdout.splitlines():
        if not line.startswith("SN"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        name = parts[1].strip().rstrip(":")
        value = parts[2].strip()
        stats[_sanitize(name)] = value
    return stats


def write_stats(stats: Dict[str, str], output: str) -> None:
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(["metric", "value"])
        for key in sorted(stats.keys()):
            writer.writerow([key, stats[key]])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute read statistics for a BAM region or BED.")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--region", help="Region string (e.g. chrM)")
    parser.add_argument("--bed", help="BED file with regions")
    parser.add_argument("--threads", type=int, default=1, help="Threads for samtools stats")
    parser.add_argument("--exclude-flags", type=int, default=None, help="Bitmask of SAM flags to exclude")
    parser.add_argument("--output", required=True, help="Output TSV path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    stats = _run_samtools_stats(args.bam, args.region, args.bed, args.threads, args.exclude_flags)
    write_stats(stats, args.output)


if __name__ == "__main__":
    main()
