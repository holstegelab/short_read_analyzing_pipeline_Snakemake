#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def load_summary(path):
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader)
    if len(rows) < 2:
        return None, None
    header = rows[0]
    record = rows[1]
    return header, record


def main():
    parser = argparse.ArgumentParser(description="Merge per-sample Kraken summaries into a single TSV")
    parser.add_argument("--summaries", help="Comma-separated list of kraken_summary.tsv files")
    parser.add_argument("--manifest", help="Path to a file with one kraken_summary.tsv path per line")
    parser.add_argument("--output", required=True, help="Output TSV path")
    args = parser.parse_args()

    if not args.summaries and not args.manifest:
        raise ValueError("Provide --summaries or --manifest")
    if args.summaries and args.manifest:
        raise ValueError("Provide either --summaries or --manifest, not both")

    if args.manifest:
        with open(args.manifest, "r") as handle:
            summary_paths = [Path(line.strip()) for line in handle if line.strip()]
    else:
        summary_paths = [Path(p) for p in args.summaries.split(",") if p]
    if not summary_paths:
        raise ValueError("No summary files provided")

    header = None
    rows = []

    for path in sorted(summary_paths):
        h, record = load_summary(path)
        if h is None:
            continue
        if header is None:
            header = h
        elif header != h:
            raise ValueError(f"Header mismatch in {path}")
        rows.append(record)

    if header is None:
        raise ValueError("No valid summary records found")

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


if __name__ == "__main__":
    main()
