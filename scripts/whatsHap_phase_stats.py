#!/usr/bin/env python3
import argparse
import os
import re
import sys
import csv
import math
import yaml
import traceback

CHR_AUTOSOMES = {f"chr{i}" for i in range(1, 23)}

SECTION_RE = re.compile(r"^-+\s*Chromosome\s*(.+?)\s*-+\s*$")
NUM_RE = re.compile(r"([0-9]+)")
FLOAT_RE = re.compile(r"([0-9]+(?:\.[0-9]+)?)")

# Metrics we aggregate per region group
METRIC_KEYS = [
    "variants",
    "het",
    "phased",
    "blocks",
    "block_sizes_sum_variants",
    "block_lengths_sum_bp",
    "block_longest_bp",
    "block_shortest_bp",
]


def init_metrics():
    return {
        "variants": 0,
        "het": 0,
        "phased": 0,
        "blocks": 0,
        "block_sizes_sum_variants": 0,
        "block_lengths_sum_bp": 0,
        "block_longest_bp": 0,
        "block_shortest_bp": None,  # None -> will set to min positive seen
    }


def add_into(acc, part):
    for k in METRIC_KEYS:
        if k in {"block_longest_bp"}:
            acc[k] = max(acc.get(k, 0) or 0, part.get(k, 0) or 0)
        elif k in {"block_shortest_bp"}:
            v = part.get(k)
            if v in (None, 0):
                continue
            if acc.get(k) in (None, 0):
                acc[k] = v
            else:
                acc[k] = min(acc[k], v)
        else:
            acc[k] = (acc.get(k, 0) or 0) + (part.get(k, 0) or 0)


def parse_stats_file(path):
    # returns dict: chrom -> metrics
    current = None
    by_chrom = {}
    # temp state for which block section we're in
    in_block_sizes = False
    in_block_lengths = False
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            s = line.strip()
            m = SECTION_RE.match(s)
            if m:
                current = m.group(1).strip()
                by_chrom.setdefault(current, init_metrics())
                in_block_sizes = False
                in_block_lengths = False
                continue

            if s.startswith("Block sizes (no. of variants)"):
                in_block_sizes = True
                in_block_lengths = False
                continue
            if s.startswith("Block lengths (basepairs)"):
                in_block_lengths = True
                in_block_sizes = False
                continue

            if current is None:
                continue

            # Counts section
            if s.startswith("Variants in VCF:"):
                m = NUM_RE.search(s)
                if m:
                    by_chrom[current]["variants"] += int(m.group(1))
                continue
            if s.startswith("Heterozygous:"):
                m = NUM_RE.search(s)
                if m:
                    by_chrom[current]["het"] += int(m.group(1))
                continue
            if s.startswith("Phased:"):
                m = NUM_RE.search(s)
                if m:
                    by_chrom[current]["phased"] += int(m.group(1))
                continue
            if s.startswith("Blocks:"):
                m = NUM_RE.search(s)
                if m:
                    by_chrom[current]["blocks"] += int(m.group(1))
                continue

            # Block sizes
            if in_block_sizes:
                if s.startswith("Sum of sizes:"):
                    m = NUM_RE.search(s)
                    if m:
                        by_chrom[current]["block_sizes_sum_variants"] += int(m.group(1))
                continue

            # Block lengths
            if in_block_lengths:
                if s.startswith("Sum of lengths:"):
                    m = NUM_RE.search(s)
                    if m:
                        by_chrom[current]["block_lengths_sum_bp"] += int(m.group(1))
                    continue
                if s.startswith("Longest block:"):
                    m = NUM_RE.search(s)
                    if m:
                        val = int(m.group(1))
                        by_chrom[current]["block_longest_bp"] = max(by_chrom[current]["block_longest_bp"] or 0, val)
                    continue
                if s.startswith("Shortest block:"):
                    m = NUM_RE.search(s)
                    if m:
                        val = int(m.group(1))
                        prev = by_chrom[current]["block_shortest_bp"]
                        if prev in (None, 0) and val > 0:
                            by_chrom[current]["block_shortest_bp"] = val
                        elif val > 0:
                            by_chrom[current]["block_shortest_bp"] = min(prev, val)
                    continue

    return by_chrom


def region_group(chrom_name):
    chrom_name = chrom_name.strip()
    # Skip aggregated section; we compute 'all' ourselves by summing groups
    if chrom_name == "ALL chromosomes (aggregated)":
        return None
    if chrom_name == "chrM":
        return "chrm"
    if chrom_name in {"chrX"}:
        return "chrx"
    # tolerate names like chr22
    if chrom_name in CHR_AUTOSOMES:
        return "auto"
    # unknown -> ignore
    return None


def safe_div(n, d):
    if d in (None, 0):
        return float("nan")
    return n / d


def format_number(x):
    if x is None:
        return "NaN"
    if isinstance(x, float):
        if math.isnan(x) or math.isinf(x):
            return "NaN"
        return f"{x:.6f}"
    return str(x)


def main():
    ap = argparse.ArgumentParser(description="Aggregate WhatsHap phase stats across files and regions")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--inputs", nargs="+", help="List of WhatsHap stats files", required=True)
    ap.add_argument("--sex-yaml", default=None, help="Path to {KMER}/{sample}.result.yaml")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    if not args.inputs:
        print("ERROR: No WhatsHap stats inputs provided", file=sys.stderr)
        sys.exit(2)

    grouped = {k: init_metrics() for k in ("auto", "chrx", "chrm", "all")}

    matched = 0
    for path in args.inputs:
        if not os.path.exists(path):
            continue
        matched += 1
        by_chrom = parse_stats_file(path)
        for chrom, metrics in by_chrom.items():
            grp = region_group(chrom)
            if grp is None:
                continue
            # accumulate per-group
            add_into(grouped[grp], metrics)
            # compute overall on the fly by summing non-aggregated groups
            add_into(grouped["all"], metrics)

    # fill None shortest to 0
    for grp in grouped.values():
        if grp["block_shortest_bp"] is None:
            grp["block_shortest_bp"] = 0

    if matched == 0:
        print(f"WARNING: No WhatsHap stats files found for sample {args.sample}. Inputs provided: {len(args.inputs)}", file=sys.stderr)

    sex = ""
    if args.sex_yaml and os.path.exists(args.sex_yaml):
        try:
            with open(args.sex_yaml, "r") as yh:
                y = yaml.safe_load(yh) or {}
                sex = str(y.get("sex", ""))
        except Exception:
            print(f"ERROR: Failed to read sex YAML {args.sex_yaml}", file=sys.stderr)
            traceback.print_exc()
            sex = ""

    header = [
        "sample",
        "sex",
        # autosomes
        "phase_auto_variants",
        "phase_auto_het",
        "phase_auto_phased",
        "phase_auto_het_phased_pct",
        "phase_auto_blocks",
        "phase_auto_block_sizes_sum_variants",
        "phase_auto_block_lengths_sum_bp",
        "phase_auto_block_longest_bp",
        "phase_auto_block_shortest_bp",
        # chrX
        "phase_chrx_variants",
        "phase_chrx_het",
        "phase_chrx_phased",
        "phase_chrx_het_phased_pct",
        "phase_chrx_blocks",
        "phase_chrx_block_sizes_sum_variants",
        "phase_chrx_block_lengths_sum_bp",
        "phase_chrx_block_longest_bp",
        "phase_chrx_block_shortest_bp",
        # chrM
        "phase_chrm_variants",
        "phase_chrm_het",
        "phase_chrm_phased",
        "phase_chrm_het_phased_pct",
        "phase_chrm_blocks",
        "phase_chrm_block_sizes_sum_variants",
        "phase_chrm_block_lengths_sum_bp",
        "phase_chrm_block_longest_bp",
        "phase_chrm_block_shortest_bp",
        # all (auto + X + M)
        "phase_all_variants",
        "phase_all_het",
        "phase_all_phased",
        "phase_all_het_phased_pct",
        "phase_all_blocks",
        "phase_all_block_sizes_sum_variants",
        "phase_all_block_lengths_sum_bp",
        "phase_all_block_longest_bp",
        "phase_all_block_shortest_bp",
    ]

    def row_for(grp):
        g = grouped[grp]
        return [
            g["variants"],
            g["het"],
            g["phased"],
            safe_div(g["phased"], g["het"]) * 100.0 if g["het"] else float("nan"),
            g["blocks"],
            g["block_sizes_sum_variants"],
            g["block_lengths_sum_bp"],
            g["block_longest_bp"],
            g["block_shortest_bp"],
        ]

    row = [args.sample, sex] + row_for("auto") + row_for("chrx") + row_for("chrm") + row_for("all")

    with open(args.output, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(header)
        w.writerow([format_number(v) for v in row])


if __name__ == "__main__":
    main()
