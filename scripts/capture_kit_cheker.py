#!/usr/bin/env python

import argparse
import json
import pybedtools


def load_mosdepth_coverage(coverage_file):
    """
    Loads a mosdepth .regions.bed file and converts it into a dictionary
    with (chrom, start, end) -> mean coverage.
    """
    coverage_intervals = {}
    with open(coverage_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            chrom, start, end, mean_coverage = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
            coverage_intervals[(chrom, start, end)] = mean_coverage
    return coverage_intervals


def compute_kit_score(kit_data, coverage_intervals, min_cov):
    """
    Computes a score for a capture kit based on how well the sample coverage
    overlaps with the expected target regions.
    """
    total_target_bases = kit_data["total_length"]
    if total_target_bases == 0:
        return 0.0

    covered_bases = 0
    for chrom, start, end in kit_data["merged_intervals"]:
        for (cov_chrom, cov_start, cov_end), mean_cov in coverage_intervals.items():
            if chrom != cov_chrom:
                continue

            # Compute overlap
            overlap_start = max(start, cov_start)
            overlap_end = min(end, cov_end)
            overlap = max(0, overlap_end - overlap_start)

            if overlap > 0 and mean_cov >= min_cov:
                covered_bases += overlap

    return covered_bases / total_target_bases


def main():
    parser = argparse.ArgumentParser(description="Infer capture kit based on sample coverage")
    parser.add_argument("--coverage", required=True, help="Path to mosdepth .regions.bed file")
    parser.add_argument("--metadata_capture", required=True, help="Expected capture kit from metadata")
    parser.add_argument("--output", required=True, help="Output file to write inferred capture kit")
    parser.add_argument("--kit_data", default="capture_kits_precomputed.json",
                        help="Precomputed capture kit data (default: capture_kits_precomputed.json)")
    parser.add_argument("--min_cov", type=float, default=20, help="Minimum coverage threshold (default: 20)")
    args = parser.parse_args()

    # Load coverage data
    coverage_intervals = load_mosdepth_coverage(args.coverage)

    # Load precomputed capture kit data
    with open(args.kit_data, "r") as f:
        precomputed_kits = json.load(f)

    kit_scores = {}
    for kit_name, kit_data in precomputed_kits.items():
        score = compute_kit_score(kit_data, coverage_intervals, args.min_cov)
        kit_scores[kit_name] = score

    # Determine best-matching kit
    best_kit = max(kit_scores, key=kit_scores.get)
    best_score = kit_scores[best_kit]
    discrepancy = (best_kit.lower() != args.metadata_capture.lower())

    # Write results
    with open(args.output, "w") as out_f:
        out_f.write(f"metadata_capture: {args.metadata_capture}\n")
        out_f.write(f"inferred_capture_kit: {best_kit}\n")
        out_f.write(f"inferred_score: {best_score:.3f}\n")
        out_f.write(f"discrepancy: {'YES' if discrepancy else 'NO'}\n")
        out_f.write("\nAll kit scores:\n")
        for kit, score in sorted(kit_scores.items(), key=lambda x: x[1], reverse=True):
            out_f.write(f"{kit}: {score:.3f}\n")

    print(f"Inferred capture kit: {best_kit} (score: {best_score:.3f})")
    if discrepancy:
        print(f"Discrepancy detected: metadata_capture is '{args.metadata_capture}'.")


if __name__ == "__main__":
    main()
