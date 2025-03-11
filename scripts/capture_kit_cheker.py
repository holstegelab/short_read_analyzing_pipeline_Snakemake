#!/usr/bin/env python

import argparse
import glob
import os
import sys
import json
import pybedtools


def compute_score_for_kit(kit_data, sample_cov_bed, min_cov):
    """
    Computes the score for a given capture kit.
    Score is defined as the fraction of the precomputed total targeted bases
    that have overlapping sample coverage with a mean coverage >= min_cov.
    """
    # Create a BedTool from precomputed merged intervals
    intervals = [pybedtools.Interval(chrom, start, end) for chrom, start, end in kit_data["merged_intervals"]]
    kit_bt = pybedtools.BedTool(intervals)

    total_length = kit_data["total_length"]
    if total_length == 0:
        return 0.0

    intersection = kit_bt.intersect(sample_cov_bed, wa=True, wb=True)
    covered_length = 0
    for entry in intersection:
        # entry: [kit_chrom, kit_start, kit_end, sample_chrom, sample_start, sample_end, mean_coverage]
        try:
            kit_start = int(entry[1])
            kit_end = int(entry[2])
            sample_start = int(entry[4])
            sample_end = int(entry[5])
            cov = float(entry[6])
        except (IndexError, ValueError):
            continue

        overlap = min(kit_end, sample_end) - max(kit_start, sample_start)
        if overlap > 0 and cov >= min_cov:
            covered_length += overlap

    score = covered_length / total_length
    return score


def main():
    parser = argparse.ArgumentParser(
        description="Infer capture kit based on sample coverage and precomputed kit data."
    )
    parser.add_argument("--coverage", required=True,
                        help="Path to sample coverage file (mosdepth .regions.bed)")
    parser.add_argument("--metadata_capture", required=True,
                        help="Capture kit specified in metadata")
    parser.add_argument("--output", required=True,
                        help="Output file to write inferred capture kit information")
    parser.add_argument("--kit_data", default="capture_kits_precomputed.json",
                        help="JSON file with precomputed capture kit data (default: capture_kits_precomputed.json)")
    parser.add_argument("--kit_dir", default="capture_kits",
                        help="Directory containing capture kit BED files (for backward compatibility)")
    parser.add_argument("--min_cov", type=float, default=20,
                        help="Minimum coverage threshold (default: 20)")
    args = parser.parse_args()

    try:
        sample_cov = pybedtools.BedTool(args.coverage)
    except Exception as e:
        sys.exit("Error loading coverage file: " + str(e))

    # Load precomputed kit data
    try:
        with open(args.kit_data, "r") as f:
            precomputed_kits = json.load(f)
    except Exception as e:
        sys.exit("Error loading precomputed kit data: " + str(e))

    kit_scores = {}
    for kit_name, kit_data in precomputed_kits.items():
        score = compute_score_for_kit(kit_data, sample_cov, args.min_cov)
        kit_scores[kit_name] = score

    best_kit = max(kit_scores, key=kit_scores.get)
    best_score = kit_scores[best_kit]
    discrepancy = (best_kit.lower() != args.metadata_capture.lower())

    try:
        with open(args.output, "w") as out_f:
            out_f.write("metadata_capture: {}\n".format(args.metadata_capture))
            out_f.write("inferred_capture_kit: {}\n".format(best_kit))
            out_f.write("inferred_score: {:.3f}\n".format(best_score))
            out_f.write("discrepancy: {}\n".format("YES" if discrepancy else "NO"))
            out_f.write("\nAll kit scores:\n")
            for kit, score in sorted(kit_scores.items(), key=lambda x: x[1], reverse=True):
                out_f.write("{}: {:.3f}\n".format(kit, score))
    except Exception as e:
        sys.exit("Error writing output: " + str(e))

    print("Inferred capture kit: {} (score: {:.3f})".format(best_kit, best_score))
    if discrepancy:
        print("Discrepancy detected: metadata_capture is '{}'.".format(args.metadata_capture))


if __name__ == "__main__":
    main()
