#!/usr/bin/env python

import argparse
import glob
import os
import sys
import pybedtools

def compute_score_for_kit(kit_bed_file, sample_cov_bed, min_cov):
    """
    Computes the score for a given capture kit BED file.
    Score is defined as the fraction of the total bases in the kit intervals
    that have overlapping sample coverage with a mean coverage >= min_cov.
    """
    kit_intervals = pybedtools.BedTool(kit_bed_file)
    total_length = sum(int(interval.end) - int(interval.start) for interval in kit_intervals)
    if total_length == 0:
        return 0.0

    # Intersect kit intervals with the sample coverage intervals.
    # Assumes kit BED file has at least 3 columns (chrom, start, end)
    # and the sample coverage file (mosdepth .regions.bed) has 4 columns:
    # [chrom, start, end, mean_coverage]
    intersection = kit_intervals.intersect(sample_cov_bed, wa=True, wb=True)

    covered_length = 0
    for entry in intersection:
        # entry fields:
        # [0]: kit chrom, [1]: kit start, [2]: kit end,
        # [3]: sample chrom, [4]: sample start, [5]: sample end, [6]: mean coverage
        try:
            kit_start = int(entry[1])
            kit_end   = int(entry[2])
            sample_start = int(entry[4])
            sample_end   = int(entry[5])
            cov = float(entry[6])
        except (IndexError, ValueError):
            continue

        # Calculate the overlap between kit interval and sample interval
        overlap = min(kit_end, sample_end) - max(kit_start, sample_start)
        if overlap > 0 and cov >= min_cov:
            covered_length += overlap

    score = covered_length / total_length
    return score

def main():
    parser = argparse.ArgumentParser(
        description="Infer capture kit based on sample coverage and known kit BED files."
    )
    parser.add_argument("--coverage", required=True,
                        help="Path to sample coverage file (mosdepth .regions.bed)")
    parser.add_argument("--metadata_capture", required=True,
                        help="Capture kit specified in metadata")
    parser.add_argument("--output", required=True,
                        help="Output file to write inferred capture kit information")
    parser.add_argument("--kit_dir", default="capture_kits",
                        help="Directory containing capture kit BED files (default: capture_kits)")
    parser.add_argument("--min_cov", type=float, default=20,
                        help="Minimum coverage threshold (default: 20)")
    args = parser.parse_args()

    # Load the sample coverage file as a BedTool object
    try:
        sample_cov = pybedtools.BedTool(args.coverage)
    except Exception as e:
        sys.exit("Error loading coverage file: " + str(e))

    # Get list of capture kit BED files from the kit_dir
    kit_files = glob.glob(os.path.join(args.kit_dir, "*.bed"))
    if not kit_files:
        sys.exit("No kit BED files found in directory: " + args.kit_dir)

    kit_scores = {}
    for kit_file in kit_files:
        kit_name = os.path.basename(kit_file).rsplit(".bed", 1)[0]
        score = compute_score_for_kit(kit_file, sample_cov, args.min_cov)
        kit_scores[kit_name] = score

    # Determine the best matching kit (highest score)
    best_kit = max(kit_scores, key=kit_scores.get)
    best_score = kit_scores[best_kit]

    # Flag a discrepancy if the inferred kit and metadata kit do not match (case-insensitive)
    discrepancy = (best_kit.lower() != args.metadata_capture.lower())

    # Write results to the output file
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

    # Also print summary to stdout
    print("Inferred capture kit: {} (score: {:.3f})".format(best_kit, best_score))
    if discrepancy:
        print("Discrepancy detected: metadata_capture is '{}'.".format(args.metadata_capture))

if __name__ == "__main__":
    main()
