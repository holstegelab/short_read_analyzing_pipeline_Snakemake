#!/usr/bin/env python

import argparse
import json
import sys
import polars as pl
from concurrent.futures import ThreadPoolExecutor


def load_sample_coverage(coverage_file, min_cov):
    """
    Load the sample coverage file (e.g., from mosdepth) using Polars.
    Expects a tab-separated file with four columns: chrom, start, end, coverage.
    Only bins with coverage >= min_cov are retained.
    """
    try:
        df = pl.read_csv(
            coverage_file,
            separator="\t",
            has_header=False,
            new_columns=["chrom", "start", "end", "coverage"]
        )
    except Exception as e:
        sys.exit(f"Error reading sample coverage file {coverage_file}: {e}")

    # Cast columns to correct types and filter on coverage threshold
    df = df.with_columns([
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
        pl.col("coverage").cast(pl.Float64)
    ]).filter(pl.col("coverage") >= min_cov)

    return df


def compute_score_for_kit(kit_name, kit_data, sample_df):
    """
    For a given capture kit, create a DataFrame from the precomputed intervals
    (which are assumed to be the bins used for coverage). Then, perform an inner join
    with the sample coverage DataFrame (on chrom, start, end). The score is the fraction
    of expected bins that are covered in the sample.
    """
    try:
        kit_df = pl.DataFrame(kit_data["intervals"], schema=["chrom", "start", "end"])
    except Exception as e:
        return kit_name, 0.0

    # Join on all three columns; this assumes the bins in the sample file match those in the kit data.
    joined = kit_df.join(sample_df, on=["chrom", "start", "end"], how="inner")
    covered_bins = joined.height
    total_bins = kit_data.get("total_bins", 1)
    score = covered_bins / total_bins
    return kit_name, score


def infer_capture_kit(coverage_file, kit_json, metadata_capture, min_cov, output_file):
    # Load sample coverage data
    sample_df = load_sample_coverage(coverage_file, min_cov)

    # Load precomputed capture kit data
    try:
        with open(kit_json, "r") as f:
            precomputed_kits = json.load(f)
    except Exception as e:
        sys.exit(f"Error loading precomputed kit data from {kit_json}: {e}")

    kit_scores = {}
    # Process each kit in parallel
    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(compute_score_for_kit, kit_name, kit_data, sample_df)
            for kit_name, kit_data in precomputed_kits.items()
        ]
        for future in futures:
            kit_name, score = future.result()
            kit_scores[kit_name] = score

    best_kit = max(kit_scores, key=kit_scores.get)
    best_score = kit_scores[best_kit]
    discrepancy = (best_kit.lower() != metadata_capture.lower())

    # Write the results to the output file
    try:
        with open(output_file, "w") as out_f:
            out_f.write(f"metadata_capture: {metadata_capture}\n")
            out_f.write(f"inferred_capture_kit: {best_kit}\n")
            out_f.write(f"inferred_score: {best_score:.3f}\n")
            out_f.write(f"discrepancy: {'YES' if discrepancy else 'NO'}\n")
            out_f.write("\nAll kit scores:\n")
            for kit, score in sorted(kit_scores.items(), key=lambda x: x[1], reverse=True):
                out_f.write(f"{kit}: {score:.3f}\n")
    except Exception as e:
        sys.exit(f"Error writing output to {output_file}: {e}")

    print(f"Inferred capture kit: {best_kit} (score: {best_score:.3f})")
    if discrepancy:
        print(f"Discrepancy detected: metadata_capture is '{metadata_capture}'.")


def main():
    parser = argparse.ArgumentParser(
        description="Infer capture kit based on sample coverage and precomputed kit data."
    )
    parser.add_argument("--coverage", required=True, help="Sample coverage file (from mosdepth) with bins")
    parser.add_argument("--metadata_capture", required=True, help="Capture kit specified in metadata")
    parser.add_argument("--output", required=True, help="Output file for inferred capture kit results")
    parser.add_argument("--kit_json", default="capture_kits_precomputed.json", help="Precomputed kit JSON file")
    parser.add_argument("--min_cov", type=float, default=20, help="Minimum coverage threshold (default: 20)")
    args = parser.parse_args()

    infer_capture_kit(args.coverage, args.kit_json, args.metadata_capture, args.min_cov, args.output)


if __name__ == "__main__":
    main()