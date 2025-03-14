import os
import json
import glob
import pybedtools
import polars as pl
import argparse

pj = os.path.join
RESOURCES = '/gpfs/work3/0/qtholstg/hg38_res_v2/'
INTERVALS_DIR = pj(RESOURCES,'intervals')
PRECOMPUTEED_BED = pj(INTERVALS_DIR, 'precomputed_kits.json')


def precompute_from_bed(bed_file):
    """
    Read a capture kit BED file and return the total number of bins and
    a list of intervals. The BED file is assumed to already contain the
    bins used for coverage (e.g. 100-bp chunks in exomes, larger bins elsewhere).
    """
    try:
        # Read BED file without header, expecting three columns: chrom, start, end
        df = pl.read_csv(
            bed_file,
            sep="\t",
            has_header=False,
            new_columns=["chrom", "start", "end"]
        )
    except Exception as e:
        raise Exception(f"Error reading {bed_file}: {e}")

    # Sort by chromosome and start coordinate (optional if your BED is already sorted)
    df = df.sort(["chrom", "start"])

    # Total number of bins equals the number of rows
    total_bins = df.height

    # Convert each row to a tuple: (chrom, start, end)
    intervals = [(row[0], int(row[1]), int(row[2])) for row in df.rows()]
    return total_bins, intervals

def main():
    parser = argparse.ArgumentParser(
        description="Precompute capture kit data from BED files used for coverage."
    )
    parser.add_argument("--kit_dir", default=INTERVALS_DIR,
                        help="Directory containing capture kit BED files (default: capture_kits)")
    parser.add_argument("--output", default=PRECOMPUTEED_BED,
                        help="Output JSON file for precomputed kit data")
    args = parser.parse_args()

    kit_files = glob.glob(os.path.join(args.kit_dir, "*.bed"))
    if not kit_files:
        raise Exception(f"No BED files found in directory {args.kit_dir}")

    precomputed = {}
    for bed_file in kit_files:
        kit_name = os.path.basename(bed_file).rsplit(".bed", 1)[0]
        total_bins, intervals = precompute_from_bed(bed_file)
        precomputed[kit_name] = {
            "total_bins": total_bins,
            "intervals": intervals
        }
        print(f"Precomputed {kit_name}: {total_bins} bins")

    with open(args.output, "w") as out_f:
        json.dump(precomputed, out_f, indent=4)
    print("Precomputed capture kit data saved to", args.output)

if __name__ == "__main__":
    main()

