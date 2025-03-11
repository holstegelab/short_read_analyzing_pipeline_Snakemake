import os
import json
import glob
import pybedtools

pj = os.path.join
RESOURCES = '/gpfs/work3/0/qtholstg/hg38_res_v2/'
INTERVALS_DIR = pj(RESOURCES,'intervals')
PRECOMPUTEED_BED = pj(INTERVALS_DIR, 'precomputed_kits.json')

# !/usr/bin/env python



def precompute_capture_kits(kit_dir = INTERVALS_DIR, output_file = PRECOMPUTEED_BED):
    """
    Precomputes merged target regions and total target base counts for each capture kit.
    """
    precomputed = {}

    for bed_path in glob.glob(os.path.join(kit_dir, "*.bed")):
        kit_name = os.path.basename(bed_path).replace(".bed", "")

        # Load and merge overlapping intervals
        kit_bed = pybedtools.BedTool(bed_path).merge()
        total_target_bases = sum(int(iv.end) - int(iv.start) for iv in kit_bed)

        # Store merged intervals as a list of tuples for fast lookups
        precomputed[kit_name] = {
            "total_length": total_target_bases,
            "merged_intervals": [(iv.chrom, int(iv.start), int(iv.end)) for iv in kit_bed]
        }
        print(f"Precomputed {kit_name}: {total_target_bases} bases")

    # Save precomputed data to JSON
    with open(output_file, "w") as f:
        json.dump(precomputed, f, indent=4)

    print(f"Precomputed kit data saved to {output_file}")


if __name__ == "__main__":
    precompute_capture_kits()
