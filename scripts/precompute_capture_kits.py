#!/usr/bin/env python

import glob
import os
import json
import pybedtools
import "../constants"

def merge_and_compute_total(bed_file):
    # Load the kit intervals and merge overlapping intervals
    kit_bt = pybedtools.BedTool(bed_file)
    merged = kit_bt.merge()
    total_length = sum(int(interval.end) - int(interval.start) for interval in merged)
    return total_length, merged


def main():
    kit_dir = constants.INTERVALS_DIR
    output_file = constants.PRECOMPUTEED_BED
    kit_files = glob.glob(os.path.join(kit_dir, "*.bed"))

    precomputed = {}
    for kit_file in kit_files:
        kit_name = os.path.basename(kit_file).rsplit(".bed", 1)[0]
        total_length, merged = merge_and_compute_total(kit_file)
        precomputed[kit_name] = {
            "total_length": total_length,
            # Optionally, save the merged intervals as a list of tuples (chrom, start, end)
            "merged_intervals": [(iv.chrom, int(iv.start), int(iv.end)) for iv in merged]
        }
        print(f"Precomputed {kit_name}: {total_length} bases")

    with open(output_file, "w") as out_f:
        json.dump(precomputed, out_f, indent=4)
    print("Precomputed kit data saved to", output_file)


if __name__ == "__main__":
    main()
