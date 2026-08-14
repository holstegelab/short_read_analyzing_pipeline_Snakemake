#!/usr/bin/env python3
"""Leaf-only aggregation of existing short-read cohort statistics."""

import argparse
import csv
import math
import os
import sys
from pathlib import Path


PIPELINE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPELINE_ROOT))

import read_stats
import utils
from common import (
    DEEPVARIANT,
    KMER,
    KRAKEN,
    SAMPLEFILE_TO_SAMPLES,
    SAMPLEINFO,
    SAMPLEINFODIR,
    STAT,
    WINDOWS_ANNOTATED,
    level0_regions,
    level1_regions,
    pj,
)


def cohort_samples(samplefile):
    key = os.path.basename(samplefile)
    if key.endswith(".tsv"):
        key = key[:-4]
    if key not in SAMPLEFILE_TO_SAMPLES:
        raise KeyError(f"Unknown samplefile: {samplefile}")
    return key, sorted(SAMPLEFILE_TO_SAMPLES[key])


def require_existing(paths, label):
    missing = [
        str(path)
        for path in paths
        if not os.path.isfile(path) or os.path.getsize(path) == 0
    ]
    if missing:
        raise FileNotFoundError(
            f"{label}: {len(missing)} existing non-empty input(s) required; "
            f"first_missing={missing[:5]}"
        )


def sample_metadata(sample):
    info = SAMPLEINFO[sample]
    if "readgroups" in info:
        return info
    path = pj(SAMPLEINFODIR, f"{sample}.dat")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Readgroup metadata missing: {path}")
    stored = utils.load(path)
    result = info.copy()
    result["readgroups"] = stored["readgroups"]
    result["alternative_names"] = info.get("alternative_names", set()).union(
        stored.get("alternative_names", set())
    )
    return result


def aggregate_bam(samplefile, output):
    _, samples = cohort_samples(samplefile)

    def paths(pattern):
        return [pj(STAT, pattern.format(sample=sample)) for sample in samples]

    stats = paths("{sample}.samtools.stat")
    exome_stats = paths("{sample}.samtools.exome.stat")
    vpca2 = [pj(STAT, "contam", f"{sample}.verifybamid.pca2.selfSM") for sample in samples]
    bam_extra_all = paths("{sample}.bam_all.tsv")
    bam_extra_exome = paths("{sample}.bam_exome.tsv")
    pre_adapter = paths("{sample}.pre_adapter_summary_metrics")
    bait_bias = paths("{sample}.bait_bias_summary_metrics")
    hs_stats = paths("{sample}.hs_metrics")
    chrm_stats = paths("{sample}.chrM_read_stats.tsv")
    numt_stats = paths("{sample}.numt_read_stats.tsv")
    require_existing(
        stats
        + exome_stats
        + vpca2
        + bam_extra_all
        + bam_extra_exome
        + pre_adapter
        + bait_bias
        + hs_stats
        + chrm_stats
        + numt_stats,
        "bam-quality aggregation",
    )
    header, data = read_stats.combine_quality_stats(
        samples,
        stats,
        exome_stats,
        vpca2,
        bam_extra_all,
        bam_extra_exome,
        pre_adapter,
        bait_bias,
        hs_stats,
        chrm_stats,
        numt_stats,
    )
    read_stats.write_tsv(output, header, data)


def aggregate_bam_rg(samplefile, output):
    _, samples = cohort_samples(samplefile)
    sample_readgroups = []
    for sample in samples:
        for readgroup in sample_metadata(sample)["readgroups"]:
            sample_readgroups.append((sample, readgroup["info"]["ID"]))
    sample_readgroups.sort()

    aremoval = [pj(STAT, f"{sample}.{rg}.adapter_removal.log") for sample, rg in sample_readgroups]
    aidentify = [pj(STAT, f"{sample}.{rg}.fastq.adapters") for sample, rg in sample_readgroups]
    mergestats = [pj(STAT, f"{sample}.{rg}.merge_stats.tsv") for sample, rg in sample_readgroups]
    dragmap_stats = [pj(STAT, f"{sample}.{rg}.dragmap.log") for sample, rg in sample_readgroups]
    dechimer_stats = [pj(STAT, f"{sample}.{rg}.dechimer_stats.tsv") for sample, rg in sample_readgroups]
    require_existing(
        aremoval + aidentify + mergestats + dragmap_stats + dechimer_stats,
        "readgroup-quality aggregation",
    )
    header, data = read_stats.combine_rg_quality_stats(
        sample_readgroups,
        aremoval,
        aidentify,
        mergestats,
        dragmap_stats,
        dechimer_stats,
    )
    read_stats.write_tsv(output, header, data)


def aggregate_oxo(samplefile, output):
    _, samples = cohort_samples(samplefile)
    pre_adapter = [pj(STAT, f"{sample}.pre_adapter_detail_metrics") for sample in samples]
    bait_bias = [pj(STAT, f"{sample}.bait_bias_detail_metrics") for sample in samples]
    require_existing(pre_adapter + bait_bias, "oxo-quality aggregation")
    header, data = read_stats.combine_oxo_stats(samples, pre_adapter, bait_bias)
    read_stats.write_tsv(output, header, data)


def aggregate_sex(samplefile, output):
    key, samples = cohort_samples(samplefile)
    kmer_stats = [pj(KMER, f"{sample}.result.yaml") for sample in samples]
    require_existing(kmer_stats, "sex aggregation")
    reported = [SAMPLEFILE_TO_SAMPLES[key][sample]["sex"] for sample in samples]
    header, data = read_stats.combine_sex_stats(samples, kmer_stats, reported)
    read_stats.write_tsv(output, header, data)


def aggregate_coverage(samplefile, output):
    key, samples = cohort_samples(samplefile)
    bam_table = f"{key}.bam_quality.tab"
    mapped = [pj(STAT, f"{sample}.samtools.stat") for sample in samples]
    coverage = [pj(STAT, "cov", f"{sample}.regions.bed.gz") for sample in samples]
    require_existing([bam_table] + mapped + coverage, "coverage aggregation")
    read_stats.write_coverage_to_hdf5(
        WINDOWS_ANNOTATED,
        samples,
        mapped,
        coverage,
        output,
    )


def aggregate_kraken(samplefile, output):
    _, samples = cohort_samples(samplefile)
    summaries = [Path(pj(KRAKEN, f"{sample}.kraken_summary.tsv")) for sample in samples]
    require_existing(summaries, "Kraken aggregation")
    header = None
    rows = []
    for path in sorted(summaries):
        with path.open(newline="") as handle:
            records = list(csv.reader(handle, delimiter="\t"))
        if len(records) < 2:
            raise ValueError(f"No Kraken data row in {path}")
        if header is None:
            header = records[0]
        elif header != records[0]:
            raise ValueError(f"Kraken header mismatch in {path}")
        rows.append(records[1])
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def aggregate_phase(samplefile, output):
    _, samples = cohort_samples(samplefile)
    paths = [pj(STAT, f"{sample}.phase_stats.tsv") for sample in samples]
    require_existing(paths, "phase aggregation")
    header = None
    rows = []
    for path in paths:
        with open(path) as handle:
            reader = csv.reader(handle, delimiter="\t")
            current_header = next(reader, None)
            record = next(reader, None)
        if current_header is None or record is None:
            raise ValueError(f"No phase data row in {path}")
        if header is None:
            header = current_header
        elif header != current_header:
            raise ValueError(f"Phase header mismatch in {path}")
        rows.append(record)
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def aggregate_deepvariant(samplefile, output):
    _, samples = cohort_samples(samplefile)
    summaries = []
    gvcfs = []
    for sample in samples:
        if "wgs" in SAMPLEINFO[sample]["sample_type"]:
            regions = level1_regions
            for region in regions:
                summaries.append(
                    pj(STAT, "deepvariant_bcftools", f"{sample}.{region}.summary.tsv")
                )
                gvcfs.append(
                    pj(
                        DEEPVARIANT,
                        "gVCF",
                        "exome_extract",
                        region,
                        f"{sample}.{region}.wg.vcf.gz",
                    )
                )
        else:
            regions = level0_regions
            for region in regions:
                summaries.append(
                    pj(STAT, "deepvariant_bcftools", f"{sample}.{region}.summary.tsv")
                )
                gvcfs.append(
                    pj(DEEPVARIANT, "gVCF", region, f"{sample}.{region}.wg.vcf.gz")
                )
    require_existing(summaries, "DeepVariant summary aggregation")
    require_existing(gvcfs, "DeepVariant gVCF guard")

    aggregated = {}
    field_order = []

    def remember_field(field):
        if field not in field_order:
            field_order.append(field)

    def parse_number(value):
        if value is None:
            return None
        text = str(value).strip()
        if text == "" or text.lower() == "nan":
            return None
        try:
            if any(ch in text for ch in [".", "e", "E"]):
                return float(text)
            return int(text)
        except ValueError:
            return None

    def region_to_slice(region):
        if not region:
            return ""
        if region.startswith("A"):
            return "A"
        if region.startswith("F"):
            return "F"
        if region.startswith("X"):
            return "XH" if region.endswith("H") else "X"
        if region.startswith("Y"):
            return "YH" if region.endswith("H") else "Y"
        return region

    for path in sorted(summaries):
        basename = os.path.basename(path)
        parts = basename.split(".")
        sample_name = parts[0] if parts else ""
        region_name = parts[1] if len(parts) > 1 else ""
        with open(path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for field in reader.fieldnames or []:
                if field not in {"sample", "region"}:
                    remember_field(field)
            for record in reader:
                sample = record.get("sample", sample_name)
                region_tok = record.get("region", region_name)
                slice_name = region_to_slice(region_tok)
                agg = aggregated.setdefault(
                    (sample, slice_name), {"sums": {}, "strings": {}}
                )
                for field, value in record.items():
                    if field in {"sample", "region"}:
                        continue
                    remember_field(field)
                    numeric = parse_number(value)
                    if numeric is None:
                        if value not in ("", None):
                            agg["strings"].setdefault(field, value)
                        continue
                    if field == "number_of_samples":
                        previous = agg["sums"].get(field)
                        agg["sums"][field] = (
                            numeric if previous is None else max(previous, numeric)
                        )
                    else:
                        agg["sums"][field] = agg["sums"].get(field, 0.0) + numeric

    if not aggregated:
        raise ValueError("No DeepVariant summary records found")

    derived_fields = [
        "tstv_ratio",
        "tstv_ratio_first_alt",
        "indel_repeat_ratio",
        "indel_repeat_consistent_frac",
        "indel_repeat_inconsistent_frac",
        "avg_length_consistent_del",
        "avg_length_inconsistent_del",
        "avg_length_consistent_ins",
        "avg_length_inconsistent_ins",
    ]

    def format_number(value):
        if isinstance(value, float):
            if math.isfinite(value) and abs(value - round(value)) < 1e-9:
                return str(int(round(value)))
            return f"{value:.6f}"
        return str(value)

    def safe_ratio(numerator, denominator):
        if denominator in (None, 0) or numerator is None:
            return "NaN"
        return f"{numerator / denominator:.6f}"

    def safe_avg(total, count):
        if count in (None, 0) or total is None:
            return "NaN"
        return f"{total / count:.6f}"

    header = ["sample", "slice"] + field_order + derived_fields
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        for sample, slice_name in sorted(aggregated):
            sums = aggregated[(sample, slice_name)]["sums"]
            strings = aggregated[(sample, slice_name)]["strings"]
            total_indels = sums.get("number_of_indels")
            consistent = sums.get("indel_repeat_consistent")
            inconsistent = sums.get("indel_repeat_inconsistent")
            derived = {
                "tstv_ratio": safe_ratio(sums.get("tstv_ts"), sums.get("tstv_tv")),
                "tstv_ratio_first_alt": safe_ratio(
                    sums.get("tstv_ts_first_alt"), sums.get("tstv_tv_first_alt")
                ),
                "indel_repeat_ratio": safe_ratio(
                    consistent, (consistent or 0) + (inconsistent or 0)
                ),
                "indel_repeat_consistent_frac": safe_ratio(consistent, total_indels),
                "indel_repeat_inconsistent_frac": safe_ratio(
                    inconsistent, total_indels
                ),
            }
            for key, avg_key in [
                ("sum_length_consistent_del", "avg_length_consistent_del"),
                ("sum_length_inconsistent_del", "avg_length_inconsistent_del"),
                ("sum_length_consistent_ins", "avg_length_consistent_ins"),
                ("sum_length_inconsistent_ins", "avg_length_inconsistent_ins"),
            ]:
                derived[avg_key] = safe_avg(
                    sums.get(key), sums.get(f"{key}_count")
                )
            row = [sample, slice_name]
            for field in field_order:
                if field in sums:
                    row.append(format_number(sums[field]))
                elif field in strings:
                    row.append(strings[field])
                else:
                    row.append("")
            row.extend(derived.get(field, "NaN") for field in derived_fields)
            writer.writerow(row)


AGGREGATORS = {
    "bam": aggregate_bam,
    "bam-rg": aggregate_bam_rg,
    "oxo": aggregate_oxo,
    "sex": aggregate_sex,
    "coverage": aggregate_coverage,
    "kraken": aggregate_kraken,
    "deepvariant": aggregate_deepvariant,
    "phase": aggregate_phase,
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=sorted(AGGREGATORS), required=True)
    parser.add_argument("--samplefile", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    output = os.path.abspath(args.output)
    temporary = output + ".tmp"
    os.makedirs(os.path.dirname(output), exist_ok=True)
    try:
        if os.path.exists(temporary):
            os.unlink(temporary)
        AGGREGATORS[args.kind](args.samplefile, temporary)
        if not os.path.isfile(temporary) or os.path.getsize(temporary) == 0:
            raise ValueError(f"Aggregation produced no data: {temporary}")
        os.replace(temporary, output)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


if __name__ == "__main__":
    main()
