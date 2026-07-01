#!/usr/bin/env python3
"""Summarize GATK gCNV postprocess outputs.

The default analysis uses genotyped segment VCFs because they are small and
contain final CNV calls. Interval VCF and denoised copy-ratio scans are optional
because those files are large for WGS runs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


SEGMENT_RE = re.compile(
    r"COHORT_(?P<cohort>\d+)_SAMPLE_(?P<sample>.+)_(?P<index>\d+)\.vcf\.gz$"
)


@dataclass
class SegmentCall:
    cohort: str
    sample: str
    sample_index: int
    contig: str
    start: int
    end: int
    length: int
    record_id: str
    alt: str
    filter: str
    qual: str
    svtype: str
    cn: str
    gt: str
    np: int | None
    qs: int | None
    qa: int | None
    qss: int | None
    qse: int | None
    passes_filter: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze GATK gCNV final outputs from a GATK_gCNV directory. "
            "Writes TSV summaries and a short text report."
        )
    )
    parser.add_argument(
        "--gatk-gcnv",
        type=Path,
        default=Path("GATK_gCNV"),
        help="Path to the GATK_gCNV output directory.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("gcnv_analysis"),
        help="Output directory for analysis tables.",
    )
    parser.add_argument(
        "--min-qs",
        type=int,
        default=50,
        help="Minimum segment QS for high-confidence calls.",
    )
    parser.add_argument(
        "--min-np",
        type=int,
        default=2,
        help="Minimum segment NP/bin count for high-confidence calls.",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=1000,
        help="Minimum segment length in bp for high-confidence calls.",
    )
    parser.add_argument(
        "--exclude-contigs",
        default="chrX,chrY,chrM",
        help="Comma-separated contigs to exclude from high-confidence summaries.",
    )
    parser.add_argument(
        "--min-cohort-size",
        type=int,
        default=30,
        help="Flag assigned cohorts smaller than this size.",
    )
    parser.add_argument(
        "--merge-gap",
        type=int,
        default=1000,
        help="Merge high-confidence CNV calls into recurrent loci when separated by at most this gap.",
    )
    parser.add_argument(
        "--min-recurrent-samples",
        type=int,
        default=2,
        help="Minimum unique samples required in recurrent_loci.tsv.",
    )
    parser.add_argument(
        "--include-intervals",
        action="store_true",
        help="Also scan genotyped interval VCFs. This can be slow and I/O heavy.",
    )
    parser.add_argument(
        "--max-interval-files",
        type=int,
        default=0,
        help="Limit interval VCF files scanned when --include-intervals is set; 0 means no limit.",
    )
    parser.add_argument(
        "--include-copy-ratio",
        action="store_true",
        help="Also scan denoised copy-ratio TSVs. This can be very slow for WGS.",
    )
    parser.add_argument(
        "--max-copy-ratio-files",
        type=int,
        default=0,
        help="Limit copy-ratio TSV files scanned when --include-copy-ratio is set; 0 means no limit.",
    )
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def int_or_none(value: str) -> int | None:
    if value in {"", "."}:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def float_or_none(value: str) -> float | None:
    if value in {"", "."}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_info(info_text: str) -> dict[str, str]:
    result = {}
    for item in info_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        else:
            key, value = item, ""
        result[key] = value
    return result


def parse_sample_format(format_text: str, sample_text: str) -> dict[str, str]:
    keys = format_text.split(":")
    values = sample_text.split(":")
    return dict(zip(keys, values))


def sample_from_segment_path(path: Path) -> tuple[str, str, int]:
    match = SEGMENT_RE.match(path.name)
    if not match:
        raise ValueError(f"Unexpected gCNV filename: {path}")
    return match.group("cohort"), match.group("sample"), int(match.group("index"))


def svtype_from_record(alt: str, cn: str, info: dict[str, str]) -> str:
    if "SVTYPE" in info:
        return info["SVTYPE"]
    cn_int = int_or_none(cn)
    if cn_int is not None:
        if cn_int < 2:
            return "DEL"
        if cn_int > 2:
            return "DUP"
    if alt in {"", "."}:
        return "REF"
    return alt.strip("<>")


def is_non_ref_call(call: SegmentCall) -> bool:
    cn_int = int_or_none(call.cn)
    if cn_int is not None and cn_int != 2:
        return True
    return call.alt not in {"", "."}


def is_pass_filter(filter_value: str) -> bool:
    return filter_value in {"", ".", "PASS"}


def passes_analysis_filter(call: SegmentCall, args: argparse.Namespace, excluded_contigs: set[str]) -> bool:
    if not is_non_ref_call(call):
        return False
    if call.contig in excluded_contigs:
        return False
    if not is_pass_filter(call.filter):
        return False
    if call.qs is None or call.qs < args.min_qs:
        return False
    if call.np is None or call.np < args.min_np:
        return False
    if call.length < args.min_length:
        return False
    return True


def iter_segment_calls(path: Path) -> Iterable[SegmentCall]:
    cohort, sample, sample_index = sample_from_segment_path(path)
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, record_id, _ref, alt, qual, filt, info_text, fmt, sample_text = fields[:10]
            start = int(pos)
            info = parse_info(info_text)
            end = int(info.get("END", start))
            values = parse_sample_format(fmt, sample_text)
            cn = values.get("CN", "")
            yield SegmentCall(
                cohort=cohort,
                sample=sample,
                sample_index=sample_index,
                contig=chrom,
                start=start,
                end=end,
                length=max(0, end - start + 1),
                record_id=record_id,
                alt=alt,
                filter=filt,
                qual=qual,
                svtype=svtype_from_record(alt, cn, info),
                cn=cn,
                gt=values.get("GT", ""),
                np=int_or_none(values.get("NP", "")),
                qs=int_or_none(values.get("QS", "")),
                qa=int_or_none(values.get("QA", "")),
                qss=int_or_none(values.get("QSS", "")),
                qse=int_or_none(values.get("QSE", "")),
                passes_filter=False,
            )


def discover_segment_vcfs(gatk_gcnv: Path) -> list[Path]:
    return sorted(gatk_gcnv.glob("GENOTYPED_CALLS_segments_*/*.vcf.gz"))


def discover_interval_vcfs(gatk_gcnv: Path, limit: int) -> list[Path]:
    paths = sorted(gatk_gcnv.glob("GENOTYPED_CALLS_intervals_*/*.vcf.gz"))
    return paths if limit <= 0 else paths[:limit]


def discover_copy_ratio_tsvs(gatk_gcnv: Path, limit: int) -> list[Path]:
    paths = sorted(gatk_gcnv.glob("GENOTYPED_CALLS_denoised_copy_ratio_*/*.tsv"))
    return paths if limit <= 0 else paths[:limit]


def read_cohort_report(gatk_gcnv: Path) -> tuple[dict[str, str], Counter, list[dict[str, str]]]:
    report_dir = gatk_gcnv / "cohort_report"
    sample_to_cohort: dict[str, str] = {}
    cohort_sizes: Counter = Counter()
    unassigned_rows: list[dict[str, str]] = []

    sample_files = sorted(report_dir.glob("*.samples_per_cohort.tsv"))
    if sample_files:
        with sample_files[0].open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sample_to_cohort[row["sample"]] = row["cohort"]
                cohort_sizes[row["cohort"]] += 1

    unassigned_files = sorted(report_dir.glob("*.unassigned_samples.tsv"))
    if unassigned_files:
        with unassigned_files[0].open(newline="") as handle:
            unassigned_rows = list(csv.DictReader(handle, delimiter="\t"))

    return sample_to_cohort, cohort_sizes, unassigned_rows


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def mean_or_blank(values: list[int | float]) -> str:
    if not values:
        return ""
    return f"{statistics.fmean(values):.4f}"


def median(values: list[int]) -> float:
    return float(statistics.median(values)) if values else 0.0


def mad(values: list[int], med: float) -> float:
    if not values:
        return 0.0
    return float(statistics.median([abs(value - med) for value in values]))


def summarize_segments(args: argparse.Namespace, excluded_contigs: set[str]):
    segment_paths = discover_segment_vcfs(args.gatk_gcnv)
    sample_rows: list[dict[str, object]] = []
    call_rows: list[dict[str, object]] = []
    high_conf_calls: list[SegmentCall] = []
    cohort_counts = Counter()

    for path in segment_paths:
        cohort, sample, sample_index = sample_from_segment_path(path)
        cohort_counts[cohort] += 1
        counters = Counter()
        bp = Counter()
        np_sum = Counter()
        qs_values: list[int] = []
        high_conf_qs_values: list[int] = []

        for call in iter_segment_calls(path):
            non_ref = is_non_ref_call(call)
            high_conf = passes_analysis_filter(call, args, excluded_contigs)
            call.passes_filter = high_conf

            counters["segment_records_total"] += 1
            counters[f"filter_{call.filter or '.'}"] += 1
            if call.contig in excluded_contigs:
                counters["excluded_contig_records"] += 1
            if non_ref:
                counters["cnv_segments_raw"] += 1
                counters[f"cnv_{call.svtype.lower()}_raw"] += 1
                bp["cnv_bp_raw"] += call.length
                if call.np is not None:
                    np_sum["cnv_np_raw"] += call.np
                if call.qs is not None:
                    qs_values.append(call.qs)

                call_rows.append(
                    {
                        "cohort": call.cohort,
                        "sample": call.sample,
                        "sample_index": call.sample_index,
                        "contig": call.contig,
                        "start": call.start,
                        "end": call.end,
                        "length": call.length,
                        "svtype": call.svtype,
                        "cn": call.cn,
                        "gt": call.gt,
                        "np": call.np if call.np is not None else "",
                        "qs": call.qs if call.qs is not None else "",
                        "qa": call.qa if call.qa is not None else "",
                        "qss": call.qss if call.qss is not None else "",
                        "qse": call.qse if call.qse is not None else "",
                        "filter": call.filter,
                        "qual": call.qual,
                        "record_id": call.record_id,
                        "passes_analysis_filter": int(high_conf),
                    }
                )

            if high_conf:
                high_conf_calls.append(call)
                counters["cnv_segments_high_conf"] += 1
                counters[f"cnv_{call.svtype.lower()}_high_conf"] += 1
                bp["cnv_bp_high_conf"] += call.length
                if call.np is not None:
                    np_sum["cnv_np_high_conf"] += call.np
                if call.qs is not None:
                    high_conf_qs_values.append(call.qs)

                cn_int = int_or_none(call.cn)
                if cn_int == 0:
                    counters["hom_del_high_conf"] += 1
                elif cn_int == 1:
                    counters["het_del_high_conf"] += 1
                elif cn_int is not None and cn_int > 2:
                    counters["dup_high_conf"] += 1

        sample_rows.append(
            {
                "cohort": cohort,
                "sample": sample,
                "sample_index": sample_index,
                "segment_records_total": counters["segment_records_total"],
                "cnv_segments_raw": counters["cnv_segments_raw"],
                "cnv_del_raw": counters["cnv_del_raw"],
                "cnv_dup_raw": counters["cnv_dup_raw"],
                "cnv_bp_raw": bp["cnv_bp_raw"],
                "cnv_np_raw": np_sum["cnv_np_raw"],
                "cnv_qs_mean_raw": mean_or_blank(qs_values),
                "cnv_segments_high_conf": counters["cnv_segments_high_conf"],
                "cnv_del_high_conf": counters["cnv_del_high_conf"],
                "cnv_dup_high_conf": counters["cnv_dup_high_conf"],
                "hom_del_high_conf": counters["hom_del_high_conf"],
                "het_del_high_conf": counters["het_del_high_conf"],
                "dup_high_conf": counters["dup_high_conf"],
                "cnv_bp_high_conf": bp["cnv_bp_high_conf"],
                "cnv_np_high_conf": np_sum["cnv_np_high_conf"],
                "cnv_qs_mean_high_conf": mean_or_blank(high_conf_qs_values),
                "excluded_contig_records": counters["excluded_contig_records"],
            }
        )

    return segment_paths, cohort_counts, sample_rows, call_rows, high_conf_calls


def contig_key(contig: str) -> tuple[int, object]:
    if contig.startswith("chr"):
        suffix = contig[3:]
        if suffix.isdigit():
            return (0, int(suffix))
        order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
        if suffix in order:
            return (0, order[suffix])
    return (1, contig)


def recurrent_loci(calls: list[SegmentCall], merge_gap: int, min_samples: int) -> list[dict[str, object]]:
    grouped = sorted(calls, key=lambda c: (contig_key(c.contig), c.svtype, c.start, c.end))
    loci: list[dict[str, object]] = []
    current: dict[str, object] | None = None

    def finish(locus: dict[str, object] | None) -> None:
        if not locus:
            return
        samples = locus.pop("_samples")
        cohorts = locus.pop("_cohorts")
        qs_values = locus.pop("_qs_values")
        np_values = locus.pop("_np_values")
        sample_count = len(samples)
        if sample_count < min_samples:
            return
        locus["sample_count"] = sample_count
        locus["samples"] = ",".join(sorted(samples))
        locus["cohorts"] = ",".join(sorted(cohorts))
        locus["mean_qs"] = mean_or_blank(qs_values)
        locus["mean_np"] = mean_or_blank(np_values)
        loci.append(locus)

    for call in grouped:
        if current is None:
            current = {
                "contig": call.contig,
                "start": call.start,
                "end": call.end,
                "svtype": call.svtype,
                "call_count": 1,
                "_samples": {call.sample},
                "_cohorts": {call.cohort},
                "_qs_values": [call.qs] if call.qs is not None else [],
                "_np_values": [call.np] if call.np is not None else [],
            }
            continue

        same_locus = (
            current["contig"] == call.contig
            and current["svtype"] == call.svtype
            and call.start <= int(current["end"]) + merge_gap
        )
        if not same_locus:
            finish(current)
            current = {
                "contig": call.contig,
                "start": call.start,
                "end": call.end,
                "svtype": call.svtype,
                "call_count": 1,
                "_samples": {call.sample},
                "_cohorts": {call.cohort},
                "_qs_values": [call.qs] if call.qs is not None else [],
                "_np_values": [call.np] if call.np is not None else [],
            }
            continue

        current["end"] = max(int(current["end"]), call.end)
        current["call_count"] = int(current["call_count"]) + 1
        current["_samples"].add(call.sample)
        current["_cohorts"].add(call.cohort)
        if call.qs is not None:
            current["_qs_values"].append(call.qs)
        if call.np is not None:
            current["_np_values"].append(call.np)

    finish(current)
    for locus in loci:
        locus["length"] = int(locus["end"]) - int(locus["start"]) + 1
    return sorted(loci, key=lambda row: (-int(row["sample_count"]), contig_key(str(row["contig"])), int(row["start"])))


def write_qc_flags(
    out_dir: Path,
    sample_rows: list[dict[str, object]],
    cohort_report_sizes: Counter,
    segment_cohort_counts: Counter,
    unassigned_rows: list[dict[str, str]],
    args: argparse.Namespace,
) -> list[dict[str, object]]:
    flags: list[dict[str, object]] = []
    burdens = [int(row["cnv_segments_high_conf"]) for row in sample_rows]
    med = median(burdens)
    burden_mad = mad(burdens, med)
    high_threshold = med + 6 * burden_mad if burden_mad > 0 else med * 2
    low_threshold = max(0, med - 6 * burden_mad) if burden_mad > 0 else 0

    for cohort, size in sorted(cohort_report_sizes.items()):
        if size < args.min_cohort_size:
            flags.append(
                {
                    "scope": "cohort",
                    "cohort": cohort,
                    "sample": "",
                    "flag": "TINY_COHORT",
                    "value": size,
                    "detail": f"Assigned cohort size is below {args.min_cohort_size}.",
                }
            )

    for cohort, size in sorted(segment_cohort_counts.items()):
        if cohort_report_sizes and size != cohort_report_sizes.get(cohort, 0):
            flags.append(
                {
                    "scope": "cohort",
                    "cohort": cohort,
                    "sample": "",
                    "flag": "SEGMENT_COUNT_MISMATCH",
                    "value": size,
                    "detail": f"Segment VCF sample count differs from cohort report count {cohort_report_sizes.get(cohort, 0)}.",
                }
            )

    for row in sample_rows:
        burden = int(row["cnv_segments_high_conf"])
        if burden > high_threshold:
            flags.append(
                {
                    "scope": "sample",
                    "cohort": row["cohort"],
                    "sample": row["sample"],
                    "flag": "HIGH_CNV_BURDEN",
                    "value": burden,
                    "detail": f"Above robust threshold {high_threshold:.2f}; median={med:.2f}, MAD={burden_mad:.2f}.",
                }
            )
        if burden < low_threshold:
            flags.append(
                {
                    "scope": "sample",
                    "cohort": row["cohort"],
                    "sample": row["sample"],
                    "flag": "LOW_CNV_BURDEN",
                    "value": burden,
                    "detail": f"Below robust threshold {low_threshold:.2f}; median={med:.2f}, MAD={burden_mad:.2f}.",
                }
            )

    if unassigned_rows:
        flags.append(
            {
                "scope": "run",
                "cohort": "",
                "sample": "",
                "flag": "UNASSIGNED_SAMPLES",
                "value": len(unassigned_rows),
                "detail": "Samples marked unassigned/noise in cohort report.",
            }
        )

    write_tsv(
        out_dir / "qc_flags.tsv",
        ["scope", "cohort", "sample", "flag", "value", "detail"],
        flags,
    )
    return flags


def summarize_interval_vcfs(paths: list[Path]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in paths:
        cohort, sample, sample_index = sample_from_segment_path(path)
        counters = Counter()
        cnq_values: list[int] = []
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    continue
                values = parse_sample_format(fields[8], fields[9])
                cn = values.get("CN", "")
                counters["intervals_total"] += 1
                counters[f"cn_{cn}"] += 1
                cnq = int_or_none(values.get("CNQ", ""))
                if cnq is not None:
                    cnq_values.append(cnq)
        rows.append(
            {
                "cohort": cohort,
                "sample": sample,
                "sample_index": sample_index,
                "intervals_total": counters["intervals_total"],
                "cn0_intervals": counters["cn_0"],
                "cn1_intervals": counters["cn_1"],
                "cn2_intervals": counters["cn_2"],
                "cn3_intervals": counters["cn_3"],
                "cn4_intervals": counters["cn_4"],
                "cnq_mean": mean_or_blank(cnq_values),
            }
        )
    return rows


def summarize_copy_ratios(paths: list[Path]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in paths:
        match = re.match(
            r"COHORT_(?P<cohort>\d+)_SAMPLE_(?P<sample>.+)_(?P<index>\d+)\.tsv$",
            path.name,
        )
        if not match:
            continue
        n = 0
        mean = 0.0
        m2 = 0.0
        low = 0
        high = 0
        with path.open() as handle:
            reader = csv.DictReader((line for line in handle if not line.startswith("@")), delimiter="\t")
            for row in reader:
                value = float_or_none(row.get("LINEAR_COPY_RATIO", ""))
                if value is None or not math.isfinite(value):
                    continue
                n += 1
                delta = value - mean
                mean += delta / n
                m2 += delta * (value - mean)
                if value < 1.5:
                    low += 1
                elif value > 2.5:
                    high += 1
        stdev = math.sqrt(m2 / (n - 1)) if n > 1 else 0.0
        rows.append(
            {
                "cohort": match.group("cohort"),
                "sample": match.group("sample"),
                "sample_index": match.group("index"),
                "intervals": n,
                "mean_linear_copy_ratio": f"{mean:.6f}" if n else "",
                "stdev_linear_copy_ratio": f"{stdev:.6f}" if n else "",
                "fraction_lt_1_5": f"{low / n:.6f}" if n else "",
                "fraction_gt_2_5": f"{high / n:.6f}" if n else "",
            }
        )
    return rows


def write_summary_text(
    out_dir: Path,
    args: argparse.Namespace,
    segment_paths: list[Path],
    cohort_report_sizes: Counter,
    segment_cohort_counts: Counter,
    unassigned_rows: list[dict[str, str]],
    sample_rows: list[dict[str, object]],
    flags: list[dict[str, object]],
    recurrent_rows: list[dict[str, object]],
) -> None:
    high_conf_counts = [int(row["cnv_segments_high_conf"]) for row in sample_rows]
    raw_counts = [int(row["cnv_segments_raw"]) for row in sample_rows]
    top_samples = sorted(sample_rows, key=lambda row: int(row["cnv_segments_high_conf"]), reverse=True)[:10]
    with (out_dir / "analysis_summary.txt").open("w") as handle:
        print(f"GATK_gCNV directory: {args.gatk_gcnv}", file=handle)
        print(f"Segment VCFs scanned: {len(segment_paths)}", file=handle)
        print(
            "High-confidence criteria: "
            f"FILTER in ./PASS, QS>={args.min_qs}, NP>={args.min_np}, "
            f"length>={args.min_length}, excluded_contigs={args.exclude_contigs}",
            file=handle,
        )
        print("", file=handle)
        print("Cohort sizes from cohort report:", file=handle)
        for cohort, count in sorted(cohort_report_sizes.items()):
            print(f"  cohort {cohort}: {count}", file=handle)
        print("Cohort sizes from segment VCFs:", file=handle)
        for cohort, count in sorted(segment_cohort_counts.items()):
            print(f"  cohort {cohort}: {count}", file=handle)
        print(f"Unassigned samples in report: {len(unassigned_rows)}", file=handle)
        print("", file=handle)
        if high_conf_counts:
            print(
                "High-confidence CNV segments per sample: "
                f"min={min(high_conf_counts)}, median={statistics.median(high_conf_counts)}, "
                f"max={max(high_conf_counts)}",
                file=handle,
            )
        if raw_counts:
            print(
                "Raw non-reference CNV segments per sample: "
                f"min={min(raw_counts)}, median={statistics.median(raw_counts)}, max={max(raw_counts)}",
                file=handle,
            )
        print("", file=handle)
        print("Top samples by high-confidence CNV burden:", file=handle)
        for row in top_samples:
            print(
                f"  cohort {row['cohort']} {row['sample']}: "
                f"{row['cnv_segments_high_conf']} high-confidence "
                f"({row['cnv_segments_raw']} raw)",
                file=handle,
            )
        print("", file=handle)
        print(f"QC flags: {len(flags)}", file=handle)
        for flag in flags[:20]:
            print(
                f"  {flag['flag']} {flag['cohort']} {flag['sample']}: "
                f"{flag['value']} {flag['detail']}",
                file=handle,
            )
        print("", file=handle)
        print(f"Recurrent loci with >= {args.min_recurrent_samples} samples: {len(recurrent_rows)}", file=handle)
        for row in recurrent_rows[:20]:
            print(
                f"  {row['contig']}:{row['start']}-{row['end']} {row['svtype']} "
                f"samples={row['sample_count']} calls={row['call_count']}",
                file=handle,
            )


def main() -> None:
    args = parse_args()
    args.gatk_gcnv = args.gatk_gcnv.resolve()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    excluded_contigs = {contig.strip() for contig in args.exclude_contigs.split(",") if contig.strip()}
    cohort_sample_map, cohort_report_sizes, unassigned_rows = read_cohort_report(args.gatk_gcnv)
    segment_paths, segment_cohort_counts, sample_rows, call_rows, high_conf_calls = summarize_segments(args, excluded_contigs)

    write_tsv(
        args.out_dir / "sample_segment_summary.tsv",
        [
            "cohort",
            "sample",
            "sample_index",
            "segment_records_total",
            "cnv_segments_raw",
            "cnv_del_raw",
            "cnv_dup_raw",
            "cnv_bp_raw",
            "cnv_np_raw",
            "cnv_qs_mean_raw",
            "cnv_segments_high_conf",
            "cnv_del_high_conf",
            "cnv_dup_high_conf",
            "hom_del_high_conf",
            "het_del_high_conf",
            "dup_high_conf",
            "cnv_bp_high_conf",
            "cnv_np_high_conf",
            "cnv_qs_mean_high_conf",
            "excluded_contig_records",
        ],
        sorted(sample_rows, key=lambda row: (str(row["cohort"]), int(row["sample_index"]))),
    )

    write_tsv(
        args.out_dir / "segment_calls.tsv",
        [
            "cohort",
            "sample",
            "sample_index",
            "contig",
            "start",
            "end",
            "length",
            "svtype",
            "cn",
            "gt",
            "np",
            "qs",
            "qa",
            "qss",
            "qse",
            "filter",
            "qual",
            "record_id",
            "passes_analysis_filter",
        ],
        call_rows,
    )

    recurrent_rows = recurrent_loci(high_conf_calls, args.merge_gap, args.min_recurrent_samples)
    write_tsv(
        args.out_dir / "recurrent_loci.tsv",
        [
            "contig",
            "start",
            "end",
            "length",
            "svtype",
            "sample_count",
            "call_count",
            "cohorts",
            "mean_qs",
            "mean_np",
            "samples",
        ],
        recurrent_rows,
    )

    cohort_rows = []
    all_cohorts = sorted(set(cohort_report_sizes) | set(segment_cohort_counts))
    for cohort in all_cohorts:
        cohort_rows.append(
            {
                "cohort": cohort,
                "assigned_samples_report": cohort_report_sizes.get(cohort, 0),
                "segment_vcf_samples": segment_cohort_counts.get(cohort, 0),
            }
        )
    write_tsv(
        args.out_dir / "cohort_summary.tsv",
        ["cohort", "assigned_samples_report", "segment_vcf_samples"],
        cohort_rows,
    )
    if unassigned_rows:
        write_tsv(
            args.out_dir / "unassigned_samples.tsv",
            ["samplefile", "sample", "reason"],
            unassigned_rows,
        )

    samples_with_segment_vcf = {str(row["sample"]) for row in sample_rows}
    missing_segment_rows = [
        {"sample": sample, "cohort": cohort}
        for sample, cohort in sorted(cohort_sample_map.items())
        if sample not in samples_with_segment_vcf
    ]
    if missing_segment_rows:
        write_tsv(args.out_dir / "missing_segment_outputs.tsv", ["sample", "cohort"], missing_segment_rows)

    flags = write_qc_flags(
        args.out_dir,
        sample_rows,
        cohort_report_sizes,
        segment_cohort_counts,
        unassigned_rows,
        args,
    )

    if args.include_intervals:
        interval_rows = summarize_interval_vcfs(discover_interval_vcfs(args.gatk_gcnv, args.max_interval_files))
        write_tsv(
            args.out_dir / "interval_summary.tsv",
            [
                "cohort",
                "sample",
                "sample_index",
                "intervals_total",
                "cn0_intervals",
                "cn1_intervals",
                "cn2_intervals",
                "cn3_intervals",
                "cn4_intervals",
                "cnq_mean",
            ],
            interval_rows,
        )

    if args.include_copy_ratio:
        copy_ratio_rows = summarize_copy_ratios(
            discover_copy_ratio_tsvs(args.gatk_gcnv, args.max_copy_ratio_files)
        )
        write_tsv(
            args.out_dir / "copy_ratio_summary.tsv",
            [
                "cohort",
                "sample",
                "sample_index",
                "intervals",
                "mean_linear_copy_ratio",
                "stdev_linear_copy_ratio",
                "fraction_lt_1_5",
                "fraction_gt_2_5",
            ],
            copy_ratio_rows,
        )

    write_summary_text(
        args.out_dir,
        args,
        segment_paths,
        cohort_report_sizes,
        segment_cohort_counts,
        unassigned_rows,
        sample_rows,
        flags,
        recurrent_rows,
    )

    print(f"Wrote gCNV analysis to {args.out_dir}")


if __name__ == "__main__":
    main()
