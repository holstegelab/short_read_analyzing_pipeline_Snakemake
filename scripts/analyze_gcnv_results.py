#!/usr/bin/env python3
"""Summarize GATK gCNV postprocess outputs.

The default analysis uses genotyped segment VCFs because they are small and
contain final CNV calls. Interval VCF and denoised copy-ratio scans are optional
because those files are large for WGS runs.

The script supports both the original gCNV layout and the WGS autosomes layout
from gCNV_gatk_wgs_autosomes.smk, where postprocess outputs live in
GATK_gCNV/wgs_autosomes_all_samples.
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
WORKFLOW_SUBDIRS = {
    "standard": Path("."),
    "wgs-autosomes": Path("wgs_autosomes_all_samples"),
}


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
    svlen: int | None
    ac_orig: str
    af_orig: str
    an_orig: str
    passes_filter: bool


@dataclass
class AnalysisResult:
    workflow: str
    gatk_gcnv: Path
    out_dir: Path
    segment_paths: list[Path]
    cohort_report_sizes: Counter
    segment_cohort_counts: Counter
    unassigned_rows: list[dict[str, str]]
    sample_to_samplefile: dict[str, str]
    sample_rows: list[dict[str, object]]
    call_rows: list[dict[str, object]]
    raw_calls: list[SegmentCall]
    high_conf_calls: list[SegmentCall]
    recurrent_rows: list[dict[str, object]]
    flags: list[dict[str, object]]


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
        help=(
            "Path to a GATK_gCNV output directory. For WGS autosomes this can be "
            "either GATK_gCNV or GATK_gCNV/wgs_autosomes_all_samples."
        ),
    )
    parser.add_argument(
        "--workflow",
        choices=["auto", "standard", "wgs-autosomes", "both"],
        default="auto",
        help=(
            "Output layout to analyze. 'both' analyzes the standard layout and "
            "wgs_autosomes_all_samples, then writes comparison tables."
        ),
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
    parser.add_argument(
        "--annotation-bed",
        type=Path,
        action="append",
        default=[],
        help=(
            "Optional BED file with gene/region/CNV annotations. May be provided "
            "multiple times. Column 4 is used as the annotation name when present."
        ),
    )
    return parser.parse_args()


def has_gcnv_outputs(path: Path) -> bool:
    return (
        any(path.glob("GENOTYPED_CALLS_segments_*/*.vcf.gz"))
        or any(path.glob("GENOTYPED_CALLS_intervals_*/*.vcf.gz"))
        or (path / "cohort_report").is_dir()
    )


def normalize_gatk_gcnv_path(path: Path) -> Path:
    resolved = path.resolve()
    if has_gcnv_outputs(resolved):
        return resolved
    child = resolved / "GATK_gCNV"
    if has_gcnv_outputs(child) or child.exists():
        return child.resolve()
    return resolved


def resolve_gatk_gcnv_dir(path: Path, workflow: str) -> tuple[Path, str]:
    resolved = normalize_gatk_gcnv_path(path)
    autosomes_subdir = WORKFLOW_SUBDIRS["wgs-autosomes"]

    if workflow == "wgs-autosomes":
        if resolved.name == autosomes_subdir.name:
            return resolved, workflow
        return (resolved / autosomes_subdir).resolve(), workflow

    if workflow == "standard":
        return resolved, workflow

    direct_has_outputs = has_gcnv_outputs(resolved)
    autosomes_dir = resolved / autosomes_subdir
    autosomes_has_outputs = has_gcnv_outputs(autosomes_dir)

    if autosomes_has_outputs and not direct_has_outputs:
        return autosomes_dir.resolve(), "wgs-autosomes"
    if resolved.name == autosomes_subdir.name:
        return resolved, "wgs-autosomes"
    return resolved, "standard"


def resolve_run_dirs(path: Path, workflow: str) -> list[tuple[str, Path]]:
    if workflow != "both":
        run_dir, resolved_workflow = resolve_gatk_gcnv_dir(path, workflow)
        return [(resolved_workflow, run_dir)]

    base = normalize_gatk_gcnv_path(path)
    autosomes_subdir = WORKFLOW_SUBDIRS["wgs-autosomes"]
    if base.name == autosomes_subdir.name:
        standard_dir = base.parent
        autosomes_dir = base
    else:
        standard_dir = base
        autosomes_dir = base / autosomes_subdir
    return [
        ("standard", standard_dir.resolve()),
        ("wgs-autosomes", autosomes_dir.resolve()),
    ]


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


def int_from_list_or_none(value: str) -> int | None:
    for item in value.split(","):
        parsed = int_or_none(item)
        if parsed is not None:
            return parsed
    return None


def float_or_none(value: str) -> float | None:
    if value in {"", "."}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def float_values(value: str) -> list[float]:
    values = []
    for item in value.split(","):
        parsed = float_or_none(item)
        if parsed is not None and math.isfinite(parsed):
            values.append(parsed)
    return values


def int_values(value: str) -> list[int]:
    values = []
    for item in value.split(","):
        parsed = int_or_none(item)
        if parsed is not None:
            values.append(parsed)
    return values


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
                svlen=int_from_list_or_none(info.get("SVLEN", "")),
                ac_orig=info.get("AC_Orig", ""),
                af_orig=info.get("AF_Orig", ""),
                an_orig=info.get("AN_Orig", ""),
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


def read_contig_lengths(vcf_paths: list[Path]) -> dict[str, int]:
    contig_lengths: dict[str, int] = {}
    for path in vcf_paths:
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if line.startswith("#CHROM"):
                    break
                if not line.startswith("##contig=<"):
                    continue
                contig_match = re.search(r"ID=([^,>]+)", line)
                length_match = re.search(r"length=(\d+)", line)
                if contig_match and length_match:
                    contig_lengths.setdefault(contig_match.group(1), int(length_match.group(1)))
        if contig_lengths:
            break
    return contig_lengths


def add_metric_value(
    store: defaultdict[tuple[str, str], list[int | float]],
    scope: str,
    metric: str,
    value: int | float | None,
) -> None:
    if value is not None:
        store[(scope, metric)].append(value)


def metric_rows(
    base: dict[str, object],
    values_by_metric: dict[tuple[str, str], list[int | float]],
) -> list[dict[str, object]]:
    rows = []
    for (scope, metric), values in sorted(values_by_metric.items()):
        row = dict(base)
        row["scope"] = scope
        row["metric"] = metric
        row.update(numeric_summary(values))
        rows.append(row)
    return rows


def grouped_calls(calls: list[SegmentCall], keys: tuple[str, ...]) -> dict[tuple[object, ...], list[SegmentCall]]:
    groups: dict[tuple[object, ...], list[SegmentCall]] = defaultdict(list)
    for call in calls:
        key = []
        for item in keys:
            key.append(getattr(call, item))
        groups[tuple(key)].append(call)
    return groups


def summarize_svtypes(raw_calls: list[SegmentCall], high_conf_calls: list[SegmentCall]) -> list[dict[str, object]]:
    rows = []
    for scope, calls in [("raw", raw_calls), ("high_conf", high_conf_calls)]:
        for (cohort, svtype), group in sorted(grouped_calls(calls, ("cohort", "svtype")).items()):
            cn_counts = Counter(call.cn or "." for call in group)
            row = {
                "scope": scope,
                "cohort": cohort,
                "svtype": svtype,
                "call_count": len(group),
                "sample_count": len({call.sample for call in group}),
                "cnv_bp": sum(call.length for call in group),
                "cn0_calls": cn_counts["0"],
                "cn1_calls": cn_counts["1"],
                "cn2_calls": cn_counts["2"],
                "cn3_calls": cn_counts["3"],
                "cn4_calls": cn_counts["4"],
                "cn5_plus_calls": sum(count for cn, count in cn_counts.items() if int_or_none(cn) is not None and int(cn) >= 5),
            }
            row.update(prefixed_summary("length", [call.length for call in group]))
            row.update(prefixed_summary("qs", [call.qs for call in group if call.qs is not None]))
            rows.append(row)
    return rows


def summarize_contigs(
    raw_calls: list[SegmentCall],
    high_conf_calls: list[SegmentCall],
    contig_lengths: dict[str, int],
    excluded_contigs: set[str],
) -> list[dict[str, object]]:
    rows = []
    for scope, calls in [("raw", raw_calls), ("high_conf", high_conf_calls)]:
        for (cohort, contig, svtype), group in sorted(
            grouped_calls(calls, ("cohort", "contig", "svtype")).items(),
            key=lambda item: (str(item[0][0]), contig_key(str(item[0][1])), str(item[0][2])),
        ):
            contig_length = contig_lengths.get(str(contig), 0)
            mb = contig_length / 1_000_000 if contig_length else 0.0
            cnv_bp = sum(call.length for call in group)
            row = {
                "scope": scope,
                "cohort": cohort,
                "contig": contig,
                "svtype": svtype,
                "excluded_from_high_conf": int(str(contig) in excluded_contigs),
                "contig_length": contig_length or "",
                "call_count": len(group),
                "sample_count": len({call.sample for call in group}),
                "cnv_bp": cnv_bp,
                "call_rate_per_mb": format_number(len(group) / mb) if mb else "",
                "bp_per_mb": format_number(cnv_bp / mb) if mb else "",
                "mean_length": mean_or_blank([call.length for call in group]),
                "median_length": median_or_blank([call.length for call in group]),
                "mean_qs": mean_or_blank([call.qs for call in group if call.qs is not None]),
                "median_qs": median_or_blank([call.qs for call in group if call.qs is not None]),
                "mean_np": mean_or_blank([call.np for call in group if call.np is not None]),
                "median_np": median_or_blank([call.np for call in group if call.np is not None]),
            }
            rows.append(row)
    return rows


def summarize_size_distributions(raw_calls: list[SegmentCall], high_conf_calls: list[SegmentCall]) -> list[dict[str, object]]:
    rows = []
    for scope, calls in [("raw", raw_calls), ("high_conf", high_conf_calls)]:
        for cohort_label, cohort_calls in [("all", calls)] + sorted(grouped_calls(calls, ("cohort",)).items()):
            if isinstance(cohort_label, tuple):
                cohort = cohort_label[0]
                group_calls = cohort_calls
            else:
                cohort = cohort_label
                group_calls = cohort_calls
            for svtype, svtype_calls in sorted(grouped_calls(group_calls, ("svtype",)).items()):
                row = {
                    "scope": scope,
                    "cohort": cohort,
                    "svtype": svtype[0],
                }
                row.update(numeric_summary([call.length for call in svtype_calls]))
                rows.append(row)
    return rows


def summarize_allele_frequencies(raw_calls: list[SegmentCall], high_conf_calls: list[SegmentCall]) -> list[dict[str, object]]:
    rows = []
    bins = [
        ("af_lt_0_001", 0.0, 0.001),
        ("af_0_001_0_01", 0.001, 0.01),
        ("af_0_01_0_05", 0.01, 0.05),
        ("af_0_05_0_10", 0.05, 0.10),
        ("af_ge_0_10", 0.10, math.inf),
    ]
    for scope, calls in [("raw", raw_calls), ("high_conf", high_conf_calls)]:
        for (cohort, svtype), group in sorted(grouped_calls(calls, ("cohort", "svtype")).items()):
            af_values = [value for call in group for value in float_values(call.af_orig)]
            ac_values = [value for call in group for value in int_values(call.ac_orig)]
            an_values = [value for call in group for value in int_values(call.an_orig)]
            row = {
                "scope": scope,
                "cohort": cohort,
                "svtype": svtype,
                "call_count": len(group),
                "calls_with_af": len([call for call in group if float_values(call.af_orig)]),
                "calls_with_ac": len([call for call in group if int_values(call.ac_orig)]),
                "calls_with_an": len([call for call in group if int_values(call.an_orig)]),
            }
            for name, lower, upper in bins:
                row[name] = sum(1 for value in af_values if lower <= value < upper)
            row.update(prefixed_summary("af_orig", af_values))
            row.update(prefixed_summary("ac_orig", ac_values))
            row.update(prefixed_summary("an_orig", an_values))
            rows.append(row)
    return rows


def summarize_group_burdens(
    sample_rows: list[dict[str, object]],
    sample_to_samplefile: dict[str, str],
) -> list[dict[str, object]]:
    metrics = [
        "cnv_segments_raw",
        "cnv_del_raw",
        "cnv_dup_raw",
        "cnv_bp_raw",
        "cnv_segments_high_conf",
        "cnv_del_high_conf",
        "cnv_dup_high_conf",
        "cnv_bp_high_conf",
    ]
    grouped_rows: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    grouped_rows[("run", "all")] = list(sample_rows)
    for row in sample_rows:
        grouped_rows[("cohort", str(row["cohort"]))].append(row)
        samplefile = sample_to_samplefile.get(str(row["sample"]), "")
        if samplefile:
            grouped_rows[("samplefile", samplefile)].append(row)

    rows = []
    for (group_type, group_value), rows_in_group in sorted(grouped_rows.items()):
        for metric in metrics:
            row = {
                "group_type": group_type,
                "group_value": group_value,
                "sample_count": len(rows_in_group),
                "metric": metric,
            }
            row.update(numeric_summary(int(row_value.get(metric, 0)) for row_value in rows_in_group))
            rows.append(row)
    return rows


def read_cohort_report(gatk_gcnv: Path) -> tuple[dict[str, str], dict[str, str], Counter, list[dict[str, str]]]:
    report_dir = gatk_gcnv / "cohort_report"
    sample_to_cohort: dict[str, str] = {}
    sample_to_samplefile: dict[str, str] = {}
    cohort_sizes: Counter = Counter()
    unassigned_rows: list[dict[str, str]] = []

    sample_files = sorted(report_dir.glob("*.samples_per_cohort.tsv"))
    for sample_file in sample_files:
        with sample_file.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sample = row.get("sample", "")
                cohort = row.get("cohort", "")
                if not sample or not cohort:
                    continue
                if sample in sample_to_cohort:
                    continue
                sample_to_cohort[sample] = cohort
                sample_to_samplefile[sample] = row.get("samplefile", "")
                cohort_sizes[cohort] += 1

    unassigned_files = sorted(report_dir.glob("*.unassigned_samples.tsv"))
    for unassigned_file in unassigned_files:
        with unassigned_file.open(newline="") as handle:
            unassigned_rows.extend(csv.DictReader(handle, delimiter="\t"))

    return sample_to_cohort, sample_to_samplefile, cohort_sizes, unassigned_rows


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


STAT_COLUMNS = [
    "count",
    "min",
    "q10",
    "q25",
    "median",
    "q75",
    "q90",
    "max",
    "mean",
    "stdev",
    "mad",
]


def format_number(value: int | float | None) -> str:
    if value is None:
        return ""
    if not math.isfinite(float(value)):
        return ""
    value_float = float(value)
    if value_float.is_integer():
        return str(int(value_float))
    return f"{value_float:.4f}"


def percentile(sorted_values: list[float], fraction: float) -> float:
    if not sorted_values:
        return 0.0
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = fraction * (len(sorted_values) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return sorted_values[lower]
    weight = position - lower
    return sorted_values[lower] * (1.0 - weight) + sorted_values[upper] * weight


def numeric_summary(values: Iterable[int | float]) -> dict[str, object]:
    numeric_values = []
    for value in values:
        try:
            parsed = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(parsed):
            numeric_values.append(parsed)
    if not numeric_values:
        return {column: "" for column in STAT_COLUMNS}

    sorted_values = sorted(numeric_values)
    med = percentile(sorted_values, 0.5)
    deviations = sorted(abs(value - med) for value in sorted_values)
    stdev = statistics.stdev(sorted_values) if len(sorted_values) > 1 else 0.0
    return {
        "count": len(sorted_values),
        "min": format_number(sorted_values[0]),
        "q10": format_number(percentile(sorted_values, 0.10)),
        "q25": format_number(percentile(sorted_values, 0.25)),
        "median": format_number(med),
        "q75": format_number(percentile(sorted_values, 0.75)),
        "q90": format_number(percentile(sorted_values, 0.90)),
        "max": format_number(sorted_values[-1]),
        "mean": format_number(statistics.fmean(sorted_values)),
        "stdev": format_number(stdev),
        "mad": format_number(percentile(deviations, 0.5)),
    }


def prefixed_summary(prefix: str, values: Iterable[int | float]) -> dict[str, object]:
    return {f"{prefix}_{key}": value for key, value in numeric_summary(values).items()}


def median_or_blank(values: list[int | float]) -> str:
    if not values:
        return ""
    return format_number(statistics.median(values))


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
    raw_calls: list[SegmentCall] = []
    high_conf_calls: list[SegmentCall] = []
    filter_rows: list[dict[str, object]] = []
    quality_rows: list[dict[str, object]] = []
    cohort_counts = Counter()

    for path in segment_paths:
        cohort, sample, sample_index = sample_from_segment_path(path)
        cohort_counts[cohort] += 1
        counters = Counter()
        bp = Counter()
        np_sum = Counter()
        qs_values: list[int] = []
        high_conf_qs_values: list[int] = []
        filter_counts = Counter()
        raw_filter_counts = Counter()
        high_conf_filter_counts = Counter()
        sample_metric_values: defaultdict[tuple[str, str], list[int | float]] = defaultdict(list)

        for call in iter_segment_calls(path):
            non_ref = is_non_ref_call(call)
            high_conf = passes_analysis_filter(call, args, excluded_contigs)
            call.passes_filter = high_conf

            counters["segment_records_total"] += 1
            counters[f"filter_{call.filter or '.'}"] += 1
            filter_value = call.filter or "."
            filter_counts[filter_value] += 1
            if call.contig in excluded_contigs:
                counters["excluded_contig_records"] += 1
            if non_ref:
                raw_calls.append(call)
                raw_filter_counts[filter_value] += 1
                counters["cnv_segments_raw"] += 1
                counters[f"cnv_{call.svtype.lower()}_raw"] += 1
                bp["cnv_bp_raw"] += call.length
                if call.np is not None:
                    np_sum["cnv_np_raw"] += call.np
                if call.qs is not None:
                    qs_values.append(call.qs)
                add_metric_value(sample_metric_values, "raw", "length", call.length)
                add_metric_value(sample_metric_values, "raw", "np", call.np)
                add_metric_value(sample_metric_values, "raw", "qs", call.qs)
                add_metric_value(sample_metric_values, "raw", "qa", call.qa)
                add_metric_value(sample_metric_values, "raw", "qss", call.qss)
                add_metric_value(sample_metric_values, "raw", "qse", call.qse)

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
                        "svlen": call.svlen if call.svlen is not None else "",
                        "ac_orig": call.ac_orig,
                        "af_orig": call.af_orig,
                        "an_orig": call.an_orig,
                        "filter": call.filter,
                        "qual": call.qual,
                        "record_id": call.record_id,
                        "passes_analysis_filter": int(high_conf),
                    }
                )

            if high_conf:
                high_conf_filter_counts[filter_value] += 1
                high_conf_calls.append(call)
                counters["cnv_segments_high_conf"] += 1
                counters[f"cnv_{call.svtype.lower()}_high_conf"] += 1
                bp["cnv_bp_high_conf"] += call.length
                if call.np is not None:
                    np_sum["cnv_np_high_conf"] += call.np
                if call.qs is not None:
                    high_conf_qs_values.append(call.qs)
                add_metric_value(sample_metric_values, "high_conf", "length", call.length)
                add_metric_value(sample_metric_values, "high_conf", "np", call.np)
                add_metric_value(sample_metric_values, "high_conf", "qs", call.qs)
                add_metric_value(sample_metric_values, "high_conf", "qa", call.qa)
                add_metric_value(sample_metric_values, "high_conf", "qss", call.qss)
                add_metric_value(sample_metric_values, "high_conf", "qse", call.qse)

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
                "cnv_length_median_raw": median_or_blank(sample_metric_values[("raw", "length")]),
                "cnv_np_median_raw": median_or_blank(sample_metric_values[("raw", "np")]),
                "cnv_qs_median_raw": median_or_blank(sample_metric_values[("raw", "qs")]),
                "cnv_segments_high_conf": counters["cnv_segments_high_conf"],
                "cnv_del_high_conf": counters["cnv_del_high_conf"],
                "cnv_dup_high_conf": counters["cnv_dup_high_conf"],
                "hom_del_high_conf": counters["hom_del_high_conf"],
                "het_del_high_conf": counters["het_del_high_conf"],
                "dup_high_conf": counters["dup_high_conf"],
                "cnv_bp_high_conf": bp["cnv_bp_high_conf"],
                "cnv_np_high_conf": np_sum["cnv_np_high_conf"],
                "cnv_qs_mean_high_conf": mean_or_blank(high_conf_qs_values),
                "cnv_length_median_high_conf": median_or_blank(sample_metric_values[("high_conf", "length")]),
                "cnv_np_median_high_conf": median_or_blank(sample_metric_values[("high_conf", "np")]),
                "cnv_qs_median_high_conf": median_or_blank(sample_metric_values[("high_conf", "qs")]),
                "excluded_contig_records": counters["excluded_contig_records"],
            }
        )
        for filter_value in sorted(filter_counts):
            filter_rows.append(
                {
                    "cohort": cohort,
                    "sample": sample,
                    "sample_index": sample_index,
                    "filter": filter_value,
                    "segment_records_total": filter_counts[filter_value],
                    "raw_nonref_records": raw_filter_counts[filter_value],
                    "high_conf_records": high_conf_filter_counts[filter_value],
                }
            )
        quality_rows.extend(
            metric_rows(
                {"cohort": cohort, "sample": sample, "sample_index": sample_index},
                sample_metric_values,
            )
        )

    return segment_paths, cohort_counts, sample_rows, call_rows, raw_calls, high_conf_calls, filter_rows, quality_rows


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


def empty_interval_stats() -> dict[str, object]:
    return {
        "counter": Counter(),
        "cnq_values": [],
        "cnlp_best_values": [],
        "cnlp_second_best_values": [],
    }


def add_interval_observation(stats: dict[str, object], cn: str, cnq: int | None, cnlp: list[int]) -> None:
    counter: Counter = stats["counter"]  # type: ignore[assignment]
    counter["intervals_total"] += 1
    counter[f"cn_{cn}"] += 1
    cn_int = int_or_none(cn)
    if cn_int is not None:
        if cn_int < 2:
            counter["deletion_intervals"] += 1
        elif cn_int == 2:
            counter["diploid_intervals"] += 1
        elif cn_int > 2:
            counter["duplication_intervals"] += 1
        if cn_int != 2:
            counter["non_diploid_intervals"] += 1
        if cn_int >= 5:
            counter["cn5_plus_intervals"] += 1
    if cnq is not None:
        stats["cnq_values"].append(cnq)  # type: ignore[index,union-attr]
    if cnlp:
        sorted_cnlp = sorted(cnlp)
        stats["cnlp_best_values"].append(sorted_cnlp[0])  # type: ignore[index,union-attr]
        if len(sorted_cnlp) > 1:
            stats["cnlp_second_best_values"].append(sorted_cnlp[1])  # type: ignore[index,union-attr]


def interval_summary_row(base: dict[str, object], stats: dict[str, object]) -> dict[str, object]:
    counter: Counter = stats["counter"]  # type: ignore[assignment]
    row = {
        **base,
        "intervals_total": counter["intervals_total"],
        "cn0_intervals": counter["cn_0"],
        "cn1_intervals": counter["cn_1"],
        "cn2_intervals": counter["cn_2"],
        "cn3_intervals": counter["cn_3"],
        "cn4_intervals": counter["cn_4"],
        "cn5_plus_intervals": counter["cn5_plus_intervals"],
        "deletion_intervals": counter["deletion_intervals"],
        "diploid_intervals": counter["diploid_intervals"],
        "duplication_intervals": counter["duplication_intervals"],
        "non_diploid_intervals": counter["non_diploid_intervals"],
        "cnq_mean": mean_or_blank(stats["cnq_values"]),  # type: ignore[arg-type]
        "cnq_median": median_or_blank(stats["cnq_values"]),  # type: ignore[arg-type]
        "cnq_min": numeric_summary(stats["cnq_values"])["min"],  # type: ignore[arg-type]
        "cnq_q10": numeric_summary(stats["cnq_values"])["q10"],  # type: ignore[arg-type]
        "cnq_q90": numeric_summary(stats["cnq_values"])["q90"],  # type: ignore[arg-type]
        "cnq_max": numeric_summary(stats["cnq_values"])["max"],  # type: ignore[arg-type]
        "cnlp_best_mean": mean_or_blank(stats["cnlp_best_values"]),  # type: ignore[arg-type]
        "cnlp_second_best_mean": mean_or_blank(stats["cnlp_second_best_values"]),  # type: ignore[arg-type]
    }
    return row


def summarize_interval_vcfs(paths: list[Path]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    contig_rows: list[dict[str, object]] = []
    for path in paths:
        cohort, sample, sample_index = sample_from_segment_path(path)
        sample_stats = empty_interval_stats()
        contig_stats: defaultdict[str, dict[str, object]] = defaultdict(empty_interval_stats)
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    continue
                values = parse_sample_format(fields[8], fields[9])
                cn = values.get("CN", "")
                cnq = int_or_none(values.get("CNQ", ""))
                cnlp = int_values(values.get("CNLP", ""))
                add_interval_observation(sample_stats, cn, cnq, cnlp)
                add_interval_observation(contig_stats[fields[0]], cn, cnq, cnlp)
        rows.append(
            interval_summary_row(
                {"cohort": cohort, "sample": sample, "sample_index": sample_index},
                sample_stats,
            )
        )
        for contig, stats in sorted(contig_stats.items(), key=lambda item: contig_key(item[0])):
            contig_rows.append(
                interval_summary_row(
                    {
                        "cohort": cohort,
                        "sample": sample,
                        "sample_index": sample_index,
                        "contig": contig,
                    },
                    stats,
                )
            )
    return rows, contig_rows


def empty_copy_ratio_stats() -> dict[str, object]:
    return {
        "values": [],
        "rows_total": 0,
        "nonfinite": 0,
        "zero": 0,
        "near_zero": 0,
        "lt_1_0": 0,
        "lt_1_5": 0,
        "gt_2_5": 0,
        "gt_3_0": 0,
    }


def add_copy_ratio_observation(stats: dict[str, object], value: float | None) -> None:
    stats["rows_total"] = int(stats["rows_total"]) + 1
    if value is None or not math.isfinite(value):
        stats["nonfinite"] = int(stats["nonfinite"]) + 1
        return
    stats["values"].append(value)  # type: ignore[index,union-attr]
    if value == 0:
        stats["zero"] = int(stats["zero"]) + 1
    if value < 0.1:
        stats["near_zero"] = int(stats["near_zero"]) + 1
    if value < 1.0:
        stats["lt_1_0"] = int(stats["lt_1_0"]) + 1
    if value < 1.5:
        stats["lt_1_5"] = int(stats["lt_1_5"]) + 1
    if value > 2.5:
        stats["gt_2_5"] = int(stats["gt_2_5"]) + 1
    if value > 3.0:
        stats["gt_3_0"] = int(stats["gt_3_0"]) + 1


def copy_ratio_summary_row(base: dict[str, object], stats: dict[str, object]) -> dict[str, object]:
    values: list[float] = stats["values"]  # type: ignore[assignment]
    n = len(values)
    summary = numeric_summary(values)
    row = {
        **base,
        "rows_total": stats["rows_total"],
        "intervals": n,
        "nonfinite_intervals": stats["nonfinite"],
        "mean_linear_copy_ratio": summary["mean"],
        "stdev_linear_copy_ratio": summary["stdev"],
        "median_linear_copy_ratio": summary["median"],
        "mad_linear_copy_ratio": summary["mad"],
        "min_linear_copy_ratio": summary["min"],
        "q10_linear_copy_ratio": summary["q10"],
        "q25_linear_copy_ratio": summary["q25"],
        "q75_linear_copy_ratio": summary["q75"],
        "q90_linear_copy_ratio": summary["q90"],
        "max_linear_copy_ratio": summary["max"],
        "fraction_zero": f"{int(stats['zero']) / n:.6f}" if n else "",
        "fraction_lt_0_1": f"{int(stats['near_zero']) / n:.6f}" if n else "",
        "fraction_lt_1_0": f"{int(stats['lt_1_0']) / n:.6f}" if n else "",
        "fraction_lt_1_5": f"{int(stats['lt_1_5']) / n:.6f}" if n else "",
        "fraction_gt_2_5": f"{int(stats['gt_2_5']) / n:.6f}" if n else "",
        "fraction_gt_3_0": f"{int(stats['gt_3_0']) / n:.6f}" if n else "",
    }
    return row


def update_segment_copy_ratio_stats(
    contig_call_lists: dict[str, list[tuple[int, SegmentCall]]],
    active_indices: dict[str, int],
    stats: list[dict[str, object]],
    contig: str,
    start: int,
    end: int,
    value: float,
) -> None:
    calls = contig_call_lists.get(contig)
    if not calls:
        return
    active_idx = active_indices.get(contig, 0)
    while active_idx < len(calls) and calls[active_idx][1].end < start:
        active_idx += 1
    active_indices[contig] = active_idx

    idx = active_idx
    while idx < len(calls) and calls[idx][1].start <= end:
        call_index, call = calls[idx]
        overlap = max(0, min(call.end, end) - max(call.start, start) + 1)
        if overlap:
            stats[call_index]["overlap_bp"] = int(stats[call_index]["overlap_bp"]) + overlap
            stats[call_index]["intervals"] = int(stats[call_index]["intervals"]) + 1
            stats[call_index]["weighted_sum"] = float(stats[call_index]["weighted_sum"]) + value * overlap
        idx += 1


def summarize_copy_ratios(
    paths: list[Path],
    raw_calls: list[SegmentCall],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    contig_rows: list[dict[str, object]] = []
    segment_rows: list[dict[str, object]] = []
    calls_by_sample: defaultdict[str, list[SegmentCall]] = defaultdict(list)
    for call in raw_calls:
        calls_by_sample[call.sample].append(call)

    for path in paths:
        match = re.match(
            r"COHORT_(?P<cohort>\d+)_SAMPLE_(?P<sample>.+)_(?P<index>\d+)\.tsv$",
            path.name,
        )
        if not match:
            continue
        sample = match.group("sample")
        sample_stats = empty_copy_ratio_stats()
        contig_stats: defaultdict[str, dict[str, object]] = defaultdict(empty_copy_ratio_stats)
        sample_calls = sorted(calls_by_sample.get(sample, []), key=lambda call: (contig_key(call.contig), call.start, call.end))
        contig_call_lists: dict[str, list[tuple[int, SegmentCall]]] = defaultdict(list)
        for call_index, call in enumerate(sample_calls):
            contig_call_lists[call.contig].append((call_index, call))
        active_indices: dict[str, int] = {}
        segment_stats = [
            {"overlap_bp": 0, "intervals": 0, "weighted_sum": 0.0}
            for _ in sample_calls
        ]

        with path.open() as handle:
            reader = csv.DictReader((line for line in handle if not line.startswith("@")), delimiter="\t")
            for row in reader:
                value = float_or_none(row.get("LINEAR_COPY_RATIO", ""))
                contig = row.get("CONTIG", "")
                add_copy_ratio_observation(sample_stats, value)
                if contig:
                    add_copy_ratio_observation(contig_stats[contig], value)
                if value is None or not math.isfinite(value):
                    continue
                start = int_or_none(row.get("START", ""))
                end = int_or_none(row.get("END", ""))
                if contig and start is not None and end is not None:
                    update_segment_copy_ratio_stats(
                        contig_call_lists,
                        active_indices,
                        segment_stats,
                        contig,
                        start,
                        end,
                        value,
                    )

        rows.append(
            copy_ratio_summary_row(
                {
                    "cohort": match.group("cohort"),
                    "sample": sample,
                    "sample_index": match.group("index"),
                },
                sample_stats,
            )
        )
        for contig, stats in sorted(contig_stats.items(), key=lambda item: contig_key(item[0])):
            contig_rows.append(
                copy_ratio_summary_row(
                    {
                        "cohort": match.group("cohort"),
                        "sample": sample,
                        "sample_index": match.group("index"),
                        "contig": contig,
                    },
                    stats,
                )
            )
        for call, stats in zip(sample_calls, segment_stats):
            overlap_bp = int(stats["overlap_bp"])
            if not overlap_bp:
                continue
            segment_rows.append(
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
                    "passes_analysis_filter": int(call.passes_filter),
                    "copy_ratio_intervals": stats["intervals"],
                    "copy_ratio_overlap_bp": overlap_bp,
                    "mean_linear_copy_ratio": format_number(float(stats["weighted_sum"]) / overlap_bp),
                }
            )
    return rows, contig_rows, segment_rows


@dataclass
class AnnotationRegion:
    source: str
    name: str
    contig: str
    start: int
    end: int


def read_annotation_bed(paths: list[Path]) -> list[AnnotationRegion]:
    regions = []
    for path in paths:
        source = path.stem
        with path.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                bed_start = int_or_none(fields[1])
                bed_end = int_or_none(fields[2])
                if bed_start is None or bed_end is None or bed_end <= bed_start:
                    continue
                name = fields[3] if len(fields) > 3 and fields[3] else f"{fields[0]}:{bed_start + 1}-{bed_end}"
                regions.append(
                    AnnotationRegion(
                        source=source,
                        name=name,
                        contig=fields[0],
                        start=bed_start + 1,
                        end=bed_end,
                    )
                )
    return sorted(regions, key=lambda region: (contig_key(region.contig), region.start, region.end, region.name))


def summarize_annotation_overlaps(
    annotation_paths: list[Path],
    raw_calls: list[SegmentCall],
    high_conf_calls: list[SegmentCall],
) -> list[dict[str, object]]:
    regions = read_annotation_bed(annotation_paths)
    if not regions:
        return []

    rows = []
    for scope, calls in [("raw", raw_calls), ("high_conf", high_conf_calls)]:
        calls_by_contig: defaultdict[str, list[SegmentCall]] = defaultdict(list)
        for call in calls:
            calls_by_contig[call.contig].append(call)
        for contig_calls in calls_by_contig.values():
            contig_calls.sort(key=lambda call: (call.start, call.end))

        active_indices: defaultdict[str, int] = defaultdict(int)
        for region in regions:
            contig_calls = calls_by_contig.get(region.contig, [])
            active_idx = active_indices[region.contig]
            while active_idx < len(contig_calls) and contig_calls[active_idx].end < region.start:
                active_idx += 1
            active_indices[region.contig] = active_idx

            overlaps_by_svtype: defaultdict[str, list[tuple[SegmentCall, int]]] = defaultdict(list)
            idx = active_idx
            while idx < len(contig_calls) and contig_calls[idx].start <= region.end:
                call = contig_calls[idx]
                overlap = max(0, min(call.end, region.end) - max(call.start, region.start) + 1)
                if overlap:
                    overlaps_by_svtype[call.svtype].append((call, overlap))
                idx += 1

            for svtype, overlaps in sorted(overlaps_by_svtype.items()):
                calls_for_region = [call for call, _ in overlaps]
                overlap_bp = sum(overlap for _, overlap in overlaps)
                rows.append(
                    {
                        "scope": scope,
                        "annotation_source": region.source,
                        "annotation_name": region.name,
                        "contig": region.contig,
                        "start": region.start,
                        "end": region.end,
                        "svtype": svtype,
                        "call_count": len(overlaps),
                        "sample_count": len({call.sample for call in calls_for_region}),
                        "overlap_bp": overlap_bp,
                        "mean_overlap_fraction": mean_or_blank(
                            [overlap / call.length for call, overlap in overlaps if call.length]
                        ),
                        "samples": ",".join(sorted({call.sample for call in calls_for_region})),
                    }
                )
    return rows


def write_pca_summaries(gatk_gcnv: Path, out_dir: Path) -> None:
    coordinate_files = sorted((gatk_gcnv / "cohort_report").glob("*.pca_2d_coordinates.tsv"))
    coordinate_rows: list[dict[str, str]] = []
    for coordinate_file in coordinate_files:
        with coordinate_file.open(newline="") as handle:
            coordinate_rows.extend(csv.DictReader(handle, delimiter="\t"))
    if not coordinate_rows:
        return

    assigned_by_cohort: defaultdict[str, list[tuple[float, float]]] = defaultdict(list)
    for row in coordinate_rows:
        pc1 = float_or_none(row.get("pc1", ""))
        pc2 = float_or_none(row.get("pc2", ""))
        cohort = row.get("cohort", "")
        if row.get("assignment") == "assigned" and cohort and pc1 is not None and pc2 is not None:
            assigned_by_cohort[cohort].append((pc1, pc2))

    centroids = {
        cohort: (
            statistics.fmean(pc1 for pc1, _ in points),
            statistics.fmean(pc2 for _, pc2 in points),
        )
        for cohort, points in assigned_by_cohort.items()
        if points
    }

    sample_distance_rows = []
    grouped_values: defaultdict[tuple[str, str], dict[tuple[str, str], list[int | float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for row in coordinate_rows:
        pc1 = float_or_none(row.get("pc1", ""))
        pc2 = float_or_none(row.get("pc2", ""))
        cohort = row.get("cohort", "")
        assignment = row.get("assignment", "")
        distance = None
        if pc1 is not None and pc2 is not None and cohort in centroids:
            centroid_pc1, centroid_pc2 = centroids[cohort]
            distance = math.sqrt((pc1 - centroid_pc1) ** 2 + (pc2 - centroid_pc2) ** 2)

        sample_distance_rows.append(
            {
                "samplefile": row.get("samplefile", ""),
                "sample": row.get("sample", ""),
                "assignment": assignment,
                "cohort": cohort,
                "pc1": pc1 if pc1 is not None else "",
                "pc2": pc2 if pc2 is not None else "",
                "distance_to_cohort_centroid": format_number(distance),
            }
        )

        for group_key in [("assignment", assignment or "unknown"), ("cohort", cohort or "unassigned")]:
            if pc1 is not None:
                grouped_values[group_key][("pca", "pc1")].append(pc1)
            if pc2 is not None:
                grouped_values[group_key][("pca", "pc2")].append(pc2)
            if distance is not None:
                grouped_values[group_key][("pca", "distance_to_cohort_centroid")].append(distance)

    summary_rows = []
    for (group_type, group_value), values_by_metric in sorted(grouped_values.items()):
        group_sample_count = sum(
            1
            for row in coordinate_rows
            if (
                (group_type == "assignment" and (row.get("assignment", "") or "unknown") == group_value)
                or (group_type == "cohort" and (row.get("cohort", "") or "unassigned") == group_value)
            )
        )
        summary_rows.extend(
            metric_rows(
                {
                    "group_type": group_type,
                    "group_value": group_value,
                    "sample_count": group_sample_count,
                },
                values_by_metric,
            )
        )

    write_tsv(
        out_dir / "sample_pca_distances.tsv",
        ["samplefile", "sample", "assignment", "cohort", "pc1", "pc2", "distance_to_cohort_centroid"],
        sample_distance_rows,
    )
    write_tsv(
        out_dir / "cohort_pca_summary.tsv",
        ["group_type", "group_value", "sample_count", "scope", "metric"] + STAT_COLUMNS,
        summary_rows,
    )


def write_summary_text(
    out_dir: Path,
    args: argparse.Namespace,
    workflow: str,
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
        print(f"Workflow layout: {workflow}", file=handle)
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


def analyze_one_run(args: argparse.Namespace, workflow: str, gatk_gcnv: Path, out_dir: Path) -> AnalysisResult:
    run_args = argparse.Namespace(**vars(args))
    run_args.gatk_gcnv = gatk_gcnv
    run_args.out_dir = out_dir
    if not has_gcnv_outputs(run_args.gatk_gcnv):
        raise SystemExit(f"No gCNV outputs found for {workflow} in {run_args.gatk_gcnv}.")
    run_args.out_dir.mkdir(parents=True, exist_ok=True)
    excluded_contigs = {contig.strip() for contig in args.exclude_contigs.split(",") if contig.strip()}
    cohort_sample_map, sample_to_samplefile, cohort_report_sizes, unassigned_rows = read_cohort_report(run_args.gatk_gcnv)
    (
        segment_paths,
        segment_cohort_counts,
        sample_rows,
        call_rows,
        raw_calls,
        high_conf_calls,
        filter_rows,
        quality_rows,
    ) = summarize_segments(
        run_args,
        excluded_contigs,
    )
    contig_lengths = read_contig_lengths(segment_paths)

    write_tsv(
        run_args.out_dir / "sample_segment_summary.tsv",
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
            "cnv_length_median_raw",
            "cnv_np_median_raw",
            "cnv_qs_median_raw",
            "cnv_segments_high_conf",
            "cnv_del_high_conf",
            "cnv_dup_high_conf",
            "hom_del_high_conf",
            "het_del_high_conf",
            "dup_high_conf",
            "cnv_bp_high_conf",
            "cnv_np_high_conf",
            "cnv_qs_mean_high_conf",
            "cnv_length_median_high_conf",
            "cnv_np_median_high_conf",
            "cnv_qs_median_high_conf",
            "excluded_contig_records",
        ],
        sorted(sample_rows, key=lambda row: (str(row["cohort"]), int(row["sample_index"]))),
    )

    write_tsv(
        run_args.out_dir / "segment_calls.tsv",
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
            "svlen",
            "ac_orig",
            "af_orig",
            "an_orig",
            "filter",
            "qual",
            "record_id",
            "passes_analysis_filter",
        ],
        call_rows,
    )

    write_tsv(
        run_args.out_dir / "sample_quality_summary.tsv",
        ["cohort", "sample", "sample_index", "scope", "metric"] + STAT_COLUMNS,
        quality_rows,
    )

    write_tsv(
        run_args.out_dir / "filter_summary.tsv",
        [
            "cohort",
            "sample",
            "sample_index",
            "filter",
            "segment_records_total",
            "raw_nonref_records",
            "high_conf_records",
        ],
        filter_rows,
    )

    write_tsv(
        run_args.out_dir / "contig_summary.tsv",
        [
            "scope",
            "cohort",
            "contig",
            "svtype",
            "excluded_from_high_conf",
            "contig_length",
            "call_count",
            "sample_count",
            "cnv_bp",
            "call_rate_per_mb",
            "bp_per_mb",
            "mean_length",
            "median_length",
            "mean_qs",
            "median_qs",
            "mean_np",
            "median_np",
        ],
        summarize_contigs(raw_calls, high_conf_calls, contig_lengths, excluded_contigs),
    )

    write_tsv(
        run_args.out_dir / "size_distribution.tsv",
        ["scope", "cohort", "svtype"] + STAT_COLUMNS,
        summarize_size_distributions(raw_calls, high_conf_calls),
    )

    write_tsv(
        run_args.out_dir / "svtype_summary.tsv",
        [
            "scope",
            "cohort",
            "svtype",
            "call_count",
            "sample_count",
            "cnv_bp",
            "cn0_calls",
            "cn1_calls",
            "cn2_calls",
            "cn3_calls",
            "cn4_calls",
            "cn5_plus_calls",
        ]
        + [f"length_{column}" for column in STAT_COLUMNS]
        + [f"qs_{column}" for column in STAT_COLUMNS],
        summarize_svtypes(raw_calls, high_conf_calls),
    )

    write_tsv(
        run_args.out_dir / "allele_frequency_summary.tsv",
        [
            "scope",
            "cohort",
            "svtype",
            "call_count",
            "calls_with_af",
            "calls_with_ac",
            "calls_with_an",
            "af_lt_0_001",
            "af_0_001_0_01",
            "af_0_01_0_05",
            "af_0_05_0_10",
            "af_ge_0_10",
        ]
        + [f"af_orig_{column}" for column in STAT_COLUMNS]
        + [f"ac_orig_{column}" for column in STAT_COLUMNS]
        + [f"an_orig_{column}" for column in STAT_COLUMNS],
        summarize_allele_frequencies(raw_calls, high_conf_calls),
    )

    write_tsv(
        run_args.out_dir / "group_burden_summary.tsv",
        ["group_type", "group_value", "sample_count", "metric"] + STAT_COLUMNS,
        summarize_group_burdens(sample_rows, sample_to_samplefile),
    )

    if run_args.annotation_bed:
        write_tsv(
            run_args.out_dir / "annotation_overlap_summary.tsv",
            [
                "scope",
                "annotation_source",
                "annotation_name",
                "contig",
                "start",
                "end",
                "svtype",
                "call_count",
                "sample_count",
                "overlap_bp",
                "mean_overlap_fraction",
                "samples",
            ],
            summarize_annotation_overlaps(run_args.annotation_bed, raw_calls, high_conf_calls),
        )

    write_pca_summaries(run_args.gatk_gcnv, run_args.out_dir)

    recurrent_rows = recurrent_loci(high_conf_calls, run_args.merge_gap, run_args.min_recurrent_samples)
    write_tsv(
        run_args.out_dir / "recurrent_loci.tsv",
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
        run_args.out_dir / "cohort_summary.tsv",
        ["cohort", "assigned_samples_report", "segment_vcf_samples"],
        cohort_rows,
    )
    if unassigned_rows:
        write_tsv(
            run_args.out_dir / "unassigned_samples.tsv",
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
        write_tsv(run_args.out_dir / "missing_segment_outputs.tsv", ["sample", "cohort"], missing_segment_rows)

    flags = write_qc_flags(
        run_args.out_dir,
        sample_rows,
        cohort_report_sizes,
        segment_cohort_counts,
        unassigned_rows,
        run_args,
    )

    if run_args.include_intervals:
        interval_rows, interval_contig_rows = summarize_interval_vcfs(
            discover_interval_vcfs(run_args.gatk_gcnv, run_args.max_interval_files)
        )
        write_tsv(
            run_args.out_dir / "interval_summary.tsv",
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
                "cn5_plus_intervals",
                "deletion_intervals",
                "diploid_intervals",
                "duplication_intervals",
                "non_diploid_intervals",
                "cnq_mean",
                "cnq_median",
                "cnq_min",
                "cnq_q10",
                "cnq_q90",
                "cnq_max",
                "cnlp_best_mean",
                "cnlp_second_best_mean",
            ],
            interval_rows,
        )
        write_tsv(
            run_args.out_dir / "interval_contig_summary.tsv",
            [
                "cohort",
                "sample",
                "sample_index",
                "contig",
                "intervals_total",
                "cn0_intervals",
                "cn1_intervals",
                "cn2_intervals",
                "cn3_intervals",
                "cn4_intervals",
                "cn5_plus_intervals",
                "deletion_intervals",
                "diploid_intervals",
                "duplication_intervals",
                "non_diploid_intervals",
                "cnq_mean",
                "cnq_median",
                "cnq_min",
                "cnq_q10",
                "cnq_q90",
                "cnq_max",
                "cnlp_best_mean",
                "cnlp_second_best_mean",
            ],
            interval_contig_rows,
        )

    if run_args.include_copy_ratio:
        copy_ratio_rows, copy_ratio_contig_rows, segment_copy_ratio_rows = summarize_copy_ratios(
            discover_copy_ratio_tsvs(run_args.gatk_gcnv, run_args.max_copy_ratio_files),
            raw_calls,
        )
        copy_ratio_fields = [
            "cohort",
            "sample",
            "sample_index",
            "rows_total",
            "intervals",
            "nonfinite_intervals",
            "mean_linear_copy_ratio",
            "stdev_linear_copy_ratio",
            "median_linear_copy_ratio",
            "mad_linear_copy_ratio",
            "min_linear_copy_ratio",
            "q10_linear_copy_ratio",
            "q25_linear_copy_ratio",
            "q75_linear_copy_ratio",
            "q90_linear_copy_ratio",
            "max_linear_copy_ratio",
            "fraction_zero",
            "fraction_lt_0_1",
            "fraction_lt_1_0",
            "fraction_lt_1_5",
            "fraction_gt_2_5",
            "fraction_gt_3_0",
        ]
        write_tsv(
            run_args.out_dir / "copy_ratio_summary.tsv",
            copy_ratio_fields,
            copy_ratio_rows,
        )
        write_tsv(
            run_args.out_dir / "copy_ratio_contig_summary.tsv",
            copy_ratio_fields[:3] + ["contig"] + copy_ratio_fields[3:],
            copy_ratio_contig_rows,
        )
        write_tsv(
            run_args.out_dir / "segment_copy_ratio_summary.tsv",
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
                "passes_analysis_filter",
                "copy_ratio_intervals",
                "copy_ratio_overlap_bp",
                "mean_linear_copy_ratio",
            ],
            segment_copy_ratio_rows,
        )

    write_summary_text(
        run_args.out_dir,
        run_args,
        workflow,
        segment_paths,
        cohort_report_sizes,
        segment_cohort_counts,
        unassigned_rows,
        sample_rows,
        flags,
        recurrent_rows,
    )

    return AnalysisResult(
        workflow=workflow,
        gatk_gcnv=run_args.gatk_gcnv,
        out_dir=run_args.out_dir,
        segment_paths=segment_paths,
        cohort_report_sizes=cohort_report_sizes,
        segment_cohort_counts=segment_cohort_counts,
        unassigned_rows=unassigned_rows,
        sample_to_samplefile=sample_to_samplefile,
        sample_rows=sample_rows,
        call_rows=call_rows,
        raw_calls=raw_calls,
        high_conf_calls=high_conf_calls,
        recurrent_rows=recurrent_rows,
        flags=flags,
    )


def workflow_prefix(workflow: str) -> str:
    return workflow.replace("-", "_")


def int_row_value(row: dict[str, object], key: str) -> int:
    value = row.get(key, 0)
    if value in {"", None}:
        return 0
    return int(value)


def write_sample_comparison(out_dir: Path, left: AnalysisResult, right: AnalysisResult) -> None:
    left_prefix = workflow_prefix(left.workflow)
    right_prefix = workflow_prefix(right.workflow)
    left_by_sample = {str(row["sample"]): row for row in left.sample_rows}
    right_by_sample = {str(row["sample"]): row for row in right.sample_rows}
    rows = []

    for sample in sorted(set(left_by_sample) | set(right_by_sample)):
        left_row = left_by_sample.get(sample, {})
        right_row = right_by_sample.get(sample, {})
        left_high = int_row_value(left_row, "cnv_segments_high_conf")
        right_high = int_row_value(right_row, "cnv_segments_high_conf")
        left_raw = int_row_value(left_row, "cnv_segments_raw")
        right_raw = int_row_value(right_row, "cnv_segments_raw")
        left_bp = int_row_value(left_row, "cnv_bp_high_conf")
        right_bp = int_row_value(right_row, "cnv_bp_high_conf")
        if left_row and right_row:
            status = "both"
        elif left_row:
            status = f"{left.workflow}_only"
        else:
            status = f"{right.workflow}_only"
        rows.append(
            {
                "sample": sample,
                "status": status,
                f"{left_prefix}_cohort": left_row.get("cohort", ""),
                f"{right_prefix}_cohort": right_row.get("cohort", ""),
                f"{left_prefix}_raw_cnv_segments": left_raw if left_row else "",
                f"{right_prefix}_raw_cnv_segments": right_raw if right_row else "",
                "delta_raw_cnv_segments": right_raw - left_raw if left_row and right_row else "",
                f"{left_prefix}_high_conf_cnv_segments": left_high if left_row else "",
                f"{right_prefix}_high_conf_cnv_segments": right_high if right_row else "",
                "delta_high_conf_cnv_segments": right_high - left_high if left_row and right_row else "",
                f"{left_prefix}_high_conf_cnv_bp": left_bp if left_row else "",
                f"{right_prefix}_high_conf_cnv_bp": right_bp if right_row else "",
                "delta_high_conf_cnv_bp": right_bp - left_bp if left_row and right_row else "",
            }
        )

    write_tsv(
        out_dir / "comparison_sample_summary.tsv",
        [
            "sample",
            "status",
            f"{left_prefix}_cohort",
            f"{right_prefix}_cohort",
            f"{left_prefix}_raw_cnv_segments",
            f"{right_prefix}_raw_cnv_segments",
            "delta_raw_cnv_segments",
            f"{left_prefix}_high_conf_cnv_segments",
            f"{right_prefix}_high_conf_cnv_segments",
            "delta_high_conf_cnv_segments",
            f"{left_prefix}_high_conf_cnv_bp",
            f"{right_prefix}_high_conf_cnv_bp",
            "delta_high_conf_cnv_bp",
        ],
        rows,
    )


def reciprocal_overlap(left: SegmentCall, right: SegmentCall) -> tuple[int, float, float]:
    overlap = max(0, min(left.end, right.end) - max(left.start, right.start) + 1)
    if overlap == 0:
        return 0, 0.0, 0.0
    return overlap, overlap / left.length if left.length else 0.0, overlap / right.length if right.length else 0.0


def write_call_overlap_comparison(
    out_dir: Path,
    left: AnalysisResult,
    right: AnalysisResult,
    min_reciprocal_overlap: float = 0.5,
) -> Counter:
    right_groups: dict[tuple[str, str, str], list[tuple[int, SegmentCall]]] = defaultdict(list)
    for idx, call in enumerate(right.high_conf_calls):
        right_groups[(call.sample, call.contig, call.svtype)].append((idx, call))

    rows = []
    matched_right: set[int] = set()
    counts = Counter()

    for left_call in left.high_conf_calls:
        best_idx = None
        best_call = None
        best_overlap = 0
        best_left_fraction = 0.0
        best_right_fraction = 0.0
        for right_idx, right_call in right_groups.get((left_call.sample, left_call.contig, left_call.svtype), []):
            if right_idx in matched_right:
                continue
            overlap, left_fraction, right_fraction = reciprocal_overlap(left_call, right_call)
            if left_fraction < min_reciprocal_overlap or right_fraction < min_reciprocal_overlap:
                continue
            if overlap > best_overlap:
                best_idx = right_idx
                best_call = right_call
                best_overlap = overlap
                best_left_fraction = left_fraction
                best_right_fraction = right_fraction

        if best_call is None or best_idx is None:
            counts["left_only"] += 1
            rows.append(
                {
                    "status": f"{left.workflow}_only",
                    "sample": left_call.sample,
                    "contig": left_call.contig,
                    "svtype": left_call.svtype,
                    f"{workflow_prefix(left.workflow)}_start": left_call.start,
                    f"{workflow_prefix(left.workflow)}_end": left_call.end,
                    f"{workflow_prefix(left.workflow)}_cn": left_call.cn,
                    f"{workflow_prefix(left.workflow)}_qs": left_call.qs if left_call.qs is not None else "",
                    f"{workflow_prefix(right.workflow)}_start": "",
                    f"{workflow_prefix(right.workflow)}_end": "",
                    f"{workflow_prefix(right.workflow)}_cn": "",
                    f"{workflow_prefix(right.workflow)}_qs": "",
                    "overlap_bp": "",
                    "left_overlap_fraction": "",
                    "right_overlap_fraction": "",
                }
            )
            continue

        matched_right.add(best_idx)
        counts["matched"] += 1
        rows.append(
            {
                "status": "matched",
                "sample": left_call.sample,
                "contig": left_call.contig,
                "svtype": left_call.svtype,
                f"{workflow_prefix(left.workflow)}_start": left_call.start,
                f"{workflow_prefix(left.workflow)}_end": left_call.end,
                f"{workflow_prefix(left.workflow)}_cn": left_call.cn,
                f"{workflow_prefix(left.workflow)}_qs": left_call.qs if left_call.qs is not None else "",
                f"{workflow_prefix(right.workflow)}_start": best_call.start,
                f"{workflow_prefix(right.workflow)}_end": best_call.end,
                f"{workflow_prefix(right.workflow)}_cn": best_call.cn,
                f"{workflow_prefix(right.workflow)}_qs": best_call.qs if best_call.qs is not None else "",
                "overlap_bp": best_overlap,
                "left_overlap_fraction": f"{best_left_fraction:.4f}",
                "right_overlap_fraction": f"{best_right_fraction:.4f}",
            }
        )

    for right_idx, right_call in enumerate(right.high_conf_calls):
        if right_idx in matched_right:
            continue
        counts["right_only"] += 1
        rows.append(
            {
                "status": f"{right.workflow}_only",
                "sample": right_call.sample,
                "contig": right_call.contig,
                "svtype": right_call.svtype,
                f"{workflow_prefix(left.workflow)}_start": "",
                f"{workflow_prefix(left.workflow)}_end": "",
                f"{workflow_prefix(left.workflow)}_cn": "",
                f"{workflow_prefix(left.workflow)}_qs": "",
                f"{workflow_prefix(right.workflow)}_start": right_call.start,
                f"{workflow_prefix(right.workflow)}_end": right_call.end,
                f"{workflow_prefix(right.workflow)}_cn": right_call.cn,
                f"{workflow_prefix(right.workflow)}_qs": right_call.qs if right_call.qs is not None else "",
                "overlap_bp": "",
                "left_overlap_fraction": "",
                "right_overlap_fraction": "",
            }
        )

    write_tsv(
        out_dir / "comparison_high_conf_call_overlap.tsv",
        [
            "status",
            "sample",
            "contig",
            "svtype",
            f"{workflow_prefix(left.workflow)}_start",
            f"{workflow_prefix(left.workflow)}_end",
            f"{workflow_prefix(left.workflow)}_cn",
            f"{workflow_prefix(left.workflow)}_qs",
            f"{workflow_prefix(right.workflow)}_start",
            f"{workflow_prefix(right.workflow)}_end",
            f"{workflow_prefix(right.workflow)}_cn",
            f"{workflow_prefix(right.workflow)}_qs",
            "overlap_bp",
            "left_overlap_fraction",
            "right_overlap_fraction",
        ],
        rows,
    )
    return counts


def write_combined_overview(out_dir: Path, results: list[AnalysisResult], overlap_counts: Counter | None = None) -> None:
    with (out_dir / "combined_overview.txt").open("w") as handle:
        print("Combined gCNV overview", file=handle)
        print("", file=handle)
        for result in results:
            high_conf_total = sum(int_row_value(row, "cnv_segments_high_conf") for row in result.sample_rows)
            raw_total = sum(int_row_value(row, "cnv_segments_raw") for row in result.sample_rows)
            samples = len(result.sample_rows)
            cohorts = len(result.segment_cohort_counts)
            print(f"{result.workflow}:", file=handle)
            print(f"  directory: {result.gatk_gcnv}", file=handle)
            print(f"  samples with segment VCFs: {samples}", file=handle)
            print(f"  cohorts with segment VCFs: {cohorts}", file=handle)
            print(f"  raw non-reference CNV segments: {raw_total}", file=handle)
            print(f"  high-confidence CNV segments: {high_conf_total}", file=handle)
            print(f"  recurrent loci: {len(result.recurrent_rows)}", file=handle)
            print(f"  QC flags: {len(result.flags)}", file=handle)
            print("", file=handle)

        if len(results) == 2:
            left, right = results
            left_samples = {str(row["sample"]) for row in left.sample_rows}
            right_samples = {str(row["sample"]) for row in right.sample_rows}
            shared_samples = left_samples & right_samples
            print("Comparison:", file=handle)
            print(f"  shared samples: {len(shared_samples)}", file=handle)
            print(f"  {left.workflow}-only samples: {len(left_samples - right_samples)}", file=handle)
            print(f"  {right.workflow}-only samples: {len(right_samples - left_samples)}", file=handle)
            if overlap_counts is not None:
                print(f"  matched high-confidence calls: {overlap_counts['matched']}", file=handle)
                print(f"  {left.workflow}-only high-confidence calls: {overlap_counts['left_only']}", file=handle)
                print(f"  {right.workflow}-only high-confidence calls: {overlap_counts['right_only']}", file=handle)


def write_comparison_outputs(out_dir: Path, results: list[AnalysisResult]) -> None:
    if len(results) != 2:
        write_combined_overview(out_dir, results)
        return
    left, right = results
    write_sample_comparison(out_dir, left, right)
    overlap_counts = write_call_overlap_comparison(out_dir, left, right)
    write_combined_overview(out_dir, results, overlap_counts)


def main() -> None:
    args = parse_args()
    run_dirs = resolve_run_dirs(args.gatk_gcnv, args.workflow)
    multi_run = len(run_dirs) > 1
    results = []
    for workflow, gatk_gcnv in run_dirs:
        out_dir = args.out_dir / workflow if multi_run else args.out_dir
        results.append(analyze_one_run(args, workflow, gatk_gcnv, out_dir))

    if multi_run:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        write_comparison_outputs(args.out_dir, results)

    print(f"Wrote gCNV analysis to {args.out_dir}")


if __name__ == "__main__":
    main()
