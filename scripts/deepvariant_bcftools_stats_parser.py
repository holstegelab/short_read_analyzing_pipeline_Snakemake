#!/usr/bin/env python3

import argparse
import csv
import re
import subprocess
from pathlib import Path

SECTION_PARSERS = {}


def sanitize_key(raw_key):
    key = raw_key.rstrip(':').lower()
    key = re.sub(r"[^0-9a-z]+", "_", key)
    return re.sub(r"_+", "_", key).strip("_")


def section(name):
    def decorator(func):
        SECTION_PARSERS[name] = func
        return func

    return decorator


@section("SN")
def parse_summary(tokens, state):
    _, _, key, value = tokens
    state.setdefault("summary", {})[sanitize_key(key)] = value


@section("TSTV")
def parse_tstv(tokens, state):
    _, _, ts, tv, ratio, ts1, tv1, ratio1 = tokens
    state.setdefault("tstv", []).append(
        {
            "ts": ts,
            "tv": tv,
            "ts_tv": ratio,
            "ts_first_alt": ts1,
            "tv_first_alt": tv1,
            "ts_tv_first_alt": ratio1,
        }
    )


@section("ICS")
def parse_ics(tokens, state):
    _, _, repeat_consistent, repeat_inconsistent, not_applicable, ratio = tokens
    state.setdefault("ics", []).append(
        {
            "repeat_consistent": repeat_consistent,
            "repeat_inconsistent": repeat_inconsistent,
            "not_applicable": not_applicable,
            "consistent_ratio": ratio,
        }
    )


@section("ICL")
def parse_icl(tokens, state):
    _, _, length, rc_del, ri_del, cons_ins, incons_ins, ratio = tokens
    state.setdefault("icl", []).append(
        {
            "length": length,
            "repeat_consistent_del": rc_del,
            "repeat_inconsistent_del": ri_del,
            "consistent_ins": cons_ins,
            "inconsistent_ins": incons_ins,
            "consistent_ratio": ratio,
        }
    )


@section("SiS")
def parse_sis(tokens, state):
    _, _, allele_count, snps, transitions, transversions, indels, repeat_consistent, repeat_inconsistent, not_applicable = tokens
    state.setdefault("sis", []).append(
        {
            "allele_count": allele_count,
            "snps": snps,
            "transitions": transitions,
            "transversions": transversions,
            "indels": indels,
            "repeat_consistent": repeat_consistent,
            "repeat_inconsistent": repeat_inconsistent,
            "not_applicable": not_applicable,
        }
    )


@section("AF")
def parse_af(tokens, state):
    _, _, allele_freq, snps, transitions, transversions, indels, repeat_consistent, repeat_inconsistent, not_applicable = tokens
    state.setdefault("af", []).append(
        {
            "allele_freq": allele_freq,
            "snps": snps,
            "transitions": transitions,
            "transversions": transversions,
            "indels": indels,
            "repeat_consistent": repeat_consistent,
            "repeat_inconsistent": repeat_inconsistent,
            "not_applicable": not_applicable,
        }
    )


@section("QUAL")
def parse_qual(tokens, state):
    _, _, qual, snps, transitions, transversions, indels = tokens
    state.setdefault("qual", []).append(
        {
            "qual": qual,
            "snps": snps,
            "transitions": transitions,
            "transversions": transversions,
            "indels": indels,
        }
    )


@section("ST")
def parse_st(tokens, state):
    _, _, mut_type, count = tokens
    state.setdefault("st", []).append({"substitution_type": mut_type, "count": count})


@section("IDD")
def parse_idd(tokens, state):
    _, _, length, n_sites, n_genotypes, mean_vaf = tokens
    state.setdefault("idd", []).append(
        {
            "length": length,
            "n_sites": n_sites,
            "n_genotypes": n_genotypes,
            "mean_vaf": mean_vaf,
        }
    )


SECTION_HEADERS = {
    "summary": [
        "number_of_samples",
        "number_of_records",
        "number_of_no_alts",
        "number_of_snps",
        "number_of_mnps",
        "number_of_indels",
        "number_of_others",
        "number_of_multiallelic_sites",
        "number_of_multiallelic_snp_sites",
    ],
    "tstv": [
        "ts",
        "tv",
        "ts_tv",
        "ts_first_alt",
        "tv_first_alt",
        "ts_tv_first_alt",
    ],
    "ics": [
        "repeat_consistent",
        "repeat_inconsistent",
        "not_applicable",
        "consistent_ratio",
    ],
    "icl": [
        "length",
        "repeat_consistent_del",
        "repeat_inconsistent_del",
        "consistent_ins",
        "inconsistent_ins",
        "consistent_ratio",
    ],
    "sis": [
        "allele_count",
        "snps",
        "transitions",
        "transversions",
        "indels",
        "repeat_consistent",
        "repeat_inconsistent",
        "not_applicable",
    ],
    "af": [
        "allele_freq",
        "snps",
        "transitions",
        "transversions",
        "indels",
        "repeat_consistent",
        "repeat_inconsistent",
        "not_applicable",
    ],
    "qual": [
        "qual",
        "snps",
        "transitions",
        "transversions",
        "indels",
    ],
    "st": [
        "substitution_type",
        "count",
    ],
    "idd": [
        "length",
        "n_sites",
        "n_genotypes",
        "mean_vaf",
    ],
}


def parse_bcftools_stats(stats_path):
    state = {}
    pattern = re.compile(r"^([A-Z]+)\t")
    with Path(stats_path).open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            match = pattern.match(line)
            if not match:
                continue
            prefix = match.group(1)
            tokens = line.strip().split("\t")
            parser = SECTION_PARSERS.get(prefix)
            if parser:
                parser(tokens, state)
    return state


def _weighted_percentiles_from_qual(qual_rows, key="snps", probs=(0.025, 0.25, 0.5, 0.75, 0.975)):
    data = []
    for row in qual_rows:
        try:
            q = float(row.get("qual", "nan"))
        except Exception:
            continue
        try:
            c = int(row.get(key, 0))
        except Exception:
            c = 0
        if c <= 0:
            continue
        data.append((q, c))
    if not data:
        return [float("nan")] * len(probs)
    data.sort(key=lambda x: x[0])
    total = sum(c for _, c in data)
    if total <= 0:
        return [float("nan")] * len(probs)
    results = []
    for p in probs:
        target = p * total
        cum = 0
        val = data[-1][0]
        for q, c in data:
            cum += c
            if cum >= target:
                val = q
                break
        results.append(val)
    return results


def _format_pctiles(values):
    if not values:
        return "NaN;NaN;NaN;NaN;NaN"
    out = []
    for v in values:
        try:
            f = float(v)
        except Exception:
            out.append("NaN")
            continue
        if f != f:
            out.append("NaN")
            continue
        if f.is_integer():
            out.append(str(int(f)))
        else:
            out.append(f"{f:.2f}")
    while len(out) < 5:
        out.append("NaN")
    return ";".join(out)


def summarize_metrics(state):
    summary = state.get("summary", {})
    tstv = state.get("tstv", [{}])[0]
    consistent = int(state.get("ics", [{}])[0].get("repeat_consistent", "0"))
    inconsistent = int(state.get("ics", [{}])[0].get("repeat_inconsistent", "0"))
    metrics = {
        "number_of_samples": summary.get("number_of_samples", "NaN"),
        "number_of_records": summary.get("number_of_records", "NaN"),
        "number_of_no_alts": summary.get("number_of_no_alts", "NaN"),
        "number_of_snps": summary.get("number_of_snps", "NaN"),
        "number_of_indels": summary.get("number_of_indels", "NaN"),
        "number_of_multiallelic_sites": summary.get("number_of_multiallelic_sites", "NaN"),
        "number_of_multiallelic_snp_sites": summary.get("number_of_multiallelic_snp_sites", "NaN"),
        "number_of_mnps": summary.get("number_of_mnps", "NaN"),
        "tstv_ts": tstv.get("ts", "NaN"),
        "tstv_tv": tstv.get("tv", "NaN"),
        "tstv_ts_first_alt": tstv.get("ts_first_alt", "NaN"),
        "tstv_tv_first_alt": tstv.get("tv_first_alt", "NaN"),
        "indel_repeat_consistent": consistent,
        "indel_repeat_inconsistent": inconsistent,
    }

    metrics.update(summarize_substitutions(state.get("st", [])))
    metrics.update(summarize_quality_thresholds(state.get("qual", [])))
    metrics.update(summarize_indel_lengths(state.get("icl", [])))

    # Add percentiles from QUAL section (as proxy for GQ percentiles), per SNPs and INDELs
    qual_rows = state.get("qual", [])
    if qual_rows:
        metrics["gq_snps_percentiles"] = _format_pctiles(_weighted_percentiles_from_qual(qual_rows, key="snps"))
        metrics["gq_indels_percentiles"] = _format_pctiles(_weighted_percentiles_from_qual(qual_rows, key="indels"))
    else:
        # empty distributions
        metrics["gq_snps_percentiles"] = _format_pctiles([])
        metrics["gq_indels_percentiles"] = _format_pctiles([])

    return metrics


def summarize_substitutions(rows):
    metrics = {}
    for row in rows:
        subtype = row.get("substitution_type", "").strip()
        try:
            count = int(row.get("count", 0))
        except ValueError:
            count = 0
        if not subtype:
            continue
        key = re.sub(r"[^A-Za-z0-9]+", "_", subtype).strip("_")
        metrics[f"substitution_{key.lower()}"] = count
    return metrics


def summarize_quality_thresholds(rows):
    thresholds = [0, 10, 20, 30, 40, 50]
    metrics = {}
    parsed = []
    for row in rows:
        try:
            qual = float(row.get("qual", "nan"))
        except ValueError:
            continue
        parsed.append(
            (
                qual,
                int(row.get("snps", 0)),
                int(row.get("transitions", 0)),
                int(row.get("transversions", 0)),
                int(row.get("indels", 0)),
            )
        )

    for threshold in thresholds:
        snps = sum(s for q, s, _, _, _ in parsed if q > threshold)
        transitions = sum(tr for q, _, tr, _, _ in parsed if q > threshold)
        transversions = sum(tv for q, _, _, tv, _ in parsed if q > threshold)
        indels = sum(ind for q, _, _, _, ind in parsed if q > threshold)

        prefix = f"qual_gt{threshold}"
        metrics[f"{prefix}_snps"] = snps
        metrics[f"{prefix}_transitions"] = transitions
        metrics[f"{prefix}_transversions"] = transversions
        metrics[f"{prefix}_indels"] = indels

    return metrics


def summarize_indel_lengths(rows):
    categories = {
        "repeat_consistent_del": "sum_length_consistent_del",
        "repeat_inconsistent_del": "sum_length_inconsistent_del",
        "consistent_ins": "sum_length_consistent_ins",
        "inconsistent_ins": "sum_length_inconsistent_ins",
    }

    totals = {key: 0 for key in categories}
    counts = {key: 0 for key in categories}

    for row in rows:
        try:
            length = float(row.get("length", "nan"))
        except ValueError:
            continue
        for source_key in categories:
            try:
                value = int(row.get(source_key, 0))
            except ValueError:
                value = 0
            totals[source_key] += length * value
            counts[source_key] += value

    metrics = {}
    for source_key, metric_key in categories.items():
        metrics[metric_key] = totals[source_key]
        metrics[f"{metric_key}_count"] = counts[source_key]

    return metrics


def write_section(output_dir, name, rows):
    if not rows:
        return None
    output_path = Path(output_dir) / f"{name}.tsv"
    fieldnames = SECTION_HEADERS.get(name)
    if fieldnames is None:
        fieldnames = sorted({k for row in rows for k in row.keys()})
    with output_path.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            cleaned = {k: row.get(k, "") for k in fieldnames}
            writer.writerow(cleaned)
    return output_path


def write_summary(output_path, metrics, append=False):
    out_path = Path(output_path)
    if append:
        exists = out_path.exists()
        size = out_path.stat().st_size if exists else 0
        if exists and size > 0:
            # Reuse existing header to guarantee consistent column order
            with out_path.open("r", newline="") as inp:
                header_line = inp.readline().rstrip("\n")
            header_fields = header_line.split("\t") if header_line else list(metrics.keys())
            with out_path.open("a", newline="") as out:
                writer = csv.DictWriter(out, fieldnames=header_fields, delimiter="\t")
                row = {k: metrics.get(k, "") for k in header_fields}
                writer.writerow(row)
        else:
            # First write: create header from current metrics
            fieldnames = list(metrics.keys())
            with out_path.open("w", newline="") as out:
                writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()
                writer.writerow(metrics)
    else:
        fieldnames = list(metrics.keys())
        with out_path.open("w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerow(metrics)


def main():
    parser = argparse.ArgumentParser(description="Extract metrics from bcftools stats output")
    parser.add_argument("stats", help="Input bcftools stats file")
    parser.add_argument("summary", help="Output summary TSV file")
    parser.add_argument(
        "--section-dir",
        help="Optional directory to write per-section TSV tables",
    )
    parser.add_argument("--sample", help="Sample name to include in summary row")
    parser.add_argument("--region", help="Region label to include in summary row")
    parser.add_argument("--append", action="store_true")
    args = parser.parse_args()

    state = parse_bcftools_stats(args.stats)
    metrics = summarize_metrics(state)
    if args.sample or args.region:
        ordered = {}
        if args.sample:
            ordered["sample"] = args.sample
        if args.region:
            ordered["region"] = args.region
        for k, v in metrics.items():
            ordered[k] = v
        metrics = ordered
    write_summary(args.summary, metrics, append=args.append)

    if args.section_dir:
        out_dir = Path(args.section_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        for key, rows in state.items():
            if key in SECTION_HEADERS:
                write_section(out_dir, key, rows if isinstance(rows, list) else [rows])


if __name__ == "__main__":
    main()
