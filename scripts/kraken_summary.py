#!/usr/bin/env python3
import argparse
import csv
import sys
from collections import defaultdict

CATEGORY_TAXIDS = {
    "viruses": {10239},
    "dna_viruses": {2731341, 2731342, 2732004},
    "herpesvirales": {548681},
    "herpesviridae": {10292},
    "bacteria": {2},
    "fungi": {4751},
    "human": {9606},
}

TARGET_RANK_CODES = {
    "D": "domain",
    "K": "kingdom",
    "P": "phylum",
    "C": "class",
    "O": "order",
    "F": "family",
    "G": "genus",
    "S": "species",
}


def parse_report(report_path):
    rows = []
    with open(report_path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            extra_columns = len(row) - 6
            percent = float(row[0])
            clade_reads = int(row[1])
            taxon_reads = int(row[2])
            rank = row[3 + extra_columns]
            taxid = int(row[4 + extra_columns])
            name = row[5 + extra_columns]
            stripped_name = name.lstrip()
            depth = len(name) - len(stripped_name)
            rows.append(
                {
                    "percent": percent,
                    "clade_reads": clade_reads,
                    "taxon_reads": taxon_reads,
                    "rank": rank,
                    "rank_code": rank[:1],
                    "taxid": taxid,
                    "name": stripped_name,
                    "depth": depth,
                }
            )
    return rows


def parse_merge_stats(paths):
    totals = defaultdict(int)
    for path in paths:
        if not path:
            continue
        with open(path) as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("adapter_length"):
                    break
                parts = stripped.split()
                if len(parts) < 2:
                    continue
                key, value = parts[0], parts[1]
                try:
                    totals[key] += int(value)
                except ValueError:
                    continue
    return totals


def parse_bracken(path, min_fraction=0.01):
    if not path:
        return []
    entries = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            try:
                fraction = float(row.get("fraction_total_reads", "0"))
            except ValueError:
                continue
            if fraction < min_fraction:
                continue
            if row.get("taxonomy_lvl") != "S":
                continue
            name = (row.get("name") or "").strip()
            if not name:
                continue
            entries.append((name, fraction))
    return entries


def summarize(rows, sample_id, merge_totals):
    if not rows:
        raise ValueError("Kraken report is empty")

    root_classified_reads = None
    unclassified_reads = 0
    unclassified_percent = None
    category_counts = defaultdict(int)
    best_by_rank = {}
    best_species = None
    best_genus = None

    for entry in rows:
        taxid = entry["taxid"]
        rank_code = entry["rank_code"]
        clade_reads = entry["clade_reads"]
        percent = entry["percent"]

        if taxid == 1:
            root_classified_reads = clade_reads
        if taxid == 0:
            unclassified_reads = clade_reads
            unclassified_percent = percent
            continue

        for label, taxids in CATEGORY_TAXIDS.items():
            if taxid in taxids:
                category_counts[label] += clade_reads

        if rank_code in TARGET_RANK_CODES:
            current_best = best_by_rank.get(rank_code)
            if current_best is None or clade_reads > current_best["clade_reads"]:
                best_by_rank[rank_code] = entry

        if rank_code == "S":
            if best_species is None or clade_reads > best_species["clade_reads"]:
                best_species = entry
        elif rank_code == "G":
            if best_genus is None or clade_reads > best_genus["clade_reads"]:
                best_genus = entry

    if root_classified_reads is None:
        raise ValueError("Root node (taxid 1) not found in Kraken report")

    total_reads = root_classified_reads + unclassified_reads

    if unclassified_reads > total_reads and unclassified_percent is not None:
        recomputed = int(round(total_reads * (unclassified_percent / 100.0)))
        print(
            f"[kraken_summary] Warning: unclassified_reads ({unclassified_reads}) > total_reads ({total_reads}). "
            f"Recomputing unclassified_reads from percent ({unclassified_percent}%) -> {recomputed}.",
            file=sys.stderr,
            flush=True,
        )
        unclassified_reads = recomputed

    classified_reads = root_classified_reads
    summary = {
        "sample_id": sample_id,
        "total_reads": total_reads,
        "classified_reads": classified_reads,
        "unclassified_reads": unclassified_reads,
    }

    summary["merge_fragments"] = merge_totals.get("fragments", 0)
    summary["merge_primary_reads"] = merge_totals.get("primary_reads", 0)
    summary["merge_chrEBV_fragments"] = merge_totals.get("chrEBV_mapped_fragments", 0)

    def pct(value):
        if total_reads == 0:
            return 0.0
        return round(100.0 * value / total_reads, 6)

    summary["classified_pct"] = pct(classified_reads)

    for label in ("viruses", "dna_viruses", "herpesvirales", "herpesviridae", "bacteria", "fungi", "human"):
        reads = category_counts.get(label, 0)
        summary[f"{label}_reads"] = reads
        summary[f"{label}_pct"] = pct(reads)

    def add_taxon(prefix, entry):
        if entry is None:
            summary[f"{prefix}_name"] = ""
            summary[f"{prefix}_taxid"] = ""
            summary[f"{prefix}_reads"] = 0
            summary[f"{prefix}_pct"] = 0.0
            summary[f"{prefix}_rank"] = ""
            return
        summary[f"{prefix}_name"] = entry["name"].strip()
        summary[f"{prefix}_taxid"] = entry["taxid"]
        summary[f"{prefix}_reads"] = entry["clade_reads"]
        summary[f"{prefix}_pct"] = pct(entry["clade_reads"])
        summary[f"{prefix}_rank"] = entry["rank"]

    add_taxon("top_genus", best_genus)
    add_taxon("top_species", best_species)

    return summary


def write_summary(output_path, summary):
    field_order = [
        "sample_id",
        "total_reads",
        "classified_reads",
        "unclassified_reads",
        "classified_pct",
        "merge_fragments",
        "merge_primary_reads",
        "merge_chrEBV_fragments",
        "viruses_reads",
        "viruses_pct",
        "dna_viruses_reads",
        "dna_viruses_pct",
        "herpesvirales_reads",
        "herpesvirales_pct",
        "herpesviridae_reads",
        "herpesviridae_pct",
        "bacteria_reads",
        "bacteria_pct",
        "fungi_reads",
        "fungi_pct",
        "human_reads",
        "human_pct",
        "top_genus_name",
        "top_genus_taxid",
        "top_genus_rank",
        "top_genus_reads",
        "top_genus_pct",
        "top_species_name",
        "top_species_taxid",
        "top_species_rank",
        "top_species_reads",
        "top_species_pct",
        "bracken_species_ge1pct",
    ]

    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(field_order)
        writer.writerow([summary[field] for field in field_order])


def main():
    parser = argparse.ArgumentParser(description="Summarize Kraken2 report into key taxon metrics")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--report", required=True, help="Path to Kraken report.tsv")
    parser.add_argument("--merge-stats", default="", help="Comma-separated list of merge_stats.tsv files for the sample")
    parser.add_argument("--bracken", default="", help="Optional Bracken report TSV for species abundances")
    parser.add_argument("--output", required=True, help="Path to write the summary TSV")
    args = parser.parse_args()

    rows = parse_report(args.report)
    merge_paths = [path for path in args.merge_stats.split(",") if path]
    merge_totals = parse_merge_stats(merge_paths)
    summary = summarize(rows, args.sample, merge_totals)
    bracken_entries = parse_bracken(args.bracken) if args.bracken else []
    summary["bracken_species_ge1pct"] = ";".join(
        f"{name}:{fraction * 100:.2f}%" for name, fraction in bracken_entries
    )
    write_summary(args.output, summary)


if __name__ == "__main__":
    main()
