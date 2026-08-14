#!/usr/bin/env python3
"""Update CRAM reference and GiB filesize fields from an rclone size listing."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path


def load_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for lineno, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not row:
                continue
            if len(row) != 2:
                raise ValueError(f"{path}:{lineno}: expected <bytes><tab><path>")
            size_text, remote_path = row
            normalized = remote_path.strip().lstrip("/")
            if normalized in sizes:
                raise ValueError(f"{path}:{lineno}: duplicate remote path")
            sizes[normalized] = int(size_text)
    return sizes


def update_config(value: str, size_gib: float) -> str:
    entries = []
    found = False
    for entry in value.split(","):
        entry = entry.strip()
        if not entry:
            continue
        if entry.startswith("filesize="):
            entries.append(f"filesize={size_gib:.6f}")
            found = True
        else:
            entries.append(entry)
    if not found:
        entries.append(f"filesize={size_gib:.6f}")
    return ",".join(entries)


def update_listing(listing: Path, sizes_path: Path, reference: Path) -> tuple[int, int]:
    sizes = load_sizes(sizes_path)
    rows = []
    used = set()

    with listing.open("r", encoding="utf-8", newline="") as handle:
        for lineno, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(row) != 9:
                raise ValueError(f"{listing}:{lineno}: expected 9 fields, got {len(row)}")
            remote_path = row[6].strip().lstrip("/")
            if remote_path not in sizes:
                raise ValueError(f"{listing}:{lineno}: remote CRAM missing from size listing")
            used.add(remote_path)
            row[7] = str(reference)
            row[8] = update_config(row[8], sizes[remote_path] / (1024 ** 3))
            rows.append(row)

    mode = listing.stat().st_mode
    fd, temp_name = tempfile.mkstemp(prefix=f".{listing.name}.", suffix=".tmp", dir=listing.parent)
    temp_path = Path(temp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerows(rows)
        os.chmod(temp_path, mode)
        os.replace(temp_path, listing)
    except Exception:
        temp_path.unlink(missing_ok=True)
        raise

    return len(rows), len(set(sizes) - used)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("listing", type=Path)
    parser.add_argument("sizes", type=Path)
    parser.add_argument("--reference", type=Path, required=True)
    args = parser.parse_args()

    listing = args.listing.expanduser().resolve()
    sizes = args.sizes.expanduser().resolve()
    reference = args.reference.expanduser().resolve()
    if not reference.is_file():
        raise FileNotFoundError(f"reference FASTA does not exist: {reference}")

    updated, unused = update_listing(listing, sizes, reference)
    print(f"updated_rows={updated} unused_remote_crams={unused}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
