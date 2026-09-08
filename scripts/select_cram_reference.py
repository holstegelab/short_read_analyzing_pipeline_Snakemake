#!/usr/bin/env python3
"""Select a known CRAM decode reference from the header sequence dictionary."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


HG19_CHR1_M5 = "1b22b98cdeb4a9304cb5d48026a85128"
HG19_CHRY_M5 = "1e86411d73e6f00a10590f976be01623"
HG19_B37_CHRY_M5 = "1fa3474750af0948bdf97d5a0ee52e51"
HG38_CHR1_M5 = "6aef897c3d6ff0c78aff06ac189178dd"
HG38_CHRY_M5 = "ce3e31103314a704255f3cd90369ecce"


def parse_sq_md5(header: str) -> dict[str, str]:
    """Return ``SN`` to lowercase ``M5`` values from SAM ``@SQ`` lines."""
    result: dict[str, str] = {}
    for line in header.splitlines():
        if not line.startswith("@SQ\t"):
            continue
        fields = {}
        for item in line.split("\t")[1:]:
            key, separator, value = item.partition(":")
            if separator:
                fields[key] = value
        name = fields.get("SN")
        checksum = fields.get("M5")
        if name and checksum:
            result[name] = checksum.lower()
    return result


def choose_reference(
    header: str,
    *,
    primary_reference: str,
    hg19_reference: str,
    hg19_b37_chry_reference: str,
    hg38_reference: str,
) -> tuple[str, str]:
    """Choose a reference and return it together with the selection reason."""
    checksums = parse_sq_md5(header)
    chr1 = checksums.get("chr1")
    chry = checksums.get("chrY")

    if chr1 == HG19_CHR1_M5 and chry == HG19_B37_CHRY_M5:
        return hg19_b37_chry_reference, "hg19+b37-chrY M5 signature"
    if chr1 == HG19_CHR1_M5 and chry in {None, HG19_CHRY_M5}:
        return hg19_reference, "hg19 M5 signature"
    if chr1 == HG38_CHR1_M5 and chry in {None, HG38_CHRY_M5}:
        return hg38_reference, "hg38 M5 signature"
    return primary_reference, "unknown M5 signature; keeping configured reference"


def read_cram_header(cram: str, samtools: str) -> str:
    process = subprocess.run(
        [samtools, "view", "-H", cram],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if process.returncode != 0:
        detail = process.stderr.strip() or "no error text"
        raise RuntimeError(
            f"failed to read CRAM header for reference selection: {detail}"
        )
    return process.stdout


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cram", required=True)
    parser.add_argument("--primary-reference", required=True)
    parser.add_argument("--hg19-reference", required=True)
    parser.add_argument("--hg19-b37-chry-reference", required=True)
    parser.add_argument("--hg38-reference", required=True)
    parser.add_argument("--samtools", default="samtools")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    header = read_cram_header(args.cram, args.samtools)
    selected, reason = choose_reference(
        header,
        primary_reference=args.primary_reference,
        hg19_reference=args.hg19_reference,
        hg19_b37_chry_reference=args.hg19_b37_chry_reference,
        hg38_reference=args.hg38_reference,
    )
    if not Path(selected).is_file():
        raise FileNotFoundError(
            f"selected CRAM reference does not exist ({reason}): {selected}"
        )
    print(selected)
    print(f"[split reference] {reason}: {selected}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
