#!/usr/bin/env python3
"""Recover read groups from Illumina QNAMEs in an RG-less CRAM.

Discovery runs in the read-group checkpoint.  The split step streams SAM from
the original CRAM to ``samtools split``; it never writes a full rescued copy.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path


_ILLUMINA_QNAME = re.compile(
    r"^(?P<instrument>[A-Za-z0-9_.-]+):(?P<run>[0-9]+):"
    r"(?P<flowcell>[A-Za-z0-9_.-]+):(?P<lane>[0-9]+):"
    r"[0-9]+:-?[0-9]+:-?[0-9]+(?:[:][A-Za-z0-9+._-]+)?$"
)
_SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")


def qname_group(qname: str) -> tuple[str, str, str, str]:
    """Return instrument, run, flowcell and lane; reject ambiguous names."""
    match = _ILLUMINA_QNAME.fullmatch(qname)
    if not match:
        raise ValueError(
            f"Cannot recover flowcell/lane from Illumina QNAME {qname!r}"
        )
    return tuple(match.group(name) for name in (
        "instrument", "run", "flowcell", "lane"
    ))


def readgroup_id(sample: str, group: tuple[str, str, str, str]) -> str:
    if not _SAFE_ID.fullmatch(sample):
        raise ValueError(f"Sample name is unsafe for a read-group ID: {sample!r}")
    instrument, run, flowcell, lane = group
    return f"{sample}_{instrument}_{run}_{flowcell}_L{lane}"


def _samtools_lines(command: list[str]):
    """Yield SAM lines while retaining the real samtools error on failure."""
    with tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as errors:
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=errors,
            text=True, encoding="utf-8",
        )
        assert process.stdout is not None
        try:
            yield from process.stdout
        except BaseException:
            process.kill()
            raise
        finally:
            process.stdout.close()
            status = process.wait()
        if status:
            errors.seek(0)
            raise RuntimeError(
                f"samtools exited {status}: {errors.read().strip()}"
            )


def discover_readgroups(
    input_cram: str, reference: str, samtools: str = "samtools"
) -> dict[tuple[str, str, str, str], int]:
    """Inspect every record; a first-read sample cannot prove lane purity."""
    counts: Counter[tuple[str, str, str, str]] = Counter()
    for line in _samtools_lines(
        [samtools, "view", "--reference", reference, input_cram]
    ):
        qname = line.split("\t", 1)[0]
        counts[qname_group(qname)] += 1
    if not counts:
        raise ValueError(f"RG-less CRAM contains no reads: {input_cram}")
    return dict(sorted(counts.items()))


def _rg_line(group: dict) -> str:
    info = group["info"]
    fields = ("ID", "SM", "LB", "PL")
    if not all(info.get(field) for field in fields):
        raise ValueError(f"Incomplete rescued read group: {info!r}")
    if any("\t" in str(info[field]) or "\n" in str(info[field]) for field in fields):
        raise ValueError(f"Unsafe rescued read-group field: {info!r}")
    return "@RG\t" + "\t".join(f"{field}:{info[field]}" for field in fields) + "\n"


def split_rescued_cram(
    input_cram: str, output_dir: str, sample: str, groups: list[dict],
    reference: str, output_fmt: str, extension: str, threads: int = 0,
    samtools: str = "samtools",
) -> None:
    """Add RG tags while streaming directly to per-RG CRAMs."""
    if not groups:
        raise ValueError("No rescued read groups were supplied")
    if extension not in {"cram", "ucram"}:
        raise ValueError(f"Unexpected CRAM extension: {extension!r}")
    expected = {}
    for group in groups:
        key = tuple(group["source_group"])
        if key in expected or group["info"]["ID"] in expected.values():
            raise ValueError("Duplicate rescued read group")
        expected[key] = group["info"]["ID"]

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    # Hold incomplete files under a private name until the whole stream has
    # been validated.  Snakemake's done marker is created by the caller.
    with tempfile.TemporaryDirectory(prefix=".rg-rescue-", dir=output_path) as work:
        filename_format = str(Path(work) / f"{sample}.%!.{extension}")
        split_command = [
            samtools, "split", "-@", str(threads), "--output-fmt", output_fmt,
            "--reference", reference, "-", "-f", filename_format,
        ]
        observed: Counter[tuple[str, str, str, str]] = Counter()
        with tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as errors:
            splitter = subprocess.Popen(
                split_command, stdin=subprocess.PIPE, stderr=errors,
                text=True, encoding="utf-8",
            )
            assert splitter.stdin is not None
            header_written = False
            try:
                for line in _samtools_lines(
                    [samtools, "view", "-h", "--reference", reference, input_cram]
                ):
                    if line.startswith("@") and not header_written:
                        if line.startswith("@RG\t"):
                            raise ValueError("Source CRAM gained @RG entries after discovery")
                        splitter.stdin.write(line)
                        continue
                    if not header_written:
                        for group in groups:
                            splitter.stdin.write(_rg_line(group))
                        header_written = True
                    if line.startswith("@"):
                        raise ValueError("Header line found after alignment records")
                    qname, separator, _ = line.partition("\t")
                    if not separator:
                        raise ValueError("Malformed SAM alignment during RG rescue")
                    key = qname_group(qname)
                    if key not in expected:
                        raise ValueError(f"Undiscovered flowcell/lane in {qname!r}")
                    if "\tRG:" in line:
                        raise ValueError(f"Source read already has an RG tag: {qname!r}")
                    observed[key] += 1
                    splitter.stdin.write(line.rstrip("\n") + f"\tRG:Z:{expected[key]}\n")
                if not header_written:
                    raise ValueError("Source CRAM has no alignment records")
                for group in groups:
                    key = tuple(group["source_group"])
                    if observed[key] != group["source_read_count"]:
                        raise ValueError(
                            f"Read count changed for {key}: discovered "
                            f"{group['source_read_count']}, observed {observed[key]}"
                        )
            except BaseException:
                splitter.kill()
                raise
            finally:
                try:
                    splitter.stdin.close()
                except BrokenPipeError:
                    pass
                status = splitter.wait()
            if status:
                errors.seek(0)
                raise RuntimeError(
                    f"samtools split exited {status}: {errors.read().strip()}"
                )
        for group in groups:
            name = f"{sample}.{group['info']['ID']}.{extension}"
            source = Path(work) / name
            if not source.is_file():
                raise RuntimeError(f"samtools split did not produce {name}")
        for group in groups:
            name = f"{sample}.{group['info']['ID']}.{extension}"
            os.replace(Path(work) / name, output_path / name)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("split", choices=["split"])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--groups-json", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--output-fmt", required=True)
    parser.add_argument("--extension", choices=("cram", "ucram"), required=True)
    parser.add_argument("--threads", type=int, default=0)
    args = parser.parse_args()
    split_rescued_cram(
        args.input, args.output_dir, args.sample, json.loads(args.groups_json),
        args.reference, args.output_fmt, args.extension, args.threads,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
