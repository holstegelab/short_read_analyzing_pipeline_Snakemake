#!/usr/bin/env python3
"""Extract one level-2 region from existing DeepVariant gVCFs and tar it."""

import argparse
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path


PIPELINE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPELINE_ROOT))

from common import (
    DEEPVARIANT,
    SAMPLEFILE_TO_SAMPLES,
    convert_to_level0,
    convert_to_level1,
    node_ssd_base,
    pj,
    region_to_file,
)


def cohort(samplefile):
    key = os.path.basename(samplefile)
    if key.endswith(".tsv"):
        key = key[:-4]
    if key not in SAMPLEFILE_TO_SAMPLES:
        raise KeyError(f"Unknown samplefile: {samplefile}")
    return key, SAMPLEFILE_TO_SAMPLES[key]


def require_source(path):
    index = Path(str(path) + ".tbi")
    missing = [
        str(candidate)
        for candidate in (path, index)
        if not candidate.is_file() or candidate.stat().st_size == 0
    ]
    if missing:
        raise FileNotFoundError(f"Existing gVCF and index required: {missing}")


def source_path(sample, sinfo, parent, kind):
    if kind == "wgs":
        if "wgs" not in sinfo["sample_type"]:
            return None
        return Path(
            pj(
                DEEPVARIANT,
                "gVCF",
                parent,
                f"{sample}.{parent}.wg.vcf.gz",
            )
        )
    if "wgs" in sinfo["sample_type"]:
        return Path(
            pj(
                DEEPVARIANT,
                "gVCF",
                "exome_extract",
                parent,
                f"{sample}.{parent}.wg.vcf.gz",
            )
        )
    source_region = convert_to_level0(parent)
    return Path(
        pj(
            DEEPVARIANT,
            "gVCF",
            source_region,
            f"{sample}.{source_region}.wg.vcf.gz",
        )
    )


def package(samplefile, region, kind, output):
    key, sampleinfo = cohort(samplefile)
    parent = convert_to_level1(region)
    interval = Path(region_to_file(region, wgs=(kind == "wgs"), extension="bed"))
    if not interval.is_file() or interval.stat().st_size == 0:
        raise FileNotFoundError(f"Non-empty interval file required: {interval}")

    selected = [
        (sample, source_path(sample, sinfo, parent, kind))
        for sample, sinfo in sorted(sampleinfo.items())
    ]
    selected = [(sample, path) for sample, path in selected if path is not None]
    if not selected:
        raise ValueError(f"No {kind.upper()} samples for {key} region {region}")
    for _, path in selected:
        require_source(path)

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp_root = Path(node_ssd_base(str(output.parent)))
    tmp_root.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(
            prefix=f"dv_existing_{kind}_{key}_{region}_",
            dir=str(tmp_root),
        )
    )
    partial = output.with_name(f".{output.name}.partial.{os.getpid()}")

    try:
        staged = []
        for sample, source in selected:
            stem = f"{sample}.{region}.dv.{kind}.g.vcf.gz"
            staged_gvcf = stage / stem
            subprocess.run(
                [
                    "bcftools",
                    "view",
                    "-R",
                    str(interval),
                    str(source),
                    "-O",
                    "z",
                    "-o",
                    str(staged_gvcf),
                ],
                check=True,
            )
            subprocess.run(
                ["bcftools", "index", "-f", "-t", str(staged_gvcf)],
                check=True,
            )
            staged_index = Path(str(staged_gvcf) + ".tbi")
            if (
                not staged_gvcf.is_file()
                or staged_gvcf.stat().st_size == 0
                or not staged_index.is_file()
                or staged_index.stat().st_size == 0
            ):
                raise FileNotFoundError(
                    f"Extraction failed for {sample} {kind} region {region}"
                )
            staged.extend((staged_gvcf, staged_index))

        with tarfile.open(partial, "w") as handle:
            for path in staged:
                handle.add(path, arcname=path.name)

        expected = {path.name for path in staged}
        with tarfile.open(partial, "r") as handle:
            observed = set(handle.getnames())
        missing = sorted(expected - observed)
        if missing:
            raise ValueError(
                f"{kind.upper()} tarball missing {len(missing)} entries: {missing[:5]}"
            )
        if partial.stat().st_size == 0:
            raise ValueError(f"Empty tarball produced: {partial}")
        os.replace(partial, output)
    finally:
        shutil.rmtree(stage, ignore_errors=True)
        try:
            partial.unlink()
        except FileNotFoundError:
            pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", required=True, choices=("wgs", "wes"))
    parser.add_argument("--samplefile", required=True)
    parser.add_argument("--region", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    package(args.samplefile, args.region, args.kind, args.output)


if __name__ == "__main__":
    main()
