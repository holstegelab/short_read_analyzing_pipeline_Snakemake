#!/usr/bin/env python3
"""Extract level-2 gVCF regions and package them inside the VCF Conda env."""

import os
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path


def writable_scratch(output_parent: Path) -> Path:
    user = os.environ.get("USER", "")
    job_id = os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_JOBID")
    candidates = []
    if user and job_id:
        candidates.append(Path(f"/scratch-node/{user}.{job_id}"))
    slurm_tmp = os.environ.get("SLURM_TMPDIR")
    if slurm_tmp:
        candidates.append(Path(slurm_tmp))
    candidates.append(output_parent)
    for candidate in candidates:
        if candidate.is_dir() and os.access(candidate, os.W_OK):
            return candidate
    return output_parent


samples = [str(value) for value in snakemake.params.samples]
region = str(snakemake.params.region)
samplefile = str(snakemake.params.samplefile)
dataset = str(snakemake.params.dataset).lower()
interval = Path(str(snakemake.params.interval))
sources_and_indexes = [Path(str(value)) for value in snakemake.input.gvcfs]
output_tar = Path(str(snakemake.output.tar))

if dataset not in {"wes", "wgs"}:
    raise ValueError(f"Unsupported DeepVariant level-2 dataset: {dataset}")
if len(sources_and_indexes) != 2 * len(samples):
    raise ValueError(
        f"Expected a gVCF and index for each of {len(samples)} samples, "
        f"got {len(sources_and_indexes)} inputs"
    )
if not interval.is_file():
    raise FileNotFoundError(f"Level-2 interval file missing: {interval}")

bcftools = shutil.which("bcftools")
if not bcftools:
    raise FileNotFoundError("bcftools is not available in the activated VCF environment")

output_tar.parent.mkdir(parents=True, exist_ok=True)
scratch = writable_scratch(output_tar.parent)
workdir = Path(
    tempfile.mkdtemp(
        prefix=f"dv_l2_{dataset}_{samplefile}_{region}_",
        dir=str(scratch),
    )
)
publish_tmp = output_tar.parent / f".{output_tar.name}.{os.getpid()}.tmp"

try:
    staged_paths = []
    for sample, source_path, source_index in zip(
        samples,
        sources_and_indexes[0::2],
        sources_and_indexes[1::2],
    ):
        if not source_path.is_file():
            raise FileNotFoundError(f"Source gVCF missing: {source_path}")
        if not source_index.is_file():
            raise FileNotFoundError(f"Source gVCF index missing: {source_index}")

        basename = f"{sample}.{region}.dv.{dataset}.g.vcf.gz"
        staged_gvcf = workdir / basename
        staged_tbi = Path(f"{staged_gvcf}.tbi")
        subprocess.run(
            [
                bcftools,
                "view",
                "-R",
                str(interval),
                str(source_path),
                "-O",
                "z",
                "-o",
                str(staged_gvcf),
            ],
            check=True,
        )
        subprocess.run(
            [bcftools, "index", "-f", "-t", str(staged_gvcf)],
            check=True,
        )
        if not staged_gvcf.is_file() or not staged_tbi.is_file():
            raise FileNotFoundError(
                f"Extraction failed for sample {sample}, region {region}"
            )
        staged_paths.extend([staged_gvcf, staged_tbi])

    if not staged_paths:
        raise ValueError(f"No staged gVCFs for {samplefile} {region}")

    staged_tar = workdir / output_tar.name
    with tarfile.open(staged_tar, "w") as tar_handle:
        for path in staged_paths:
            tar_handle.add(path, arcname=path.name)

    expected = {path.name for path in staged_paths}
    with tarfile.open(staged_tar, "r") as tar_handle:
        members = set(tar_handle.getnames())
    missing = sorted(expected - members)
    if missing:
        raise ValueError(f"Tarball missing entries: {missing}")

    shutil.copyfile(staged_tar, publish_tmp)
    os.replace(publish_tmp, output_tar)
finally:
    publish_tmp.unlink(missing_ok=True)
    shutil.rmtree(workdir, ignore_errors=True)
