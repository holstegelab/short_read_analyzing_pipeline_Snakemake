#!/usr/bin/env python3
"""Extract level-2 gVCF regions and package them inside the VCF Conda env."""

import math
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from pipeline_runtime import assigned_scratch


def writable_scratch(output_parent: Path) -> Path:
    return assigned_scratch(shared_fallback=os.fspath(output_parent))


def log(message):
    print(f"[DeepVariant level2] {message}", file=sys.stderr, flush=True)


def effective_workers(smk):
    """Respect both Snakemake threads and an older, already-granted lease."""
    limits = [getattr(smk, "threads", 1)]
    resources = getattr(smk, "resources", None)
    reserved = getattr(resources, "n", None)
    if reserved is not None:
        limits.append(reserved)
    granted = os.environ.get("ZSLURM_LEASE_MAX_CORES")
    if granted is not None:
        limits.append(granted)
    counts = []
    for value in limits:
        try:
            cores = float(value)
            counts.append(max(1, math.floor(cores)) if math.isfinite(cores) else 1)
        except (TypeError, ValueError):
            counts.append(1)
    return min(counts)


def release_extraction_cores(smk):
    """Best-effort CPU-only shrink once all extraction workers have stopped."""
    command = getattr(smk.params, "lease_command", None)
    if not command or not all(os.environ.get(key) for key in (
        "ZSLURM_LEASE_SOCKET", "ZSLURM_LEASE_TOKEN", "ZSLURM_JOB_ID"
    )):
        return
    try:
        result = subprocess.run(
            [sys.executable, str(command), "--json", "set", "--cores", "1",
             "--phase", "level2_tar_publish"],
            capture_output=True, text=True, timeout=30,
        )
        if result.returncode:
            log("CPU lease shrink unavailable; continuing with existing reservation")
        else:
            log("Released extraction cores; tar/publish phase holds one core")
    except (OSError, subprocess.TimeoutExpired):
        log("CPU lease shrink unavailable; continuing with existing reservation")


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
if len(set(samples)) != len(samples):
    raise ValueError("Duplicate sample names would overwrite extraction outputs")
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

def extract_sample(task):
    sample, source_path, source_index = task
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
    return [staged_gvcf, staged_tbi]


try:
    tasks = list(zip(samples, sources_and_indexes[0::2], sources_and_indexes[1::2]))
    workers = min(effective_workers(snakemake), max(1, len(tasks)))
    log(f"Extracting {len(tasks)} samples with {workers} parallel workers; scratch={workdir}")
    ordered_results = [None] * len(tasks)
    completed = 0
    last_progress = time.monotonic()
    # Only queue one task per worker. On failure, wait for those workers before
    # cleaning scratch; never publish a partial tar or launch the rest of a batch.
    with ThreadPoolExecutor(max_workers=workers) as executor:
        next_task = min(workers, len(tasks))
        pending = {executor.submit(extract_sample, tasks[i]): i for i in range(next_task)}
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                index = pending.pop(future)
                ordered_results[index] = future.result()
                completed += 1
                if completed == len(tasks) or completed % 100 == 0 or time.monotonic() - last_progress >= 60:
                    log(f"Extracted and indexed {completed}/{len(tasks)} samples")
                    last_progress = time.monotonic()
                if next_task < len(tasks):
                    pending[executor.submit(extract_sample, tasks[next_task])] = next_task
                    next_task += 1
    staged_paths = [path for result in ordered_results for path in result]

    if not staged_paths:
        raise ValueError(f"No staged gVCFs for {samplefile} {region}")

    release_extraction_cores(snakemake)
    log(f"Packing {len(staged_paths)} files in original sample order")
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

    log("Publishing tar to shared storage")
    shutil.copyfile(staged_tar, publish_tmp)
    os.replace(publish_tmp, output_tar)
    log("Tar published successfully")
finally:
    publish_tmp.unlink(missing_ok=True)
    shutil.rmtree(workdir, ignore_errors=True)
