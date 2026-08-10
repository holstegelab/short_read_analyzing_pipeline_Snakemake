#!/usr/bin/env python3
"""Run samtools sort with temporary files on the assigned node-local SSD."""

from __future__ import annotations

import argparse
import os
import shutil
import tempfile
from pathlib import Path

from io_profile import run_profiled


def assigned_scratch(explicit: str | None = None) -> Path:
    """Resolve this job's scratch root; never borrow another job's directory."""
    if explicit:
        candidate: Path | None = Path(explicit)
    else:
        user = os.environ.get("USER")
        job_id = os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_JOBID")
        candidate = Path(f"/scratch-node/{user}.{job_id}") if user and job_id else None
        if candidate is None or not candidate.is_dir():
            slurm_tmp = os.environ.get("SLURM_TMPDIR")
            if slurm_tmp and Path(slurm_tmp).resolve().is_relative_to(Path("/scratch-node")):
                candidate = Path(slurm_tmp)
    if candidate is None or not candidate.is_dir() or not os.access(candidate, os.W_OK):
        raise RuntimeError(
            "ssd_use=required but this job has no writable assigned /scratch-node directory"
        )
    return candidate.resolve()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-bam", required=True)
    parser.add_argument("--output-bai", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--label", default="sort_bam_alignment")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--memory-mb", type=int, default=6000)
    parser.add_argument("--compression-level", type=int, default=1)
    parser.add_argument("--ssd-gb", type=float, required=True)
    parser.add_argument("--poll-interval", type=float, default=5.0)
    parser.add_argument("--scratch-base", help="Explicit scratch root for tests")
    parser.add_argument("--samtools", default="samtools")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_bam = Path(args.input)
    output_bam = Path(args.output_bam)
    output_bai = Path(args.output_bai)
    if not input_bam.is_file():
        raise FileNotFoundError(input_bam)
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    output_bai.parent.mkdir(parents=True, exist_ok=True)
    metrics = Path(args.metrics)
    metrics.parent.mkdir(parents=True, exist_ok=True)

    samtools = shutil.which(args.samtools) if os.sep not in args.samtools else args.samtools
    if not samtools or not Path(samtools).is_file():
        raise FileNotFoundError(f"samtools executable not found: {args.samtools}")

    scratch = assigned_scratch(args.scratch_base)
    parent = scratch / "aligner_sort"
    parent.mkdir(parents=True, exist_ok=True)
    job_tmp = Path(tempfile.mkdtemp(prefix=f"{output_bam.stem}.", dir=parent))
    temp_prefix = job_tmp / "sort"
    command = [
        str(samtools),
        "sort",
        "-T",
        str(temp_prefix),
        "-@",
        str(args.threads),
        "-l",
        str(args.compression_level),
        "-m",
        f"{args.memory_mb}M",
        "--write-index",
        "-o",
        f"{output_bam}##idx##{output_bai}",
        str(input_bam),
    ]
    return run_profiled(
        command,
        metrics_path=metrics,
        label=args.label,
        local_paths=[job_tmp],
        input_paths=[input_bam],
        output_paths=[output_bam, output_bai],
        cleanup_paths=[job_tmp],
        requested_ssd_gb=args.ssd_gb,
        threads=args.threads,
        memory_mb=args.memory_mb,
        poll_interval=args.poll_interval,
    )


if __name__ == "__main__":
    raise SystemExit(main())
