#!/usr/bin/env python3
"""Merge read-group BAMs and mark duplicates in one node-local SSD job."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Sequence

from io_profile import run_profiled
from pipeline_runtime import (
    _atomic_copy,
    _atomic_json,
    assigned_scratch,
    executable,
    lease_preflight,
    shrink_lease,
)


def _q(value: object) -> str:
    return shlex.quote(os.fspath(value))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-bam", nargs="+", required=True)
    parser.add_argument("--input-bai", nargs="+", required=True)
    parser.add_argument("--check-marker", nargs="+", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-bam", required=True)
    parser.add_argument("--output-bai", required=True)
    parser.add_argument("--output-stat", required=True)
    parser.add_argument("--merge-log", required=True)
    parser.add_argument("--markdup-log", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--no-dedup", type=int, choices=(0, 1), default=0)
    parser.add_argument("--optical-distance", type=int, default=2500)
    parser.add_argument("--merge-threads", type=int, default=3)
    parser.add_argument("--initial-cores", type=float, required=True)
    parser.add_argument("--initial-memory-mb", type=float, required=True)
    parser.add_argument("--markdup-cores", type=float, default=0.95)
    parser.add_argument("--markdup-memory-mb", type=float, required=True)
    parser.add_argument(
        "--lease-mode",
        choices=("required", "optional", "disabled"),
        default="required",
    )
    parser.add_argument("--lease-command", default="zslurm_lease")
    parser.add_argument("--ssd-gb", type=float, required=True)
    parser.add_argument("--scratch-base", help="Explicit scratch root for tests")
    parser.add_argument("--poll-interval", type=float, default=5.0)
    return parser.parse_args()


def _require_files(paths: Sequence[str], *, nonempty: bool) -> None:
    for raw in paths:
        path = Path(raw)
        if not path.is_file():
            raise FileNotFoundError(path)
        if nonempty and path.stat().st_size <= 0:
            raise RuntimeError(f"required input is empty: {path}")


def _load_metrics(paths: Sequence[Path]) -> list[dict]:
    result = []
    for path in paths:
        try:
            result.append(json.loads(path.read_text(encoding="utf-8")))
        except FileNotFoundError:
            pass
    return result


def _logged_command(command: Sequence[str], log_path: Path) -> list[str]:
    rendered = " ".join(_q(value) for value in command)
    return ["/usr/bin/bash", "-c", f"exec {rendered} 2> {_q(log_path)}"]


def main() -> int:
    args = parse_args()
    if not (
        len(args.input_bam)
        == len(args.input_bai)
        == len(args.check_marker)
    ):
        raise ValueError("BAM, BAI and validation-marker counts must match")
    _require_files(args.input_bam, nonempty=True)
    _require_files(args.input_bai, nonempty=True)
    _require_files(args.check_marker, nonempty=False)
    samtools = executable(args.samtools)

    scratch = assigned_scratch(args.scratch_base)
    parent = scratch / "merge_markdup_fused"
    parent.mkdir(parents=True, exist_ok=True)
    job_tmp = Path(tempfile.mkdtemp(prefix=f"{args.sample}.", dir=parent))
    phase_paths: list[Path] = []
    lease: dict = {}
    success = False
    started = time.time()
    local_merge_log = job_tmp / "merge_rgs.log"
    local_markdup_log = job_tmp / "markdup.log"

    try:
        lease = lease_preflight(
            args.lease_mode,
            args.lease_command,
            initial_cores=args.initial_cores,
            initial_memory_mb=args.initial_memory_mb,
        )

        markdup_input = Path(args.input_bam[0])
        if len(args.input_bam) > 1:
            merged = job_tmp / "merged.bam"
            merge_metrics = job_tmp / "merge.io.json"
            phase_paths.append(merge_metrics)
            merge_command = [
                samtools,
                "merge",
                "-@",
                str(args.merge_threads),
                str(merged),
                *args.input_bam,
            ]
            run_profiled(
                _logged_command(merge_command, local_merge_log),
                metrics_path=merge_metrics,
                label="merge_markdup_fused.merge",
                local_paths=[job_tmp],
                input_paths=args.input_bam,
                output_paths=[merged, local_merge_log],
                requested_ssd_gb=args.ssd_gb,
                threads=args.merge_threads,
                memory_mb=int(args.initial_memory_mb),
                poll_interval=args.poll_interval,
            )
            if not merged.is_file() or merged.stat().st_size <= 0:
                raise RuntimeError(f"merged BAM is absent or empty: {merged}")
            markdup_input = merged
            lease = shrink_lease(
                lease,
                cores=args.markdup_cores,
                memory_mb=args.markdup_memory_mb,
                phase="markdup",
            )
        else:
            local_merge_log.touch()

        final_bam = job_tmp / "markdup.bam"
        final_bai = job_tmp / "markdup.bam.bai"
        final_stat = job_tmp / "markdup.stat"
        markdup_metrics = job_tmp / "markdup.io.json"
        phase_paths.append(markdup_metrics)

        if args.no_dedup:
            # Preserve the old no-dedup contract: copy/merge the BAM, rebuild
            # its index and create an empty statistics file.
            if markdup_input.is_relative_to(job_tmp):
                os.replace(markdup_input, final_bam)
            else:
                shutil.copyfile(markdup_input, final_bam)
            final_stat.touch()
            command = [samtools, "index", str(final_bam), str(final_bai)]
            profiled_inputs = [final_bam]
        else:
            command = [
                samtools,
                "markdup",
                "-T",
                str(job_tmp / "markdup-temp"),
                "-f",
                str(final_stat),
                "-S",
                "-d",
                str(args.optical_distance),
                str(markdup_input),
                "--write-index",
                f"{final_bam}##idx##{final_bai}",
            ]
            profiled_inputs = [markdup_input]

        run_profiled(
            _logged_command(command, local_markdup_log),
            metrics_path=markdup_metrics,
            label="merge_markdup_fused.markdup",
            local_paths=[job_tmp],
            input_paths=profiled_inputs,
            output_paths=[final_bam, final_bai, final_stat, local_markdup_log],
            requested_ssd_gb=args.ssd_gb,
            threads=1,
            memory_mb=int(args.markdup_memory_mb),
            poll_interval=args.poll_interval,
        )
        for path in (final_bam, final_bai):
            if not path.is_file() or path.stat().st_size <= 0:
                raise RuntimeError(f"required fused output is absent or empty: {path}")
        if not final_stat.is_file():
            raise RuntimeError(f"required fused output is absent: {final_stat}")
        subprocess.run([samtools, "quickcheck", str(final_bam)], check=True)

        # Publish the large BAM last. Until it appears, downstream rules cannot
        # observe a partially published output set.
        for source, destination in (
            (final_bai, Path(args.output_bai)),
            (final_stat, Path(args.output_stat)),
            (final_bam, Path(args.output_bam)),
        ):
            _atomic_copy(source, destination)
        success = True
        return 0
    finally:
        for source, destination in (
            (local_merge_log, Path(args.merge_log)),
            (local_markdup_log, Path(args.markdup_log)),
        ):
            if source.is_file():
                try:
                    _atomic_copy(source, destination)
                except Exception:
                    pass
        phases = _load_metrics(phase_paths)
        shutil.rmtree(job_tmp, ignore_errors=True)
        _atomic_json(
            Path(args.metrics),
            {
                "schema_version": 1,
                "label": "merge_markdup_fused",
                "sample": args.sample,
                "readgroups": len(args.input_bam),
                "no_dedup": bool(args.no_dedup),
                "success": success,
                "started_at_epoch": started,
                "duration_seconds": round(time.time() - started, 6),
                "requested": {
                    "initial_cores": args.initial_cores,
                    "initial_memory_mb": args.initial_memory_mb,
                    "markdup_cores": args.markdup_cores,
                    "markdup_memory_mb": args.markdup_memory_mb,
                    "merge_threads": args.merge_threads,
                    "ssd_gb": args.ssd_gb,
                },
                "lease": lease,
                "phases": phases,
                "scratch_job_directory": str(job_tmp),
                "scratch_removed": not job_tmp.exists(),
            },
        )


if __name__ == "__main__":
    raise SystemExit(main())
