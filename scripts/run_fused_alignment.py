#!/usr/bin/env python3
"""Run alignment, merge/dechimer/check, and coordinate sort as one job."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Sequence

from io_profile import run_profiled


from pipeline_runtime import (
    _atomic_json, _atomic_copy, _atomic_publish, assigned_scratch,
    executable, lease_preflight, shrink_lease, acquire_lease,
)


def _q(value: object) -> str:
    return shlex.quote(os.fspath(value))


def _checker_command(args: argparse.Namespace, stats: Path, checked: Path) -> str:
    ignore = " --ignore-qual-checksum-diff" if args.ignore_qual_checksum_diff else ""
    return (
        f"{_q(sys.executable)} {_q(args.bam_stats)} -i - --threads "
        f"{args.low_tool_threads} --fastq-stats {_q(args.fastq_stats)}{ignore} "
        f"-s {_q(stats)} -c {_q(checked)} > /dev/null"
    )


def _checked_pipeline(command: str, fifo: Path, checker: str) -> list[str]:
    script = (
        "set -o pipefail\n"
        f"rm -f {_q(fifo)}\n"
        f"mkfifo {_q(fifo)}\n"
        f"{checker} < {_q(fifo)} &\n"
        "checker_pid=$!\n"
        f"{command}\n"
        "pipeline_rc=$?\n"
        "wait $checker_pid\n"
        "checker_rc=$?\n"
        f"rm -f {_q(fifo)}\n"
        "if [ $pipeline_rc -ne 0 ]; then exit $pipeline_rc; fi\n"
        "exit $checker_rc\n"
    )
    return ["/usr/bin/bash", "-c", script]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared-fastq1", required=True)
    parser.add_argument("--prepared-fastq2", required=True)
    parser.add_argument("--source-fastq1", required=True)
    parser.add_argument("--source-fastq2", required=True)
    parser.add_argument("--fastq-stats", required=True)
    parser.add_argument("--reference-dir", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--readgroup", required=True)
    parser.add_argument("--output-bam", required=True)
    parser.add_argument("--output-bai", required=True)
    parser.add_argument("--dragmap-log", required=True)
    parser.add_argument("--dechimer-stats", required=True)
    parser.add_argument("--badmap-fastq1", required=True)
    parser.add_argument("--badmap-fastq2", required=True)
    parser.add_argument("--merge-stats", required=True)
    parser.add_argument("--checked", required=True)
    parser.add_argument("--check-stats", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--runner-log-label", default="align_reads_fused")
    parser.add_argument("--dragen", default="dragen-os")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--bam-merge", required=True)
    parser.add_argument("--dechimer", required=True)
    parser.add_argument("--bam-stats", required=True)
    parser.add_argument("--dechimer-threshold", type=float, default=0.005)
    parser.add_argument("--align-threads", type=int, default=24)
    parser.add_argument("--alignment-view-threads", type=int, default=2)
    parser.add_argument("--low-tool-threads", type=int, default=1)
    parser.add_argument("--initial-cores", type=float, required=True)
    parser.add_argument("--initial-memory-mb", type=float, required=True)
    parser.add_argument("--low-cores", type=float, default=6.0)
    parser.add_argument("--low-memory-mb", type=float, required=True)
    parser.add_argument("--sort-threads", type=int, default=2)
    parser.add_argument("--sort-memory-mb", type=int, default=6000)
    parser.add_argument("--sort-compression-level", type=int, default=1)
    parser.add_argument("--lease-mode", choices=("required", "optional", "disabled"), default="required")
    parser.add_argument("--lease-command", default="zslurm_lease")
    parser.add_argument("--ssd-gb", type=float, required=True)
    parser.add_argument("--scratch-base", help="Explicit scratch root for tests")
    parser.add_argument("--poll-interval", type=float, default=5.0)
    parser.add_argument("--ignore-qual-checksum-diff", action="store_true")
    return parser.parse_args()


def _require_inputs(paths: Sequence[str]) -> None:
    for raw in paths:
        path = Path(raw)
        if not path.is_file():
            raise FileNotFoundError(path)


def _load_metrics(paths: Sequence[Path]) -> list[dict[str, Any]]:
    result = []
    for path in paths:
        try:
            with path.open("rt", encoding="utf-8") as handle:
                result.append(json.load(handle))
        except FileNotFoundError:
            pass
    return result


def main() -> int:
    args = parse_args()
    _require_inputs([
        args.prepared_fastq1, args.prepared_fastq2, args.source_fastq1,
        args.source_fastq2, args.fastq_stats,
    ])
    dragen = executable(args.dragen)
    samtools = executable(args.samtools)
    bam_merge = executable(args.bam_merge)
    dechimer = executable(args.dechimer)
    if not Path(args.bam_stats).is_file():
        raise FileNotFoundError(args.bam_stats)

    scratch = assigned_scratch(args.scratch_base)
    scratch_parent = scratch / "aligner_fused"
    scratch_parent.mkdir(parents=True, exist_ok=True)
    unique = f"{args.sample}.{args.readgroup}."
    job_tmp = Path(tempfile.mkdtemp(prefix=unique, dir=scratch_parent))
    metrics_output = Path(args.metrics)
    phase_metric_paths: list[Path] = []
    lease: dict[str, Any] = {}
    success = False
    started = time.time()
    local_dragmap_log = job_tmp / "dragmap.log"

    try:
        lease = lease_preflight(
            args.lease_mode,
            args.lease_command,
            initial_cores=args.initial_cores,
            initial_memory_mb=args.initial_memory_mb,
        )

        aligned = job_tmp / "aligned.bam"
        alignment_metrics = job_tmp / "alignment.io.json"
        phase_metric_paths.append(alignment_metrics)
        alignment = (
            "set -o pipefail; "
            f"({_q(dragen)} -r {_q(args.reference_dir)} "
            f"-1 {_q(args.prepared_fastq1)} -2 {_q(args.prepared_fastq2)} "
            f"--RGID {_q(args.readgroup)} --RGSM {_q(args.sample)} "
            f"--num-threads {args.align_threads} | "
            f"{_q(samtools)} view -@ {args.alignment_view_threads} "
            f"-o {_q(aligned)} -) 2> {_q(local_dragmap_log)}"
        )
        run_profiled(
            ["/usr/bin/bash", "-c", alignment],
            metrics_path=alignment_metrics,
            label="align_reads_fused.alignment",
            local_paths=[job_tmp],
            input_paths=[args.prepared_fastq1, args.prepared_fastq2],
            output_paths=[aligned, local_dragmap_log],
            requested_ssd_gb=args.ssd_gb,
            threads=args.align_threads,
            memory_mb=int(args.initial_memory_mb),
            poll_interval=args.poll_interval,
        )

        lease = shrink_lease(
            lease,
            cores=args.low_cores,
            memory_mb=args.low_memory_mb,
            phase="alignment_tail",
        )

        merge_stats = job_tmp / "merge_stats.tsv"
        check_stats = job_tmp / "bam_check_stats.tsv"
        checked = job_tmp / "bam_checked"
        badmap1 = job_tmp / "badmap_R1.fastq.gz"
        badmap2 = job_tmp / "badmap_R2.fastq.gz"
        merged = job_tmp / "merged.bam"
        fifo1 = job_tmp / "check-stage1.fifo"
        merge_prefix = _q(bam_merge)
        stage1 = (
            f"{_q(samtools)} view -h --threads {args.low_tool_threads} {_q(aligned)} "
            f"| {merge_prefix} -a {_q(args.source_fastq1)} -b {_q(args.source_fastq2)} "
            f"-ua {_q(badmap1)} -ub {_q(badmap2)} -s {_q(merge_stats)} "
            f"| tee {_q(fifo1)} "
            f"| {_q(samtools)} fixmate -@ {args.low_tool_threads} -u -O BAM -m - {_q(merged)}"
        )
        stage1_metrics = job_tmp / "merge_check.io.json"
        phase_metric_paths.append(stage1_metrics)
        run_profiled(
            _checked_pipeline(
                stage1, fifo1, _checker_command(args, check_stats, checked)
            ),
            metrics_path=stage1_metrics,
            label="align_reads_fused.merge_check",
            local_paths=[job_tmp],
            input_paths=[aligned, args.source_fastq1, args.source_fastq2],
            output_paths=[merged, merge_stats, check_stats, checked, badmap1, badmap2],
            requested_ssd_gb=args.ssd_gb,
            threads=int(args.low_cores),
            memory_mb=int(args.low_memory_mb),
            poll_interval=args.poll_interval,
        )
        aligned.unlink()

        ratio = 0.0
        with merge_stats.open("rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("primary_soft_clipped_bp_ratio"):
                    try:
                        ratio = float(line.split("\t", 1)[1].strip())
                    except (IndexError, ValueError):
                        ratio = 0.0
                    break

        dechimer_stats = job_tmp / "dechimer_stats.tsv"
        final_local = merged
        if ratio > args.dechimer_threshold:
            final_local = job_tmp / "dechimer.bam"
            fifo2 = job_tmp / "check-stage2.fifo"
            stage2 = (
                f"{_q(samtools)} view -h --threads {args.low_tool_threads} {_q(merged)} "
                f"| {_q(dechimer)} --min_align_length 40 --loose_ends -i - "
                f"-s {_q(dechimer_stats)} | tee {_q(fifo2)} "
                f"| {_q(samtools)} fixmate -@ {args.low_tool_threads} -u -O BAM -m - {_q(final_local)}"
            )
            stage2_metrics = job_tmp / "dechimer_check.io.json"
            phase_metric_paths.append(stage2_metrics)
            run_profiled(
                _checked_pipeline(
                    stage2, fifo2, _checker_command(args, check_stats, checked)
                ),
                metrics_path=stage2_metrics,
                label="align_reads_fused.dechimer_check",
                local_paths=[job_tmp],
                input_paths=[merged],
                output_paths=[final_local, dechimer_stats, check_stats, checked],
                requested_ssd_gb=args.ssd_gb,
                threads=int(args.low_cores),
                memory_mb=int(args.low_memory_mb),
                poll_interval=args.poll_interval,
            )
            merged.unlink()
        else:
            dechimer_stats.touch()

        required_nonempty = [final_local, merge_stats, check_stats, checked]
        for path in required_nonempty:
            if not path.is_file() or path.stat().st_size <= 0:
                raise RuntimeError(f"required fused output is absent or empty: {path}")
        for path in (badmap1, badmap2, dechimer_stats, local_dragmap_log):
            if not path.is_file():
                raise RuntimeError(f"required fused output is absent: {path}")

        sorted_local = job_tmp / "sorted.bam"
        sorted_index_local = job_tmp / "sorted.bam.bai"
        sort_metrics = job_tmp / "sort.io.json"
        phase_metric_paths.append(sort_metrics)
        run_profiled(
            [
                samtools,
                "sort",
                "-T",
                str(job_tmp / "sort"),
                "-@",
                str(args.sort_threads),
                "-l",
                str(args.sort_compression_level),
                "-m",
                f"{args.sort_memory_mb}M",
                "--write-index",
                "-o",
                f"{sorted_local}##idx##{sorted_index_local}",
                str(final_local),
            ],
            metrics_path=sort_metrics,
            label="align_reads_fused.sort",
            local_paths=[job_tmp],
            input_paths=[final_local],
            output_paths=[sorted_local, sorted_index_local],
            requested_ssd_gb=args.ssd_gb,
            threads=args.sort_threads,
            memory_mb=args.sort_memory_mb,
            poll_interval=args.poll_interval,
        )
        for path in (sorted_local, sorted_index_local):
            if not path.is_file() or path.stat().st_size <= 0:
                raise RuntimeError(f"required fused sort output is absent or empty: {path}")

        copies = [
            (local_dragmap_log, Path(args.dragmap_log)),
            (dechimer_stats, Path(args.dechimer_stats)),
            (badmap1, Path(args.badmap_fastq1)),
            (badmap2, Path(args.badmap_fastq2)),
            (merge_stats, Path(args.merge_stats)),
            (checked, Path(args.checked)),
            (check_stats, Path(args.check_stats)),
            (sorted_index_local, Path(args.output_bai)),
            # Publish the large terminal BAM last so its appearance remains
            # the final completion boundary for downstream consumers.
            (sorted_local, Path(args.output_bam)),
        ]
        for source, destination in copies:
            _atomic_copy(source, destination)
        success = True
        return 0
    finally:
        if local_dragmap_log.is_file() and not Path(args.dragmap_log).is_file():
            try:
                _atomic_copy(local_dragmap_log, Path(args.dragmap_log))
            except Exception:
                pass
        phases = _load_metrics(phase_metric_paths)
        shutil.rmtree(job_tmp, ignore_errors=True)
        payload = {
            "schema_version": 1,
            "label": args.runner_log_label,
            "sample": args.sample,
            "readgroup": args.readgroup,
            "success": success,
            "started_at_epoch": started,
            "duration_seconds": round(time.time() - started, 6),
            "requested": {
                "initial_cores": args.initial_cores,
                "initial_memory_mb": args.initial_memory_mb,
                "low_cores": args.low_cores,
                "low_memory_mb": args.low_memory_mb,
                "sort_threads": args.sort_threads,
                "sort_memory_mb": args.sort_memory_mb,
                "ssd_gb": args.ssd_gb,
            },
            "lease": lease,
            "phases": phases,
            "scratch_job_directory": str(job_tmp),
            "scratch_removed": not job_tmp.exists(),
        }
        _atomic_json(metrics_output, payload)


if __name__ == "__main__":
    raise SystemExit(main())
