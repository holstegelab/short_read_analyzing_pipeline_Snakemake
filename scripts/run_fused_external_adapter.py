#!/usr/bin/env python3
"""Extract FASTQs from one BAM/CRAM read group and remove adapters on SSD."""

from __future__ import annotations

import argparse
import bz2
import fcntl
import gzip
import json
import os
import shlex
import shutil
import sys
import tempfile
import time
from pathlib import Path

from io_profile import run_profiled
from run_fused_alignment import (
    _atomic_copy,
    _atomic_json,
    assigned_scratch,
    executable,
    lease_preflight,
    shrink_lease,
)


def q(value: object) -> str:
    return shlex.quote(str(value))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-alignment", required=True)
    parser.add_argument("--cram-options", default="")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--readgroup", required=True)
    parser.add_argument("--adapter-list", required=True)
    parser.add_argument("--fastq-stats-script", required=True)
    parser.add_argument("--remove-duplicates-script", required=True)
    parser.add_argument("--pair-rescue-script", required=True)
    parser.add_argument("--remove-duplicated-reads", type=int, choices=(0, 1), default=0)
    parser.add_argument("--attempt", type=int, default=1)
    parser.add_argument("--error-file", required=True)
    parser.add_argument("--output-raw-forward", required=True)
    parser.add_argument("--output-raw-reverse", required=True)
    parser.add_argument("--output-singletons", required=True)
    parser.add_argument("--output-forward", required=True)
    parser.add_argument("--output-reverse", required=True)
    parser.add_argument("--output-adapter-log", required=True)
    parser.add_argument("--output-fastq-stats", required=True)
    parser.add_argument("--output-adapters", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--pigz", default="pigz")
    parser.add_argument("--adapter-removal", default="AdapterRemoval")
    parser.add_argument("--initial-cores", type=float, required=True)
    parser.add_argument("--initial-memory-mb", type=float, required=True)
    parser.add_argument("--adapter-cores", type=float, default=5.0)
    parser.add_argument("--adapter-memory-mb", type=float, default=1024)
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


def quality_base(path: Path, max_records: int = 1000) -> tuple[int, int]:
    if path.name.endswith(".gz"):
        handle = gzip.open(path, "rt", encoding="utf-8", errors="replace")
    elif path.name.endswith(".bz2"):
        handle = bz2.open(path, "rt", encoding="utf-8", errors="replace")
    else:
        handle = path.open("rt", encoding="utf-8", errors="replace")
    qmin = 10**9
    qmax = -1
    seen = 0
    try:
        while seen < max_records:
            header = handle.readline()
            if not header:
                break
            handle.readline()
            handle.readline()
            quality = handle.readline()
            if not quality:
                break
            for char in quality.rstrip("\r\n"):
                value = ord(char)
                qmin = min(qmin, value)
                qmax = max(qmax, value)
            seen += 1
    finally:
        handle.close()
    if qmax < 0 or qmin > qmax:
        return 33, 42
    base = 64 if qmin >= 64 else 33
    return base, max(0, qmax - base)


def record_rescue(error_file: Path, sample: str) -> None:
    error_file.parent.mkdir(parents=True, exist_ok=True)
    with error_file.open("a+", encoding="utf-8") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        handle.seek(0)
        existing = {line.strip() for line in handle if line.strip()}
        if sample not in existing:
            handle.write(sample + "\n")
            handle.flush()
        fcntl.flock(handle, fcntl.LOCK_UN)


def load_metrics(paths: list[Path]) -> list[dict]:
    result = []
    for path in paths:
        try:
            result.append(json.loads(path.read_text(encoding="utf-8")))
        except FileNotFoundError:
            pass
    return result


def main() -> int:
    args = parse_args()
    alignment = Path(args.input_alignment)
    if not alignment.is_file():
        raise FileNotFoundError(alignment)
    for raw in (
        args.adapter_list,
        args.fastq_stats_script,
        args.remove_duplicates_script,
        args.pair_rescue_script,
    ):
        if not Path(raw).is_file():
            raise FileNotFoundError(raw)
    samtools = executable(args.samtools)
    pigz = executable(args.pigz)
    adapter_removal = executable(args.adapter_removal)

    scratch = assigned_scratch(args.scratch_base)
    parent = scratch / "external_adapter_fused"
    parent.mkdir(parents=True, exist_ok=True)
    job_tmp = Path(
        tempfile.mkdtemp(prefix=f"{args.sample}.{args.readgroup}.", dir=parent)
    )
    phase_paths: list[Path] = []
    lease: dict = {}
    success = False
    started = time.time()
    raw1 = job_tmp / "raw_R1.fastq.gz"
    raw2 = job_tmp / "raw_R2.fastq.gz"
    singletons = job_tmp / "singletons.fastq.gz"

    try:
        lease = lease_preflight(
            args.lease_mode,
            args.lease_command,
            initial_cores=args.initial_cores,
            initial_memory_mb=args.initial_memory_mb,
        )
        cram_tokens = " ".join(q(token) for token in shlex.split(args.cram_options))
        if cram_tokens:
            cram_tokens += " "
        extract_script = job_tmp / "extract.sh"
        extract_script.write_text(
            "set -euo pipefail\n"
            f"{q(samtools)} view -@ 2 -u -h {cram_tokens}{q(alignment)} "
            f"| {q(samtools)} reset -@ 2 --output-fmt BAM,level=0 --no-PG --no-RG --keep-tag OQ "
            f"| {q(samtools)} sort -T {q(job_tmp / 'name-sort')} -@ 2 -u -n -m 6000M "
            f"| {q(samtools)} fastq -O -N -@ 2 -0 /dev/null -1 {q(raw1)} -2 {q(raw2)} -s {q(singletons)}\n",
            encoding="utf-8",
        )
        extract_metrics = job_tmp / "extract.io.json"
        phase_paths.append(extract_metrics)
        run_profiled(
            ["/usr/bin/bash", str(extract_script)],
            metrics_path=extract_metrics,
            label="external_adapter_fused.extract",
            local_paths=[job_tmp],
            input_paths=[alignment],
            output_paths=[raw1, raw2, singletons],
            requested_ssd_gb=args.ssd_gb,
            threads=5,
            memory_mb=int(args.initial_memory_mb),
            poll_interval=args.poll_interval,
        )
        for path in (raw1, raw2, singletons):
            if not path.is_file():
                raise RuntimeError(f"FASTQ extraction output is absent: {path}")

        lease = shrink_lease(
            lease, cores=args.adapter_cores, memory_mb=args.adapter_memory_mb
        )
        base1, max1 = quality_base(raw1)
        base2, max2 = quality_base(raw2)
        selected_base = 64 if base1 == 64 and base2 == 64 else 33
        quality_base_flag = f"--qualitybase {selected_base}" if selected_base == 64 else ""
        quality_output_flag = "--qualitybase-output 33" if selected_base == 64 else ""
        quality_max = 62

        prepared1 = job_tmp / "cut_1.fq.gz"
        prepared2 = job_tmp / "cut_2.fq.gz"
        adapter_log = job_tmp / "adapter_removal.log"
        fastq_stats = job_tmp / "fastq.stats.tsv"
        adapters = job_tmp / "adapters.txt"
        if args.attempt >= 2:
            record_rescue(Path(args.error_file), args.sample)
            source = (
                f"{q(sys.executable)} {q(args.pair_rescue_script)} --r1 {q(raw1)} --r2 {q(raw2)} "
                f"--decompress external --threads 5 --buffer-size 2048 --quiet"
            )
        else:
            source = (
                f"paste <({q(pigz)} -cd {q(raw1)} | paste - - - -) "
                f"<({q(pigz)} -cd {q(raw2)} | paste - - - -) | tr '\\t' '\\n'"
            )
        dedup = (
            f"| {q(sys.executable)} {q(args.remove_duplicates_script)} --input - --output - --quiet"
            if args.remove_duplicated_reads
            else ""
        )
        adapter_script = job_tmp / "adapter.sh"
        adapter_script.write_text(
            "set -euo pipefail\n"
            f"{source} {dedup} "
            f"| tee >({q(sys.executable)} {q(args.fastq_stats_script)} --interleaved --input - -s {q(fastq_stats)}) "
            f"| tee >({q(adapter_removal)} --identify-adapters --adapter-list {q(args.adapter_list)} "
            f"--interleaved-input --file1 /dev/stdin --threads 4 {quality_base_flag} {quality_output_flag} > {q(adapters)}) "
            f"| {q(adapter_removal)} --adapter-list {q(args.adapter_list)} --interleaved-input --file1 /dev/stdin "
            f"--gzip --gzip-level 1 --output1 {q(prepared1)} --output2 {q(prepared2)} "
            f"--settings {q(adapter_log)} --minlength 40 --singleton /dev/null --discarded /dev/null "
            f"--threads 4 {quality_base_flag} {quality_output_flag} --qualitymax {quality_max}\n",
            encoding="utf-8",
        )
        adapter_metrics = job_tmp / "adapter.io.json"
        phase_paths.append(adapter_metrics)
        run_profiled(
            ["/usr/bin/bash", str(adapter_script)],
            metrics_path=adapter_metrics,
            label="external_adapter_fused.adapter_removal",
            local_paths=[job_tmp],
            input_paths=[raw1, raw2],
            output_paths=[prepared1, prepared2, adapter_log, fastq_stats, adapters],
            requested_ssd_gb=args.ssd_gb,
            threads=5,
            memory_mb=int(args.adapter_memory_mb),
            poll_interval=args.poll_interval,
        )
        for path in (prepared1, prepared2, adapter_log, fastq_stats, adapters):
            if not path.is_file() or path.stat().st_size <= 0:
                raise RuntimeError(f"required adapter output is absent or empty: {path}")
        for source_path, destination in (
            (raw1, Path(args.output_raw_forward)),
            (raw2, Path(args.output_raw_reverse)),
            (singletons, Path(args.output_singletons)),
            (prepared1, Path(args.output_forward)),
            (prepared2, Path(args.output_reverse)),
            (fastq_stats, Path(args.output_fastq_stats)),
            (adapters, Path(args.output_adapters)),
            (adapter_log, Path(args.output_adapter_log)),
        ):
            _atomic_copy(source_path, destination)
        success = True
        return 0
    finally:
        phases = load_metrics(phase_paths)
        shutil.rmtree(job_tmp, ignore_errors=True)
        _atomic_json(
            Path(args.metrics),
            {
                "schema_version": 1,
                "label": "external_adapter_fused",
                "sample": args.sample,
                "readgroup": args.readgroup,
                "attempt": args.attempt,
                "success": success,
                "started_at_epoch": started,
                "duration_seconds": round(time.time() - started, 6),
                "requested": {
                    "initial_cores": args.initial_cores,
                    "initial_memory_mb": args.initial_memory_mb,
                    "adapter_cores": args.adapter_cores,
                    "adapter_memory_mb": args.adapter_memory_mb,
                    "ssd_gb": args.ssd_gb,
                },
                "observed_fastq_encoding": {
                    "r1_base": base1 if 'base1' in locals() else None,
                    "r1_max": max1 if 'max1' in locals() else None,
                    "r2_base": base2 if 'base2' in locals() else None,
                    "r2_max": max2 if 'max2' in locals() else None,
                },
                "lease": lease,
                "phases": phases,
                "scratch_job_directory": str(job_tmp),
                "scratch_removed": not job_tmp.exists(),
            },
        )


if __name__ == "__main__":
    raise SystemExit(main())
