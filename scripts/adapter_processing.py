"""Single adapter-processing implementation for native FASTQ and BAM/CRAM.

Keep identification and trimming concurrent, with the same tool arguments,
quality encoding detection, optional deduplication and retry rescue in both
routes. Extraction and resource leases belong to the calling rule/runner.
"""
from __future__ import annotations

import bz2
import fcntl
import gzip
import shlex
import sys
from pathlib import Path

from io_profile import run_profiled


def q(value: object) -> str:
    return shlex.quote(str(value))


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


def prepare_adapters(args, raw1: Path, raw2: Path, job_tmp: Path,
                     phase_paths: list[Path], *, pigz: str,
                     adapter_removal: str, label: str):
    """Run the shared phase; return prepared paths and detected encodings."""
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
        label=label,
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
    return (prepared1, prepared2, adapter_log, fastq_stats, adapters), {
        "r1_base": base1, "r1_max": max1,
        "r2_base": base2, "r2_max": max2,
    }
