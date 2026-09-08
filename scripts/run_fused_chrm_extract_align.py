#!/usr/bin/env python3
"""Extract chrM/NUMT reads and create four realigned BAMs in job scratch."""

from __future__ import annotations

import argparse
import json
import shlex
import shutil
import tempfile
import time
from pathlib import Path

from io_profile import run_profiled
from run_fused_alignment import (
    _atomic_json,
    _atomic_publish,
    assigned_scratch,
    executable,
)


def q(value: object) -> str:
    return shlex.quote(str(value))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-bam", required=True)
    parser.add_argument("--numts-bed", required=True)
    parser.add_argument("--original-reference", required=True)
    parser.add_argument("--shifted-reference", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-bam-chrm", required=True)
    parser.add_argument("--output-bai-chrm", required=True)
    parser.add_argument("--output-bam-shifted-chrm", required=True)
    parser.add_argument("--output-bai-shifted-chrm", required=True)
    parser.add_argument("--output-bam-numts", required=True)
    parser.add_argument("--output-bai-numts", required=True)
    parser.add_argument("--output-bam-shifted-numts", required=True)
    parser.add_argument("--output-bai-shifted-numts", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--bwa", default="bwa")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--memory-mb", type=int, default=2000)
    parser.add_argument("--ssd-gb", type=float, required=True)
    parser.add_argument("--scratch-base", help="Explicit scratch root for tests")
    parser.add_argument(
        "--shared-scratch-base",
        help="Workflow-local shared fallback when this job has no assigned SSD",
    )
    parser.add_argument("--poll-interval", type=float, default=5.0)
    return parser.parse_args()


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
    inputs = [
        Path(args.input_bam),
        Path(args.numts_bed),
        Path(args.original_reference),
        Path(args.shifted_reference),
    ]
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.threads < 1:
        raise ValueError("--threads must be positive")

    samtools = executable(args.samtools)
    bwa = executable(args.bwa)
    scratch = assigned_scratch(
        args.scratch_base,
        shared_fallback=args.shared_scratch_base,
    )
    parent = scratch / "chrm_extract_align_fused"
    parent.mkdir(parents=True, exist_ok=True)
    job_tmp = Path(tempfile.mkdtemp(prefix=f"{args.sample}.", dir=parent))
    phase_paths: list[Path] = []
    success = False
    started = time.time()

    chrm_fq1 = job_tmp / "chrM.R1.fastq.gz"
    chrm_fq2 = job_tmp / "chrM.R2.fastq.gz"
    numt_fq1 = job_tmp / "NUMTs.R1.fastq.gz"
    numt_fq2 = job_tmp / "NUMTs.R2.fastq.gz"
    local_outputs = {
        "bam_chrm": job_tmp / "chrM_orig.reads.bam",
        "bai_chrm": job_tmp / "chrM_orig.reads.bai",
        "bam_shifted_chrm": job_tmp / "chrM_shifted.reads.bam",
        "bai_shifted_chrm": job_tmp / "chrM_shifted.reads.bai",
        "bam_numts": job_tmp / "NUMTs.realign.bam",
        "bai_numts": job_tmp / "NUMTs.realign.bai",
        "bam_shifted_numts": job_tmp / "NUMTs_shifted.reads.bam",
        "bai_shifted_numts": job_tmp / "NUMTs_shifted.reads.bai",
    }

    try:
        tmpbed = job_tmp / "NUMTs_plus_chrM.bed"
        extract_script = job_tmp / "extract.sh"
        extract_script.write_text(
            "set -euo pipefail\n"
            f"{q(samtools)} view -@ {args.threads} -b -o {q(job_tmp / 'chrM.raw.bam')} "
            f"{q(args.input_bam)} chrM\n"
            f"{q(samtools)} sort -n -@ {args.threads} -T {q(job_tmp / 'chrM.sorttmp')} "
            f"-o {q(job_tmp / 'chrM.sorted.bam')} {q(job_tmp / 'chrM.raw.bam')}\n"
            f"{q(samtools)} collate -@ {args.threads} -o {q(job_tmp / 'chrM.collated.bam')} "
            f"{q(job_tmp / 'chrM.sorted.bam')}\n"
            f"{q(samtools)} fastq -O -N -@ {args.threads} -0 /dev/null -s /dev/null "
            f"-1 {q(chrm_fq1)} -2 {q(chrm_fq2)} {q(job_tmp / 'chrM.collated.bam')}\n"
            f"cp -- {q(args.numts_bed)} {q(tmpbed)}\n"
            f"printf 'chrM\\t1\\t999999999\\n' >> {q(tmpbed)}\n"
            f"{q(samtools)} view -@ {args.threads} -b -L {q(tmpbed)} "
            f"-o {q(job_tmp / 'NUMTs.raw.bam')} {q(args.input_bam)}\n"
            f"{q(samtools)} sort -n -@ {args.threads} -T {q(job_tmp / 'NUMTs.sorttmp')} "
            f"-o {q(job_tmp / 'NUMTs.sorted.bam')} {q(job_tmp / 'NUMTs.raw.bam')}\n"
            f"{q(samtools)} collate -@ {args.threads} -o {q(job_tmp / 'NUMTs.collated.bam')} "
            f"{q(job_tmp / 'NUMTs.sorted.bam')}\n"
            f"{q(samtools)} fastq -O -N -@ {args.threads} -0 /dev/null -s /dev/null "
            f"-1 {q(numt_fq1)} -2 {q(numt_fq2)} {q(job_tmp / 'NUMTs.collated.bam')}\n",
            encoding="utf-8",
        )
        extract_metrics = job_tmp / "extract.io.json"
        phase_paths.append(extract_metrics)
        run_profiled(
            ["/usr/bin/bash", str(extract_script)],
            metrics_path=extract_metrics,
            label="chrm_extract_align_fused.extract",
            local_paths=[job_tmp],
            input_paths=[args.input_bam, args.numts_bed],
            output_paths=[chrm_fq1, chrm_fq2, numt_fq1, numt_fq2],
            requested_ssd_gb=args.ssd_gb,
            threads=args.threads,
            memory_mb=args.memory_mb,
            poll_interval=args.poll_interval,
        )
        for path in (chrm_fq1, chrm_fq2, numt_fq1, numt_fq2):
            if not path.is_file() or path.stat().st_size <= 0:
                raise RuntimeError(f"required extraction output is absent or empty: {path}")

        rg = f"@RG\\tID:{args.sample}\\tSM:{args.sample}"
        alignments = (
            ("chrM_orig", args.original_reference, chrm_fq1, chrm_fq2, "bam_chrm", "bai_chrm"),
            ("chrM_shifted", args.shifted_reference, chrm_fq1, chrm_fq2, "bam_shifted_chrm", "bai_shifted_chrm"),
            ("NUMTs_orig", args.original_reference, numt_fq1, numt_fq2, "bam_numts", "bai_numts"),
            ("NUMTs_shifted", args.shifted_reference, numt_fq1, numt_fq2, "bam_shifted_numts", "bai_shifted_numts"),
        )
        commands = ["set -euo pipefail"]
        for label, reference, fq1, fq2, bam_key, bai_key in alignments:
            bam = local_outputs[bam_key]
            bai = local_outputs[bai_key]
            commands.append(
                f"{q(bwa)} mem -t {args.threads} -R {q(rg)} {q(reference)} {q(fq1)} {q(fq2)} "
                f"| {q(samtools)} sort -T {q(job_tmp / (label + '.sorttmp'))} -O bam "
                f"-@ {args.threads} -o {q(bam)}"
            )
            commands.append(
                f"{q(samtools)} index -@ {args.threads} -o {q(bai)} {q(bam)}"
            )
        align_script = job_tmp / "align.sh"
        align_script.write_text("\n".join(commands) + "\n", encoding="utf-8")
        align_metrics = job_tmp / "align.io.json"
        phase_paths.append(align_metrics)
        run_profiled(
            ["/usr/bin/bash", str(align_script)],
            metrics_path=align_metrics,
            label="chrm_extract_align_fused.align",
            local_paths=[job_tmp],
            input_paths=[chrm_fq1, chrm_fq2, numt_fq1, numt_fq2],
            output_paths=list(local_outputs.values()),
            requested_ssd_gb=args.ssd_gb,
            threads=args.threads * 2,
            memory_mb=args.memory_mb,
            poll_interval=args.poll_interval,
        )

        destinations = {
            "bam_chrm": args.output_bam_chrm,
            "bai_chrm": args.output_bai_chrm,
            "bam_shifted_chrm": args.output_bam_shifted_chrm,
            "bai_shifted_chrm": args.output_bai_shifted_chrm,
            "bam_numts": args.output_bam_numts,
            "bai_numts": args.output_bai_numts,
            "bam_shifted_numts": args.output_bam_shifted_numts,
            "bai_shifted_numts": args.output_bai_shifted_numts,
        }
        for key, source in local_outputs.items():
            if not source.is_file() or source.stat().st_size <= 0:
                raise RuntimeError(f"required alignment output is absent or empty: {source}")
            _atomic_publish(source, Path(destinations[key]))
        success = True
        return 0
    finally:
        phases = load_metrics(phase_paths)
        shutil.rmtree(job_tmp, ignore_errors=True)
        _atomic_json(
            Path(args.metrics),
            {
                "schema_version": 1,
                "label": "chrm_extract_align_fused",
                "sample": args.sample,
                "success": success,
                "started_at_epoch": started,
                "duration_seconds": round(time.time() - started, 6),
                "requested": {
                    "threads_per_tool": args.threads,
                    "memory_mb": args.memory_mb,
                    "ssd_gb": args.ssd_gb,
                },
                "phases": phases,
                "scratch_job_directory": str(job_tmp),
                "scratch_removed": not job_tmp.exists(),
            },
        )


if __name__ == "__main__":
    raise SystemExit(main())
