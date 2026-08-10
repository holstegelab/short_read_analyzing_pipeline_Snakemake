#!/usr/bin/env python3
"""Run the chrM/NUMT Mutect, filter, and BP-resolution tail on assigned SSD."""

from __future__ import annotations

import argparse
import json
import shlex
import shutil
import tempfile
import time
from pathlib import Path

from io_profile import run_profiled
from run_fused_alignment import _atomic_copy, _atomic_json, assigned_scratch, executable


def q(value: object) -> str:
    return shlex.quote(str(value))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam-chrm", required=True)
    parser.add_argument("--bai-chrm", required=True)
    parser.add_argument("--bam-shifted-chrm", required=True)
    parser.add_argument("--bai-shifted-chrm", required=True)
    parser.add_argument("--bam-numts", required=True)
    parser.add_argument("--bai-numts", required=True)
    parser.add_argument("--bam-shifted-numts", required=True)
    parser.add_argument("--bai-shifted-numts", required=True)
    parser.add_argument("--original-reference", required=True)
    parser.add_argument("--shifted-reference", required=True)
    parser.add_argument("--chain", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-chrm-gvcf", required=True)
    parser.add_argument("--output-chrm-tbi", required=True)
    parser.add_argument("--output-numt-gvcf", required=True)
    parser.add_argument("--output-numt-tbi", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--gatk", default="gatk")
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--tabix", default="tabix")
    parser.add_argument("--memory-mb", type=int, default=5000)
    parser.add_argument("--ssd-gb", type=float, required=True)
    parser.add_argument("--scratch-base", help="Explicit scratch root for tests")
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


def gatk_command(gatk: str, cpus: int, tool: str, options: str) -> str:
    java = (
        f"-XX:ActiveProcessorCount={cpus} "
        f"-XX:ParallelGCThreads={cpus} -XX:ConcGCThreads={cpus}"
    )
    return f"{q(gatk)} --java-options {q(java)} {tool} {options}"


def main() -> int:
    args = parse_args()
    source_paths = [
        Path(args.bam_chrm),
        Path(args.bai_chrm),
        Path(args.bam_shifted_chrm),
        Path(args.bai_shifted_chrm),
        Path(args.bam_numts),
        Path(args.bai_numts),
        Path(args.bam_shifted_numts),
        Path(args.bai_shifted_numts),
        Path(args.original_reference),
        Path(args.shifted_reference),
        Path(args.chain),
    ]
    for path in source_paths:
        if not path.is_file():
            raise FileNotFoundError(path)

    gatk = executable(args.gatk)
    bcftools = executable(args.bcftools)
    tabix = executable(args.tabix)
    scratch = assigned_scratch(args.scratch_base)
    parent = scratch / "chrm_tail_fused"
    parent.mkdir(parents=True, exist_ok=True)
    job_tmp = Path(tempfile.mkdtemp(prefix=f"{args.sample}.", dir=parent))
    phase_paths: list[Path] = []
    success = False
    started = time.time()

    local_inputs = {
        "bam_chrm": job_tmp / "chrM_orig.reads.bam",
        "bai_chrm": job_tmp / "chrM_orig.reads.bai",
        "bam_shifted_chrm": job_tmp / "chrM_shifted.reads.bam",
        "bai_shifted_chrm": job_tmp / "chrM_shifted.reads.bai",
        "bam_numts": job_tmp / "NUMTs.realign.bam",
        "bai_numts": job_tmp / "NUMTs.realign.bai",
        "bam_shifted_numts": job_tmp / "NUMTs_shifted.reads.bam",
        "bai_shifted_numts": job_tmp / "NUMTs_shifted.reads.bai",
    }
    source_bams = source_paths[:8]

    try:
        stage_script = job_tmp / "stage.sh"
        stage_script.write_text(
            "set -euo pipefail\n"
            + "\n".join(
                f"cp -- {q(source)} {q(destination)}"
                for source, destination in zip(source_bams, local_inputs.values())
            )
            + "\n",
            encoding="utf-8",
        )
        stage_metrics = job_tmp / "stage.io.json"
        phase_paths.append(stage_metrics)
        run_profiled(
            ["/usr/bin/bash", str(stage_script)],
            metrics_path=stage_metrics,
            label="chrm_tail_fused.stage",
            local_paths=[job_tmp],
            input_paths=source_bams,
            output_paths=list(local_inputs.values()),
            requested_ssd_gb=args.ssd_gb,
            threads=1,
            memory_mb=args.memory_mb,
            poll_interval=args.poll_interval,
        )
        for source, staged in zip(source_bams, local_inputs.values()):
            if staged.stat().st_size != source.stat().st_size:
                raise RuntimeError(f"staged input size mismatch: {source} -> {staged}")

        original = q(args.original_reference)
        shifted = q(args.shifted_reference)
        chain = q(args.chain)
        chrm_orig = job_tmp / "chrM_orig.vcf.gz"
        chrm_shift = job_tmp / "chrM_shifted.vcf.gz"
        chrm_back = job_tmp / "chrM_shifted_back.vcf.gz"
        chrm_merged_stats = job_tmp / "chrM_merged.vcf.gz.stats"
        chrm_merged = job_tmp / "chrM_merged.vcf.gz"
        chrm_filtered = job_tmp / "chrM_filtered.vcf.gz"
        numt_orig = job_tmp / "NUMT_orig.vcf.gz"
        numt_shift = job_tmp / "NUMT_shifted.vcf.gz"
        numt_back = job_tmp / "NUMT_shifted_back.vcf.gz"
        numt_merged_stats = job_tmp / "NUMT_merged.vcf.gz.stats"
        numt_merged = job_tmp / "NUMT_merged.vcf.gz"
        numt_filtered = job_tmp / "NUMT_filtered.vcf.gz"

        standard = ["set -euo pipefail"]
        for bam_key, reference, interval, output in (
            ("bam_chrm", original, "chrM:4142-12425", chrm_orig),
            ("bam_shifted_chrm", shifted, "chrM:4142-12426", chrm_shift),
            ("bam_numts", original, "chrM:4142-12425", numt_orig),
            ("bam_shifted_numts", shifted, "chrM:4142-12426", numt_shift),
        ):
            standard.append(
                gatk_command(
                    gatk,
                    1,
                    "Mutect2",
                    f"-R {reference} -L {interval} --mitochondria-mode "
                    f"-I {q(local_inputs[bam_key])} -O {q(output)}",
                )
            )
            standard.append(f"{q(tabix)} -f -p vcf {q(output)}")
        for shifted_vcf, back_vcf in ((chrm_shift, chrm_back), (numt_shift, numt_back)):
            standard.append(
                gatk_command(
                    gatk,
                    1,
                    "LiftoverVcf",
                    f"-I {q(shifted_vcf)} -O {q(back_vcf)} -C {chain} "
                    f"-R {original} --REJECT /dev/null",
                )
            )
            standard.append(f"{q(tabix)} -f -p vcf {q(back_vcf)}")
        for orig, shifted_vcf, back_vcf, merged_stats, merged, filtered in (
            (chrm_orig, chrm_shift, chrm_back, chrm_merged_stats, chrm_merged, chrm_filtered),
            (numt_orig, numt_shift, numt_back, numt_merged_stats, numt_merged, numt_filtered),
        ):
            standard.extend(
                [
                    gatk_command(
                        gatk,
                        4,
                        "MergeMutectStats",
                        f"--stats {q(str(orig) + '.stats')} "
                        f"--stats {q(str(shifted_vcf) + '.stats')} -O {q(merged_stats)}",
                    ),
                    gatk_command(
                        gatk,
                        4,
                        "MergeVcfs",
                        f"-I {q(back_vcf)} -I {q(orig)} -O {q(merged)}",
                    ),
                    f"{q(tabix)} -f -p vcf {q(merged)}",
                    gatk_command(
                        gatk,
                        4,
                        "FilterMutectCalls",
                        f"-OVI true -V {q(merged)} -R {original} "
                        f"--mitochondria-mode True -O {q(filtered)}",
                    ),
                    f"{q(tabix)} -f -p vcf {q(filtered)}",
                ]
            )
        standard_script = job_tmp / "mutect_filter.sh"
        standard_script.write_text("\n".join(standard) + "\n", encoding="utf-8")
        standard_metrics = job_tmp / "mutect_filter.io.json"
        phase_paths.append(standard_metrics)
        run_profiled(
            ["/usr/bin/bash", str(standard_script)],
            metrics_path=standard_metrics,
            label="chrm_tail_fused.mutect_filter",
            local_paths=[job_tmp],
            input_paths=list(local_inputs.values()),
            output_paths=[chrm_filtered, Path(str(chrm_filtered) + ".tbi"), numt_filtered, Path(str(numt_filtered) + ".tbi")],
            requested_ssd_gb=args.ssd_gb,
            threads=4,
            memory_mb=args.memory_mb,
            poll_interval=args.poll_interval,
        )

        chrm_bp = job_tmp / "chrM_orig_BP.g.vcf.gz"
        chrm_shift_bp = job_tmp / "chrM_shifted_BP.g.vcf.gz"
        chrm_back_bp = job_tmp / "chrM_shifted_back_BP.g.vcf.gz"
        chrm_merged_bp = job_tmp / "chrM_merged_BP.g.vcf.gz"
        chrm_norm = job_tmp / "chrM_merged_BP_norm.g.vcf.gz"
        chrm_final = job_tmp / "chrM_merged_BP_annotated.g.vcf.gz"
        numt_bp = job_tmp / "NUMT_orig_BP.g.vcf.gz"
        numt_shift_bp = job_tmp / "NUMT_shifted_BP.g.vcf.gz"
        numt_back_bp = job_tmp / "NUMT_shifted_back_BP.g.vcf.gz"
        numt_merged_bp = job_tmp / "NUMT_merged_BP.g.vcf.gz"
        numt_norm = job_tmp / "NUMT_merged_BP_norm.g.vcf.gz"
        numt_final = job_tmp / "NUMT_merged_BP_annotated.g.vcf.gz"

        bp = ["set -euo pipefail"]
        for bam_key, reference, interval, output in (
            ("bam_chrm", original, "chrM:4142-12425", chrm_bp),
            ("bam_shifted_chrm", shifted, "chrM:4142-12426", chrm_shift_bp),
            ("bam_numts", original, "chrM:4142-12425", numt_bp),
            ("bam_shifted_numts", shifted, "chrM:4142-12426", numt_shift_bp),
        ):
            bp.append(
                gatk_command(
                    gatk,
                    1,
                    "Mutect2",
                    f"-ERC BP_RESOLUTION -R {reference} -L {interval} "
                    f"--mitochondria-mode -I {q(local_inputs[bam_key])} -O {q(output)}",
                )
            )
        for shifted_vcf, back_vcf in ((chrm_shift_bp, chrm_back_bp), (numt_shift_bp, numt_back_bp)):
            bp.append(
                gatk_command(
                    gatk,
                    1,
                    "LiftoverVcf",
                    f"-I {q(shifted_vcf)} -O {q(back_vcf)} -C {chain} "
                    f"-R {original} --REJECT /dev/null",
                )
            )
        for orig, back, merged, norm, annotation, final in (
            (chrm_bp, chrm_back_bp, chrm_merged_bp, chrm_norm, chrm_filtered, chrm_final),
            (numt_bp, numt_back_bp, numt_merged_bp, numt_norm, numt_filtered, numt_final),
        ):
            bp.extend(
                [
                    gatk_command(gatk, 1, "MergeVcfs", f"-I {q(back)} -I {q(orig)} -O {q(merged)}"),
                    f"{q(bcftools)} norm -d exact -o {q(norm)} -O z {q(merged)}",
                    f"{q(tabix)} -f -p vcf {q(norm)}",
                    f"{q(bcftools)} annotate -a {q(annotation)} -c FILTER -O z -o {q(final)} {q(norm)}",
                    f"{q(tabix)} -f -p vcf {q(final)}",
                ]
            )
        bp_script = job_tmp / "bp_resolution.sh"
        bp_script.write_text("\n".join(bp) + "\n", encoding="utf-8")
        bp_metrics = job_tmp / "bp_resolution.io.json"
        phase_paths.append(bp_metrics)
        final_local = [chrm_final, Path(str(chrm_final) + ".tbi"), numt_final, Path(str(numt_final) + ".tbi")]
        run_profiled(
            ["/usr/bin/bash", str(bp_script)],
            metrics_path=bp_metrics,
            label="chrm_tail_fused.bp_resolution",
            local_paths=[job_tmp],
            input_paths=[chrm_filtered, numt_filtered],
            output_paths=final_local,
            requested_ssd_gb=args.ssd_gb,
            threads=1,
            memory_mb=args.memory_mb,
            poll_interval=args.poll_interval,
        )

        destinations = [
            Path(args.output_chrm_gvcf),
            Path(args.output_chrm_tbi),
            Path(args.output_numt_gvcf),
            Path(args.output_numt_tbi),
        ]
        for source, destination in zip(final_local, destinations):
            if not source.is_file() or source.stat().st_size <= 0:
                raise RuntimeError(f"required final output is absent or empty: {source}")
            _atomic_copy(source, destination)
        success = True
        return 0
    finally:
        phases = load_metrics(phase_paths)
        shutil.rmtree(job_tmp, ignore_errors=True)
        _atomic_json(
            Path(args.metrics),
            {
                "schema_version": 1,
                "label": "chrm_tail_fused",
                "sample": args.sample,
                "success": success,
                "started_at_epoch": started,
                "duration_seconds": round(time.time() - started, 6),
                "requested": {"memory_mb": args.memory_mb, "ssd_gb": args.ssd_gb},
                "phases": phases,
                "scratch_job_directory": str(job_tmp),
                "scratch_removed": not job_tmp.exists(),
            },
        )


if __name__ == "__main__":
    raise SystemExit(main())
