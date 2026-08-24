#!/usr/bin/env python3
"""Run regional DeepVariant and its Whatshap/merge tail in one SSD job."""

from __future__ import annotations

import argparse
import json
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
    parser.add_argument("--sample", required=True)
    parser.add_argument("--region", required=True)
    parser.add_argument("--bed", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--bai", required=True)
    parser.add_argument("--validated-sex", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--deepvariant-reference", required=True)
    parser.add_argument("--deepvariant-runner", required=True)
    parser.add_argument("--model-type", required=True)
    parser.add_argument("--haploid-contigs", required=True)
    parser.add_argument("--ploidy", type=int, choices=(1, 2), required=True)
    parser.add_argument("--skip-sex", type=int, choices=(0, 1), required=True)
    parser.add_argument("--interval-bed", required=True)
    parser.add_argument("--capture-auto-bed", required=True)
    parser.add_argument("--capture-x-bed", required=True)
    parser.add_argument("--capture-y-bed", required=True)
    parser.add_argument("--merge-script", required=True)
    parser.add_argument("--stats-parser", required=True)
    parser.add_argument("--output-vcf", required=True)
    parser.add_argument("--output-vcf-tbi", required=True)
    parser.add_argument("--output-wstats", required=True)
    parser.add_argument("--output-merge-stats", required=True)
    parser.add_argument("--output-bcftools-stats", required=True)
    parser.add_argument("--output-bcftools-summary", required=True)
    parser.add_argument("--output-tmp-gvcf", required=True)
    parser.add_argument("--output-gvcf", required=True)
    parser.add_argument("--output-gvcf-tbi", required=True)
    parser.add_argument("--output-exome-gvcf", required=True)
    parser.add_argument("--output-exome-gvcf-tbi", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--whatshap", default="whatshap")
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--bgzip", default="bgzip")
    parser.add_argument("--tabix", default="tabix")
    parser.add_argument("--num-shards", type=int, default=8)
    parser.add_argument("--initial-cores", type=float, required=True)
    parser.add_argument("--initial-memory-mb", type=float, required=True)
    parser.add_argument("--low-cores", type=float, default=1.0)
    parser.add_argument("--low-memory-mb", type=float, default=3200)
    parser.add_argument("--attempt", type=int, default=1)
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


def load_metrics(paths: list[Path]) -> list[dict]:
    result = []
    for path in paths:
        try:
            result.append(json.loads(path.read_text(encoding="utf-8")))
        except FileNotFoundError:
            pass
    return result


def empty_vcf_commands(bgzip: str, tabix: str, target: Path) -> str:
    return (
        "printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n' "
        f"| {q(bgzip)} -c > {q(target)}\n"
        f"{q(tabix)} -f -p vcf {q(target)}\n"
    )


def stats_commands(args: argparse.Namespace, bcftools: str, phased: Path, stats: Path, summary: Path) -> str:
    parser = q(Path(args.stats_parser))
    if args.region == "F":
        commands = [f"rm -f {q(summary)}"]
        for region, bed in (
            ("A", args.capture_auto_bed),
            ("X", args.capture_x_bed),
            ("Y", args.capture_y_bed),
        ):
            commands.extend(
                [
                    f"{q(bcftools)} stats -R {q(bed)} -F {q(args.reference)} {q(phased)} > {q(stats)}",
                    f"{q(sys.executable)} {parser} {q(stats)} {q(summary)} --sample {q(args.sample)} --region {region} --append",
                ]
            )
        return "\n".join(commands) + "\n"
    capture = (
        args.capture_auto_bed
        if args.region.startswith(("A", "F"))
        else args.capture_x_bed
        if args.region.startswith("X")
        else args.capture_y_bed
    )
    return (
        f"{q(bcftools)} stats -R {q(capture)} -F {q(args.reference)} {q(phased)} > {q(stats)}\n"
        f"{q(sys.executable)} {parser} {q(stats)} {q(summary)}\n"
    )


def main() -> int:
    args = parse_args()
    for raw in (
        args.bed,
        args.bam,
        args.bai,
        args.validated_sex,
        args.reference,
        args.deepvariant_reference,
        args.interval_bed,
        args.capture_auto_bed,
        args.capture_x_bed,
        args.capture_y_bed,
        args.merge_script,
        args.stats_parser,
    ):
        if not Path(raw).is_file():
            raise FileNotFoundError(raw)
    deepvariant = executable(args.deepvariant_runner)
    whatshap = executable(args.whatshap)
    bcftools = executable(args.bcftools)
    bgzip = executable(args.bgzip)
    tabix = executable(args.tabix)

    scratch = assigned_scratch(args.scratch_base)
    parent = scratch / "deepvariant_phasing_fused"
    parent.mkdir(parents=True, exist_ok=True)
    job_tmp = Path(
        tempfile.mkdtemp(prefix=f"{args.sample}.{args.region}.", dir=parent)
    )
    phase_paths: list[Path] = []
    lease: dict = {}
    success = False
    started = time.time()

    raw_vcf = job_tmp / "raw.vcf.gz"
    raw_gvcf = job_tmp / "raw.g.vcf.gz"
    raw_vcf_tbi = Path(str(raw_vcf) + ".tbi")
    raw_gvcf_tbi = Path(str(raw_gvcf) + ".tbi")

    try:
        lease = lease_preflight(
            args.lease_mode,
            args.lease_command,
            initial_cores=args.initial_cores,
            initial_memory_mb=args.initial_memory_mb,
        )
        intermediate = job_tmp / "intermediate"
        intermediate.mkdir()
        dv_script = job_tmp / "run-deepvariant.sh"
        if args.skip_sex:
            dv_body = (
                "set -euo pipefail\n"
                + empty_vcf_commands(bgzip, tabix, raw_vcf)
                + empty_vcf_commands(bgzip, tabix, raw_gvcf)
            )
        else:
            config_string = (
                "config_string=inter_op_parallelism_threads: "
                f"{args.num_shards} intra_op_parallelism_threads: {args.num_shards} "
                f"device_count: {{ key: 'CPU' value: {args.num_shards} }}"
            )
            dv_body = (
                "set -euo pipefail\n"
                # DeepVariant 1.9 logs its complete child environment.  The
                # lease credentials are only needed by this parent runner,
                # after DeepVariant exits, so do not expose them in that log.
                "unset ZSLURM_LEASE_TOKEN ZSLURM_LEASE_SOCKET\n"
                "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1\n"
                f"export TF_NUM_INTRAOP_THREADS={args.num_shards} TF_NUM_INTEROP_THREADS={args.num_shards}\n"
                f"{q(deepvariant)} "
                f"--make_examples_extra_args {q(f'normalize_reads=true,regions={args.bed},small_model_call_multiallelics=false')} "
                f"--call_variants_extra_args {q(config_string)} "
                f"--num_shards={args.num_shards} --model_type={q(args.model_type)} "
                f"--ref={q(args.deepvariant_reference)} --reads={q(args.bam)} "
                f"--output_vcf={q(raw_vcf)} --output_gvcf={q(raw_gvcf)} "
                f"--haploid_contigs {q(args.haploid_contigs)} "
                f"--intermediate_results_dir {q(intermediate)} "
                f"--postprocess_cpus {args.num_shards}\n"
            )
        dv_script.write_text(dv_body, encoding="utf-8")
        dv_metrics = job_tmp / "deepvariant.io.json"
        phase_paths.append(dv_metrics)
        run_profiled(
            ["/usr/bin/bash", str(dv_script)],
            metrics_path=dv_metrics,
            label="deepvariant_phasing_fused.deepvariant",
            local_paths=[job_tmp],
            input_paths=[args.bam, args.bai, args.bed],
            output_paths=[raw_vcf, raw_vcf_tbi, raw_gvcf, raw_gvcf_tbi],
            requested_ssd_gb=args.ssd_gb,
            threads=args.num_shards,
            memory_mb=int(args.initial_memory_mb),
            poll_interval=args.poll_interval,
        )
        for path in (raw_vcf, raw_vcf_tbi, raw_gvcf, raw_gvcf_tbi):
            if not path.is_file() or path.stat().st_size <= 0:
                raise RuntimeError(f"DeepVariant did not create required output: {path}")

        lease = shrink_lease(
            lease,
            cores=args.low_cores,
            memory_mb=args.low_memory_mb,
            phase="phasing_merge",
        )

        phased = job_tmp / "phased.vcf.gz"
        phased_tbi = Path(str(phased) + ".tbi")
        wstats = job_tmp / "whatshap.stats"
        merge_stats = job_tmp / "merge.stats"
        bcftools_stats = job_tmp / "bcftools.stats.txt"
        bcftools_summary = job_tmp / "bcftools.summary.tsv"
        tmp_gvcf = job_tmp / "merged.g.vcf"
        final_gvcf = job_tmp / "merged.g.vcf.gz"
        final_gvcf_tbi = Path(str(final_gvcf) + ".tbi")
        exome_gvcf = job_tmp / "exome.g.vcf.gz"
        exome_gvcf_tbi = Path(str(exome_gvcf) + ".tbi")

        phase_script = job_tmp / "run-phasing.sh"
        body = "set -euo pipefail\n"
        if args.ploidy == 2 and not args.skip_sex:
            body += (
                f"{q(whatshap)} phase --ignore-read-groups --reference {q(args.reference)} {q(raw_vcf)} {q(args.bam)} -o {q(phased)}\n"
                f"{q(bcftools)} index -f -t {q(phased)}\n"
                f"{q(whatshap)} stats {q(phased)} > {q(wstats)}\n"
            )
            body += stats_commands(args, bcftools, phased, bcftools_stats, bcftools_summary)
            body += (
                f"{q(sys.executable)} {q(args.merge_script)} {q(raw_gvcf)} {q(phased)} {q(tmp_gvcf)} {q(merge_stats)}\n"
                f"{q(bcftools)} view {q(tmp_gvcf)} -o {q(final_gvcf)}\n"
                f"{q(bcftools)} index --tbi {q(final_gvcf)}\n"
            )
        else:
            body += (
                f"cp {q(raw_vcf)} {q(phased)}\n"
                f"cp {q(raw_vcf_tbi)} {q(phased_tbi)}\n"
                f"cp {q(raw_gvcf)} {q(final_gvcf)}\n"
                f"cp {q(raw_gvcf_tbi)} {q(final_gvcf_tbi)}\n"
                f"touch {q(tmp_gvcf)} {q(wstats)} {q(merge_stats)}\n"
            )
            body += stats_commands(args, bcftools, phased, bcftools_stats, bcftools_summary)
        if args.skip_sex:
            body += empty_vcf_commands(bgzip, tabix, exome_gvcf)
        else:
            body += (
                f"{q(bcftools)} view -R {q(args.interval_bed)} {q(final_gvcf)} -O z -o {q(exome_gvcf)}\n"
                f"{q(bcftools)} index -f -t {q(exome_gvcf)}\n"
            )
        phase_script.write_text(body, encoding="utf-8")
        phase_metrics = job_tmp / "phasing.io.json"
        phase_paths.append(phase_metrics)
        run_profiled(
            ["/usr/bin/bash", str(phase_script)],
            metrics_path=phase_metrics,
            label="deepvariant_phasing_fused.phasing_merge",
            local_paths=[job_tmp],
            input_paths=[raw_vcf, raw_gvcf, args.bam, args.bai],
            output_paths=[
                phased,
                phased_tbi,
                wstats,
                merge_stats,
                bcftools_stats,
                bcftools_summary,
                tmp_gvcf,
                final_gvcf,
                final_gvcf_tbi,
                exome_gvcf,
                exome_gvcf_tbi,
            ],
            requested_ssd_gb=args.ssd_gb,
            threads=1,
            memory_mb=int(args.low_memory_mb),
            poll_interval=args.poll_interval,
        )
        for path in (
            phased,
            phased_tbi,
            bcftools_summary,
            final_gvcf,
            final_gvcf_tbi,
            exome_gvcf,
            exome_gvcf_tbi,
        ):
            if not path.is_file() or path.stat().st_size <= 0:
                raise RuntimeError(f"required phasing output is absent or empty: {path}")
        for path in (wstats, merge_stats, bcftools_stats, tmp_gvcf):
            if not path.is_file():
                raise RuntimeError(f"required phasing output is absent: {path}")

        for source, destination in (
            (phased_tbi, Path(args.output_vcf_tbi)),
            (phased, Path(args.output_vcf)),
            (wstats, Path(args.output_wstats)),
            (merge_stats, Path(args.output_merge_stats)),
            (bcftools_stats, Path(args.output_bcftools_stats)),
            (bcftools_summary, Path(args.output_bcftools_summary)),
            (tmp_gvcf, Path(args.output_tmp_gvcf)),
            (final_gvcf_tbi, Path(args.output_gvcf_tbi)),
            (final_gvcf, Path(args.output_gvcf)),
            (exome_gvcf_tbi, Path(args.output_exome_gvcf_tbi)),
            (exome_gvcf, Path(args.output_exome_gvcf)),
        ):
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
                "label": "deepvariant_phasing_fused",
                "sample": args.sample,
                "region": args.region,
                "attempt": args.attempt,
                "success": success,
                "started_at_epoch": started,
                "duration_seconds": round(time.time() - started, 6),
                "requested": {
                    "initial_cores": args.initial_cores,
                    "initial_memory_mb": args.initial_memory_mb,
                    "low_cores": args.low_cores,
                    "low_memory_mb": args.low_memory_mb,
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
