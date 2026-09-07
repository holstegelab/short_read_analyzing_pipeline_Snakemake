#!/usr/bin/env python3
"""Stage one markdup BAM and run its independent QC consumers from SSD."""

from __future__ import annotations

import argparse
import json
import shlex
import shutil
import tempfile
import time
from pathlib import Path

from io_profile import run_profiled
from pipeline_runtime import (
    _atomic_copy,
    _atomic_json,
    assigned_scratch,
    executable,
    lease_preflight,
)


QC_RELEASE_RESOURCES = {
    # Relative shares of the 6-core/7500-MB median-sized reservation.  The
    # tools may briefly use more; the reservation only controls ZSlurm packing.
    "verifybamid": (1.0, 350.0),
    "hs_metrics": (1.0, 2650.0),
    "artifact_oxog": (1.0, 2650.0),
    "samtools_stats": (1.0, 350.0),
    "bamstats": (1.0, 550.0),
    "mosdepth": (1.0, 950.0),
}


def q(value: object) -> str:
    return shlex.quote(str(value))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--bai", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--svd-prefix", required=True)
    parser.add_argument("--hs-interval", required=True)
    parser.add_argument("--targets-interval", required=True)
    parser.add_argument("--artifact-interval", required=True)
    parser.add_argument("--dbsnp", required=True)
    parser.add_argument("--capture-bed", required=True)
    parser.add_argument("--windows-bed", required=True)
    parser.add_argument("--bamstats-script", required=True)
    parser.add_argument("--output-selfsm", required=True)
    parser.add_argument("--output-ancestry", required=True)
    parser.add_argument("--output-hs-metrics", required=True)
    parser.add_argument("--output-bait-summary", required=True)
    parser.add_argument("--output-pre-adapter-summary", required=True)
    parser.add_argument("--output-bait-detail", required=True)
    parser.add_argument("--output-pre-adapter-detail", required=True)
    parser.add_argument("--output-error-summary", required=True)
    parser.add_argument("--output-oxog", required=True)
    parser.add_argument("--output-samtools-genome", required=True)
    parser.add_argument("--output-samtools-exome", required=True)
    parser.add_argument("--output-bamstats-all", required=True)
    parser.add_argument("--output-bamstats-exome", required=True)
    parser.add_argument("--output-coverage-regions", required=True)
    parser.add_argument("--output-coverage-csi", required=True)
    parser.add_argument("--output-coverage-global-dist", required=True)
    parser.add_argument("--output-coverage-summary", required=True)
    parser.add_argument("--output-coverage-region-dist", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--gatk", default="gatk")
    parser.add_argument("--verifybamid", default="verifybamid2")
    parser.add_argument("--mosdepth", default="mosdepth")
    parser.add_argument("--pypy", default="pypy")
    parser.add_argument("--cores", type=int, default=6)
    parser.add_argument("--memory-mb", type=int, default=7500)
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


def write_parallel_qc_script(
    path: Path,
    tasks: list[dict],
    lease: dict,
    lease_mode: str,
    job_tmp: Path,
    sample: str,
) -> dict[str, Path]:
    """Write parallel task wrappers that release their lease share on exit."""

    lines = ["set -uo pipefail", "pids=()"]
    release_paths: dict[str, Path] = {}
    for task in tasks:
        name = task["name"]
        task_script = job_tmp / f"qc.{name}.sh"
        task_script.write_text(
            "set -euo pipefail\n" + task["command"] + "\n",
            encoding="utf-8",
        )
        release_path = job_tmp / f"lease_release.{name}.json"
        release_paths[name] = release_path
        lines.extend(
            [
                "(",
                "  set +e",
                f"  /usr/bin/bash {q(task_script)}",
                "  task_rc=$?",
                "  release_rc=0",
            ]
        )
        if lease.get("available"):
            release_tmp = Path(str(release_path) + ".tmp")
            release_id = f"bam-qc:{sample}:{name}"
            release_command = " ".join(
                [
                    q(lease["command"]),
                    "--json",
                    "release",
                    "--cores",
                    q(task["cores"]),
                    "--mem-mb",
                    q(task["memory_mb"]),
                    "--release-id",
                    q(release_id),
                ]
            )
            lines.extend(
                [
                    f"  {release_command} > {q(release_tmp)}",
                    "  release_rc=$?",
                    "  if [ \"$release_rc\" -eq 0 ]; then",
                    f"    mv -f -- {q(release_tmp)} {q(release_path)}",
                    "  else",
                    f"    rm -f -- {q(release_tmp)}",
                    f"    echo 'lease release failed for {name}' >&2",
                    "  fi",
                ]
            )
        lines.extend(
            [
                "  if [ \"$task_rc\" -ne 0 ]; then exit \"$task_rc\"; fi",
            ]
        )
        if lease.get("available") and lease_mode == "required":
            lines.append(
                "  if [ \"$release_rc\" -ne 0 ]; then exit \"$release_rc\"; fi"
            )
        lines.extend(["  exit 0", ") &", 'pids+=("$!")'])

    lines.extend(
        [
            "status=0",
            'for pid in "${pids[@]}"; do',
            '  if ! wait "$pid"; then status=1; fi',
            "done",
            'exit "$status"',
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return release_paths


def load_task_releases(
    tasks: list[dict], release_paths: dict[str, Path], lease: dict
) -> list[dict]:
    results = []
    for task in tasks:
        result = {
            "name": task["name"],
            "requested_release_cores": task["cores"],
            "requested_release_mem_mb": task["memory_mb"],
            "performed": False,
        }
        path = release_paths.get(task["name"])
        if path is not None and path.is_file():
            try:
                response = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                result["reason"] = f"invalid release response: {exc}"
            else:
                result["performed"] = bool(response.get("ok"))
                result["response"] = response
        elif not lease.get("available"):
            result["reason"] = "lease unavailable"
        else:
            result["reason"] = "release command did not produce a response"
        results.append(result)
    return results


def main() -> int:
    args = parse_args()
    files = [
        Path(args.bam),
        Path(args.bai),
        Path(args.reference),
        Path(args.hs_interval),
        Path(args.targets_interval),
        Path(args.artifact_interval),
        Path(args.dbsnp),
        Path(args.capture_bed),
        Path(args.windows_bed),
        Path(args.bamstats_script),
    ]
    for path in files:
        if not path.is_file():
            raise FileNotFoundError(path)
    svd = Path(args.svd_prefix)
    if not any(svd.parent.glob(svd.name + "*")):
        raise FileNotFoundError(f"no VerifyBamID SVD files found for prefix {svd}")
    if args.cores < 1 or args.memory_mb < 1:
        raise ValueError("cores and memory must be positive")

    samtools = executable(args.samtools)
    gatk = executable(args.gatk)
    verifybamid = executable(args.verifybamid)
    mosdepth = executable(args.mosdepth)
    pypy = executable(args.pypy)
    scratch = assigned_scratch(args.scratch_base)
    parent = scratch / "bam_qc_fused"
    parent.mkdir(parents=True, exist_ok=True)
    job_tmp = Path(tempfile.mkdtemp(prefix=f"{args.sample}.", dir=parent))
    phase_paths: list[Path] = []
    task_release_paths: dict[str, Path] = {}
    lease: dict = {}
    tasks: list[dict] = []
    success = False
    started = time.time()
    local_bam = job_tmp / "markdup.bam"
    local_bai = job_tmp / "markdup.bam.bai"

    local = {
        "selfsm": job_tmp / "verifybamid.pca2.selfSM",
        "ancestry": job_tmp / "verifybamid.pca2.Ancestry",
        "hs": job_tmp / "hs_metrics",
        "bait_summary": job_tmp / "artifact.bait_bias_summary_metrics",
        "pre_summary": job_tmp / "artifact.pre_adapter_summary_metrics",
        "bait_detail": job_tmp / "artifact.bait_bias_detail_metrics",
        "pre_detail": job_tmp / "artifact.pre_adapter_detail_metrics",
        "error_summary": job_tmp / "artifact.error_summary_metrics",
        "oxog": job_tmp / "OXOG",
        "samtools_genome": job_tmp / "samtools.stat",
        "samtools_exome": job_tmp / "samtools.exome.stat",
        "bamstats_all": job_tmp / "bam_all.tsv",
        "bamstats_exome": job_tmp / "bam_exome.tsv",
        "coverage_regions": job_tmp / "coverage.regions.bed.gz",
        "coverage_csi": job_tmp / "coverage.regions.bed.gz.csi",
        "coverage_global_dist": job_tmp / "coverage.mosdepth.global.dist.txt",
        "coverage_summary": job_tmp / "coverage.mosdepth.summary.txt",
        "coverage_region_dist": job_tmp / "coverage.mosdepth.region.dist.txt",
    }

    try:
        lease = lease_preflight(
            args.lease_mode,
            args.lease_command,
            initial_cores=args.cores,
            initial_memory_mb=args.memory_mb,
        )

        stage_script = job_tmp / "stage.sh"
        stage_script.write_text(
            "set -euo pipefail\n"
            f"cp -- {q(args.bam)} {q(local_bam)}\n"
            f"cp -- {q(args.bai)} {q(local_bai)}\n",
            encoding="utf-8",
        )
        stage_metrics = job_tmp / "stage.io.json"
        phase_paths.append(stage_metrics)
        run_profiled(
            ["/usr/bin/bash", str(stage_script)],
            metrics_path=stage_metrics,
            label="bam_qc_fused.stage",
            local_paths=[job_tmp],
            input_paths=[args.bam, args.bai],
            output_paths=[local_bam, local_bai],
            requested_ssd_gb=args.ssd_gb,
            threads=1,
            memory_mb=args.memory_mb,
            poll_interval=args.poll_interval,
        )
        if local_bam.stat().st_size != Path(args.bam).stat().st_size:
            raise RuntimeError("staged BAM size does not match source")
        if local_bai.stat().st_size != Path(args.bai).stat().st_size:
            raise RuntimeError("staged BAM index size does not match source")

        java = "-Xmx3500M -XX:ActiveProcessorCount=2 -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2"
        artifact_prefix = job_tmp / "artifact"
        coverage_prefix = job_tmp / "coverage"
        task_commands = {
            "verifybamid": (
                f"{q(verifybamid)} --BamFile {q(local_bam)} --SVDPrefix {q(args.svd_prefix)} "
                f"--Reference {q(args.reference)} --DisableSanityCheck --NumThread 2 "
                f"--Output {q(job_tmp / 'verifybamid.pca2')}"
            ),
            "hs_metrics": (
                f"mkdir -p {q(job_tmp / 'hs_tmp')}\n"
                f"{q(gatk)} --java-options {q(java)} CollectHsMetrics --TMP_DIR {q(job_tmp / 'hs_tmp')} "
                f"-I {q(local_bam)} -R {q(args.reference)} -BI {q(args.hs_interval)} "
                f"-TI {q(args.targets_interval)} -Q 10 -MQ 10 -O {q(local['hs'])}"
            ),
            "artifact_oxog": (
                f"mkdir -p {q(job_tmp / 'artifact_tmp')}\n"
                f"{q(gatk)} --java-options {q(java)} CollectSequencingArtifactMetrics "
                f"--TMP_DIR {q(job_tmp / 'artifact_tmp')} -I {q(local_bam)} -O {q(artifact_prefix)} "
                f"-R {q(args.reference)} --DB_SNP {q(args.dbsnp)} --INTERVALS {q(args.artifact_interval)}\n"
                f"{q(gatk)} --java-options {q(java)} CollectOxoGMetrics -I {q(local_bam)} "
                f"-O {q(local['oxog'])} -R {q(args.reference)} --INTERVALS {q(args.artifact_interval)}"
            ),
            "samtools_stats": (
                f"{q(samtools)} stat -@ 2 -r {q(args.reference)} -d -p {q(local_bam)} > {q(local['samtools_genome'])}\n"
                f"{q(samtools)} stat -@ 2 -t {q(args.capture_bed)} -d -p -r {q(args.reference)} "
                f"{q(local_bam)} > {q(local['samtools_exome'])}"
            ),
            "bamstats": (
                f"{q(samtools)} view -s 0.05 -h {q(local_bam)} --threads 1 "
                f"| {q(pypy)} {q(args.bamstats_script)} stats > {q(local['bamstats_all'])}\n"
                f"{q(samtools)} view -s 0.05 -h {q(local_bam)} --threads 1 -L {q(args.capture_bed)} "
                f"| {q(pypy)} {q(args.bamstats_script)} stats > {q(local['bamstats_exome'])}"
            ),
            "mosdepth": (
                f"{q(mosdepth)} --threads 2 -b {q(args.windows_bed)} --no-per-base "
                f"{q(coverage_prefix)} {q(local_bam)}"
            ),
        }
        tasks = [
            {
                "name": name,
                "command": command,
                "cores": QC_RELEASE_RESOURCES[name][0],
                "memory_mb": QC_RELEASE_RESOURCES[name][1],
            }
            for name, command in task_commands.items()
        ]
        release_cores = sum(task["cores"] for task in tasks)
        release_memory_mb = sum(task["memory_mb"] for task in tasks)
        if release_cores > args.cores + 1e-9:
            raise ValueError(
                f"QC task releases require {release_cores} cores, but only "
                f"{args.cores} were reserved"
            )
        if release_memory_mb > args.memory_mb + 1e-9:
            raise ValueError(
                f"QC task releases require {release_memory_mb} MB, but only "
                f"{args.memory_mb} MB were reserved"
            )
        qc_script = job_tmp / "qc.sh"
        task_release_paths = write_parallel_qc_script(
            qc_script,
            tasks,
            lease,
            args.lease_mode,
            job_tmp,
            args.sample,
        )
        qc_metrics = job_tmp / "qc.io.json"
        phase_paths.append(qc_metrics)
        run_profiled(
            ["/usr/bin/bash", str(qc_script)],
            metrics_path=qc_metrics,
            label="bam_qc_fused.qc",
            local_paths=[job_tmp],
            input_paths=[local_bam, local_bai],
            output_paths=list(local.values()),
            requested_ssd_gb=args.ssd_gb,
            threads=args.cores,
            memory_mb=args.memory_mb,
            poll_interval=args.poll_interval,
        )

        destinations = {
            "selfsm": args.output_selfsm,
            "ancestry": args.output_ancestry,
            "hs": args.output_hs_metrics,
            "bait_summary": args.output_bait_summary,
            "pre_summary": args.output_pre_adapter_summary,
            "bait_detail": args.output_bait_detail,
            "pre_detail": args.output_pre_adapter_detail,
            "error_summary": args.output_error_summary,
            "oxog": args.output_oxog,
            "samtools_genome": args.output_samtools_genome,
            "samtools_exome": args.output_samtools_exome,
            "bamstats_all": args.output_bamstats_all,
            "bamstats_exome": args.output_bamstats_exome,
            "coverage_regions": args.output_coverage_regions,
            "coverage_csi": args.output_coverage_csi,
            "coverage_global_dist": args.output_coverage_global_dist,
            "coverage_summary": args.output_coverage_summary,
            "coverage_region_dist": args.output_coverage_region_dist,
        }
        nonempty = {
            "pre_summary",
            "bait_detail",
            "pre_detail",
            "error_summary",
            "samtools_genome",
            "samtools_exome",
            "bamstats_all",
            "bamstats_exome",
        }
        for key, source in local.items():
            if not source.is_file():
                raise RuntimeError(f"required QC output is absent: {source}")
            if key in nonempty and source.stat().st_size <= 0:
                raise RuntimeError(f"required QC output is empty: {source}")
            _atomic_copy(source, Path(destinations[key]))
        success = True
        return 0
    finally:
        phases = load_metrics(phase_paths)
        task_releases = load_task_releases(tasks, task_release_paths, lease)
        shutil.rmtree(job_tmp, ignore_errors=True)
        _atomic_json(
            Path(args.metrics),
            {
                "schema_version": 1,
                "label": "bam_qc_fused",
                "sample": args.sample,
                "success": success,
                "started_at_epoch": started,
                "duration_seconds": round(time.time() - started, 6),
                "requested": {
                    "cores": args.cores,
                    "memory_mb": args.memory_mb,
                    "ssd_gb": args.ssd_gb,
                },
                "parallel_consumers": len(tasks) if "tasks" in locals() else 0,
                "lease": lease,
                "task_releases": task_releases,
                "phases": phases,
                "scratch_job_directory": str(job_tmp),
                "scratch_removed": not job_tmp.exists(),
            },
        )


if __name__ == "__main__":
    raise SystemExit(main())
