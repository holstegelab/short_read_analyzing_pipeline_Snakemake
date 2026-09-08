#!/usr/bin/env python3
"""Build a sample KMC database and derive sex statistics in one SSD job."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq", nargs="+", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output-yaml", required=True)
    parser.add_argument("--output-chry", required=True)
    parser.add_argument("--output-chrx", required=True)
    parser.add_argument("--output-chrm", required=True)
    parser.add_argument("--output-auto", required=True)
    parser.add_argument("--kmer-chry", required=True)
    parser.add_argument("--kmer-chrx", required=True)
    parser.add_argument("--kmer-chrm", required=True)
    parser.add_argument("--kmer-auto", required=True)
    parser.add_argument("--process-sex", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--kmc", default="kmc")
    parser.add_argument("--kmc-tools", default="kmc_tools")
    parser.add_argument("--initial-cores", type=float, required=True)
    parser.add_argument("--initial-memory-mb", type=float, required=True)
    parser.add_argument("--kmc-threads", type=int, default=2)
    parser.add_argument("--low-cores", type=float, default=0.5)
    parser.add_argument("--low-memory-mb", type=float, default=3000)
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


def require_kmc_prefix(prefix: str) -> None:
    for suffix in (".kmc_pre", ".kmc_suf"):
        path = Path(prefix + suffix)
        if not path.is_file():
            raise FileNotFoundError(path)


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
    for raw in args.fastq:
        if not Path(raw).is_file():
            raise FileNotFoundError(raw)
    for prefix in (args.kmer_chry, args.kmer_chrx, args.kmer_chrm, args.kmer_auto):
        require_kmc_prefix(prefix)
    process_sex = Path(args.process_sex)
    if not process_sex.is_file():
        raise FileNotFoundError(process_sex)
    kmc = executable(args.kmc)
    kmc_tools = executable(args.kmc_tools)

    scratch = assigned_scratch(args.scratch_base)
    parent = scratch / "kmer_sex_fused"
    parent.mkdir(parents=True, exist_ok=True)
    job_tmp = Path(tempfile.mkdtemp(prefix=f"{args.sample}.", dir=parent))
    phase_paths: list[Path] = []
    lease: dict = {}
    success = False
    started = time.time()

    try:
        lease = lease_preflight(
            args.lease_mode,
            args.lease_command,
            initial_cores=args.initial_cores,
            initial_memory_mb=args.initial_memory_mb,
        )

        input_list = job_tmp / "inputs.lst"
        input_list.write_text(
            "".join(f"{Path(path).resolve()}\n" for path in args.fastq),
            encoding="utf-8",
        )
        database = job_tmp / args.sample
        kmc_tmp = job_tmp / "kmc_tmp"
        kmc_tmp.mkdir()
        kmc_metrics = job_tmp / "kmc.io.json"
        phase_paths.append(kmc_metrics)
        run_profiled(
            [
                kmc,
                "-fq",
                "-k32",
                "-cs8192",
                "-sf12",
                "-sp12",
                "-sr1",
                "-m36",
                f"@{input_list}",
                str(database),
                str(kmc_tmp),
            ],
            metrics_path=kmc_metrics,
            label="kmer_sex_fused.kmc",
            local_paths=[job_tmp],
            input_paths=args.fastq,
            output_paths=[Path(str(database) + ".kmc_pre"), Path(str(database) + ".kmc_suf")],
            requested_ssd_gb=args.ssd_gb,
            threads=args.kmc_threads,
            memory_mb=int(args.initial_memory_mb),
            poll_interval=args.poll_interval,
        )
        require_kmc_prefix(str(database))

        lease = shrink_lease(
            lease,
            cores=args.low_cores,
            memory_mb=args.low_memory_mb,
            phase="sex_check",
        )

        chry = job_tmp / "chry.tsv"
        chrx = job_tmp / "chrx.tsv"
        chrm = job_tmp / "chrm.tsv"
        auto = job_tmp / "auto.tsv"
        yaml_output = job_tmp / "result.yaml"
        sex_script = job_tmp / "run-sex.sh"
        sex_script.write_text(
            "set -euo pipefail\n"
            f"{shlex_quote(kmc_tools)} -t1 simple {shlex_quote(database)} {shlex_quote(args.kmer_chry)} -cx1 intersect {shlex_quote(str(chry) + '.tmp')} -ocleft\n"
            f"{shlex_quote(kmc_tools)} -t1 simple {shlex_quote(str(chry) + '.tmp')} {shlex_quote(args.kmer_chry)} union {shlex_quote(chry)} -ocsum\n"
            f"{shlex_quote(kmc_tools)} -t1 simple {shlex_quote(database)} {shlex_quote(args.kmer_chrx)} -cx1 intersect {shlex_quote(chrx)} -ocleft\n"
            f"{shlex_quote(kmc_tools)} -t1 simple {shlex_quote(database)} {shlex_quote(args.kmer_chrm)} -cx1 intersect {shlex_quote(chrm)} -ocleft\n"
            f"{shlex_quote(kmc_tools)} -t1 simple {shlex_quote(database)} {shlex_quote(args.kmer_auto)} -cx1 intersect {shlex_quote(auto)} -ocleft\n"
            f"{shlex_quote(kmc_tools)} -t1 transform {shlex_quote(chry)} dump {shlex_quote(chry)}\n"
            f"{shlex_quote(kmc_tools)} -t1 transform {shlex_quote(chrx)} dump {shlex_quote(chrx)}\n"
            f"{shlex_quote(kmc_tools)} -t1 transform {shlex_quote(chrm)} dump {shlex_quote(chrm)}\n"
            f"{shlex_quote(kmc_tools)} -t1 transform {shlex_quote(auto)} dump {shlex_quote(auto)}\n"
            f"rm -f {shlex_quote(str(chry) + '.kmc_pre')} {shlex_quote(str(chry) + '.kmc_suf')} {shlex_quote(str(chry) + '.tmp.kmc_pre')} {shlex_quote(str(chry) + '.tmp.kmc_suf')}\n"
            f"rm -f {shlex_quote(str(chrx) + '.kmc_pre')} {shlex_quote(str(chrx) + '.kmc_suf')} {shlex_quote(str(chrm) + '.kmc_pre')} {shlex_quote(str(chrm) + '.kmc_suf')} {shlex_quote(str(auto) + '.kmc_pre')} {shlex_quote(str(auto) + '.kmc_suf')}\n"
            f"{shlex_quote(sys.executable)} {shlex_quote(process_sex)} {shlex_quote(auto)} {shlex_quote(chry)} {shlex_quote(chrx)} {shlex_quote(chrm)} {shlex_quote(yaml_output)}\n",
            encoding="utf-8",
        )
        sex_metrics = job_tmp / "sex.io.json"
        phase_paths.append(sex_metrics)
        run_profiled(
            ["/usr/bin/bash", str(sex_script)],
            metrics_path=sex_metrics,
            label="kmer_sex_fused.sex",
            local_paths=[job_tmp],
            input_paths=[Path(str(database) + ".kmc_pre"), Path(str(database) + ".kmc_suf")],
            output_paths=[yaml_output, chry, chrx, chrm, auto],
            requested_ssd_gb=args.ssd_gb,
            threads=1,
            memory_mb=int(args.low_memory_mb),
            poll_interval=args.poll_interval,
        )

        for path in (yaml_output, chry, chrx, chrm, auto):
            if not path.is_file() or path.stat().st_size <= 0:
                raise RuntimeError(f"required fused KMC/sex output is absent or empty: {path}")
        for source, destination in (
            (chry, Path(args.output_chry)),
            (chrx, Path(args.output_chrx)),
            (chrm, Path(args.output_chrm)),
            (auto, Path(args.output_auto)),
            (yaml_output, Path(args.output_yaml)),
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
                "label": "kmer_sex_fused",
                "sample": args.sample,
                "success": success,
                "started_at_epoch": started,
                "duration_seconds": round(time.time() - started, 6),
                "requested": {
                    "initial_cores": args.initial_cores,
                    "initial_memory_mb": args.initial_memory_mb,
                    "kmc_threads": args.kmc_threads,
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


def shlex_quote(value: object) -> str:
    import shlex

    return shlex.quote(str(value))


if __name__ == "__main__":
    raise SystemExit(main())
