#!/usr/bin/env python3
"""Create, encrypt and upload one CRAM without GPFS payload intermediates."""

from __future__ import annotations

import argparse
import json
import os
import posixpath
import shlex
import shutil
import subprocess
import sys
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
    parser.add_argument("--input-bam", required=True)
    parser.add_argument("--input-bai", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--private-key", required=True)
    parser.add_argument("--recipient-keys", nargs="+", required=True)
    parser.add_argument("--upload-script", required=True)
    parser.add_argument("--upload-config", required=True)
    parser.add_argument("--upload-remote", required=True)
    parser.add_argument("--upload-directory", required=True)
    parser.add_argument("--ada", required=True)
    parser.add_argument("--output-copied", required=True)
    parser.add_argument("--output-checksum", required=True)
    parser.add_argument("--cram-log", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--cram-threads", type=int, default=2)
    parser.add_argument("--initial-cores", type=float, required=True)
    parser.add_argument("--initial-memory-mb", type=float, required=True)
    parser.add_argument("--encrypt-cores", type=float, default=0.45)
    parser.add_argument("--encrypt-memory-mb", type=float, default=200)
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


def _require_files(paths: Sequence[str]) -> None:
    for raw in paths:
        path = Path(raw)
        if not path.is_file():
            raise FileNotFoundError(path)
        if path.stat().st_size <= 0:
            raise RuntimeError(f"required input is empty: {path}")


def _load_metrics(paths: Sequence[Path]) -> list[dict]:
    result = []
    for path in paths:
        try:
            result.append(json.loads(path.read_text(encoding="utf-8")))
        except FileNotFoundError:
            pass
    return result


def _shell_with_redirects(
    command: Sequence[str], *, stdin: Path | None = None,
    stdout: Path | None = None, stderr: Path | None = None,
) -> list[str]:
    rendered = "exec " + " ".join(_q(value) for value in command)
    if stdin is not None:
        rendered += f" < {_q(stdin)}"
    if stdout is not None:
        rendered += f" > {_q(stdout)}"
    if stderr is not None:
        rendered += f" 2> {_q(stderr)}"
    return ["/usr/bin/bash", "-c", rendered]


def main() -> int:
    args = parse_args()
    _require_files(
        [
            args.input_bam,
            args.input_bai,
            args.reference,
            args.private_key,
            args.upload_script,
            args.upload_config,
            args.ada,
        ]
        + args.recipient_keys
    )
    samtools = executable(args.samtools)

    scratch = assigned_scratch(args.scratch_base)
    parent = scratch / "cram_encrypt_fused"
    parent.mkdir(parents=True, exist_ok=True)
    job_tmp = Path(tempfile.mkdtemp(prefix=f"{args.sample}.", dir=parent))
    phase_paths: list[Path] = []
    lease: dict = {}
    success = False
    started = time.time()
    local_cram_log = job_tmp / "mCRAM.log"

    cram_name = f"{args.sample}.mapped_hg38.cram"
    local_cram = job_tmp / cram_name
    local_crai = job_tmp / (cram_name + ".crai")
    local_encrypted = job_tmp / (cram_name + ".c4gh")
    local_cram_checksum = job_tmp / (local_encrypted.name + ".ADLER32")
    local_crai_checksum = job_tmp / (local_crai.name + ".ADLER32")
    local_copied = job_tmp / (cram_name + ".copied")

    try:
        lease = lease_preflight(
            args.lease_mode,
            args.lease_command,
            initial_cores=args.initial_cores,
            initial_memory_mb=args.initial_memory_mb,
        )

        cram_metrics = job_tmp / "cram.io.json"
        phase_paths.append(cram_metrics)
        cram_command = [
            samtools,
            "view",
            "--output-fmt",
            "cram,version=3.1,archive",
            "--reference",
            args.reference,
            "-@",
            str(args.cram_threads),
            "--write-index",
            "-o",
            f"{local_cram}##idx##{local_crai}",
            args.input_bam,
        ]
        run_profiled(
            _shell_with_redirects(cram_command, stderr=local_cram_log),
            metrics_path=cram_metrics,
            label="cram_encrypt_fused.cram",
            local_paths=[job_tmp],
            input_paths=[args.input_bam],
            output_paths=[local_cram, local_crai, local_cram_log],
            requested_ssd_gb=args.ssd_gb,
            threads=args.cram_threads,
            memory_mb=int(args.initial_memory_mb),
            poll_interval=args.poll_interval,
        )
        for path in (local_cram, local_crai):
            if not path.is_file() or path.stat().st_size <= 0:
                raise RuntimeError(f"required local CRAM output is absent or empty: {path}")
        subprocess.run([samtools, "quickcheck", str(local_cram)], check=True)

        lease = shrink_lease(
            lease,
            cores=args.encrypt_cores,
            memory_mb=args.encrypt_memory_mb,
            phase="cram_encrypt",
        )
        encryption_metrics = job_tmp / "encrypt.io.json"
        phase_paths.append(encryption_metrics)
        encryption_command = [
            sys.executable,
            "-m",
            "crypt4gh",
            "encrypt",
            "--sk",
            args.private_key,
        ]
        for public_key in args.recipient_keys:
            encryption_command.extend(("--recipient_pk", public_key))
        run_profiled(
            _shell_with_redirects(
                encryption_command,
                stdin=local_cram,
                stdout=local_encrypted,
            ),
            metrics_path=encryption_metrics,
            label="cram_encrypt_fused.encrypt",
            local_paths=[job_tmp],
            input_paths=[local_cram],
            output_paths=[local_encrypted],
            requested_ssd_gb=args.ssd_gb,
            threads=1,
            memory_mb=int(args.encrypt_memory_mb),
            poll_interval=args.poll_interval,
        )
        if not local_encrypted.is_file() or local_encrypted.stat().st_size <= 0:
            raise RuntimeError(f"encrypted CRAM is absent or empty: {local_encrypted}")

        upload_directory = "/" + args.upload_directory.strip("/")
        for label, source, checksum in (
            ("cram", local_encrypted, local_cram_checksum),
            ("crai", local_crai, local_crai_checksum),
        ):
            upload_metrics = job_tmp / f"upload_{label}.io.json"
            phase_paths.append(upload_metrics)
            upload_command = [
                sys.executable,
                args.upload_script,
                "upload",
                "--config",
                args.upload_config,
                "--remote",
                args.upload_remote,
                "--source",
                str(source),
                "--destination",
                posixpath.join(upload_directory, source.name),
                "--checksum-output",
                str(checksum),
                "--ada",
                args.ada,
            ]
            run_profiled(
                upload_command,
                metrics_path=upload_metrics,
                label=f"cram_encrypt_fused.upload_{label}",
                local_paths=[job_tmp],
                input_paths=[source],
                output_paths=[checksum],
                requested_ssd_gb=args.ssd_gb,
                threads=1,
                memory_mb=int(args.encrypt_memory_mb),
                poll_interval=args.poll_interval,
            )
            if not checksum.is_file() or checksum.stat().st_size <= 0:
                raise RuntimeError(
                    f"verified {label} upload produced no checksum: {checksum}"
                )

        # Publish only tiny durable receipts, checksum first and completion
        # marker last. Plaintext CRAM, CRAI and encrypted CRAM never leave the
        # assigned SSD; absence of .copied makes every partial remote upload
        # safe to overwrite on retry.
        _atomic_copy(local_cram_checksum, Path(args.output_checksum))
        local_copied.touch()
        _atomic_copy(local_copied, Path(args.output_copied))
        success = True
        return 0
    finally:
        if local_cram_log.is_file():
            try:
                _atomic_copy(local_cram_log, Path(args.cram_log))
            except Exception:
                pass
        phases = _load_metrics(phase_paths)
        shutil.rmtree(job_tmp, ignore_errors=True)
        _atomic_json(
            Path(args.metrics),
            {
                "schema_version": 1,
                "label": "cram_encrypt_fused",
                "sample": args.sample,
                "success": success,
                "started_at_epoch": started,
                "duration_seconds": round(time.time() - started, 6),
                "requested": {
                    "initial_cores": args.initial_cores,
                    "initial_memory_mb": args.initial_memory_mb,
                    "encrypt_cores": args.encrypt_cores,
                    "encrypt_memory_mb": args.encrypt_memory_mb,
                    "cram_threads": args.cram_threads,
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
