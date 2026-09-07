#!/usr/bin/env python3
"""Extract FASTQs from one BAM/CRAM read group and remove adapters on SSD."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import sys
import tempfile
import time
from pathlib import Path

from io_profile import run_profiled
from adapter_processing import prepare_adapters
from pipeline_runtime import (
    _atomic_copy,
    _atomic_json,
    assigned_scratch,
    executable,
    lease_preflight,
    shrink_lease,
)
from select_cram_reference import choose_reference, read_cram_header


DEFAULT_HG19_REFERENCE = "/gpfs/work3/0/qtholstg/hg38_res_v2/cram_refs/hg19.fa"
DEFAULT_HG19_B37_CHRY_REFERENCE = (
    "/gpfs/work3/0/qtholstg/marc/genome/hg19_b37chrY.fa"
)
DEFAULT_HG38_REFERENCE = (
    "/gpfs/work3/0/qtholstg/hg38_res_v2/cram_refs/"
    "GRCh38_full_analysis_set_plus_decoy_hla.fa"
)


def q(value: object) -> str:
    return shlex.quote(str(value))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-alignment", required=True)
    parser.add_argument("--cram-options", default="")
    # Defaults keep already-queued commands from a running Snakemake process
    # compatible. New workflow parses pass these paths explicitly.
    parser.add_argument("--hg19-reference", default=DEFAULT_HG19_REFERENCE)
    parser.add_argument(
        "--hg19-b37-chry-reference",
        default=DEFAULT_HG19_B37_CHRY_REFERENCE,
    )
    parser.add_argument("--hg38-reference", default=DEFAULT_HG38_REFERENCE)
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


def _configured_reference(tokens: list[str]) -> str | None:
    for index, token in enumerate(tokens):
        if token in {"--reference", "-T"}:
            if index + 1 >= len(tokens):
                raise ValueError(f"missing value after CRAM option {token}")
            return tokens[index + 1]
        if token.startswith("--reference="):
            return token.partition("=")[2]
        if token.startswith("-T") and token != "-T":
            return token[2:]
    return None


def _replace_reference(tokens: list[str], selected: str) -> list[str]:
    result = list(tokens)
    for index, token in enumerate(result):
        if token in {"--reference", "-T"}:
            if index + 1 >= len(result):
                raise ValueError(f"missing value after CRAM option {token}")
            result[index + 1] = selected
            return result
        if token.startswith("--reference="):
            result[index] = f"--reference={selected}"
            return result
        if token.startswith("-T") and token != "-T":
            result[index] = f"-T{selected}"
            return result
    result.extend(("--reference", selected))
    return result


def resolve_cram_options(
    alignment: Path,
    raw_options: str,
    *,
    samtools: str,
    hg19_reference: str,
    hg19_b37_chry_reference: str,
    hg38_reference: str,
) -> tuple[list[str], dict]:
    """Select the decode FASTA from a CRAM's M5 dictionary.

    The input/output contract stays unchanged. Only the reference argument used
    by the extraction subprocess can be corrected.
    """
    tokens = shlex.split(raw_options)
    configured = _configured_reference(tokens)
    selection = {
        "is_cram": alignment.suffix.lower() == ".cram",
        "configured_reference": configured,
        "selected_reference": configured,
        "reason": "not a CRAM; keeping configured options",
        "changed": False,
    }
    if not selection["is_cram"]:
        return tokens, selection
    if configured is None:
        raise ValueError("CRAM input has no configured --reference option")

    header = read_cram_header(str(alignment), samtools)
    selected, reason = choose_reference(
        header,
        primary_reference=configured,
        hg19_reference=hg19_reference,
        hg19_b37_chry_reference=hg19_b37_chry_reference,
        hg38_reference=hg38_reference,
    )
    if not Path(selected).is_file():
        raise FileNotFoundError(
            f"selected CRAM reference does not exist ({reason}): {selected}"
        )
    selection.update(
        {
            "selected_reference": selected,
            "reason": reason,
            "changed": selected != configured,
        }
    )
    return _replace_reference(tokens, selected), selection


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
    reference_selection: dict = {}
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
        cram_options, reference_selection = resolve_cram_options(
            alignment,
            args.cram_options,
            samtools=samtools,
            hg19_reference=args.hg19_reference,
            hg19_b37_chry_reference=args.hg19_b37_chry_reference,
            hg38_reference=args.hg38_reference,
        )
        print(
            "[external adapter reference] "
            f"{reference_selection['reason']}; "
            f"configured={reference_selection['configured_reference']}; "
            f"selected={reference_selection['selected_reference']}",
            file=sys.stderr,
        )
        cram_tokens = " ".join(q(token) for token in cram_options)
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
            lease,
            cores=args.adapter_cores,
            memory_mb=args.adapter_memory_mb,
            phase="adapter_removal",
        )
        prepared, encoding = prepare_adapters(
            args, raw1, raw2, job_tmp, phase_paths, pigz=pigz,
            adapter_removal=adapter_removal,
            label="external_adapter_fused.adapter_removal",
        )
        prepared1, prepared2, adapter_log, fastq_stats, adapters = prepared
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
                "observed_fastq_encoding": encoding if "encoding" in locals() else {},
                "reference_selection": reference_selection,
                "lease": lease,
                "phases": phases,
                "scratch_job_directory": str(job_tmp),
                "scratch_removed": not job_tmp.exists(),
            },
        )


if __name__ == "__main__":
    raise SystemExit(main())
