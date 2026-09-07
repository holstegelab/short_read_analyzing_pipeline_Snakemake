#!/usr/bin/env python3
"""Native paired-FASTQ wrapper around the shared adapter-processing phase."""
from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

from adapter_processing import prepare_adapters
from pipeline_runtime import _atomic_json, _atomic_publish, executable


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "input-forward", "input-reverse", "sample", "readgroup", "adapter-list",
        "fastq-stats-script", "remove-duplicates-script", "pair-rescue-script",
        "error-file", "output-forward", "output-reverse", "output-adapter-log",
        "output-fastq-stats", "output-adapters", "metrics",
    ):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--remove-duplicated-reads", type=int, choices=(0, 1), default=0)
    parser.add_argument("--attempt", type=int, default=1)
    parser.add_argument("--adapter-memory-mb", type=float, default=512)
    parser.add_argument("--pigz", default="pigz")
    parser.add_argument("--adapter-removal", default="AdapterRemoval")
    parser.add_argument("--poll-interval", type=float, default=5.0)
    args = parser.parse_args()
    # This route never required an SSD or a phase-changing scheduler lease.
    args.ssd_gb = 0
    return args


def main():
    args = parse_args()
    raw1, raw2 = Path(args.input_forward), Path(args.input_reverse)
    for path in (raw1, raw2, Path(args.adapter_list), Path(args.fastq_stats_script),
                 Path(args.remove_duplicates_script), Path(args.pair_rescue_script)):
        if not path.is_file():
            raise FileNotFoundError(path)
    pigz = executable(args.pigz)
    adapter_removal = executable(args.adapter_removal)
    parent = Path(args.output_forward).parent
    parent.mkdir(parents=True, exist_ok=True)
    phases = []
    encoding = {}
    success = False
    # Use the output filesystem: no extra SSD requirement or cross-disk copy.
    with tempfile.TemporaryDirectory(prefix=".adapter-", dir=parent) as temporary:
        job_tmp = Path(temporary)
        try:
            prepared, encoding = prepare_adapters(
                args, raw1, raw2, job_tmp, phases, pigz=pigz,
                adapter_removal=adapter_removal, label="adapter_removal",
            )
            for source, destination in zip(prepared, (
                args.output_forward, args.output_reverse, args.output_adapter_log,
                args.output_fastq_stats, args.output_adapters,
            )):
                _atomic_publish(source, Path(destination))
            success = True
        finally:
            _atomic_json(Path(args.metrics), {
                "schema_version": 1, "label": "adapter_removal",
                "sample": args.sample, "readgroup": args.readgroup,
                "attempt": args.attempt, "success": success,
                "observed_fastq_encoding": encoding,
                "phases": [json.loads(p.read_text()) for p in phases if p.is_file()],
            })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
