#!/usr/bin/env python3
"""Import/CLI entry point for the shared large-write quota gate.

Snakemake shell rules run this before their existing commands. Lease changes
belong to the scheduler job, so the following command keeps the restored lease.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from fastq_quota import DEFAULT_LEASE_COMMAND, wait_for_storage_quota


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", nargs="+", required=True)
    parser.add_argument("--quota-mode", choices=("auto", "required", "disabled"), default="auto")
    parser.add_argument("--cores", type=float, required=True)
    parser.add_argument("--memory-mb", type=float, required=True)
    parser.add_argument("--lease-command", default=DEFAULT_LEASE_COMMAND)
    parser.add_argument("--phase", default="storage_work")
    args = parser.parse_args(argv)
    wait_for_storage_quota(
        [Path(path) for path in args.destination], mode=args.quota_mode,
        resume_cores=args.cores, resume_memory_mb=args.memory_mb,
        lease_command=args.lease_command, resume_phase=args.phase,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
