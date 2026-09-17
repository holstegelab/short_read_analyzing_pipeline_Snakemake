#!/usr/bin/env python3
"""Print the scheduler-assigned scratch directory for shell-based rules."""

from __future__ import annotations

import argparse
import os

from pipeline_runtime import assigned_scratch


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fallback",
        help="shared writable fallback for rules whose ssd_use is possible",
    )
    args = parser.parse_args()
    print(os.fspath(assigned_scratch(shared_fallback=args.fallback)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
