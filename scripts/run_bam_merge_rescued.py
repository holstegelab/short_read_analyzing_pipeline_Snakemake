#!/usr/bin/env python3
"""Restore reads from the same rescued FASTQ stream used for adapter statistics.

SAM remains on stdin/stdout. The rescued, untrimmed FASTQs reach bam_merge
through a separate inherited pipe, with no extra scratch copy or recompression.
Both the producer and merger must succeed.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def merge_rescued(r1: str, r2: str, command: list[str], *, remove_duplicates: bool = False) -> int:
    scripts = Path(__file__).resolve().parent
    processes: list[subprocess.Popen] = []
    try:
        source = subprocess.Popen(
            [sys.executable, str(scripts / "fastq_pair_rescue.py"),
             "--r1", r1, "--r2", r2, "--decompress", "external",
             "--threads", "1", "--buffer-size", "2048", "--quiet", "--summary"],
            stdout=subprocess.PIPE,
        )
        processes.append(source)
        if remove_duplicates:
            dedup = subprocess.Popen(
                [sys.executable, str(scripts / "remove_interleaved_duplicates.py"),
                 "--input", "-", "--output", "-", "--quiet"],
                stdin=source.stdout, stdout=subprocess.PIPE,
            )
            processes.append(dedup)
            source.stdout.close()
            source = dedup
        pipe = source.stdout
        try:
            result = subprocess.run(
                [*command, "--interleaved-fastq", f"/dev/fd/{pipe.fileno()}"],
                pass_fds=(pipe.fileno(),), check=False,
            )
        finally:
            # Closing our copy also gives the producer a broken pipe when the
            # merger exits early, allowing its decompressor cleanup to run.
            pipe.close()
        if result.returncode:
            return result.returncode
        for process in reversed(processes):
            code = process.wait()
            if code:
                print(f"[source rescue] producer exited with status {code}", file=sys.stderr)
                return code if code > 0 else 1
        return 0
    finally:
        for process in reversed(processes):
            if process.stdout is not None:
                process.stdout.close()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--remove-duplicates", action="store_true")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("a bam_merge command is required after --")
    return merge_rescued(args.r1, args.r2, command, remove_duplicates=args.remove_duplicates)


if __name__ == "__main__":
    raise SystemExit(main())
