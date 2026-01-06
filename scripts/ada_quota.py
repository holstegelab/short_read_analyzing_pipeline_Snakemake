#!/usr/bin/env python3

import argparse
import json
import subprocess
import sys
from pathlib import Path


def _default_ada_path() -> str:
    repo_root = Path(__file__).resolve().parents[1]
    return str(repo_root / "ada")


def _default_tokenfile() -> str:
    cands = [
        Path("~/macaroons/agh_full_snellius.conf").expanduser(),
        Path("~/macaroon/agh_full_snellius.conf").expanduser(),
    ]
    for p in cands:
        if p.exists():
            return str(p)
    return str(cands[0])


def _tb(x: int) -> float:
    return float(x) / 1e12


def main() -> None:
    p = argparse.ArgumentParser(description="Show ada --space quota summary in TB, including pinned.")
    p.add_argument(
        "--space",
        default="agh_rwtapepools",
        help="Poolgroup name for ada --space (default: agh_rwtapepools)",
    )
    p.add_argument("--tokenfile", default=_default_tokenfile(), help="Token file to pass to ada")
    p.add_argument("--api", default=None, help="API base URL (passed to ada)")
    p.add_argument("--ada", default=_default_ada_path(), help="Path to ada executable")
    args = p.parse_args()

    cmd = [args.ada]
    if args.tokenfile:
        cmd += ["--tokenfile", str(Path(args.tokenfile).expanduser())]
    if args.api:
        cmd += ["--api", args.api]
    cmd += ["--space", args.space]

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise subprocess.CalledProcessError(proc.returncode, cmd, output=proc.stdout, stderr=proc.stderr)

    data = json.loads(proc.stdout)
    if not isinstance(data, dict):
        raise TypeError(f"Unexpected JSON top-level type: {type(data)}")

    for k in ("total", "free", "precious", "removable"):
        if k not in data:
            raise KeyError(f"Missing key in ada --space output: {k}")
        if not isinstance(data[k], int):
            raise TypeError(f"Expected integer bytes for {k}, got {type(data[k])}")

    total = data["total"]
    free = data["free"]
    precious = data["precious"]
    removable = data["removable"]
    pinned = total - free - precious - removable

    print(f"space {args.space} (TB, 1 TB = 10^12 bytes)")
    print(f"total\t{_tb(total):.3f}")
    print(f"free\t{_tb(free):.3f}")
    print(f"precious\t{_tb(precious):.3f}")
    print(f"removable\t{_tb(removable):.3f}")
    print(f"pinned\t{_tb(pinned):.3f}")


if __name__ == "__main__":
    main()
