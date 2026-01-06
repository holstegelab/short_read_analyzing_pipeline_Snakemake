#!/usr/bin/env python3

import argparse
import json
import re
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import stat


_PIN_KEY_RE = re.compile(
    r"(?i)(?:\bpin\b|sticky).*(?:expir|expire|until|valid|life)|(?:expir|expire|until|valid).*(?:\bpin\b|sticky)"
)


def _parse_epoch(value: object) -> datetime | None:
    if isinstance(value, bool):
        return None

    if isinstance(value, (int, float)):
        v = float(value)
        if v <= 0:
            return None

        if v >= 1e17:
            ts = v / 1e9
        elif v >= 1e14:
            ts = v / 1e6
        elif v >= 1e11:
            ts = v / 1e3
        elif v >= 1e9:
            ts = v
        else:
            return None

        try:
            return datetime.fromtimestamp(ts, tz=timezone.utc)
        except (OverflowError, OSError, ValueError):
            return None
        return None

    if isinstance(value, str):
        s = value.strip()
        if s.isdigit():
            return _parse_epoch(int(s))
        if s.endswith("Z"):
            s = s[:-1] + "+00:00"
        try:
            dt = datetime.fromisoformat(s)
        except ValueError:
            return None
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return dt.astimezone(timezone.utc)

    return None


def _walk_for_pin_times(obj: object) -> list[datetime]:
    out: list[datetime] = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(k, str) and _PIN_KEY_RE.search(k):
                dt = _parse_epoch(v)
                if dt is not None:
                    out.append(dt)
            out.extend(_walk_for_pin_times(v))
    elif isinstance(obj, list):
        for v in obj:
            out.extend(_walk_for_pin_times(v))
    return out


def _format_timedelta_seconds(seconds: int) -> str:
    if seconds < 0:
        return "expired"

    days, rem = divmod(seconds, 86_400)
    hours, rem = divmod(rem, 3600)
    minutes, sec = divmod(rem, 60)

    if days > 0:
        return f"{days}d{hours:02}h"
    if hours > 0:
        return f"{hours}h{minutes:02}m"
    if minutes > 0:
        return f"{minutes}m{sec:02}s"
    return f"{sec}s"


def _extract_pin_display(entry: dict) -> str:
    for k in ("pinLifetime", "pin_lifetime", "stickyLifetime", "sticky_lifetime"):
        if k in entry and isinstance(entry[k], (str, int, float)):
            return str(entry[k])

    times = _walk_for_pin_times(entry)
    if not times:
        return "-"

    now = datetime.now(tz=timezone.utc)
    future = sorted(t for t in times if t >= now)
    if future:
        seconds_left = int((future[0] - now).total_seconds())
        return _format_timedelta_seconds(seconds_left)

    latest = max(times)
    seconds_left = int((latest - now).total_seconds())
    return _format_timedelta_seconds(seconds_left)


def _format_mtime_ms(mtime_ms: object) -> str:
    dt = _parse_epoch(mtime_ms)
    if dt is None:
        return "-"
    return dt.strftime("%Y-%m-%d %H:%M UTC")


def _filemode(file_type: str | None, mode_value: object) -> str:
    if not isinstance(mode_value, int):
        mode_value = 0

    if file_type == "DIR":
        return stat.filemode(stat.S_IFDIR | mode_value)
    if file_type == "LINK":
        return stat.filemode(stat.S_IFLNK | mode_value)
    return stat.filemode(stat.S_IFREG | mode_value)


@dataclass(frozen=True)
class _Row:
    mode: str
    nlink: str
    owner: str
    group: str
    size: str
    mtime: str
    pin: str
    qos: str
    locality: str
    name: str


def _render_table(rows: list[_Row]) -> None:
    cols = [
        [r.mode for r in rows],
        [r.nlink for r in rows],
        [r.owner for r in rows],
        [r.group for r in rows],
        [r.size for r in rows],
        [r.mtime for r in rows],
        [r.pin for r in rows],
        [r.qos for r in rows],
        [r.locality for r in rows],
    ]
    widths = [max((len(x) for x in col), default=0) for col in cols]

    for r in rows:
        left = (
            f"{r.mode:<{widths[0]}} "
            f"{r.nlink:>{widths[1]}} "
            f"{r.owner:>{widths[2]}} "
            f"{r.group:>{widths[3]}} "
            f"{r.size:>{widths[4]}} "
            f"{r.mtime:<{widths[5]}} "
            f"{r.pin:>{widths[6]}} "
            f"{r.qos:<{widths[7]}} "
            f"{r.locality:<{widths[8]}}"
        )
        print(f"{left} {r.name}")


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


def main() -> None:
    p = argparse.ArgumentParser(
        description="Call ./ada --stat and print an ls -l style listing with pin lifetime.",
    )
    p.add_argument("path", help="dCache path to list, e.g. /uploads")
    p.add_argument("--tokenfile", default=_default_tokenfile(), help="Token file to pass to ada")
    p.add_argument("--api", default=None, help="API base URL (passed to ada)")
    p.add_argument("--ada", default=_default_ada_path(), help="Path to ada executable")

    args = p.parse_args()

    cmd = [args.ada]
    if args.tokenfile:
        cmd += ["--tokenfile", str(Path(args.tokenfile).expanduser())]
    if args.api:
        cmd += ["--api", args.api]
    cmd += ["--stat", args.path]

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        raise subprocess.CalledProcessError(
            proc.returncode, cmd, output=proc.stdout, stderr=proc.stderr
        )

    data = json.loads(proc.stdout)

    entries: list[dict]
    if isinstance(data, dict) and isinstance(data.get("children"), list):
        entries = [e for e in data["children"] if isinstance(e, dict)]
    elif isinstance(data, dict):
        entries = [data]
    else:
        raise TypeError(f"Unexpected JSON top-level type: {type(data)}")

    rows: list[_Row] = []
    for e in entries:
        file_type = e.get("fileType")
        name = e.get("fileName")
        if not isinstance(name, str):
            name = Path(args.path).name or args.path

        if file_type == "DIR" and not name.endswith("/"):
            name = f"{name}/"

        rows.append(
            _Row(
                mode=_filemode(file_type if isinstance(file_type, str) else None, e.get("mode")),
                nlink=str(e.get("nlink", "-")),
                owner=str(e.get("owner", "-")),
                group=str(e.get("group", "-")),
                size=str(e.get("size", "-")),
                mtime=_format_mtime_ms(e.get("mtime")),
                pin=_extract_pin_display(e),
                qos=str(e.get("currentQos", "-")),
                locality=str(e.get("fileLocality", "-")),
                name=name,
            )
        )

    rows.sort(key=lambda r: r.name)
    _render_table(rows)


if __name__ == "__main__":
    main()
