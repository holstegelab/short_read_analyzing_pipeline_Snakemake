#!/usr/bin/env python3
"""Lightweight process-tree and local-scratch profiling for pipeline jobs."""

from __future__ import annotations

import json
import os
import shutil
import signal
import socket
import subprocess
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def path_bytes(path: str | os.PathLike[str]) -> int | None:
    """Return allocated object size, or None when the path is absent."""
    target = Path(path)
    try:
        if target.is_symlink() or target.is_file():
            return target.stat().st_size
        if target.is_dir():
            total = 0
            for root, _, filenames in os.walk(target):
                for filename in filenames:
                    try:
                        total += (Path(root) / filename).stat().st_size
                    except FileNotFoundError:
                        pass
            return total
    except FileNotFoundError:
        pass
    return None


def _path_sizes(paths: Iterable[str | os.PathLike[str]]) -> list[dict[str, object]]:
    return [{"path": os.fspath(path), "bytes": path_bytes(path)} for path in paths]


def _filesystem_sample(path: str | os.PathLike[str]) -> dict[str, int | str] | None:
    target = Path(path)
    while not target.exists() and target != target.parent:
        target = target.parent
    try:
        stats = os.statvfs(target)
    except OSError:
        return None
    return {
        "path": str(target.resolve()),
        "total_bytes": stats.f_blocks * stats.f_frsize,
        "free_bytes": stats.f_bavail * stats.f_frsize,
        "used_bytes": (stats.f_blocks - stats.f_bfree) * stats.f_frsize,
    }


def _status_values(pid: int) -> tuple[int | None, dict[str, int]]:
    parent = None
    values: dict[str, int] = {}
    try:
        with open(f"/proc/{pid}/status", "rt", encoding="utf-8") as handle:
            for line in handle:
                key, _, raw = line.partition(":")
                if key == "PPid":
                    parent = int(raw.strip())
                elif key in {"VmRSS", "VmHWM"}:
                    values[key] = int(raw.split()[0]) * 1024
    except (FileNotFoundError, ProcessLookupError, PermissionError, ValueError):
        pass
    return parent, values


def _process_tree(root_pid: int) -> set[int]:
    children: dict[int, list[int]] = {}
    live: set[int] = set()
    try:
        entries = os.listdir("/proc")
    except OSError:
        return set()
    for name in entries:
        if not name.isdigit():
            continue
        pid = int(name)
        parent, _ = _status_values(pid)
        if parent is not None:
            live.add(pid)
            children.setdefault(parent, []).append(pid)
    result: set[int] = set()
    pending = [root_pid]
    while pending:
        pid = pending.pop()
        if pid in result or pid not in live:
            continue
        result.add(pid)
        pending.extend(children.get(pid, ()))
    return result


def _process_sample(root_pid: int) -> dict[str, int]:
    totals = {
        "processes": 0,
        "rss_bytes": 0,
        "high_water_rss_bytes": 0,
        "pss_bytes": 0,
        "read_bytes": 0,
        "write_bytes": 0,
        "read_chars": 0,
        "write_chars": 0,
    }
    for pid in _process_tree(root_pid):
        _, status = _status_values(pid)
        totals["processes"] += 1
        totals["rss_bytes"] += status.get("VmRSS", 0)
        totals["high_water_rss_bytes"] += status.get("VmHWM", 0)
        try:
            with open(f"/proc/{pid}/smaps_rollup", "rt", encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith("Pss:"):
                        totals["pss_bytes"] += int(line.split()[1]) * 1024
                        break
        except (FileNotFoundError, ProcessLookupError, PermissionError, ValueError):
            pass
        try:
            with open(f"/proc/{pid}/io", "rt", encoding="utf-8") as handle:
                io_values = {}
                for line in handle:
                    key, _, raw = line.partition(":")
                    io_values[key] = int(raw.strip())
            totals["read_bytes"] += io_values.get("read_bytes", 0)
            totals["write_bytes"] += io_values.get("write_bytes", 0)
            totals["read_chars"] += io_values.get("rchar", 0)
            totals["write_chars"] += io_values.get("wchar", 0)
        except (FileNotFoundError, ProcessLookupError, PermissionError, ValueError):
            pass
    return totals


def _clean_paths(paths: Iterable[str | os.PathLike[str]]) -> list[dict[str, object]]:
    results = []
    for raw_path in paths:
        path = Path(raw_path)
        result: dict[str, object] = {"path": str(path), "removed": False}
        try:
            resolved = path.resolve()
            if resolved == Path("/") or len(resolved.parts) < 3:
                raise ValueError(f"refusing unsafe cleanup path: {resolved}")
            if path.is_dir() and not path.is_symlink():
                shutil.rmtree(path)
            else:
                path.unlink(missing_ok=True)
            result["removed"] = not path.exists()
        except Exception as exc:  # cleanup diagnostics must survive the job
            result["error"] = f"{type(exc).__name__}: {exc}"
        results.append(result)
    return results


def _write_json_atomic(path: str | os.PathLike[str], payload: dict[str, object]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{target.name}.", dir=target.parent)
    try:
        with os.fdopen(fd, "wt", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, target)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def run_profiled(
    command: Sequence[str],
    *,
    metrics_path: str | os.PathLike[str],
    label: str,
    local_paths: Sequence[str | os.PathLike[str]],
    input_paths: Sequence[str | os.PathLike[str]] = (),
    output_paths: Sequence[str | os.PathLike[str]] = (),
    cleanup_paths: Sequence[str | os.PathLike[str]] = (),
    requested_ssd_gb: float | None = None,
    threads: int | None = None,
    memory_mb: int | None = None,
    poll_interval: float = 5.0,
) -> int:
    """Run *command*, sample its process tree/local paths, and write JSON metrics."""
    if not command:
        raise ValueError("command must not be empty")
    start_wall = time.time()
    start_monotonic = time.monotonic()
    started_at = _utc_now()
    local_start = _path_sizes(local_paths)
    fs_start = _filesystem_sample(local_paths[0]) if local_paths else None
    input_start = _path_sizes(input_paths)
    peaks = {
        "local_bytes": sum(item["bytes"] or 0 for item in local_start),
        "processes": 0,
        "rss_bytes": 0,
        "high_water_rss_bytes": 0,
        "pss_bytes": 0,
        "read_bytes": 0,
        "write_bytes": 0,
        "read_chars": 0,
        "write_chars": 0,
    }
    min_fs_free = fs_start["free_bytes"] if fs_start else None
    samples = 0
    return_code = None
    launch_error = None
    received_signal = None
    proc: subprocess.Popen[str] | None = None
    old_handlers: dict[int, object] = {}

    def forward_signal(signum: int, _frame: object) -> None:
        nonlocal received_signal
        received_signal = signum
        if proc is not None and proc.poll() is None:
            try:
                os.killpg(proc.pid, signum)
            except ProcessLookupError:
                pass

    try:
        for signum in (signal.SIGTERM, signal.SIGINT):
            old_handlers[signum] = signal.signal(signum, forward_signal)
        proc = subprocess.Popen(list(command), start_new_session=True, text=True)
        while True:
            process = _process_sample(proc.pid)
            local_size = sum(path_bytes(path) or 0 for path in local_paths)
            filesystem = _filesystem_sample(local_paths[0]) if local_paths else None
            samples += 1
            peaks["local_bytes"] = max(peaks["local_bytes"], local_size)
            for key, value in process.items():
                peaks[key] = max(peaks[key], value)
            if filesystem is not None:
                free = filesystem["free_bytes"]
                min_fs_free = free if min_fs_free is None else min(min_fs_free, free)
            return_code = proc.poll()
            if return_code is not None:
                break
            try:
                return_code = proc.wait(timeout=max(0.05, poll_interval))
                break
            except subprocess.TimeoutExpired:
                pass
    except Exception as exc:
        launch_error = f"{type(exc).__name__}: {exc}"
        if proc is not None and proc.poll() is None:
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            return_code = proc.wait()
    finally:
        for signum, handler in old_handlers.items():
            signal.signal(signum, handler)

    local_end_before_cleanup = _path_sizes(local_paths)
    fs_end = _filesystem_sample(local_paths[0]) if local_paths else None
    cleanup = _clean_paths(cleanup_paths)
    cleanup_ok = all(item.get("removed") for item in cleanup)
    ended_at = _utc_now()
    payload: dict[str, object] = {
        "schema_version": 1,
        "label": label,
        "hostname": socket.gethostname(),
        "pid": os.getpid(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_JOBID"),
        "zslurm_job_id": os.environ.get("ZSLURM_JOB_ID"),
        "started_at": started_at,
        "ended_at": ended_at,
        "duration_seconds": round(time.monotonic() - start_monotonic, 6),
        "wall_start_epoch": start_wall,
        "command": list(command),
        "return_code": return_code,
        "received_signal": received_signal,
        "launch_error": launch_error,
        "requested": {
            "ssd_gb": requested_ssd_gb,
            "threads": threads,
            "memory_mb": memory_mb,
        },
        "inputs_start": input_start,
        "outputs_end": _path_sizes(output_paths),
        "local_start": local_start,
        "local_end_before_cleanup": local_end_before_cleanup,
        "peaks": peaks,
        "filesystem_start": fs_start,
        "filesystem_end": fs_end,
        "filesystem_min_free_bytes": min_fs_free,
        "samples": samples,
        "cleanup": cleanup,
        "cleanup_ok": cleanup_ok,
    }
    _write_json_atomic(metrics_path, payload)
    if launch_error is not None:
        raise RuntimeError(launch_error)
    if return_code:
        raise subprocess.CalledProcessError(return_code, list(command))
    return 0
