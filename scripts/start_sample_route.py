#!/usr/bin/env python3
"""Idempotent filesystem helpers for the routed ``start_sample`` rule."""

from __future__ import annotations

import json
import os
import tempfile
import time
from pathlib import Path
from typing import Mapping, Sequence


MANIFEST_NAME = ".start_sample_manifest.json"
SUPPORTED_ROUTES = {"active", "archive", "dcache"}


class StartSampleRouteError(RuntimeError):
    pass


def sample_route(sample: Mapping[str, object]) -> str:
    route = sample.get("from_external")
    if route in (None, False, "", "active"):
        return "active"
    route = str(route).lower()
    if route not in SUPPORTED_ROUTES:
        raise StartSampleRouteError(f"unsupported sample source route: {route!r}")
    return route


def safe_relative_path(filename: str | os.PathLike[str]) -> Path:
    value = os.fspath(filename)
    if any(character in value for character in ("\0", "\n", "\r")):
        raise StartSampleRouteError("source filename contains a control character")
    normalized = Path(os.path.normpath(value.lstrip("/")))
    if str(normalized) in ("", ".") or ".." in normalized.parts:
        raise StartSampleRouteError(f"unsafe relative source filename: {value!r}")
    return normalized


def routed_relative_path(filename: str | os.PathLike[str], route: str) -> Path:
    relative = safe_relative_path(filename)
    if route == "archive" and relative.suffix == ".bz2":
        return relative.with_suffix(".gz")
    return relative


def sample_filenames(sample: Mapping[str, object]) -> list[str]:
    result: list[str] = []
    for key in ("file1", "file2"):
        values = sample.get(key, ()) or ()
        if isinstance(values, (str, os.PathLike)):
            values = [values]
        for value in values:
            if value:
                result.append(os.fspath(value))
    if not result:
        raise StartSampleRouteError("sample has no source files")
    return result


def expected_relative_files(
    sample: Mapping[str, object], route: str
) -> list[Path]:
    if route not in {"archive", "dcache"}:
        return []
    result: list[Path] = []
    for filename in sample_filenames(sample):
        # Archive sample sheets can contain an absolute path that is already
        # active; the legacy archive rule deliberately did not copy those.
        if route == "archive" and os.path.isabs(filename):
            continue
        result.append(routed_relative_path(filename, route))
    if not result:
        raise StartSampleRouteError(
            f"{route} sample has no files that need materializing"
        )
    if len(set(result)) != len(result):
        raise StartSampleRouteError("multiple source files map to one destination")
    return result


def _read_manifest(root: Path) -> dict[str, object] | None:
    path = root / MANIFEST_NAME
    try:
        with path.open("rt", encoding="utf-8") as handle:
            payload = json.load(handle)
    except FileNotFoundError:
        return None
    except (OSError, json.JSONDecodeError) as exc:
        raise StartSampleRouteError(f"invalid route manifest {path}: {exc}") from exc
    if payload.get("schema_version") != 1:
        raise StartSampleRouteError(f"unsupported route manifest schema in {path}")
    return payload


def validate_materialized(
    root: str | os.PathLike[str],
    sample: Mapping[str, object],
    route: str,
) -> list[dict[str, object]]:
    destination = Path(root)
    if not destination.is_dir():
        raise StartSampleRouteError(f"destination directory is absent: {destination}")
    expected = expected_relative_files(sample, route)
    manifest = _read_manifest(destination)
    manifest_sizes: dict[str, int] = {}
    if manifest is not None:
        if manifest.get("route") != route:
            raise StartSampleRouteError(
                f"route manifest says {manifest.get('route')!r}, expected {route!r}"
            )
        for item in manifest.get("files", []):
            if isinstance(item, dict) and "path" in item and "bytes" in item:
                manifest_sizes[str(item["path"])] = int(item["bytes"])
        if set(manifest_sizes) != {str(path) for path in expected}:
            raise StartSampleRouteError("route manifest file set does not match sample")

    records: list[dict[str, object]] = []
    for relative in expected:
        path = destination / relative
        if not path.is_file():
            raise StartSampleRouteError(f"materialized source file is absent: {path}")
        size = path.stat().st_size
        if size <= 0:
            raise StartSampleRouteError(f"materialized source file is empty: {path}")
        recorded = manifest_sizes.get(str(relative))
        if recorded is not None and recorded != size:
            raise StartSampleRouteError(
                f"materialized file size changed for {path}: {size} != {recorded}"
            )
        records.append({"path": str(relative), "bytes": size})
    return records


def write_manifest(
    root: str | os.PathLike[str],
    *,
    sample_name: str,
    route: str,
    records: Sequence[Mapping[str, object]],
) -> Path:
    destination = Path(root)
    destination.mkdir(parents=True, exist_ok=True)
    path = destination / MANIFEST_NAME
    _write_json_atomic(
        path,
        {
            "schema_version": 1,
            "sample": sample_name,
            "route": route,
            "files": [dict(record) for record in records],
            "written_at_epoch": time.time(),
        },
    )
    return path


def quarantine_incomplete(path: str | os.PathLike[str]) -> Path | None:
    target = Path(path)
    if not target.exists() and not target.is_symlink():
        return None
    stamp = time.strftime("%Y%m%dT%H%M%S", time.gmtime())
    quarantine = target.with_name(
        f"{target.name}.incomplete.{stamp}.{os.getpid()}"
    )
    os.replace(target, quarantine)
    return quarantine


def make_partial_directory(path: str | os.PathLike[str]) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    return Path(
        tempfile.mkdtemp(prefix=f".{target.name}.partial.", dir=target.parent)
    )


def promote_directory(
    partial: str | os.PathLike[str], destination: str | os.PathLike[str]
) -> None:
    source = Path(partial)
    target = Path(destination)
    if not source.is_dir():
        raise StartSampleRouteError(f"partial directory is absent: {source}")
    if target.exists() or target.is_symlink():
        raise StartSampleRouteError(f"destination already exists: {target}")
    os.replace(source, target)


def write_completion_markers(
    *,
    started: str | os.PathLike[str],
    route_ready: str | os.PathLike[str],
    route: str,
    sample_name: str,
    files: Sequence[Mapping[str, object]],
    legacy_marker: str | os.PathLike[str] | None = None,
) -> None:
    payload = {
        "schema_version": 1,
        "sample": sample_name,
        "route": route,
        "files": [dict(item) for item in files],
        "completed_at_epoch": time.time(),
    }
    _write_json_atomic(route_ready, payload)
    if legacy_marker is not None:
        _write_text_atomic(legacy_marker, "")
    # The historical marker is written last: its presence now means the route
    # validation and the universal route-ready marker both succeeded.
    _write_json_atomic(started, payload)


def _write_text_atomic(path: str | os.PathLike[str], content: str) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{target.name}.", dir=target.parent)
    temporary_path = Path(temporary)
    try:
        with os.fdopen(fd, "wt", encoding="utf-8") as handle:
            handle.write(content)
        os.replace(temporary_path, target)
    finally:
        temporary_path.unlink(missing_ok=True)


def _write_json_atomic(path: str | os.PathLike[str], payload: object) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{target.name}.", dir=target.parent)
    temporary_path = Path(temporary)
    try:
        with os.fdopen(fd, "wt", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary_path, target)
    finally:
        temporary_path.unlink(missing_ok=True)
