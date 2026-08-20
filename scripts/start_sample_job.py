#!/usr/bin/env python3
"""Execute one routed ``start_sample`` materialization job.

The Snakemake rules deliberately stay route-specific so that archive and
dCache materializations can be declared as real temporary directory outputs.
This module keeps the transfer implementation shared between those rules.
"""

from __future__ import annotations

import bz2
import gzip
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Callable, Mapping
from urllib.parse import urlsplit

from start_sample_route import (
    StartSampleRouteError,
    make_partial_directory,
    promote_directory,
    quarantine_incomplete,
    routed_relative_path,
    safe_relative_path,
    sample_filenames,
    sample_route,
    validate_materialized,
    write_completion_markers,
    write_manifest,
)


def run_start_sample_job(
    *,
    sample: Mapping[str, object],
    sample_name: str,
    expected_route: str,
    destination: str | os.PathLike[str] | None,
    started: str | os.PathLike[str],
    route_ready: str | os.PathLike[str],
    append_prefix: Callable[[str, str], str],
    dcache_source_file: Callable[
        [Mapping[str, object], str], tuple[str, str]
    ],
    transfer_script: str | os.PathLike[str],
    source_dir: str | os.PathLike[str],
    dcache_download_workers: int,
    dcache_download_lock_slots: int,
    aws_cli: str = "aws",
    s3_max_attempts: int = 6,
    s3_initial_backoff_seconds: float = 30,
    s3_max_backoff_seconds: float = 300,
) -> None:
    """Validate or materialize one sample, then publish atomic markers."""

    route = sample_route(sample)
    if route != expected_route:
        raise StartSampleRouteError(
            f"sample {sample_name} selected {route!r}, but rule expects "
            f"{expected_route!r}"
        )

    records: list[dict[str, object]] = []
    materialized = Path(destination) if destination is not None else None

    if route == "active":
        if materialized is not None:
            raise StartSampleRouteError(
                "active route must not own a temporary materialization directory"
            )
        for filename in sample_filenames(sample):
            source = Path(append_prefix(str(sample["prefix"]), filename))
            if not source.is_file():
                raise FileNotFoundError(
                    f"active source file is absent for {sample_name}: {source}"
                )
            records.append({"path": str(source), "bytes": source.stat().st_size})
    else:
        if materialized is None:
            raise StartSampleRouteError(
                f"{route} route requires a materialization directory output"
            )
        try:
            records = validate_materialized(materialized, sample, route)
            print(
                f"[start_sample] adopting validated {route} destination "
                f"for {sample_name}: {materialized}",
                flush=True,
            )
        except StartSampleRouteError as existing_error:
            if materialized.exists() or materialized.is_symlink():
                quarantine = quarantine_incomplete(materialized)
                print(
                    f"[start_sample] quarantined incomplete destination "
                    f"{quarantine}: {existing_error}",
                    file=sys.stderr,
                    flush=True,
                )
            partial = make_partial_directory(materialized)
            try:
                if route == "archive":
                    _copy_archive_sample(
                        sample=sample,
                        partial=partial,
                        append_prefix=append_prefix,
                    )
                elif route == "dcache":
                    _download_dcache_sample(
                        sample=sample,
                        sample_name=sample_name,
                        partial=partial,
                        dcache_source_file=dcache_source_file,
                        transfer_script=transfer_script,
                        source_dir=source_dir,
                        download_workers=dcache_download_workers,
                        download_lock_slots=dcache_download_lock_slots,
                    )
                elif route == "s3":
                    _download_s3_sample(
                        sample=sample,
                        sample_name=sample_name,
                        partial=partial,
                        append_prefix=append_prefix,
                        aws_cli=aws_cli,
                        max_attempts=s3_max_attempts,
                        initial_backoff_seconds=s3_initial_backoff_seconds,
                        max_backoff_seconds=s3_max_backoff_seconds,
                    )
                else:  # Protected by sample_route(), retained as a hard guard.
                    raise StartSampleRouteError(
                        f"no materializer for source route {route!r}"
                    )

                records = validate_materialized(partial, sample, route)
                write_manifest(
                    partial,
                    sample_name=sample_name,
                    route=route,
                    records=records,
                )
                promote_directory(partial, materialized)
            except Exception:
                if partial.exists() or partial.is_symlink():
                    quarantine_incomplete(partial)
                raise

            records = validate_materialized(materialized, sample, route)

        # Adopted legacy directories may not yet contain a manifest.
        write_manifest(
            materialized,
            sample_name=sample_name,
            route=route,
            records=records,
        )

        if route == "archive":
            _release_archive_sources(sample, append_prefix)

    write_completion_markers(
        started=started,
        route_ready=route_ready,
        route=route,
        sample_name=sample_name,
        files=records,
    )


def _copy_archive_sample(
    *,
    sample: Mapping[str, object],
    partial: Path,
    append_prefix: Callable[[str, str], str],
) -> None:
    for filename in sample_filenames(sample):
        if os.path.isabs(filename):
            source = Path(filename)
            if not source.is_file():
                raise FileNotFoundError(source)
            continue
        source_value = append_prefix(str(sample["prefix"]), filename)
        if ":/" in source_value:
            source_value = source_value.split(":/", 1)[1]
        source = Path(source_value)
        relative = safe_relative_path(filename)
        copied = partial / relative
        copied.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["rsync", "--size-only", "--partial", str(source), str(copied)],
            check=True,
        )
        if copied.suffix == ".bz2":
            converted = partial / routed_relative_path(filename, "archive")
            temporary_gz = converted.with_name(
                converted.name + f".tmp.{os.getpid()}"
            )
            converted.parent.mkdir(parents=True, exist_ok=True)
            with bz2.open(copied, "rb") as source_handle, gzip.open(
                temporary_gz, "wb", compresslevel=1
            ) as destination_handle:
                while chunk := source_handle.read(16 * 1024 * 1024):
                    destination_handle.write(chunk)
            os.replace(temporary_gz, converted)
            copied.unlink()


def _download_dcache_sample(
    *,
    sample: Mapping[str, object],
    sample_name: str,
    partial: Path,
    dcache_source_file: Callable[
        [Mapping[str, object], str], tuple[str, str]
    ],
    transfer_script: str | os.PathLike[str],
    source_dir: str | os.PathLike[str],
    download_workers: int,
    download_lock_slots: int,
) -> None:
    remote = sample.get("source_remote")
    config_path = sample.get("source_config")
    if not remote or not config_path:
        raise ValueError(f"Missing dCache remote/config for {sample_name}")

    rows: list[tuple[str, Path]] = []
    for filename in sample_filenames(sample):
        file_remote, remote_path = dcache_source_file(sample, filename)
        if file_remote != remote:
            raise ValueError(
                f"Mixed dCache remotes for {sample_name}: "
                f"{remote!r} and {file_remote!r}"
            )
        rows.append((remote_path, partial / safe_relative_path(filename)))

    fd, list_path = tempfile.mkstemp(
        prefix=f".{sample_name}.dcache-download-",
        suffix=".tsv",
        dir=source_dir,
        text=True,
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            for remote_path, local_path in rows:
                handle.write(f"{remote_path}\t{local_path}\n")
        subprocess.run(
            [
                sys.executable,
                str(transfer_script),
                "download",
                "--config",
                str(config_path),
                "--remote",
                str(remote),
                "--file-list",
                list_path,
                "--workers",
                str(min(len(rows), max(1, int(download_workers)))),
                "--download-lock-slots",
                str(max(1, int(download_lock_slots))),
                "--no-stage",
            ],
            check=True,
        )
    finally:
        try:
            os.unlink(list_path)
        except FileNotFoundError:
            pass


def _validated_s3_uri(value: str) -> str:
    """Validate an S3 object URI without accepting embedded authentication."""
    if any(character in value for character in ("\0", "\n", "\r")):
        raise StartSampleRouteError("S3 URI contains a control character")
    try:
        parsed = urlsplit(value)
        port = parsed.port
    except ValueError as exc:
        raise StartSampleRouteError(f"invalid S3 URI {value!r}: {exc}") from exc
    if parsed.scheme != "s3" or not parsed.netloc or not parsed.path.lstrip("/"):
        raise StartSampleRouteError(f"invalid S3 object URI: {value!r}")
    if parsed.username or parsed.password or port is not None:
        raise StartSampleRouteError(
            "S3 source URI must not contain credentials or a port"
        )
    if parsed.query or parsed.fragment:
        raise StartSampleRouteError(
            "S3 source URI must not contain a query or fragment"
        )
    if any(
        character
        not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.-"
        for character in parsed.netloc
    ):
        raise StartSampleRouteError(f"invalid S3 bucket: {parsed.netloc!r}")
    if any(part == ".." for part in parsed.path.split("/")):
        raise StartSampleRouteError(f"S3 object URI escapes its root: {value!r}")
    return f"s3://{parsed.netloc}/{parsed.path.lstrip('/')}"


def _download_s3_sample(
    *,
    sample: Mapping[str, object],
    sample_name: str,
    partial: Path,
    append_prefix: Callable[[str, str], str],
    aws_cli: str,
    max_attempts: int,
    initial_backoff_seconds: float,
    max_backoff_seconds: float,
) -> None:
    """Download requester-pays objects with bounded exponential backoff.

    The command intentionally uses the normal AWS credential chain.  Anonymous
    access is not enabled for the NIAGADS bucket, and credentials must never be
    serialized into a sample sheet or command-line argument.
    """
    max_attempts = max(1, int(max_attempts))
    initial_backoff_seconds = max(0.0, float(initial_backoff_seconds))
    max_backoff_seconds = max(0.0, float(max_backoff_seconds))
    environment = os.environ.copy()
    environment.setdefault("AWS_EC2_METADATA_DISABLED", "true")
    environment.setdefault("AWS_RETRY_MODE", "adaptive")
    environment.setdefault("AWS_MAX_ATTEMPTS", "10")
    environment.setdefault("AWS_PAGER", "")

    for filename in sample_filenames(sample):
        source_uri = _validated_s3_uri(
            append_prefix(str(sample["prefix"]), filename)
        )
        destination = partial / safe_relative_path(filename)
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(
            f".{destination.name}.s3-part-{os.getpid()}"
        )
        last_error: Exception | None = None

        for attempt in range(1, max_attempts + 1):
            temporary.unlink(missing_ok=True)
            command = [
                str(aws_cli),
                "--cli-connect-timeout",
                "30",
                "--cli-read-timeout",
                "0",
                "s3",
                "cp",
                source_uri,
                str(temporary),
                "--request-payer",
                "requester",
                "--only-show-errors",
            ]
            try:
                subprocess.run(command, check=True, env=environment)
                if not temporary.is_file() or temporary.stat().st_size <= 0:
                    raise StartSampleRouteError(
                        f"AWS CLI returned success without a non-empty object "
                        f"for {sample_name}: {source_uri}"
                    )
                os.replace(temporary, destination)
                last_error = None
                break
            except (OSError, subprocess.CalledProcessError, StartSampleRouteError) as exc:
                last_error = exc
                temporary.unlink(missing_ok=True)
                if attempt >= max_attempts:
                    break
                delay = min(
                    max_backoff_seconds,
                    initial_backoff_seconds * (2 ** (attempt - 1)),
                )
                print(
                    f"[start_sample] S3 download attempt {attempt}/"
                    f"{max_attempts} failed for {sample_name}; retrying in "
                    f"{delay:g}s: {exc}",
                    file=sys.stderr,
                    flush=True,
                )
                time.sleep(delay)

        if last_error is not None:
            raise StartSampleRouteError(
                f"S3 download failed after {max_attempts} attempts for "
                f"{sample_name}: {source_uri}"
            ) from last_error


def _release_archive_sources(
    sample: Mapping[str, object],
    append_prefix: Callable[[str, str], str],
) -> None:
    for filename in sample_filenames(sample):
        if os.path.isabs(filename):
            continue
        source_value = append_prefix(str(sample["prefix"]), filename)
        if ":/" in source_value:
            source_value = source_value.split(":/", 1)[1]
        subprocess.run(
            ["/opt/dacommands/bin/darelease", source_value],
            check=False,
        )
