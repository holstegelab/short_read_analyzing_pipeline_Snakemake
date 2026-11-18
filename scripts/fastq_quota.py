"""Shared best-effort project-quota backoff for large writes on Snellius.

No shared locks, reservations or scheduler changes: each producer checks only
its destination fileset, and staggered retries spread recovery over time.
"""
from __future__ import annotations

import os
import random
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Sequence

from pipeline_runtime import LEASE_ENV, acquire_lease, lease_preflight, shrink_lease


MMLSQUOTA = "/usr/lpp/mmfs/bin/mmlsquota"
DEFAULT_LEASE_COMMAND = str(Path(__file__).resolve().with_name("zslurm_lease_client.py"))


class QuotaReadError(RuntimeError):
    """Quota could not be established; do not treat this as free space."""


def project_quota_target(destination: Path) -> tuple[str, str] | None:
    """Use the project/mount mapping from Snellius's prjspc-quota wrapper."""
    match = re.match(r"^/gpfs/work([1-5])/[0-2]/([^/]+)(?:/|$)",
                     str(destination.resolve()))
    if not match:
        return None
    number, fileset = match.groups()
    prefix = "lstor" if number == "5" else "wstor"
    return fileset, f"{prefix}_work{number}"


def parse_project_quota(output: str, fileset: str, filesystem: str) -> dict[str, float]:
    header = None
    try:
        for line in output.splitlines():
            fields = line.split(":")
            if fields[:2] != ["mmlsquota", "fileset"]:
                continue
            if fields[2] == "HEADER":
                header = fields
                continue
            if header is None:
                continue
            row = dict(zip(header, fields))
            if (row.get("filesystemName") != filesystem
                    or row.get("name") != fileset
                    or row.get("quotaType") != "FILESET"):
                continue
            result = {}
            for label, prefix in (("storage_percent", "block"),
                                  ("inode_percent", "files")):
                values = [int(row[prefix + suffix]) for suffix in
                          ("Usage", "InDoubt", "Quota", "Limit")]
                if any(value < 0 for value in values):
                    raise ValueError("negative quota value")
                usage, in_doubt, quota, limit = values
                limits = [value for value in (quota, limit) if value > 0]
                # GPFS reports zero for an unlimited quota axis.
                result[label] = 100.0 * (usage + in_doubt) / min(limits) if limits else 0.0
            return result
    except (IndexError, KeyError, ValueError) as exc:
        raise QuotaReadError(f"invalid quota output for {fileset}@{filesystem}: {exc}") from exc
    raise QuotaReadError(f"no quota row for {fileset}@{filesystem}")


def read_project_quota(fileset: str, filesystem: str) -> dict[str, float]:
    try:
        result = subprocess.run(
            [MMLSQUOTA, "-e", "-Y", "-j", fileset, filesystem],
            capture_output=True, text=True, timeout=30, check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise QuotaReadError(f"quota query failed for {fileset}@{filesystem}: {exc}") from exc
    if result.returncode:
        raise QuotaReadError(
            f"quota query failed for {fileset}@{filesystem}: "
            f"{(result.stderr or result.stdout).strip()}"
        )
    return parse_project_quota(result.stdout, fileset, filesystem)


def wait_for_storage_quota(
    destinations: Sequence[Path],
    lease: dict[str, Any] | None = None,
    *,
    mode: str = "auto",
    threshold_percent: float = 97.0,
    resume_cores: float = 1.0,
    resume_memory_mb: float = 768.0,
    lease_command: str = DEFAULT_LEASE_COMMAND,
    resume_phase: str = "storage_work",
    wait_phase: str = "storage_quota_wait",
    log_label: str = "storage quota",
    wait_note: str = "before shared-storage writes",
) -> dict[str, Any]:
    """Check before a large write; after blockage, poll and stagger restart.

    Auto mode checks Snellius project paths only. Disabled mode performs no
    path lookup, quota command, sleep or lease adjustment (e.g. on Spider).
    SSD reservations are deliberately untouched. Quota-query errors back off
    like a full quota. This is not an atomic free-space reservation.
    """
    if mode not in {"auto", "required", "disabled"}:
        raise ValueError("quota mode must be auto, required, or disabled")
    if not 0 < threshold_percent < 100:
        raise ValueError("quota threshold must be between 0 and 100 percent")
    report: dict[str, Any] = {"enabled": False, "checks": 0, "wait_seconds": 0.0}
    if mode == "disabled":
        report["reason"] = "explicitly disabled"
        return report
    targets = set()
    for destination in destinations:
        target = project_quota_target(destination)
        if target is None:
            if mode == "required":
                raise ValueError(f"not a Snellius project destination: {destination}")
        else:
            targets.add(target)
    if not targets:
        report["reason"] = "no Snellius project destinations"
        return report
    report.update(enabled=True, threshold_percent=threshold_percent,
                  targets=[list(target) for target in sorted(targets)])
    # Workers may reparse updated rules after submission. Never reacquire more
    # than the actual original reservation advertised by their chief.
    if lease is None:
        resume_cores = min(resume_cores, float(os.environ.get(
            "ZSLURM_LEASE_MAX_CORES", resume_cores)))
        resume_memory_mb = min(resume_memory_mb, float(os.environ.get(
            "ZSLURM_LEASE_MAX_MEM_MB", resume_memory_mb)))
    waiting = False
    released = False

    def check() -> bool:
        report["checks"] += 1
        allowed = True
        for fileset, filesystem in sorted(targets):
            try:
                usage = read_project_quota(fileset, filesystem)
                print(f"[{log_label}] {fileset}@{filesystem}: "
                      f"storage={usage['storage_percent']:.4f}%, "
                      f"inodes={usage['inode_percent']:.4f}%; "
                      f"write threshold={threshold_percent:g}%", file=sys.stderr)
                if max(usage.values()) >= threshold_percent:
                    allowed = False
            except QuotaReadError as exc:
                print(f"[{log_label}] {exc}; waiting", file=sys.stderr)
                allowed = False
        return allowed

    def pause(low: float, high: float, reason: str) -> None:
        delay = random.uniform(low, high)
        print(f"[{log_label}] {reason}; waiting {delay:.1f}s {wait_note}", file=sys.stderr)
        report["wait_seconds"] += delay
        time.sleep(delay)

    allowed = check()
    while True:
        if not allowed:
            waiting = True
            if not released:
                # Discover leases lazily: healthy/disabled checks should not
                # add scheduler traffic. Local execution can wait without one.
                if lease is None:
                    lease = lease_preflight(
                        "required" if all(os.environ.get(key) for key in LEASE_ENV) else "optional",
                        lease_command, initial_cores=resume_cores,
                        initial_memory_mb=resume_memory_mb,
                    )
                shrink_lease(lease, cores=min(0.1, resume_cores),
                             memory_mb=min(1024.0, resume_memory_mb),
                             phase=wait_phase)
                released = True
            pause(240.0, 360.0, "quota unavailable")
            allowed = check()
            continue
        if not waiting:
            return report
        pause(0.0, 600.0, "quota available; staggering write restart")
        allowed = check()
        if not allowed:
            continue
        if lease is not None and lease.get("available"):
            acquire_lease(lease, cores=resume_cores, memory_mb=resume_memory_mb,
                          phase=resume_phase, wait_seconds=60.0, attempts=1)
            if not lease["acquire"]["acquired"]:
                pause(240.0, 360.0, "waiting for work CPU/RAM lease")
                allowed = check()
                continue
            released = False
            # Reacquiring can take time; do not use the pre-acquisition quota.
            allowed = check()
            if not allowed:
                continue
        return report


def wait_for_fastq_quota(
    destinations: Sequence[Path],
    lease: dict[str, Any],
    *,
    mode: str = "auto",
    threshold_percent: float = 97.0,
    copy_cores: float = 1.0,
    copy_memory_mb: float = 768.0,
) -> dict[str, Any]:
    """Keep the original SSD FASTQ-publication interface and phase labels."""
    return wait_for_storage_quota(
        destinations, lease, mode=mode, threshold_percent=threshold_percent,
        resume_cores=copy_cores, resume_memory_mb=copy_memory_mb,
        resume_phase="fastq_copy", wait_phase="fastq_quota_wait",
        log_label="FASTQ quota", wait_note="on SSD",
    )
