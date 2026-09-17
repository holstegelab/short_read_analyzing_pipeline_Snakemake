"""Shared scratch, publication and scheduler-lease helpers for pipeline runners.

This module has no alignment-specific processing. Scheduler behavior is kept
unchanged when a rule uses a different execution backend or input route.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Sequence


LEASE_ENV = (
    "ZSLURM_LEASE_SOCKET",
    "ZSLURM_LEASE_TOKEN",
    "ZSLURM_JOB_ID",
)


class LeaseError(RuntimeError):
    def __init__(
        self,
        message: str,
        *,
        response: dict[str, Any] | None = None,
        returncode: int | None = None,
    ) -> None:
        super().__init__(message)
        self.response = response
        self.returncode = returncode


def assigned_scratch(
    explicit: str | None = None,
    *,
    shared_fallback: str | None = None,
) -> Path:
    """Resolve this job's scratch root; never borrow another job's directory.

    ZSlurm exports ``ZSLURM_SCRATCH_DIR`` for the individual child job.  The
    legacy Snellius path and ``SLURM_TMPDIR`` keep direct Slurm execution
    working where the scheduler provides a dedicated variable. Generic
    ``TMPDIR`` is deliberately not guessed: on one site it can be local SSD,
    while on another it can be the node's small system filesystem. ZSlurm
    normalizes either layout through its explicit child-job variable.

    Required-SSD callers omit ``shared_fallback``.  Rules with
    ``ssd_use=possible`` explicitly supply a workflow-local shared directory.
    """
    candidates: list[Path] = []
    if explicit:
        candidates.append(Path(explicit))

    assigned = os.environ.get("ZSLURM_SCRATCH_DIR")
    if assigned:
        candidates.append(Path(assigned))

    user = os.environ.get("USER")
    job_id = os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_JOBID")
    if user and job_id:
        candidates.append(Path(f"/scratch-node/{user}.{job_id}"))

    slurm_tmp = os.environ.get("SLURM_TMPDIR")
    if slurm_tmp:
        candidates.append(Path(slurm_tmp))

    seen: set[Path] = set()
    for candidate in candidates:
        try:
            resolved = candidate.expanduser().resolve()
        except OSError:
            continue
        if resolved in seen:
            continue
        seen.add(resolved)
        if resolved.is_dir() and os.access(resolved, os.W_OK):
            return resolved

    if shared_fallback:
        fallback = Path(shared_fallback)
        fallback.mkdir(parents=True, exist_ok=True)
        if fallback.is_dir() and os.access(fallback, os.W_OK):
            return fallback.resolve()

    raise RuntimeError(
        "ssd_use=required but this job has no writable assigned scratch "
        "directory (ZSLURM_SCRATCH_DIR, /scratch-node job path, "
        "or SLURM_TMPDIR)"
    )


def executable(value: str) -> str:
    path = value if os.sep in value else shutil.which(value)
    if not path or not Path(path).is_file() or not os.access(path, os.X_OK):
        raise FileNotFoundError(f"executable not found: {value}")
    return str(Path(path).resolve())


def _atomic_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "wt", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def _atomic_copy(source: Path, destination: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{destination.name}.", dir=destination.parent)
    os.close(fd)
    try:
        shutil.copyfile(source, temporary)
        os.replace(temporary, destination)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def _atomic_publish(source: Path, destination: Path) -> None:
    """Publish a final output atomically, moving it when both paths share a FS."""
    if not source.is_file():
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_raw = tempfile.mkstemp(
        prefix=f".{destination.name}.", dir=destination.parent
    )
    os.close(fd)
    temporary = Path(temporary_raw)
    try:
        if source.stat().st_dev == destination.parent.stat().st_dev:
            os.replace(source, temporary)
        else:
            shutil.copyfile(source, temporary)
        os.replace(temporary, destination)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def _lease_request(command: str, arguments: Sequence[str]) -> dict[str, Any]:
    process = subprocess.run(
        [command, "--json", *arguments],
        check=False,
        capture_output=True,
        text=True,
    )
    response: dict[str, Any] | None = None
    if process.stdout.strip():
        try:
            parsed = json.loads(process.stdout)
        except json.JSONDecodeError:
            parsed = None
        if isinstance(parsed, dict):
            response = parsed
    if process.returncode:
        detail = (process.stderr or process.stdout).strip()
        raise LeaseError(
            f"{' '.join(arguments)} failed with exit {process.returncode}: {detail}",
            response=response,
            returncode=process.returncode,
        )
    if response is None:
        raise LeaseError(f"lease command returned invalid JSON: {process.stdout!r}")
    if not response.get("ok"):
        raise LeaseError(
            f"lease request was rejected: {response.get('status')}: "
            f"{response.get('message')}",
            response=response,
        )
    return response


def lease_preflight(
    mode: str,
    lease_command: str,
    *,
    initial_cores: float,
    initial_memory_mb: float,
) -> dict[str, Any]:
    """Verify dynamic leases before expensive work, or explicitly disable them."""
    if mode == "disabled":
        return {"mode": mode, "available": False, "reason": "explicitly disabled"}
    missing = [name for name in LEASE_ENV if not os.environ.get(name)]
    command = shutil.which(lease_command) if os.sep not in lease_command else lease_command
    if missing or not command or not Path(command).is_file():
        reason = (
            f"missing environment {missing}" if missing else
            f"lease command not found: {lease_command}"
        )
        if mode == "required":
            raise LeaseError(reason)
        print(f"[fused-alignment] lease unavailable; retaining maximum: {reason}", file=sys.stderr)
        return {"mode": mode, "available": False, "reason": reason}
    try:
        status = _lease_request(str(command), ["status"])
        if float(status["held_cores"]) + 1e-6 < float(initial_cores):
            raise LeaseError(
                f"initial lease holds {status['held_cores']} cores, expected "
                f"at least {initial_cores}"
            )
        if float(status["held_mem_mb"]) + 1e-6 < float(initial_memory_mb):
            raise LeaseError(
                f"initial lease holds {status['held_mem_mb']} MB, expected "
                f"at least {initial_memory_mb} MB"
            )
    except Exception as exc:
        if mode == "required":
            raise
        print(
            f"[fused-alignment] lease preflight failed; retaining maximum: {exc}",
            file=sys.stderr,
        )
        return {"mode": mode, "available": False, "reason": str(exc)}
    return {
        "mode": mode,
        "available": True,
        "command": str(command),
        "status_before": status,
    }


def shrink_lease(
    lease: dict[str, Any],
    *,
    cores: float,
    memory_mb: float,
    phase: str | None = None,
) -> dict[str, Any]:
    if not lease.get("available"):
        lease["shrink"] = {"performed": False, "reason": "lease unavailable"}
        return lease
    last_error: Exception | None = None
    last_response: dict[str, Any] | None = None
    # zslurm_chief samples live PSS every five seconds by default. Cover at
    # least two fresh samples so a just-exited high-memory phase cannot leave
    # the job at its maximum reservation for the entire low-memory tail.
    attempts = 11
    for attempt in range(1, attempts + 1):
        arguments = [
            "set",
            "--cores",
            str(cores),
            "--mem-mb",
            str(memory_mb),
            "--wait",
            "0",
        ]
        # One logical phase starts on the first request. Safety-floor retries
        # must not fragment reporting into a series of duplicate phases.
        if phase is not None and attempt == 1:
            arguments.extend(["--phase", phase])
        try:
            response = _lease_request(str(lease["command"]), arguments)
            last_response = response
            target_reached = (
                abs(float(response["held_cores"]) - float(cores)) <= 1e-6
                and abs(float(response["held_mem_mb"]) - float(memory_mb)) <= 1e-6
            )
            if target_reached:
                adjustment = {
                    "performed": True,
                    "target_reached": True,
                    "attempts": attempt,
                    "response": response,
                }
                lease["shrink"] = adjustment
                lease.setdefault("adjustments", []).append(adjustment)
                return lease
            last_error = LeaseError(
                "lease safety floor retained "
                f"{response['held_cores']} cores and {response['held_mem_mb']} MB"
            )
        except Exception as exc:
            last_error = exc
        if attempt < attempts:
            # A completed high-memory subprocess can remain briefly visible in
            # the chief's process-tree sample. Retry the idempotent absolute
            # target after that observation has had time to settle.
            time.sleep(1.0)
    if last_response is not None:
        adjustment = {
            "performed": True,
            "target_reached": False,
            "attempts": attempts,
            "response": last_response,
            "reason": str(last_error),
        }
        lease["shrink"] = adjustment
        lease.setdefault("adjustments", []).append(adjustment)
        # Retaining more than requested after a shrink is safe (and expected
        # while the live-usage safety floor catches up). Failing to reacquire
        # requested capacity is not safe for a required lease.
        if lease["mode"] == "required" and (
            float(last_response["held_cores"]) + 1e-6 < float(cores)
            or float(last_response["held_mem_mb"]) + 1e-6 < float(memory_mb)
        ):
            raise LeaseError(f"could not grow required lease: {last_error}")
        return lease
    try:
        status = _lease_request(str(lease["command"]), ["status"])
    except Exception:
        status = None
    if status is not None and (
        abs(float(status["held_cores"]) - float(cores)) <= 1e-6
        and abs(float(status["held_mem_mb"]) - float(memory_mb)) <= 1e-6
    ):
        adjustment = {
            "performed": True,
            "target_reached": True,
            "attempts": attempts,
            "response": status,
            "verified_after_lost_reply": True,
        }
        lease["shrink"] = adjustment
        lease.setdefault("adjustments", []).append(adjustment)
        return lease
    if lease["mode"] == "required":
        raise LeaseError(f"could not set required lease: {last_error}")
    print(
        f"[fused-alignment] lease adjustment failed: {last_error}",
        file=sys.stderr,
    )
    lease["shrink"] = {"performed": False, "reason": str(last_error)}
    return lease


def _lease_holds_at_least(
    response: dict[str, Any], *, cores: float, memory_mb: float
) -> bool:
    """Return whether a response confirms all capacity required by a phase."""

    try:
        return (
            float(response["held_cores"]) + 1e-6 >= float(cores)
            and float(response["held_mem_mb"]) + 1e-6 >= float(memory_mb)
        )
    except (KeyError, TypeError, ValueError):
        return False


def acquire_lease(
    lease: dict[str, Any],
    *,
    cores: float,
    memory_mb: float,
    phase: str | None = None,
    wait_seconds: float = 3600.0,
    attempts: int = 3,
    retry_delay_seconds: float = 1.0,
) -> dict[str, Any]:
    """Safely acquire a larger phase lease without making capacity fatal.

    ZSlurm keeps the previous lease while a FIFO growth request waits. A
    capacity timeout therefore means the caller must use its lower-resource
    fallback; it is not a reason to fail otherwise-valid work. Transient
    transport/manager errors are retried after checking ``status`` because an
    absolute request may have committed even when its reply was lost.

    Callers must inspect ``lease["acquire"]["acquired"]`` before starting the
    higher-resource phase.
    """

    if attempts < 1:
        raise ValueError("lease acquire attempts must be at least one")
    if wait_seconds < 0:
        raise ValueError("lease acquire wait_seconds cannot be negative")
    if retry_delay_seconds < 0:
        raise ValueError("lease acquire retry_delay_seconds cannot be negative")

    if not lease.get("available"):
        adjustment = {
            "kind": "acquire",
            "performed": False,
            "target_reached": False,
            "acquired": False,
            "attempts": 0,
            "requested_cores": float(cores),
            "requested_mem_mb": float(memory_mb),
            "wait_seconds": float(wait_seconds),
            "reason": "lease unavailable",
        }
        lease["acquire"] = adjustment
        lease.setdefault("adjustments", []).append(adjustment)
        return lease

    command = str(lease["command"])
    last_error: Exception | None = None
    last_response: dict[str, Any] | None = None
    attempts_used = 0
    timed_out = False

    def record_success(
        response: dict[str, Any], *, verified_after_lost_reply: bool = False
    ) -> dict[str, Any]:
        adjustment = {
            "kind": "acquire",
            "performed": True,
            "target_reached": True,
            "acquired": True,
            "attempts": attempts_used,
            "requested_cores": float(cores),
            "requested_mem_mb": float(memory_mb),
            "wait_seconds": float(wait_seconds),
            "response": response,
        }
        if verified_after_lost_reply:
            adjustment["verified_after_lost_reply"] = True
        lease["acquire"] = adjustment
        lease.setdefault("adjustments", []).append(adjustment)
        return lease

    for attempt in range(1, attempts + 1):
        attempts_used = attempt
        arguments = [
            "set",
            "--cores",
            str(cores),
            "--mem-mb",
            str(memory_mb),
            "--wait",
            str(wait_seconds),
        ]
        if phase is not None:
            arguments.extend(["--phase", phase])

        error_response: dict[str, Any] | None = None
        try:
            response = _lease_request(command, arguments)
            last_response = response
            if _lease_holds_at_least(response, cores=cores, memory_mb=memory_mb):
                return record_success(response)
            last_error = LeaseError(
                "lease response did not confirm requested capacity",
                response=response,
            )
            error_response = response
        except Exception as exc:
            last_error = exc
            if isinstance(exc, LeaseError):
                error_response = exc.response
                if error_response is not None:
                    last_response = error_response

        # A set request is absolute and may have committed even if the response
        # was lost. Verify before retrying and before selecting a fallback.
        try:
            status = _lease_request(command, ["status"])
        except Exception:
            status = None
        if status is not None:
            last_response = status
            if _lease_holds_at_least(status, cores=cores, memory_mb=memory_mb):
                return record_success(status, verified_after_lost_reply=True)

        response_status = (
            str(error_response.get("status", "")) if error_response else ""
        )
        response_code = error_response.get("code") if error_response else None
        if response_status == "timeout" or response_code == 4:
            timed_out = True
            break
        # Invalid/denied requests cannot improve through an unchanged retry.
        if response_code in {2, 3}:
            break
        if attempt < attempts:
            time.sleep(retry_delay_seconds)

    adjustment = {
        "kind": "acquire",
        "performed": True,
        "target_reached": False,
        "acquired": False,
        "attempts": attempts_used,
        "requested_cores": float(cores),
        "requested_mem_mb": float(memory_mb),
        "wait_seconds": float(wait_seconds),
        "timed_out": timed_out,
        "reason": str(last_error),
    }
    if last_response is not None:
        adjustment["response"] = last_response
    lease["acquire"] = adjustment
    lease.setdefault("adjustments", []).append(adjustment)
    print(
        "[fused-alignment] lease acquire unavailable; caller must use its "
        f"lower-resource fallback: {last_error}",
        file=sys.stderr,
    )
    return lease
