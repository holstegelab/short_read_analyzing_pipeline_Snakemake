import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
sys.path.insert(0, str(SCRIPTS))

from io_profile import path_bytes, run_profiled  # noqa: E402
import io_profile  # noqa: E402
from pipeline_runtime import assigned_scratch  # noqa: E402


def test_path_bytes_counts_nested_files(tmp_path):
    nested = tmp_path / "nested"
    nested.mkdir()
    (tmp_path / "one").write_bytes(b"123")
    (nested / "two").write_bytes(b"4567")
    assert path_bytes(tmp_path) == 7
    assert path_bytes(tmp_path / "absent") is None


def test_path_bytes_tolerates_runfiles_directory_replaced_during_walk(tmp_path, monkeypatch):
    runfiles = tmp_path / "Bazel.runfiles" / "runfiles"
    runfiles.mkdir(parents=True)
    (runfiles / "payload").write_bytes(b"temporary")
    (tmp_path / "stable").write_bytes(b"123")
    walk = io_profile.os.walk

    def changing_walk(path):
        for root, directories, filenames in walk(path):
            if Path(root) == runfiles:
                (runfiles / "payload").unlink()
                runfiles.rmdir()
                runfiles.write_bytes(b"replaced")
            yield root, directories, filenames

    monkeypatch.setattr(io_profile.os, "walk", changing_walk)
    assert path_bytes(tmp_path) == 3
    assert path_bytes(runfiles / "payload") is None


@pytest.mark.parametrize("sampler,operation", [
    ("_path_sizes", "local_paths"),
    ("_process_sample", "process_tree"),
    ("_filesystem_sample", "filesystem"),
])
def test_sampling_failure_does_not_abort_command(tmp_path, monkeypatch, sampler, operation):
    def failed_sample(*args):
        raise OSError("temporary filesystem measurement failure")

    monkeypatch.setattr(io_profile, sampler, failed_sample)
    output = tmp_path / "completed"
    metrics = tmp_path / "metrics.json"
    run_profiled(
        [sys.executable, "-c", "import time; from pathlib import Path; "
         f"time.sleep(0.15); Path({str(output)!r}).write_text('complete')"],
        metrics_path=metrics, label="measurement-failure", local_paths=[tmp_path],
        poll_interval=0.05,
    )
    payload = json.loads(metrics.read_text())
    assert output.read_text() == "complete"
    assert payload["return_code"] == 0
    assert payload["launch_error"] is None
    assert any(e["operation"] == operation for e in payload["sampling_errors"])


@pytest.mark.parametrize("exit_code", [0, 7])
def test_metrics_quota_failure_preserves_command_exit_and_cleanup(tmp_path, monkeypatch, exit_code, capsys):
    local = tmp_path / "scratch" / "job"
    local.mkdir(parents=True)

    def quota_failure(*args):
        raise OSError(122, "Disk quota exceeded")

    monkeypatch.setattr(io_profile, "_write_json_atomic", quota_failure)
    arguments = dict(
        metrics_path=tmp_path / "metrics.json", label="quota-test",
        local_paths=[local], cleanup_paths=[local], poll_interval=0.05,
    )
    command = [sys.executable, "-c", f"raise SystemExit({exit_code})"]
    if exit_code:
        with pytest.raises(subprocess.CalledProcessError) as error:
            run_profiled(command, **arguments)
        assert error.value.returncode == exit_code
    else:
        assert run_profiled(command, **arguments) == 0
    assert not local.exists()
    assert "Disk quota exceeded" in capsys.readouterr().err


def test_failed_orchestrator_stops_surviving_shards_before_cleanup(tmp_path):
    local = tmp_path / "scratch" / "job"
    local.mkdir(parents=True)
    ready = tmp_path / "shard.pid"
    child = (
        "import os, signal, time; from pathlib import Path; "
        "signal.signal(signal.SIGTERM, signal.SIG_IGN); "
        f"Path({str(ready)!r}).write_text(str(os.getpid())); "
        "time.sleep(0.3); "
        f"Path({str(local / 'late-write')!r}).parent.mkdir(parents=True, exist_ok=True); "
        f"Path({str(local / 'late-write')!r}).write_text('still running'); time.sleep(10)"
    )
    parent = (
        "import subprocess, sys, time; from pathlib import Path; "
        f"subprocess.Popen([sys.executable, '-c', {child!r}]); "
        f"ready = Path({str(ready)!r}); "
        "\nwhile not ready.exists(): time.sleep(0.01)\nraise SystemExit(7)"
    )
    try:
        with pytest.raises(subprocess.CalledProcessError) as error:
            run_profiled(
                [sys.executable, "-c", parent],
                metrics_path=tmp_path / "failed.json", label="failed-shards",
                local_paths=[local], cleanup_paths=[local], poll_interval=0.05,
            )
        assert error.value.returncode == 7
        time.sleep(0.4)
        assert not local.exists()
    finally:
        if ready.exists():
            try:
                os.kill(int(ready.read_text()), signal.SIGKILL)
            except ProcessLookupError:
                pass


def test_profile_records_failure_and_cleans_local_path(tmp_path):
    local = tmp_path / "scratch" / "job"
    local.mkdir(parents=True)
    metrics = tmp_path / "failed.json"
    command = [
        sys.executable,
        "-c",
        "from pathlib import Path; Path(r'%s').write_bytes(b'x' * 4096); raise SystemExit(7)"
        % (local / "spill"),
    ]
    with pytest.raises(subprocess.CalledProcessError) as error:
        run_profiled(
            command,
            metrics_path=metrics,
            label="failure-test",
            local_paths=[local],
            cleanup_paths=[local],
            poll_interval=0.05,
        )
    assert error.value.returncode == 7
    payload = json.loads(metrics.read_text())
    assert payload["return_code"] == 7
    assert payload["cleanup_ok"] is True
    assert not local.exists()


def test_assigned_scratch_requires_explicit_or_slurm_path(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    monkeypatch.delenv("SLURM_JOBID", raising=False)
    monkeypatch.delenv("ZSLURM_SCRATCH_DIR", raising=False)
    monkeypatch.delenv("SLURM_TMPDIR", raising=False)
    monkeypatch.delenv("TMPDIR", raising=False)
    with pytest.raises(RuntimeError, match="no writable assigned"):
        assigned_scratch()
    assert assigned_scratch(str(tmp_path)) == tmp_path.resolve()


def test_profile_records_success_and_cleans_local_path(tmp_path):
    local = tmp_path / "scratch" / "job"
    local.mkdir(parents=True)
    output = tmp_path / "published"
    metrics = tmp_path / "success.io.json"
    run_profiled(
        [sys.executable, "-c", "from pathlib import Path; import time; "
         f"Path({str(local / 'spill')!r}).write_bytes(b's' * 8192); "
         f"Path({str(output)!r}).write_bytes(b'complete'); time.sleep(0.15)"],
        metrics_path=metrics, label="success-test", local_paths=[local],
        output_paths=[output], cleanup_paths=[local], requested_ssd_gb=6,
        poll_interval=0.05,
    )
    payload = json.loads(metrics.read_text())
    assert output.read_bytes() == b"complete"
    assert payload["return_code"] == 0
    assert payload["requested"]["ssd_gb"] == 6.0
    assert payload["peaks"]["local_bytes"] >= 8192
    assert payload["cleanup_ok"] is True
    assert not local.exists()
