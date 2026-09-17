import json
import subprocess
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
sys.path.insert(0, str(SCRIPTS))

from io_profile import path_bytes, run_profiled  # noqa: E402
from pipeline_runtime import assigned_scratch  # noqa: E402


def test_path_bytes_counts_nested_files(tmp_path):
    nested = tmp_path / "nested"
    nested.mkdir()
    (tmp_path / "one").write_bytes(b"123")
    (nested / "two").write_bytes(b"4567")
    assert path_bytes(tmp_path) == 7
    assert path_bytes(tmp_path / "absent") is None


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
