import json
import subprocess
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
sys.path.insert(0, str(SCRIPTS))

from io_profile import path_bytes, run_profiled  # noqa: E402
from run_samtools_sort_ssd import assigned_scratch  # noqa: E402


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
    monkeypatch.delenv("SLURM_TMPDIR", raising=False)
    with pytest.raises(RuntimeError, match="no writable assigned"):
        assigned_scratch()
    assert assigned_scratch(str(tmp_path)) == tmp_path.resolve()


def test_sort_runner_profiles_and_removes_job_scratch(tmp_path):
    fake_samtools = tmp_path / "samtools"
    fake_samtools.write_text(
        """#!/usr/bin/env python3
import sys
import time
from pathlib import Path

args = sys.argv[1:]
assert args[0] == "sort"
temp_prefix = Path(args[args.index("-T") + 1])
temp_prefix.parent.mkdir(parents=True, exist_ok=True)
(temp_prefix.parent / "spill.tmp").write_bytes(b"s" * 8192)
output, index = args[args.index("-o") + 1].split("##idx##")
Path(output).write_bytes(b"bam")
Path(index).write_bytes(b"bai")
time.sleep(0.15)
"""
    )
    fake_samtools.chmod(0o755)
    input_bam = tmp_path / "input.bam"
    input_bam.write_bytes(b"unsorted")
    output_bam = tmp_path / "output.bam"
    output_bai = tmp_path / "output.bam.bai"
    metrics = tmp_path / "sort.io.json"
    scratch = tmp_path / "node-ssd"
    scratch.mkdir()

    subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "run_samtools_sort_ssd.py"),
            "--input",
            str(input_bam),
            "--output-bam",
            str(output_bam),
            "--output-bai",
            str(output_bai),
            "--metrics",
            str(metrics),
            "--ssd-gb",
            "6",
            "--scratch-base",
            str(scratch),
            "--samtools",
            str(fake_samtools),
            "--poll-interval",
            "0.05",
        ],
        check=True,
    )

    payload = json.loads(metrics.read_text())
    assert output_bam.read_bytes() == b"bam"
    assert output_bai.read_bytes() == b"bai"
    assert payload["return_code"] == 0
    assert payload["requested"]["ssd_gb"] == 6.0
    assert payload["peaks"]["local_bytes"] >= 8192
    assert payload["cleanup_ok"] is True
    assert list((scratch / "aligner_sort").iterdir()) == []
