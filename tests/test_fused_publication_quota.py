"""Exercise quota backoff in the actual runners without large sequencing jobs."""
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import fastq_quota as quota
import run_fused_deepvariant_phasing as dv
import run_fused_merge_markdup as markdup


def install_quota(monkeypatch, runner, scratch, outputs, *, shared, readings):
    readings = iter(readings)
    events = []
    completed_phases = []

    def target(path):
        if shared or not path.is_relative_to(scratch):
            return "test-project", "wstor_work3"
        return None

    def read(*_args):
        value = next(readings)
        events.append(("check", value))
        return {"storage_percent": value, "inode_percent": 10}

    def pause(seconds):
        events.append(("sleep", seconds))
        # No publication, and no loss of local results, during either wait.
        assert not outputs.exists()
        if completed_phases:
            assert list(scratch.rglob("*.bam")) or list(scratch.rglob("*.vcf.gz"))
        else:
            # Shared-scratch check must precede intermediate/cache writes.
            assert list(scratch.iterdir()) == []

    def shrink(lease, **kwargs):
        events.append(("shrink", kwargs))
        return lease

    def acquire(lease, **kwargs):
        events.append(("acquire", kwargs))
        lease["acquire"] = {"acquired": True}
        return lease

    def profile(command, **kwargs):
        if shared:
            assert any(e[0] == "acquire" for e in events)
        completed_phases.append(kwargs["label"])
        for path in kwargs["output_paths"]:
            Path(path).write_bytes(b"result\n")
        Path(kwargs["metrics_path"]).write_text(json.dumps({"label": kwargs["label"]}))

    monkeypatch.setattr(quota, "project_quota_target", target)
    monkeypatch.setattr(quota, "read_project_quota", read)
    monkeypatch.setattr(quota.time, "sleep", pause)
    monkeypatch.setattr(quota.random, "uniform", lambda a, b: (a + b) / 2)
    monkeypatch.setattr(quota, "shrink_lease", shrink)
    monkeypatch.setattr(quota, "acquire_lease", acquire)
    monkeypatch.setattr(runner, "lease_preflight", lambda *_a, **_k: {"available": True})
    monkeypatch.setattr(runner, "shrink_lease", shrink)
    monkeypatch.setattr(runner, "run_profiled", profile)
    monkeypatch.setattr(runner, "executable", lambda _value: "/usr/bin/true")
    return events, completed_phases


@pytest.mark.parametrize("no_dedup", [0, 1])
def test_markdup_retains_completed_bam_until_quota_recovers(tmp_path, monkeypatch, no_dedup):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    outputs = tmp_path / "outputs"
    inputs = []
    for name in ("input.bam", "input.bam.bai", "checked"):
        path = tmp_path / name
        path.write_bytes(b"input\n")
        inputs.append(str(path))
    argv = [str(markdup.__file__), "--input-bam", inputs[0], "--input-bai", inputs[1],
            "--check-marker", inputs[2], "--sample", "S1", "--initial-cores", "0.95",
            "--initial-memory-mb", "30000", "--markdup-memory-mb", "30000",
            "--no-dedup", str(no_dedup), "--ssd-gb", "12", "--scratch-base", str(scratch)]
    for flag, name in (("output-bam", "markdup.bam"), ("output-bai", "markdup.bam.bai"),
                       ("output-stat", "markdup.stat"), ("merge-log", "merge.log"),
                       ("markdup-log", "markdup.log"), ("metrics", "metrics.json")):
        argv.extend(["--" + flag, str(outputs / name)])
    monkeypatch.setattr(sys, "argv", argv)
    events, phases = install_quota(monkeypatch, markdup, scratch, outputs, shared=False,
                                   readings=[98, 96, 96, 96])
    assert markdup.main() == 0
    assert phases == ["merge_markdup_fused.markdup"]
    assert (outputs / "markdup.bam").read_bytes() == b"result\n"
    assert [e[1]["phase"] for e in events if e[0] == "shrink"] == ["markdup_quota_wait"]
    resumed = [e[1] for e in events if e[0] == "acquire"]
    assert [(e["cores"], e["memory_mb"]) for e in resumed] == [(0.95, 1024)]
    metrics = json.loads((outputs / "metrics.json").read_text())
    assert metrics["success"] and metrics["scratch_removed"]
    assert metrics["quota_wait"]["checks"] == 4
    assert metrics["quota_wait"]["wait_seconds"] == 600


@pytest.mark.parametrize("shared", [False, True])
def test_deepvariant_gates_shared_scratch_and_retains_calls_before_publication(
    tmp_path, monkeypatch, shared
):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    outputs = tmp_path / "outputs"
    input_file = tmp_path / "input"
    input_file.write_bytes(b"input\n")
    argv = [str(dv.__file__), "--sample", "S1", "--region", "A1", "--model-type", "WGS",
            "--haploid-contigs", "chrNONE", "--ploidy", "2", "--skip-sex", "0",
            "--initial-cores", "7.5", "--initial-memory-mb", "8700", "--ssd-gb", "16",
            "--scratch-base", str(scratch), "--quota-mode", "required"]
    for flag in ("bed", "bam", "bai", "validated-sex", "reference", "deepvariant-reference",
                 "interval-bed", "capture-auto-bed", "capture-x-bed", "capture-y-bed",
                 "merge-script", "stats-parser"):
        argv.extend(["--" + flag, str(input_file)])
    # deepvariant_executable() requires a prepared native runtime:
    # PREFIX/bin/run_deepvariant (executable) plus PREFIX/.deepvariant-native.ready.
    native = tmp_path / "native"
    (native / "bin").mkdir(parents=True)
    runner = native / "bin" / "run_deepvariant"
    runner.write_text("#!/bin/sh\nexit 0\n")
    runner.chmod(0o755)
    (native / ".deepvariant-native.ready").touch()
    argv.extend(["--deepvariant-runner", str(runner)])
    for name in ("vcf", "vcf-tbi", "wstats", "merge-stats", "bcftools-stats", "bcftools-summary",
                 "gvcf", "gvcf-tbi", "exome-gvcf", "exome-gvcf-tbi"):
        argv.extend(["--output-" + name, str(outputs / name)])
    argv.extend(["--metrics", str(outputs / "metrics.json")])
    monkeypatch.setattr(sys, "argv", argv)
    events, phases = install_quota(monkeypatch, dv, scratch, outputs, shared=shared,
                                   readings=[98, 96, 96, 96] * (2 if shared else 1))
    assert dv.main() == 0
    assert len(phases) == 2
    assert (outputs / "gvcf").read_bytes() == b"result\n"
    resumed = [e[1] for e in events if e[0] == "acquire"]
    expected = [(7.5, 8700), (1, 1024)] if shared else [(1, 1024)]
    assert [(e["cores"], e["memory_mb"]) for e in resumed] == expected
    metrics = json.loads((outputs / "metrics.json").read_text())
    assert metrics["success"] and metrics["scratch_removed"]
    assert metrics["scratch_quota_wait"]["enabled"] == shared
    assert metrics["quota_wait"]["checks"] == 4
