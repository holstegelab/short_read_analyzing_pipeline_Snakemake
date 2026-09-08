import os
import runpy
import subprocess
import tarfile
import threading
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/extract_and_tar_deepvariant_level2.py"


def setup_run(tmp_path, monkeypatch, workers=4, count=8, dataset="wgs"):
    for key in list(os.environ):
        if key.startswith("ZSLURM_") or key.startswith("SLURM_"):
            monkeypatch.delenv(key)
    root = tmp_path / "out"
    root.mkdir()
    bed = tmp_path / "region.bed"
    bed.write_text("chr1\t0\t100\n")
    samples = [f"sample_{i}" for i in range(count)]
    inputs = []
    for sample in samples:
        vcf = tmp_path / f"{sample}.vcf.gz"
        index = Path(str(vcf) + ".tbi")
        vcf.write_bytes(("input:" + sample).encode())
        index.write_bytes(b"input-index")
        inputs.extend([str(vcf), str(index)])
    smk = SimpleNamespace(
        threads=workers, resources=SimpleNamespace(n=workers),
        params=SimpleNamespace(samples=samples, region="A1", samplefile="cohort",
                               dataset=dataset, interval=str(bed)),
        input=SimpleNamespace(gvcfs=inputs),
        output=SimpleNamespace(tar=str(root / "cohort.tar")),
    )
    state = {"active": 0, "maximum": 0, "calls": [], "fail_sample": None, "no_index": False}
    lock = threading.Lock()

    def fake_run(args, **kwargs):
        with lock:
            state["calls"].append(args)
            state["active"] += 1
            state["maximum"] = max(state["maximum"], state["active"])
        try:
            if args[1] == "view":
                source = Path(args[4])
                sample = source.name.split(".")[0]
                if sample == state["fail_sample"]:
                    raise subprocess.CalledProcessError(7, args)
                # Deliberately finish samples out of order.
                time.sleep(0.003 * (8 - int(sample.split("_")[1]) % 8))
                Path(args[-1]).write_bytes(b"extracted:" + source.read_bytes())
            elif args[1] == "index" and not state["no_index"]:
                Path(args[-1] + ".tbi").write_bytes(b"generated-index")
            return subprocess.CompletedProcess(args, 0)
        finally:
            with lock:
                state["active"] -= 1

    monkeypatch.setattr("shutil.which", lambda name: "/fake/bcftools")
    monkeypatch.setattr(subprocess, "run", fake_run)
    return smk, state


def contents(path):
    with tarfile.open(path, "r:") as archive:
        names = archive.getnames()
        return names, [archive.extractfile(name).read() for name in names]


@pytest.mark.parametrize("dataset", ["wgs", "wes"])
def test_parallel_payloads_and_member_order_match_serial(tmp_path, monkeypatch, dataset):
    smk, state = setup_run(tmp_path, monkeypatch, workers=4, dataset=dataset)
    runpy.run_path(str(SCRIPT), init_globals={"snakemake": smk})
    parallel = contents(smk.output.tar)
    assert 2 <= state["maximum"] <= 4
    assert state["active"] == 0
    assert len(state["calls"]) == 16
    assert all("--threads" not in command for command in state["calls"])
    smk.threads = 1
    state["maximum"] = 0
    state["calls"].clear()
    runpy.run_path(str(SCRIPT), init_globals={"snakemake": smk})
    assert contents(smk.output.tar) == parallel
    assert state["maximum"] == 1
    assert len(state["calls"]) == 16  # Existing tar does NOT trigger reuse.
    assert all(Path(path).is_file() for path in smk.input.gvcfs)


@pytest.mark.parametrize("threads,reserved,granted,expected", [
    (8, 8, "1.5", 1),  # Previously queued ZSlurm job.
    (8, "1.5", None, 1),
    (8, 8, "4", 4),
    (2, 8, "8", 2),
    (8, 8, "invalid", 1),
    (8, 8, "nan", 1),
    (8, 8, "0.9", 1),
])
def test_pool_respects_reservation_and_legacy_leases(tmp_path, monkeypatch, threads, reserved, granted, expected):
    smk, state = setup_run(tmp_path, monkeypatch, workers=threads)
    smk.resources.n = reserved
    if granted is not None:
        monkeypatch.setenv("ZSLURM_LEASE_MAX_CORES", granted)
    runpy.run_path(str(SCRIPT), init_globals={"snakemake": smk})
    assert state["maximum"] == expected


def test_failure_waits_for_workers_and_never_publishes_partial_archive(tmp_path, monkeypatch):
    smk, state = setup_run(tmp_path, monkeypatch, workers=2, count=20)
    original = b"previous-result-is-not-a-reuse-shortcut"
    Path(smk.output.tar).write_bytes(original)
    state["fail_sample"] = "sample_0"
    with pytest.raises(subprocess.CalledProcessError):
        runpy.run_path(str(SCRIPT), init_globals={"snakemake": smk})
    assert state["active"] == 0
    assert sum(args[1] == "view" for args in state["calls"]) <= 2
    assert Path(smk.output.tar).read_bytes() == original
    assert list(Path(smk.output.tar).parent.iterdir()) == [Path(smk.output.tar)]


def test_missing_generated_index_fails_without_publish(tmp_path, monkeypatch):
    smk, state = setup_run(tmp_path, monkeypatch)
    state["no_index"] = True
    with pytest.raises(FileNotFoundError, match="Extraction failed"):
        runpy.run_path(str(SCRIPT), init_globals={"snakemake": smk})
    assert state["active"] == 0
    assert not Path(smk.output.tar).exists()
    assert not list(Path(smk.output.tar).parent.iterdir())


def test_duplicate_samples_rejected_before_any_work(tmp_path, monkeypatch):
    smk, state = setup_run(tmp_path, monkeypatch)
    smk.params.samples[1] = smk.params.samples[0]
    with pytest.raises(ValueError, match="Duplicate sample"):
        runpy.run_path(str(SCRIPT), init_globals={"snakemake": smk})
    assert not state["calls"]


@pytest.mark.parametrize("failure", [None, "rejected", "timeout"])
def test_pack_phase_only_releases_cpu_after_workers_stop(tmp_path, monkeypatch, failure):
    smk, state = setup_run(tmp_path, monkeypatch)
    smk.params.lease_command = "/fake/lease_client.py"
    for key in ("ZSLURM_LEASE_SOCKET", "ZSLURM_LEASE_TOKEN", "ZSLURM_JOB_ID"):
        monkeypatch.setenv(key, "test-only")
    real_fake_run = subprocess.run
    lease_calls = []

    def lease_run(args, **kwargs):
        if "/fake/lease_client.py" not in args:
            return real_fake_run(args, **kwargs)
        assert state["active"] == 0
        lease_calls.append(args)
        assert "--mem-mb" not in args
        if failure == "timeout":
            raise subprocess.TimeoutExpired(args, 30)
        return subprocess.CompletedProcess(args, int(failure == "rejected"))

    monkeypatch.setattr(subprocess, "run", lease_run)
    runpy.run_path(str(SCRIPT), init_globals={"snakemake": smk})
    assert len(lease_calls) == 1
    assert lease_calls[0][-4:] == ["--cores", "1", "--phase", "level2_tar_publish"]
    assert Path(smk.output.tar).is_file()


def test_empty_batch_still_fails_without_publishing(tmp_path, monkeypatch):
    smk, state = setup_run(tmp_path, monkeypatch, count=0)
    with pytest.raises(ValueError, match="No staged gVCFs"):
        runpy.run_path(str(SCRIPT), init_globals={"snakemake": smk})
    assert not Path(smk.output.tar).exists()


def test_both_rules_reserve_workers_without_changing_memory():
    source = (SCRIPT.parents[1] / "Deepvariant.smk").read_text()
    assert 'config.get("deepvariant_level2_workers", 8)' in source
    for dataset in ("wgs", "wes"):
        block = source.split(f"rule extract_and_tar_deepvariant_level2_{dataset}:", 1)[1].split("\nrule ", 1)[0]
        assert "threads: DEEPVARIANT_LEVEL2_WORKERS" in block
        assert "n=lambda wildcards, threads: str(threads)" in block
        assert "mem_mb=4000" in block
        assert 'lease_command=srcdir("scripts/zslurm_lease_client.py")' in block
