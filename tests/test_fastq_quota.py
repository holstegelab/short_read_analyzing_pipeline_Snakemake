import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import fastq_quota as quota


HEADER = (
    "mmlsquota:fileset:HEADER:version:reserved:reserved:filesystemName:"
    "quotaType:id:name:blockUsage:blockQuota:blockLimit:blockInDoubt:"
    "blockGrace:filesUsage:filesQuota:filesLimit:filesInDoubt:filesGrace:"
    "remarks:fid:filesetname:\n"
)


def output(*, used=900, in_doubt=80, soft=1000, hard=1100):
    return HEADER + (
        f"mmlsquota:fileset:0:1:::wstor_work3:FILESET:119:qtholstg:"
        f"{used}:{soft}:{hard}:{in_doubt}:none:100:1000:1100:10:none:host:::\n"
    )


def test_parse_includes_in_doubt_and_uses_the_lower_positive_limit():
    usage = quota.parse_project_quota(output(), "qtholstg", "wstor_work3")
    assert usage == {"storage_percent": 98.0, "inode_percent": 11.0}
    assert quota.parse_project_quota(
        output(soft=0, hard=1000), "qtholstg", "wstor_work3"
    )["storage_percent"] == 98.0
    assert quota.parse_project_quota(
        output(soft=0, hard=0), "qtholstg", "wstor_work3"
    )["storage_percent"] == 0.0


@pytest.mark.parametrize("text,fileset", [
    ("", "qtholstg"),
    (output(), "wrong-project"),
    (output(used=-1), "qtholstg"),
    (output().replace(":900:", ":invalid:"), "qtholstg"),
    (HEADER + "mmlsquota:fileset:0:1:::wstor_work3:FILESET:119:qtholstg:\n", "qtholstg"),
])
def test_bad_or_wrong_project_quota_is_not_free_space(text, fileset):
    with pytest.raises(quota.QuotaReadError):
        quota.parse_project_quota(text, fileset, "wstor_work3")


@pytest.mark.parametrize("number,prefix", [(1, "wstor"), (3, "wstor"), (5, "lstor")])
def test_only_destination_project_and_mount_are_selected(number, prefix):
    assert quota.project_quota_target(
        Path(f"/gpfs/work{number}/0/a-project/sample/fq/nonexistent.fq.gz")
    ) == ("a-project", f"{prefix}_work{number}")


def test_destination_symlinks_are_resolved(tmp_path):
    alias = tmp_path / "data"
    alias.symlink_to("/gpfs/work3/0/qtholstg/marc", target_is_directory=True)
    assert quota.project_quota_target(alias / "fq/new.fq.gz") == ("qtholstg", "wstor_work3")


def test_quota_command_does_not_query_home_or_scratch(monkeypatch):
    calls = []

    def run(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(returncode=0, stdout=output(), stderr="")

    monkeypatch.setattr(quota.subprocess, "run", run)
    assert quota.read_project_quota("qtholstg", "wstor_work3")["storage_percent"] == 98
    assert calls == [([
        quota.MMLSQUOTA, "-e", "-Y", "-j", "qtholstg", "wstor_work3",
    ], dict(capture_output=True, text=True, timeout=30, check=False))]


@pytest.mark.parametrize("failure", [
    FileNotFoundError("no mmlsquota"),
    subprocess.TimeoutExpired("mmlsquota", 30),
    SimpleNamespace(returncode=1, stdout="", stderr="quota unavailable"),
])
def test_quota_command_failures_are_reported(monkeypatch, failure):
    def run(*args, **kwargs):
        if isinstance(failure, Exception):
            raise failure
        return failure

    monkeypatch.setattr(quota.subprocess, "run", run)
    with pytest.raises(quota.QuotaReadError):
        quota.read_project_quota("qtholstg", "wstor_work3")


def setup_wait(monkeypatch, readings, *, acquire_results=(True,)):
    readings = iter(readings)
    acquire_results = iter(acquire_results)
    events = []

    def read(fileset, filesystem):
        events.append(("check", fileset, filesystem))
        value = next(readings)
        if isinstance(value, Exception):
            raise value
        return {"storage_percent": value, "inode_percent": 10.0}

    def shrink(lease, **kwargs):
        events.append(("shrink", kwargs))

    def acquire(lease, **kwargs):
        events.append(("acquire", kwargs))
        lease["acquire"] = {"acquired": next(acquire_results)}

    monkeypatch.setattr(quota, "read_project_quota", read)
    monkeypatch.setattr(quota, "shrink_lease", shrink)
    monkeypatch.setattr(quota, "acquire_lease", acquire)
    monkeypatch.setattr(quota.random, "uniform", lambda low, high: (low + high) / 2)
    monkeypatch.setattr(quota.time, "sleep", lambda seconds: events.append(("sleep", seconds)))
    return events


DESTINATIONS = [Path("/gpfs/work3/0/qtholstg/marc/fq/raw.fq.gz"),
                Path("/gpfs/work3/0/qtholstg/marc/fq/trimmed.fq.gz")]


def test_healthy_quota_copies_immediately_with_one_query_per_fileset(monkeypatch):
    events = setup_wait(monkeypatch, [96.0])
    report = quota.wait_for_fastq_quota(DESTINATIONS, {"available": True})
    assert events == [("check", "qtholstg", "wstor_work3")]
    assert report["checks"] == 1
    assert report["wait_seconds"] == 0


def test_wait_releases_cpu_ram_then_jitters_and_rechecks_before_copy(monkeypatch):
    events = setup_wait(monkeypatch, [97.0, 96.0, 96.0, 96.0])
    report = quota.wait_for_fastq_quota(DESTINATIONS, {"available": True})
    assert [event[0] for event in events] == [
        "check", "shrink", "sleep", "check", "sleep", "check", "acquire", "check",
    ]
    assert events[1][1] == dict(cores=0.1, memory_mb=768.0, phase="fastq_quota_wait")
    assert events[2] == ("sleep", 300.0)
    assert events[4] == ("sleep", 300.0)
    assert events[6][1] == dict(cores=1.0, memory_mb=768.0, phase="fastq_copy",
                              wait_seconds=60.0, attempts=1)
    assert report["checks"] == 4
    assert report["wait_seconds"] == 600.0


def test_recovery_rebound_returns_to_polling_without_acquiring_early(monkeypatch):
    events = setup_wait(monkeypatch, [98, 96, 98, 96, 96, 96])
    quota.wait_for_fastq_quota(DESTINATIONS, {"available": True})
    assert sum(event[0] == "shrink" for event in events) == 1
    assert sum(event[0] == "acquire" for event in events) == 1
    assert sum(event[0] == "check" for event in events) == 6
    assert sum(event[0] == "sleep" for event in events) == 4


def test_quota_growth_during_lease_acquisition_releases_again(monkeypatch):
    events = setup_wait(monkeypatch, [98, 96, 96, 98, 96, 96, 96],
                        acquire_results=(True, True))
    quota.wait_for_fastq_quota(DESTINATIONS, {"available": True})
    assert sum(event[0] == "shrink" for event in events) == 2
    assert sum(event[0] == "acquire" for event in events) == 2


def test_unreadable_quota_waits_instead_of_starting_copy(monkeypatch):
    events = setup_wait(monkeypatch, [quota.QuotaReadError("timeout"), 96, 96, 96])
    quota.wait_for_fastq_quota(DESTINATIONS, {"available": True})
    assert events[1][0] == "shrink"
    assert events[2] == ("sleep", 300.0)


def test_lease_timeout_does_not_copy_without_reacquiring(monkeypatch):
    events = setup_wait(monkeypatch, [98, 96, 96, 96, 96, 96],
                        acquire_results=(False, True))
    quota.wait_for_fastq_quota(DESTINATIONS, {"available": True})
    assert sum(event[0] == "acquire" for event in events) == 2
    assert events[-1][0] == "check"


def test_disabled_does_not_even_resolve_destinations(monkeypatch):
    def unexpected(*args, **kwargs):
        raise AssertionError("disabled check must have no side effects")

    monkeypatch.setattr(quota, "project_quota_target", unexpected)
    report = quota.wait_for_fastq_quota(DESTINATIONS, {}, mode="disabled")
    assert report["enabled"] is False
    assert report["reason"] == "explicitly disabled"


def test_auto_skips_spider_and_required_rejects_non_snellius_path(monkeypatch):
    events = setup_wait(monkeypatch, [])
    destinations = [Path("/project/holstegelab/data/fq/sample.fq.gz")]
    report = quota.wait_for_fastq_quota(destinations, {})
    assert report["enabled"] is False
    assert events == []
    with pytest.raises(ValueError, match="not a Snellius project destination"):
        quota.wait_for_fastq_quota(destinations, {}, mode="required")


def test_spider_profile_explicitly_disables_fastq_quota():
    profile = yaml.safe_load((REPO / "profiles/spider/config.yaml").read_text())
    assert "external_adapter_quota_mode=disabled" in profile["config"]
    assert "storage_quota_mode=disabled" in profile["config"]


def test_spider_profile_setting_reaches_snakemake(tmp_path):
    snakefile = tmp_path / "Snakefile"
    snakefile.write_text(
        'print("FASTQ_QUOTA_MODE=" + config.get("external_adapter_quota_mode", "missing"))\n'
        'print("STORAGE_QUOTA_MODE=" + config.get("storage_quota_mode", "missing"))\n'
        'rule all:\n    output: "not-created"\n    shell: "touch {output}"\n'
    )
    result = subprocess.run(
        [sys.executable, "-m", "snakemake", "--profile", str(REPO / "profiles/spider"),
         "--snakefile", str(snakefile), "--executor", "local", "--cores", "1", "--dry-run",
         "--conda-prefix", str(tmp_path / "conda"),
         "--apptainer-prefix", str(tmp_path / "apptainer")],
        cwd=tmp_path, text=True, capture_output=True, timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "FASTQ_QUOTA_MODE=disabled" in result.stdout
    assert "STORAGE_QUOTA_MODE=disabled" in result.stdout
    assert not (tmp_path / "not-created").exists()


def test_quota_gate_is_scoped_to_fastq_producer_and_precedes_publication():
    runner = (REPO / "scripts/run_fused_external_adapter.py").read_text()
    assert runner.index("quota_wait = wait_for_fastq_quota(") < runner.index(
        "for source_path, destination in publications:"
    )
    assert '"quota_wait": quota_wait' in runner
    for filename in ("run_fused_alignment.py", "run_fused_kmer_sex.py"):
        assert "wait_for_fastq_quota" not in (REPO / "scripts" / filename).read_text()
        assert "wait_for_storage_quota" not in (REPO / "scripts" / filename).read_text()


@pytest.mark.parametrize("cores,memory", [(5.0, 512.0), (8.0, 14000.0), (0.5, 512.0)])
def test_shared_gate_restores_full_work_profile_and_discovers_lease_only_when_blocked(
    monkeypatch, cores, memory
):
    for key in quota.LEASE_ENV:
        monkeypatch.setenv(key, "test")
    monkeypatch.delenv("ZSLURM_LEASE_MAX_CORES", raising=False)
    monkeypatch.delenv("ZSLURM_LEASE_MAX_MEM_MB", raising=False)
    events = setup_wait(monkeypatch, [98, 96, 96, 96])

    def preflight(mode, command, **kwargs):
        events.append(("preflight", mode, command, kwargs))
        return {"available": True}

    monkeypatch.setattr(quota, "lease_preflight", preflight)
    quota.wait_for_storage_quota(
        DESTINATIONS, resume_cores=cores, resume_memory_mb=memory,
        lease_command="/pipeline/lease", resume_phase="compute",
    )
    assert [event[0] for event in events] == [
        "check", "preflight", "shrink", "sleep", "check", "sleep", "check", "acquire", "check",
    ]
    assert events[1][1:] == ("required", "/pipeline/lease", dict(
        initial_cores=cores, initial_memory_mb=memory,
    ))
    assert events[2][1] == dict(cores=0.1, memory_mb=min(1024.0, memory),
                              phase="storage_quota_wait")
    assert events[-2][1] == dict(cores=cores, memory_mb=memory, phase="compute",
                               wait_seconds=60.0, attempts=1)


def test_shared_gate_caps_reparsed_worker_at_original_reservation(monkeypatch):
    monkeypatch.setenv("ZSLURM_LEASE_MAX_CORES", "0.9")
    monkeypatch.setenv("ZSLURM_LEASE_MAX_MEM_MB", "512")
    events = setup_wait(monkeypatch, [98, 96, 96, 96])
    calls = []

    def preflight(mode, command, **kwargs):
        calls.append(kwargs)
        return {"available": True}

    monkeypatch.setattr(quota, "lease_preflight", preflight)
    quota.wait_for_storage_quota(DESTINATIONS, resume_cores=3, resume_memory_mb=7500)
    assert calls == [dict(initial_cores=0.9, initial_memory_mb=512.0)]
    assert events[-2][1]["cores"] == 0.9
    assert events[-2][1]["memory_mb"] == 512.0


@pytest.mark.parametrize("mode,destinations,readings", [
    ("auto", DESTINATIONS, [96]),
    ("disabled", DESTINATIONS, []),
    ("auto", [Path("/project/holstegelab/fq/new.fq.gz")], []),
])
def test_shared_gate_does_not_contact_scheduler_when_not_blocked(
    monkeypatch, mode, destinations, readings
):
    events = setup_wait(monkeypatch, readings)

    def unexpected(*args, **kwargs):
        raise AssertionError("no lease discovery for healthy/disabled/non-GPFS checks")

    monkeypatch.setattr(quota, "lease_preflight", unexpected)
    quota.wait_for_storage_quota(destinations, mode=mode)
    assert all(event[0] == "check" for event in events)
