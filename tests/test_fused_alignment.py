import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "scripts" / "run_fused_alignment.py"
sys.path.insert(0, str(REPO / "scripts"))

import run_fused_alignment as fused_alignment


def _script(path, body):
    path.write_text("#!/usr/bin/env python3\n" + body)
    path.chmod(0o755)
    return path


def test_assigned_scratch_uses_shared_fallback_only_when_supplied(
    tmp_path, monkeypatch
):
    monkeypatch.setenv("USER", "zslurm-no-such-test-user")
    monkeypatch.setenv("SLURM_JOB_ID", "999999999")
    monkeypatch.delenv("SLURM_TMPDIR", raising=False)

    with pytest.raises(RuntimeError, match="ssd_use=required"):
        fused_alignment.assigned_scratch()

    fallback = tmp_path / "shared" / "scratch"
    assert fused_alignment.assigned_scratch(
        shared_fallback=str(fallback)
    ) == fallback.resolve()
    assert fallback.is_dir()


def test_atomic_publish_moves_a_shared_filesystem_output(tmp_path):
    source = tmp_path / "shared-scratch" / "result"
    source.parent.mkdir()
    source.write_bytes(b"complete")
    destination = tmp_path / "outputs" / "result"

    fused_alignment._atomic_publish(source, destination)

    assert destination.read_bytes() == b"complete"
    assert not source.exists()


def test_lease_request_preserves_structured_nonzero_response(monkeypatch):
    response = {
        "ok": False,
        "code": 4,
        "status": "timeout",
        "message": "capacity did not become available",
    }
    process = subprocess.CompletedProcess(
        args=["/lease"], returncode=4, stdout=json.dumps(response), stderr=""
    )
    monkeypatch.setattr(
        fused_alignment.subprocess, "run", lambda *args, **kwargs: process
    )

    with pytest.raises(fused_alignment.LeaseError) as error:
        fused_alignment._lease_request("/lease", ["set"])

    assert error.value.returncode == 4
    assert error.value.response == response


def test_shrink_lease_retries_a_temporary_memory_safety_floor(monkeypatch):
    responses = iter(
        [
            {"ok": True, "held_cores": 5, "held_mem_mb": 14250},
            {"ok": True, "held_cores": 5, "held_mem_mb": 1024},
        ]
    )
    calls = []

    def fake_request(command, arguments):
        calls.append((command, arguments))
        return next(responses)

    monkeypatch.setattr(fused_alignment, "_lease_request", fake_request)
    monkeypatch.setattr(fused_alignment.time, "sleep", lambda _seconds: None)
    lease = {"available": True, "mode": "required", "command": "/lease"}

    result = fused_alignment.shrink_lease(
        lease, cores=5, memory_mb=1024, phase="alignment_tail"
    )

    assert len(calls) == 2
    assert calls[0][1][-2:] == ["--phase", "alignment_tail"]
    assert "--phase" not in calls[1][1]
    assert result["shrink"]["performed"] is True
    assert result["shrink"]["target_reached"] is True
    assert result["shrink"]["attempts"] == 2
    assert result["shrink"]["response"]["held_mem_mb"] == 1024


def test_shrink_lease_records_a_persistent_safety_floor(monkeypatch):
    response = {"ok": True, "held_cores": 5, "held_mem_mb": 4096}
    monkeypatch.setattr(
        fused_alignment, "_lease_request", lambda _command, _arguments: response
    )
    monkeypatch.setattr(fused_alignment.time, "sleep", lambda _seconds: None)
    lease = {"available": True, "mode": "required", "command": "/lease"}

    result = fused_alignment.shrink_lease(lease, cores=5, memory_mb=1024)

    assert result["shrink"]["performed"] is True
    assert result["shrink"]["target_reached"] is False
    assert result["shrink"]["attempts"] == 11
    assert "safety floor" in result["shrink"]["reason"]


def test_required_lease_rejects_an_underfilled_growth_target(monkeypatch):
    response = {"ok": True, "held_cores": 2.5, "held_mem_mb": 768}
    monkeypatch.setattr(
        fused_alignment, "_lease_request", lambda _command, _arguments: response
    )
    monkeypatch.setattr(fused_alignment.time, "sleep", lambda _seconds: None)
    lease = {"available": True, "mode": "required", "command": "/lease"}

    with pytest.raises(fused_alignment.LeaseError, match="could not grow"):
        fused_alignment.shrink_lease(
            lease,
            cores=5,
            memory_mb=768,
            phase="adapter_removal",
        )


def test_acquire_lease_waits_for_and_confirms_growth(monkeypatch):
    calls = []

    def fake_request(command, arguments):
        calls.append((command, arguments))
        return {"ok": True, "held_cores": 5, "held_mem_mb": 2048}

    monkeypatch.setattr(fused_alignment, "_lease_request", fake_request)
    lease = {"available": True, "mode": "required", "command": "/lease"}

    result = fused_alignment.acquire_lease(
        lease,
        cores=5,
        memory_mb=768,
        phase="adapter_removal",
        wait_seconds=123,
    )

    assert calls == [
        (
            "/lease",
            [
                "set", "--cores", "5", "--mem-mb", "768", "--wait", "123",
                "--phase", "adapter_removal",
            ],
        )
    ]
    assert result["acquire"]["acquired"] is True
    # More memory than requested is a safe memory-floor result.
    assert result["acquire"]["response"]["held_mem_mb"] == 2048


def test_acquire_lease_capacity_timeout_is_nonfatal(monkeypatch):
    calls = []

    def fake_request(_command, arguments):
        calls.append(arguments)
        if arguments == ["status"]:
            return {"ok": True, "held_cores": 3, "held_mem_mb": 768}
        raise fused_alignment.LeaseError(
            "timed out",
            response={
                "ok": False,
                "code": 4,
                "status": "timeout",
                "message": "capacity did not become available",
            },
            returncode=4,
        )

    monkeypatch.setattr(fused_alignment, "_lease_request", fake_request)
    lease = {"available": True, "mode": "required", "command": "/lease"}

    result = fused_alignment.acquire_lease(
        lease, cores=5, memory_mb=768, wait_seconds=60
    )

    assert calls[-1] == ["status"]
    assert result["acquire"]["acquired"] is False
    assert result["acquire"]["timed_out"] is True
    assert result["acquire"]["attempts"] == 1
    assert result["acquire"]["response"]["held_cores"] == 3


def test_acquire_lease_verifies_growth_after_a_lost_reply(monkeypatch):
    calls = []

    def fake_request(_command, arguments):
        calls.append(arguments)
        if arguments == ["status"]:
            return {"ok": True, "held_cores": 5, "held_mem_mb": 768}
        raise fused_alignment.LeaseError("connection closed after request")

    monkeypatch.setattr(fused_alignment, "_lease_request", fake_request)
    lease = {"available": True, "mode": "required", "command": "/lease"}

    result = fused_alignment.acquire_lease(
        lease, cores=5, memory_mb=768, wait_seconds=60
    )

    assert calls[-1] == ["status"]
    assert result["acquire"]["acquired"] is True
    assert result["acquire"]["verified_after_lost_reply"] is True
    assert result["acquire"]["attempts"] == 1


def test_acquire_lease_retries_a_transient_manager_error(monkeypatch):
    set_attempts = 0
    sleeps = []

    def fake_request(_command, arguments):
        nonlocal set_attempts
        if arguments == ["status"]:
            return {"ok": True, "held_cores": 3, "held_mem_mb": 768}
        set_attempts += 1
        if set_attempts == 1:
            raise fused_alignment.LeaseError(
                "temporary manager error",
                response={"ok": False, "code": 5, "status": "manager-error"},
                returncode=5,
            )
        return {"ok": True, "held_cores": 5, "held_mem_mb": 768}

    monkeypatch.setattr(fused_alignment, "_lease_request", fake_request)
    monkeypatch.setattr(
        fused_alignment.time, "sleep", lambda seconds: sleeps.append(seconds)
    )
    lease = {"available": True, "mode": "required", "command": "/lease"}

    result = fused_alignment.acquire_lease(
        lease,
        cores=5,
        memory_mb=768,
        attempts=2,
        retry_delay_seconds=2,
    )

    assert set_attempts == 2
    assert sleeps == [2]
    assert result["acquire"]["acquired"] is True
    assert result["acquire"]["attempts"] == 2


def test_fused_runner_shrinks_lease_and_keeps_intermediate_local(tmp_path):
    tools = tmp_path / "tools"
    tools.mkdir()
    dragen = _script(
        tools / "dragen-os",
        "import sys\nsys.stderr.write('dragmap log\\n')\nsys.stdout.buffer.write(b'aligned-bam')\n",
    )
    samtools = _script(
        tools / "samtools",
        """import sys
from pathlib import Path
args = sys.argv[1:]
if args[0] == 'view' and '-o' in args:
    Path(args[args.index('-o') + 1]).write_bytes(sys.stdin.buffer.read())
elif args[0] == 'view':
    sys.stdout.buffer.write(Path(args[-1]).read_bytes())
elif args[0] == 'fixmate':
    Path(args[-1]).write_bytes(sys.stdin.buffer.read())
elif args[0] == 'sort':
    bam, bai = args[args.index('-o') + 1].split('##idx##')
    Path(bam).write_bytes(Path(args[-1]).read_bytes())
    Path(bai).write_bytes(b'index')
else:
    raise SystemExit('unexpected samtools args: ' + repr(args))
""",
    )
    bam_merge = _script(
        tools / "bam_merge",
        """import sys
from pathlib import Path
args = sys.argv[1:]
Path(args[args.index('-ua') + 1]).write_bytes(b'bad1')
Path(args[args.index('-ub') + 1]).write_bytes(b'bad2')
Path(args[args.index('-s') + 1]).write_text('primary_soft_clipped_bp_ratio\\t0.01\\n')
sys.stdout.buffer.write(sys.stdin.buffer.read())
""",
    )
    dechimer = _script(
        tools / "dechimer",
        """import sys
from pathlib import Path
args = sys.argv[1:]
Path(args[args.index('-s') + 1]).write_text('dechimered\\t1\\n')
sys.stdout.buffer.write(sys.stdin.buffer.read())
""",
    )
    checker = _script(
        tools / "checker.py",
        """import sys
from pathlib import Path
args = sys.argv[1:]
sys.stdin.buffer.read()
Path(args[args.index('-s') + 1]).write_text('compare_all\\t1\\n')
Path(args[args.index('-c') + 1]).write_text('Match\\n')
""",
    )
    lease = _script(
        tools / "zslurm_lease",
        """import json
import os
import sys
from pathlib import Path
args = sys.argv[1:]
command = args[1]
with Path(os.environ['FAKE_LEASE_LOG']).open('a') as handle:
    handle.write(command + '\\n')
if command == 'status':
    held_cores, held_mem = 22.75, 38000
else:
    held_cores = float(args[args.index('--cores') + 1])
    held_mem = float(args[args.index('--mem-mb') + 1])
print(json.dumps({'ok': True, 'code': 0, 'status': 'granted',
                  'held_cores': held_cores, 'held_mem_mb': held_mem,
                  'max_cores': 22.75, 'max_mem_mb': 38000, 'epoch': 1}))
""",
    )

    inputs = tmp_path / "inputs"
    inputs.mkdir()
    prepared1 = inputs / "cut1.fq.gz"
    prepared2 = inputs / "cut2.fq.gz"
    source1 = inputs / "raw1.fq.gz"
    source2 = inputs / "raw2.fq.gz"
    fastq_stats = inputs / "fastq.stats.tsv"
    for path in (prepared1, prepared2, source1, source2):
        path.write_bytes(b"reads")
    fastq_stats.write_text("compare_fastq_nrow\t1\n")

    outputs = tmp_path / "outputs"
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    lease_log = tmp_path / "lease.log"
    command = [
        sys.executable, str(RUNNER),
        "--prepared-fastq1", str(prepared1), "--prepared-fastq2", str(prepared2),
        "--source-fastq1", str(source1), "--source-fastq2", str(source2),
        "--fastq-stats", str(fastq_stats), "--reference-dir", str(tmp_path / "ref"),
        "--sample", "S1", "--readgroup", "RG1",
        "--output-bam", str(outputs / "final.bam"),
        "--output-bai", str(outputs / "final.bam.bai"),
        "--dragmap-log", str(outputs / "dragmap.log"),
        "--dechimer-stats", str(outputs / "dechimer.tsv"),
        "--badmap-fastq1", str(outputs / "bad1.fq.gz"),
        "--badmap-fastq2", str(outputs / "bad2.fq.gz"),
        "--merge-stats", str(outputs / "merge.tsv"),
        "--checked", str(outputs / "checked"),
        "--check-stats", str(outputs / "check.tsv"),
        "--metrics", str(outputs / "fused.io.json"),
        "--dragen", str(dragen), "--samtools", str(samtools),
        "--bam-merge", str(bam_merge), "--dechimer", str(dechimer),
        "--bam-stats", str(checker), "--initial-cores", "22.75",
        "--initial-memory-mb", "38000", "--low-cores", "2",
        "--low-memory-mb", "13500", "--lease-mode", "required",
        "--lease-command", str(lease), "--ssd-gb", "16",
        "--scratch-base", str(scratch), "--poll-interval", "0.05",
    ]
    environment = os.environ.copy()
    environment.update({
        "ZSLURM_LEASE_SOCKET": "fake.socket",
        "ZSLURM_LEASE_TOKEN": "fake-token",
        "ZSLURM_JOB_ID": "123",
        "FAKE_LEASE_LOG": str(lease_log),
    })
    subprocess.run(command, check=True, env=environment)

    assert (outputs / "final.bam").read_bytes() == b"aligned-bam"
    assert (outputs / "final.bam.bai").read_bytes() == b"index"
    assert lease_log.read_text().splitlines() == ["status", "set"]
    metrics = json.loads((outputs / "fused.io.json").read_text())
    assert metrics["success"] is True
    assert metrics["lease"]["shrink"]["response"]["held_cores"] == 2
    assert [phase["label"] for phase in metrics["phases"]] == [
        "align_reads_fused.alignment", "align_reads_fused.merge_check",
        "align_reads_fused.dechimer_check", "align_reads_fused.sort",
    ]
    assert metrics["scratch_removed"] is True
    assert list((scratch / "aligner_fused").iterdir()) == []
