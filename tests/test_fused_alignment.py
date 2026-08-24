import json
import os
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "scripts" / "run_fused_alignment.py"
sys.path.insert(0, str(REPO / "scripts"))

import run_fused_alignment as fused_alignment


def _script(path, body):
    path.write_text("#!/usr/bin/env python3\n" + body)
    path.chmod(0o755)
    return path


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
    assert all(call[1][-2:] == ["--phase", "alignment_tail"] for call in calls)
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
    held_cores, held_mem = 22.75, 40000
else:
    held_cores = float(args[args.index('--cores') + 1])
    held_mem = float(args[args.index('--mem-mb') + 1])
print(json.dumps({'ok': True, 'code': 0, 'status': 'granted',
                  'held_cores': held_cores, 'held_mem_mb': held_mem,
                  'max_cores': 22.75, 'max_mem_mb': 40000, 'epoch': 1}))
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
        "--initial-memory-mb", "40000", "--low-cores", "6",
        "--low-memory-mb", "15000", "--lease-mode", "required",
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
    assert metrics["lease"]["shrink"]["response"]["held_cores"] == 6
    assert [phase["label"] for phase in metrics["phases"]] == [
        "align_reads_fused.alignment", "align_reads_fused.merge_check",
        "align_reads_fused.dechimer_check", "align_reads_fused.sort",
    ]
    assert metrics["scratch_removed"] is True
    assert list((scratch / "aligner_fused").iterdir()) == []
