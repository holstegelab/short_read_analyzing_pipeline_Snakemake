import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "scripts" / "run_fused_merge_markdup.py"


def _script(path: Path, body: str) -> Path:
    path.write_text("#!/usr/bin/env python3\n" + body, encoding="utf-8")
    path.chmod(0o755)
    return path


def _samtools(*args: object, **kwargs) -> subprocess.CompletedProcess:
    return subprocess.run(["samtools", *map(str, args)], check=True, **kwargs)


def _readgroup_bam(root: Path, group: int) -> tuple[Path, Path, Path]:
    sequence = "ACGT" * 12 + "AC"
    quality = "I" * 50
    rg = f"rg{group}"
    sam = root / f"{rg}.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:queryname\n"
        "@SQ\tSN:chr1\tLN:10000\n"
        f"@RG\tID:{rg}\tSM:S1\tLB:lib{group}\tPL:ILLUMINA\n"
        + "".join(
            (
                f"pairA{group}\t99\tchr1\t101\t60\t50M\t=\t201\t150\t{sequence}\t{quality}\tRG:Z:{rg}\n",
                f"pairA{group}\t147\tchr1\t201\t60\t50M\t=\t101\t-150\t{sequence}\t{quality}\tRG:Z:{rg}\n",
                f"pairB{group}\t99\tchr1\t101\t60\t50M\t=\t201\t150\t{sequence}\t{quality}\tRG:Z:{rg}\n",
                f"pairB{group}\t147\tchr1\t201\t60\t50M\t=\t101\t-150\t{sequence}\t{quality}\tRG:Z:{rg}\n",
                f"secondary{group}\t256\tchr1\t401\t20\t50M\t*\t0\t0\t{sequence}\t{quality}\tRG:Z:{rg}\n",
                f"supplementary{group}\t2048\tchr1\t501\t20\t50M\t*\t0\t0\t{sequence}\t{quality}\tRG:Z:{rg}\n",
                f"unmapped{group}\t4\t*\t0\t0\t*\t*\t0\t0\t{sequence}\t{quality}\tRG:Z:{rg}\n",
            )
        ),
        encoding="utf-8",
    )
    raw = root / f"{rg}.raw.bam"
    fixed = root / f"{rg}.fixed.bam"
    bam = root / f"{rg}.sorted.bam"
    _samtools("view", "-b", "-o", raw, sam, capture_output=True)
    _samtools("fixmate", "-m", raw, fixed, capture_output=True)
    _samtools("sort", "-o", bam, fixed, capture_output=True)
    _samtools("index", bam, capture_output=True)
    marker = root / f"{rg}.bam_checked"
    marker.touch()
    return bam, Path(str(bam) + ".bai"), marker


def _records(path: Path) -> list[str]:
    result = _samtools("view", path, capture_output=True, text=True)
    return result.stdout.splitlines()


def _stable_header(path: Path) -> list[str]:
    result = _samtools("view", "-H", path, capture_output=True, text=True)
    # Scratch paths necessarily change @PG command lines; sequence/readgroup
    # dictionaries and every non-provenance header record must be identical.
    return [line for line in result.stdout.splitlines() if not line.startswith("@PG")]


def _flag_count(path: Path, flag: int) -> int:
    result = _samtools("view", "-c", "-f", flag, path, capture_output=True, text=True)
    return int(result.stdout.strip())


def _normalized_markdup_stat(path: Path) -> list[str]:
    # samtools records its literal command (including scratch paths); all
    # duplicate metrics below that provenance line must remain identical.
    return [line for line in path.read_text().splitlines() if not line.startswith("COMMAND:")]


@pytest.mark.skipif(not shutil.which("samtools"), reason="samtools required")
@pytest.mark.parametrize("groups,no_dedup", [(1, False), (2, False), (2, True)])
def test_fused_merge_markdup_matches_former_commands(tmp_path, groups, no_dedup):
    inputs = [_readgroup_bam(tmp_path, group) for group in range(groups)]
    bams = [item[0] for item in inputs]
    bais = [item[1] for item in inputs]
    checks = [item[2] for item in inputs]

    baseline_input = bams[0]
    if groups > 1:
        baseline_input = tmp_path / "baseline.merged.bam"
        _samtools("merge", "-@", "3", baseline_input, *bams, capture_output=True)
    baseline_bam = tmp_path / "baseline.markdup.bam"
    baseline_bai = tmp_path / "baseline.markdup.bam.bai"
    baseline_stat = tmp_path / "baseline.markdup.stat"
    if no_dedup:
        shutil.copyfile(baseline_input, baseline_bam)
        _samtools("index", baseline_bam, baseline_bai, capture_output=True)
        baseline_stat.touch()
    else:
        _samtools(
            "markdup", "-T", tmp_path / "baseline-temp", "-f", baseline_stat,
            "-S", "-d", "2500", baseline_input, "--write-index",
            f"{baseline_bam}##idx##{baseline_bai}", capture_output=True,
        )

    tools = tmp_path / "tools"
    tools.mkdir()
    lease_log = tmp_path / "lease.log"
    lease = _script(
        tools / "zslurm_lease",
        """import json
import os
import sys
args = sys.argv[1:]
command = args[1]
with open(os.environ['FAKE_LEASE_LOG'], 'a', encoding='utf-8') as handle:
    handle.write(command + '\\n')
if command == 'status':
    cores = float(os.environ['FAKE_INITIAL_CORES'])
    memory = float(os.environ['FAKE_INITIAL_MEMORY'])
else:
    cores = float(args[args.index('--cores') + 1])
    memory = float(args[args.index('--mem-mb') + 1])
print(json.dumps({'ok': True, 'held_cores': cores, 'held_mem_mb': memory}))
""",
    )
    outputs = tmp_path / "outputs"
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    initial_cores = 2.3 if groups > 1 else 0.95
    command = [
        sys.executable, str(RUNNER), "--input-bam", *map(str, bams),
        "--input-bai", *map(str, bais), "--check-marker", *map(str, checks),
        "--sample", "S1", "--output-bam", str(outputs / "markdup.bam"),
        "--output-bai", str(outputs / "markdup.bam.bai"),
        "--output-stat", str(outputs / "markdup.stat"),
        "--merge-log", str(outputs / "merge.log"),
        "--markdup-log", str(outputs / "markdup.log"),
        "--metrics", str(outputs / "metrics.json"),
        "--no-dedup", str(int(no_dedup)), "--merge-threads", "3",
        "--initial-cores", str(initial_cores), "--initial-memory-mb", "3000",
        "--markdup-cores", "0.95", "--markdup-memory-mb", "3000",
        "--lease-command", str(lease), "--ssd-gb", "12",
        "--scratch-base", str(scratch), "--poll-interval", "0.05",
    ]
    environment = os.environ.copy()
    environment.update(
        {
            "ZSLURM_LEASE_SOCKET": "fake.socket",
            "ZSLURM_LEASE_TOKEN": "fake-token",
            "ZSLURM_JOB_ID": "123",
            "FAKE_LEASE_LOG": str(lease_log),
            "FAKE_INITIAL_CORES": str(initial_cores),
            "FAKE_INITIAL_MEMORY": "3000",
        }
    )
    subprocess.run(command, check=True, env=environment)

    result_bam = outputs / "markdup.bam"
    _samtools("quickcheck", result_bam)
    assert _stable_header(result_bam) == _stable_header(baseline_bam)
    assert _records(result_bam) == _records(baseline_bam)
    assert _normalized_markdup_stat(outputs / "markdup.stat") == _normalized_markdup_stat(baseline_stat)
    for flag in (4, 256, 1024, 2048):
        assert _flag_count(result_bam, flag) == _flag_count(baseline_bam, flag)
    assert _flag_count(result_bam, 4) > 0
    assert _flag_count(result_bam, 256) > 0
    assert _flag_count(result_bam, 2048) > 0
    if not no_dedup:
        assert _flag_count(result_bam, 1024) > 0

    expected_lease = ["status", "set"] if groups > 1 else ["status"]
    assert lease_log.read_text().splitlines() == expected_lease
    metrics = json.loads((outputs / "metrics.json").read_text())
    assert metrics["success"] is True
    assert metrics["readgroups"] == groups
    assert metrics["scratch_removed"] is True
    assert list((scratch / "merge_markdup_fused").iterdir()) == []


def test_fused_merge_markdup_failure_publishes_no_bam_and_cleans_scratch(tmp_path):
    tools = tmp_path / "tools"
    tools.mkdir()
    samtools = _script(
        tools / "samtools",
        """import sys
if sys.argv[1] == 'markdup':
    raise SystemExit(17)
raise SystemExit('unexpected command')
""",
    )
    bam = tmp_path / "input.bam"
    bai = tmp_path / "input.bam.bai"
    checked = tmp_path / "input.bam_checked"
    bam.write_bytes(b"bam")
    bai.write_bytes(b"index")
    checked.touch()
    outputs = tmp_path / "outputs"
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = subprocess.run(
        [
            sys.executable, str(RUNNER), "--input-bam", str(bam),
            "--input-bai", str(bai), "--check-marker", str(checked),
            "--sample", "S1", "--output-bam", str(outputs / "markdup.bam"),
            "--output-bai", str(outputs / "markdup.bam.bai"),
            "--output-stat", str(outputs / "markdup.stat"),
            "--merge-log", str(outputs / "merge.log"),
            "--markdup-log", str(outputs / "markdup.log"),
            "--metrics", str(outputs / "metrics.json"), "--samtools", str(samtools),
            "--initial-cores", "0.95", "--initial-memory-mb", "3000",
            "--markdup-memory-mb", "3000", "--lease-mode", "disabled",
            "--ssd-gb", "12", "--scratch-base", str(scratch),
            "--poll-interval", "0.05",
        ],
        check=False,
    )
    assert result.returncode != 0
    assert not (outputs / "markdup.bam").exists()
    assert not (outputs / "markdup.bam.bai").exists()
    metrics = json.loads((outputs / "metrics.json").read_text())
    assert metrics["success"] is False
    assert metrics["scratch_removed"] is True
    assert list((scratch / "merge_markdup_fused").iterdir()) == []
