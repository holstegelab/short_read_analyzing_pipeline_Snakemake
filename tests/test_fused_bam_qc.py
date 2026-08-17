import json
import os
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "scripts" / "run_fused_bam_qc.py"
sys.path.insert(0, str(REPO / "scripts"))

from run_fused_bam_qc import write_parallel_qc_script


def _script(path, body):
    path.write_text("#!/usr/bin/env python3\n" + body)
    path.chmod(0o755)
    return path


def test_bam_qc_stages_once_and_publishes_all_legacy_outputs(tmp_path):
    tools = tmp_path / "tools"
    tools.mkdir()
    samtools = _script(
        tools / "samtools",
        """import sys
args = sys.argv[1:]
if args[0] == 'stat':
    print('SN\\traw total sequences:\\t1')
elif args[0] == 'view':
    print('@HD\\tVN:1.6')
    print('r1\\t4\\t*\\t0\\t0\\t*\\t*\\t0\\t0\\tACGT\\t!!!!')
elif args[0] == 'quickcheck':
    pass
else:
    raise SystemExit('unexpected samtools args: ' + repr(args))
""",
    )
    pypy = _script(
        tools / "pypy",
        """import sys
sys.stdin.read()
print('read_count\\t1')
""",
    )
    verify = _script(
        tools / "verifybamid2",
        """import sys
from pathlib import Path
args = sys.argv[1:]
prefix = Path(args[args.index('--Output') + 1])
Path(str(prefix) + '.selfSM').write_text('FREEMIX\\n0.01\\n')
Path(str(prefix) + '.Ancestry').write_text('ANCESTRY\\nEUR\\n')
""",
    )
    gatk = _script(
        tools / "gatk",
        """import sys
from pathlib import Path
args = sys.argv[1:]
tool = args[2] if args[:1] == ['--java-options'] else args[0]
output = Path(args[args.index('-O') + 1])
if tool == 'CollectSequencingArtifactMetrics':
    for suffix in ('bait_bias_summary_metrics', 'pre_adapter_summary_metrics', 'bait_bias_detail_metrics', 'pre_adapter_detail_metrics', 'error_summary_metrics'):
        Path(str(output) + '.' + suffix).write_text(suffix + '\\n')
else:
    output.write_text(tool + '\\n')
""",
    )
    mosdepth = _script(
        tools / "mosdepth",
        """import gzip
import sys
from pathlib import Path
prefix = Path(sys.argv[-2])
with gzip.open(str(prefix) + '.regions.bed.gz', 'wt') as handle:
    handle.write('chr1\\t0\\t1\\t1\\n')
Path(str(prefix) + '.regions.bed.gz.csi').write_bytes(b'index')
for suffix in ('mosdepth.global.dist.txt', 'mosdepth.summary.txt', 'mosdepth.region.dist.txt'):
    Path(str(prefix) + '.' + suffix).write_text(suffix + '\\n')
""",
    )
    lease_log = tmp_path / "lease.jsonl"
    lease = _script(
        tools / "zslurm_lease",
        f"""import json
import sys
args = sys.argv[1:]
if 'status' in args:
    print(json.dumps({{
        'ok': True, 'status': 'current', 'held_cores': 6,
        'max_cores': 6, 'held_mem_mb': 7500, 'max_mem_mb': 7500,
        'epoch': 0,
    }}))
elif 'release' in args:
    record = {{
        'release_id': args[args.index('--release-id') + 1],
        'cores': float(args[args.index('--cores') + 1]),
        'mem_mb': float(args[args.index('--mem-mb') + 1]),
    }}
    with open({str(lease_log)!r}, 'a', encoding='utf-8') as handle:
        handle.write(json.dumps(record, sort_keys=True) + '\\n')
    print(json.dumps({{
        'ok': True, 'status': 'released', 'held_cores': 0.1,
        'max_cores': 6, 'held_mem_mb': 1504, 'max_mem_mb': 7500,
        'epoch': 1, 'release_id': record['release_id'],
        'requested_release_cores': record['cores'],
        'requested_release_mem_mb': record['mem_mb'],
        'released_cores': record['cores'],
        'released_mem_mb': record['mem_mb'], 'duplicate': False,
    }}))
else:
    raise SystemExit('unexpected lease args: ' + repr(args))
""",
    )

    inputs = tmp_path / "inputs"
    inputs.mkdir()
    bam = inputs / "sample.bam"
    bai = inputs / "sample.bam.bai"
    bam.write_bytes(b"bam")
    bai.write_bytes(b"index")
    named = {}
    for name in ("ref", "hs", "targets", "artifact", "dbsnp", "capture", "windows", "bamstats"):
        named[name] = inputs / name
        named[name].write_text("data\n")
    svd = inputs / "svd"
    (inputs / "svd.UD").write_text("svd\n")
    outputs = tmp_path / "outputs"
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    output_flags = {
        "--output-selfsm": "selfSM",
        "--output-ancestry": "Ancestry",
        "--output-hs-metrics": "hs_metrics",
        "--output-bait-summary": "bait_summary",
        "--output-pre-adapter-summary": "pre_summary",
        "--output-bait-detail": "bait_detail",
        "--output-pre-adapter-detail": "pre_detail",
        "--output-error-summary": "error_summary",
        "--output-oxog": "OXOG",
        "--output-samtools-genome": "samtools.stat",
        "--output-samtools-exome": "samtools.exome.stat",
        "--output-bamstats-all": "bam_all.tsv",
        "--output-bamstats-exome": "bam_exome.tsv",
        "--output-coverage-regions": "coverage.regions.bed.gz",
        "--output-coverage-csi": "coverage.regions.bed.gz.csi",
        "--output-coverage-global-dist": "coverage.global.txt",
        "--output-coverage-summary": "coverage.summary.txt",
        "--output-coverage-region-dist": "coverage.region.txt",
    }
    command = [
        sys.executable, str(RUNNER), "--bam", str(bam), "--bai", str(bai),
        "--sample", "S1", "--reference", str(named["ref"]),
        "--svd-prefix", str(svd), "--hs-interval", str(named["hs"]),
        "--targets-interval", str(named["targets"]),
        "--artifact-interval", str(named["artifact"]), "--dbsnp", str(named["dbsnp"]),
        "--capture-bed", str(named["capture"]), "--windows-bed", str(named["windows"]),
        "--bamstats-script", str(named["bamstats"]),
    ]
    for flag, name in output_flags.items():
        command.extend((flag, str(outputs / name)))
    command.extend(
        (
            "--metrics", str(outputs / "metrics.json"), "--samtools", str(samtools),
            "--gatk", str(gatk), "--verifybamid", str(verify),
            "--mosdepth", str(mosdepth), "--pypy", str(pypy),
            "--cores", "6", "--memory-mb", "7500", "--ssd-gb", "32",
            "--lease-mode", "required", "--lease-command", str(lease),
            "--scratch-base", str(scratch), "--poll-interval", "0.05",
        )
    )
    test_env = os.environ.copy()
    test_env.update(
        {
            "ZSLURM_LEASE_SOCKET": str(tmp_path / "fake.sock"),
            "ZSLURM_LEASE_TOKEN": "fake-token",
            "ZSLURM_JOB_ID": "fake-job",
        }
    )
    subprocess.run(command, check=True, env=test_env)

    for name in output_flags.values():
        assert (outputs / name).is_file()
    metrics = json.loads((outputs / "metrics.json").read_text())
    assert metrics["success"] is True
    assert metrics["parallel_consumers"] == 6
    assert metrics["lease"]["available"] is True
    assert len(metrics["task_releases"]) == 6
    assert all(item["performed"] for item in metrics["task_releases"])
    assert [phase["label"] for phase in metrics["phases"]] == [
        "bam_qc_fused.stage",
        "bam_qc_fused.qc",
    ]
    assert metrics["scratch_removed"] is True
    assert list((scratch / "bam_qc_fused").iterdir()) == []
    releases = [json.loads(line) for line in lease_log.read_text().splitlines()]
    assert len(releases) == 6
    assert len({release["release_id"] for release in releases}) == 6
    assert sum(release["cores"] for release in releases) == 6
    assert sum(release["mem_mb"] for release in releases) == 7500


def test_parallel_qc_releases_resources_when_a_task_fails(tmp_path):
    release_log = tmp_path / "release_ids.txt"
    lease = _script(
        tmp_path / "zslurm_lease",
        f"""import json
import sys
args = sys.argv[1:]
release_id = args[args.index('--release-id') + 1]
with open({str(release_log)!r}, 'a', encoding='utf-8') as handle:
    handle.write(release_id + '\\n')
print(json.dumps({{
    'ok': True, 'status': 'released', 'held_cores': 1,
    'max_cores': 2, 'held_mem_mb': 100, 'max_mem_mb': 200,
    'epoch': 1, 'release_id': release_id,
    'requested_release_cores': 1, 'requested_release_mem_mb': 100,
    'released_cores': 1, 'released_mem_mb': 100, 'duplicate': False,
}}))
""",
    )
    tasks = [
        {"name": "success", "command": "sleep 0.05", "cores": 1, "memory_mb": 100},
        {"name": "failure", "command": "exit 7", "cores": 1, "memory_mb": 100},
    ]
    qc_script = tmp_path / "qc.sh"
    release_paths = write_parallel_qc_script(
        qc_script,
        tasks,
        {"available": True, "command": str(lease)},
        "required",
        tmp_path,
        "S1",
    )

    result = subprocess.run(["/usr/bin/bash", str(qc_script)], check=False)

    assert result.returncode != 0
    assert sorted(release_log.read_text().splitlines()) == [
        "bam-qc:S1:failure",
        "bam-qc:S1:success",
    ]
    assert all(path.is_file() for path in release_paths.values())
