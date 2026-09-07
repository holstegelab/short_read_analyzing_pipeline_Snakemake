import gzip
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "scripts" / "run_fused_external_adapter.py"
sys.path.insert(0, str(REPO / "scripts"))
import select_cram_reference as REFERENCE_SELECTOR

SPEC = importlib.util.spec_from_file_location("run_fused_external_adapter", RUNNER)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def _script(path, body):
    path.write_text("#!/usr/bin/env python3\n" + body)
    path.chmod(0o755)
    return path


def test_adapter_identification_and_trimming_share_one_five_core_lease():
    runner = RUNNER.read_text()

    assert runner.count("shrink_lease(") == 1
    assert "cores=args.adapter_cores" in runner
    shared = (REPO / "scripts" / "adapter_processing.py").read_text()
    assert "--identify-adapters" in shared
    assert shared.count("--threads 4") == 2


def test_runtime_cram_reference_selection_rewrites_existing_command(
    tmp_path, monkeypatch
):
    configured = tmp_path / "hg19.fa"
    hg19_b37_chry = tmp_path / "hg19_b37chrY.fa"
    hg38 = tmp_path / "hg38.fa"
    for reference in (configured, hg19_b37_chry, hg38):
        reference.write_text(">chr1\nA\n")
    alignment = tmp_path / "input.cram"
    alignment.write_bytes(b"CRAM")
    header = (
        f"@SQ\tSN:chr1\tLN:1\tM5:{REFERENCE_SELECTOR.HG19_CHR1_M5}\n"
        f"@SQ\tSN:chrY\tLN:1\tM5:{REFERENCE_SELECTOR.HG19_B37_CHRY_M5}\n"
    )
    monkeypatch.setattr(MODULE, "read_cram_header", lambda cram, samtools: header)

    tokens, selection = MODULE.resolve_cram_options(
        alignment,
        f"--reference {configured} --input-fmt-option required_fields=0x0fff",
        samtools="samtools",
        hg19_reference=str(configured),
        hg19_b37_chry_reference=str(hg19_b37_chry),
        hg38_reference=str(hg38),
    )

    assert tokens[:2] == ["--reference", str(hg19_b37_chry)]
    assert tokens[2:] == ["--input-fmt-option", "required_fields=0x0fff"]
    assert selection == {
        "is_cram": True,
        "configured_reference": str(configured),
        "selected_reference": str(hg19_b37_chry),
        "reason": "hg19+b37-chrY M5 signature",
        "changed": True,
    }


def test_external_adapter_fusion_extracts_once_and_publishes_legacy_contract(tmp_path):
    tools = tmp_path / "tools"
    tools.mkdir()
    samtools = _script(
        tools / "samtools",
        """import gzip
import sys
from pathlib import Path
args = sys.argv[1:]
if args[0] == 'view':
    sys.stdout.buffer.write(Path(args[-1]).read_bytes())
elif args[0] in ('reset', 'sort'):
    sys.stdout.buffer.write(sys.stdin.buffer.read())
elif args[0] == 'fastq':
    sys.stdin.buffer.read()
    for flag in ('-1', '-2', '-s'):
        name = b'r2' if flag == '-2' else b'r1'
        data = b'@' + name + b'\\nACGT\\n+\\n!!!!\\n'
        with gzip.open(args[args.index(flag) + 1], 'wb') as handle:
            handle.write(data)
else:
    raise SystemExit('unexpected: ' + repr(args))
""",
    )
    pigz = _script(
        tools / "pigz",
        """import gzip
import sys
with gzip.open(sys.argv[-1], 'rb') as handle:
    sys.stdout.buffer.write(handle.read())
""",
    )
    adapter_removal = _script(
        tools / "AdapterRemoval",
        """import gzip
import sys
from pathlib import Path
args = sys.argv[1:]
data = sys.stdin.buffer.read()
if '--identify-adapters' in args:
    print('adapter')
else:
    for flag in ('--output1', '--output2'):
        with gzip.open(args[args.index(flag) + 1], 'wb') as handle:
            handle.write(b'@r1\\nACGT\\n+\\n!!!!\\n')
    Path(args[args.index('--settings') + 1]).write_text('settings\\n')
""",
    )
    fastq_stats = _script(
        tools / "fastq_stats.py",
        """import sys
from pathlib import Path
sys.stdin.buffer.read()
Path(sys.argv[sys.argv.index('-s') + 1]).write_text('compare_fastq_nrow\\t1\\n')
""",
    )
    passthrough = _script(
        tools / "passthrough.py",
        "import sys\nsys.stdout.buffer.write(sys.stdin.buffer.read())\n",
    )
    rescue = _script(
        tools / "rescue.py",
        "import sys\nraise SystemExit('rescue should not run on attempt 1')\n",
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
cores = 5 if command == 'status' else float(args[args.index('--cores') + 1])
memory = 14250 if command == 'status' else float(args[args.index('--mem-mb') + 1])
print(json.dumps({'ok': True, 'held_cores': cores, 'held_mem_mb': memory}))
""",
    )

    alignment = tmp_path / "input.bam"
    alignment.write_bytes(b"bam")
    adapters = tmp_path / "adapters.txt"
    adapters.write_text("adapter\n")
    outputs = tmp_path / "outputs"
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    lease_log = tmp_path / "lease.log"
    command = [
        sys.executable, str(RUNNER),
        "--input-alignment", str(alignment), "--sample", "S1", "--readgroup", "RG1",
        "--adapter-list", str(adapters), "--fastq-stats-script", str(fastq_stats),
        "--remove-duplicates-script", str(passthrough), "--pair-rescue-script", str(rescue),
        "--error-file", str(tmp_path / "errors"),
        "--output-raw-forward", str(outputs / "raw1.fq.gz"),
        "--output-raw-reverse", str(outputs / "raw2.fq.gz"),
        "--output-singletons", str(outputs / "singletons.fq.gz"),
        "--output-forward", str(outputs / "cut1.fq.gz"),
        "--output-reverse", str(outputs / "cut2.fq.gz"),
        "--output-adapter-log", str(outputs / "adapter.log"),
        "--output-fastq-stats", str(outputs / "stats.tsv"),
        "--output-adapters", str(outputs / "adapters.txt"),
        "--metrics", str(outputs / "metrics.json"),
        "--samtools", str(samtools), "--pigz", str(pigz),
        "--adapter-removal", str(adapter_removal),
        "--initial-cores", "5", "--initial-memory-mb", "14250",
        "--lease-command", str(lease), "--ssd-gb", "32",
        "--scratch-base", str(scratch), "--poll-interval", "0.05",
    ]
    environment = os.environ.copy()
    environment.update(
        {
            "ZSLURM_LEASE_SOCKET": "fake.socket",
            "ZSLURM_LEASE_TOKEN": "fake-token",
            "ZSLURM_JOB_ID": "123",
            "FAKE_LEASE_LOG": str(lease_log),
        }
    )
    subprocess.run(command, check=True, env=environment)

    assert gzip.open(outputs / "raw1.fq.gz", "rt").read().startswith("@r1")
    assert gzip.open(outputs / "raw2.fq.gz", "rt").read().startswith("@r2")
    assert (outputs / "singletons.fq.gz").is_file()
    assert gzip.open(outputs / "cut1.fq.gz", "rt").read().startswith("@r1")
    assert lease_log.read_text().splitlines() == ["status", "set"]
    metrics = json.loads((outputs / "metrics.json").read_text())
    assert metrics["success"] is True
    assert metrics["requested"]["adapter_cores"] == 5.0
    assert [
        adjustment["response"]["held_cores"]
        for adjustment in metrics["lease"]["adjustments"]
    ] == [5.0]
    assert [phase["label"] for phase in metrics["phases"]] == [
        "external_adapter_fused.extract",
        "external_adapter_fused.adapter_removal",
    ]
    assert metrics["scratch_removed"] is True
    assert list((scratch / "external_adapter_fused").iterdir()) == []
