import json
import os
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "scripts" / "run_fused_kmer_sex.py"


def _script(path, body):
    path.write_text("#!/usr/bin/env python3\n" + body)
    path.chmod(0o755)
    return path


def test_fused_kmer_sex_keeps_database_local_and_shrinks(tmp_path):
    tools = tmp_path / "tools"
    tools.mkdir()
    kmc = _script(
        tools / "kmc",
        """import sys
from pathlib import Path
prefix = Path(sys.argv[-2])
Path(str(prefix) + '.kmc_pre').write_bytes(b'pre')
Path(str(prefix) + '.kmc_suf').write_bytes(b'suf')
""",
    )
    kmc_tools = _script(
        tools / "kmc_tools",
        """import sys
from pathlib import Path
args = sys.argv[1:]
if 'simple' in args:
    out = Path(args[-2])
    Path(str(out) + '.kmc_pre').write_bytes(b'pre')
    Path(str(out) + '.kmc_suf').write_bytes(b'suf')
elif 'transform' in args:
    Path(args[-1]).write_text('kmer\\tcount\\nAAAA\\t2\\n')
else:
    raise SystemExit('unexpected args: ' + repr(args))
""",
    )
    process_sex = _script(
        tools / "process_sex.py",
        """import sys
from pathlib import Path
Path(sys.argv[-1]).write_text('sex: F\\n')
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
cores = 2 if command == 'status' else float(args[args.index('--cores') + 1])
memory = 36000 if command == 'status' else float(args[args.index('--mem-mb') + 1])
print(json.dumps({'ok': True, 'held_cores': cores, 'held_mem_mb': memory}))
""",
    )

    fastqs = [tmp_path / "r1.fq.gz", tmp_path / "r2.fq.gz"]
    for path in fastqs:
        path.write_bytes(b"reads")
    references = []
    for name in ("y", "x", "m", "a"):
        prefix = tmp_path / f"ref-{name}"
        Path(str(prefix) + ".kmc_pre").write_bytes(b"pre")
        Path(str(prefix) + ".kmc_suf").write_bytes(b"suf")
        references.append(prefix)

    outputs = tmp_path / "outputs"
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    lease_log = tmp_path / "lease.log"
    command = [
        sys.executable,
        str(RUNNER),
        "--fastq",
        *map(str, fastqs),
        "--sample",
        "S1",
        "--output-yaml",
        str(outputs / "result.yaml"),
        "--output-chry",
        str(outputs / "chry.tsv"),
        "--output-chrx",
        str(outputs / "chrx.tsv"),
        "--output-chrm",
        str(outputs / "chrm.tsv"),
        "--output-auto",
        str(outputs / "auto.tsv"),
        "--kmer-chry",
        str(references[0]),
        "--kmer-chrx",
        str(references[1]),
        "--kmer-chrm",
        str(references[2]),
        "--kmer-auto",
        str(references[3]),
        "--process-sex",
        str(process_sex),
        "--metrics",
        str(outputs / "metrics.json"),
        "--kmc",
        str(kmc),
        "--kmc-tools",
        str(kmc_tools),
        "--initial-cores",
        "2",
        "--initial-memory-mb",
        "36000",
        "--lease-command",
        str(lease),
        "--ssd-gb",
        "32",
        "--scratch-base",
        str(scratch),
        "--poll-interval",
        "0.05",
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

    assert (outputs / "result.yaml").read_text() == "sex: F\n"
    assert lease_log.read_text().splitlines() == ["status", "set"]
    metrics = json.loads((outputs / "metrics.json").read_text())
    assert metrics["success"] is True
    assert [phase["label"] for phase in metrics["phases"]] == [
        "kmer_sex_fused.kmc",
        "kmer_sex_fused.sex",
    ]
    assert metrics["scratch_removed"] is True
    assert list((scratch / "kmer_sex_fused").iterdir()) == []
