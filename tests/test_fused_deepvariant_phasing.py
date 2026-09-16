import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "scripts" / "run_fused_deepvariant_phasing.py"


def _script(path, body):
    path.write_text("#!/usr/bin/env python3\n" + body)
    path.chmod(0o755)
    return path


@pytest.mark.parametrize("shared_fallback", [False, True])
@pytest.mark.parametrize("mode", ["success", "deepvariant_failure", "phasing_failure", "skip_sex"])
def test_fused_deepvariant_phasing_keeps_raw_calls_local(tmp_path, shared_fallback, mode):
    tools = tmp_path / "tools"
    tools.mkdir()
    deepvariant = _script(
        tools / "run_deepvariant",
        """import importlib.util
import json
import os
import sys
import tempfile
from pathlib import Path
Path(os.environ['FAKE_DV_ENV_LOG']).write_text(
    os.environ.get('ZSLURM_LEASE_TOKEN', '') + '|' +
    os.environ.get('ZSLURM_LEASE_SOCKET', '')
)
temporary = Path(tempfile.mkdtemp(prefix='Bazel.runfiles_'))
source = temporary / '__autograph_generated_file_test.py'
source.write_text('VALUE = 42')
spec = importlib.util.spec_from_file_location('auto_test', source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
assert module.VALUE == 42
assert Path(module.__cached__).is_file()
Path(os.environ['FAKE_TEMP_LOG']).write_text(json.dumps({
    'env': {key: os.environ.get(key) for key in ('TMPDIR', 'TMP', 'TEMP', 'TEMPDIR', 'XDG_CACHE_HOME', 'PYTHONPYCACHEPREFIX')},
    'tempfile': str(temporary), 'bytecode': module.__cached__,
}))
if os.environ['FAKE_MODE'] == 'deepvariant_failure':
    raise SystemExit(7)
values = {a.split('=', 1)[0]: a.split('=', 1)[1] for a in sys.argv[1:] if '=' in a}
for key in ('--output_vcf', '--output_gvcf'):
    out = Path(values[key])
    out.write_text('vcf\\n')
    Path(str(out) + '.tbi').write_bytes(b'index')
""",
    )
    whatshap = _script(
        tools / "whatshap",
        """import os
import shutil
import sys
import tempfile
from pathlib import Path
Path(os.environ['FAKE_PHASE_TEMP_LOG']).write_text(tempfile.mkdtemp(prefix='phasing_'))
if os.environ['FAKE_MODE'] == 'phasing_failure':
    raise SystemExit(8)
args = sys.argv[1:]
if args[0] == 'phase':
    shutil.copyfile(args[args.index('--reference') + 2], args[args.index('-o') + 1])
elif args[0] == 'stats':
    print('phased\\t1')
""",
    )
    bcftools = _script(
        tools / "bcftools",
        """import shutil
import sys
from pathlib import Path
args = sys.argv[1:]
if args[0] == 'index':
    Path(args[-1] + '.tbi').write_bytes(b'index')
elif args[0] == 'stats':
    print('SN\\t0\\tnumber of records:\\t1')
elif args[0] == 'view':
    out = Path(args[args.index('-o') + 1])
    source = Path(args[1]) if args[1] != '-R' else Path(args[args.index('-O') - 1])
    shutil.copyfile(source, out)
else:
    raise SystemExit('unexpected: ' + repr(args))
""",
    )
    bgzip = _script(tools / "bgzip", "import sys\nsys.stdout.buffer.write(sys.stdin.buffer.read())\n")
    tabix = _script(
        tools / "tabix",
        """import sys
from pathlib import Path
Path(sys.argv[-1] + '.tbi').write_bytes(b'index')
""",
    )
    merge = _script(
        tools / "merge.py",
        """import shutil
import sys
from pathlib import Path
shutil.copyfile(sys.argv[1], sys.argv[3])
Path(sys.argv[4]).write_text('merged\\t1\\n')
""",
    )
    stats_parser = _script(
        tools / "stats.py",
        """import sys
from pathlib import Path
Path(sys.argv[2]).write_text('sample\\tvalue\\nS1\\t1\\n')
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
cores = 8 if command == 'status' else float(args[args.index('--cores') + 1])
memory = 10000 if command == 'status' else float(args[args.index('--mem-mb') + 1])
print(json.dumps({'ok': True, 'held_cores': cores, 'held_mem_mb': memory}))
""",
    )

    inputs = tmp_path / "inputs"
    inputs.mkdir()
    named = {}
    for name in ("bed", "bam", "bai", "sex", "ref", "dvref", "auto", "x", "y"):
        named[name] = inputs / name
        named[name].write_text("data\n")
    outputs = tmp_path / "outputs"
    scratch = tmp_path / "scratch with 'quotes'"
    scratch.mkdir()
    lease_log = tmp_path / "lease.log"
    dv_env_log = tmp_path / "deepvariant.env"
    output_args = {
        "--output-vcf": outputs / "phased.vcf.gz",
        "--output-vcf-tbi": outputs / "phased.vcf.gz.tbi",
        "--output-wstats": outputs / "whatshap.stats",
        "--output-merge-stats": outputs / "merge.stats",
        "--output-bcftools-stats": outputs / "bcftools.stats",
        "--output-bcftools-summary": outputs / "summary.tsv",
        "--output-gvcf": outputs / "merged.g.vcf.gz",
        "--output-gvcf-tbi": outputs / "merged.g.vcf.gz.tbi",
        "--output-exome-gvcf": outputs / "exome.g.vcf.gz",
        "--output-exome-gvcf-tbi": outputs / "exome.g.vcf.gz.tbi",
    }
    command = [
        sys.executable,
        str(RUNNER),
        "--sample", "S1", "--region", "A0", "--bed", str(named["bed"]),
        "--bam", str(named["bam"]), "--bai", str(named["bai"]),
        "--validated-sex", str(named["sex"]), "--reference", str(named["ref"]),
        "--deepvariant-reference", str(named["dvref"]),
        "--deepvariant-runner", str(deepvariant), "--model-type", "WGS",
        "--haploid-contigs", "chrNONE", "--ploidy", "2", "--skip-sex", str(int(mode == "skip_sex")),
        "--interval-bed", str(named["bed"]),
        "--capture-auto-bed", str(named["auto"]),
        "--capture-x-bed", str(named["x"]), "--capture-y-bed", str(named["y"]),
        "--merge-script", str(merge), "--stats-parser", str(stats_parser),
    ]
    for flag, path in output_args.items():
        command.extend([flag, str(path)])
    command.extend(
        [
            "--metrics", str(outputs / "metrics.json"),
            "--whatshap", str(whatshap), "--bcftools", str(bcftools),
            "--bgzip", str(bgzip), "--tabix", str(tabix),
            "--initial-cores", "8", "--initial-memory-mb", "10000",
            "--attempt", "3",
            "--lease-command", str(lease), "--ssd-gb", "16",
            "--scratch-base", str(scratch), "--poll-interval", "0.05",
        ]
    )
    environment = os.environ.copy()
    inherited_tmp = tmp_path / "inherited shared tmp"
    inherited_tmp.mkdir()
    inherited_cache = tmp_path / "inherited shared cache"
    if shared_fallback:
        command[command.index("--scratch-base") + 1] = str(tmp_path / "absent SSD")
        command.extend(["--shared-scratch-base", str(scratch)])
    environment.update(
        {
            "ZSLURM_LEASE_SOCKET": "fake.socket",
            "ZSLURM_LEASE_TOKEN": "fake-token",
            "ZSLURM_JOB_ID": "123",
            "FAKE_LEASE_LOG": str(lease_log),
            "FAKE_DV_ENV_LOG": str(dv_env_log),
            "FAKE_MODE": mode,
            "FAKE_TEMP_LOG": str(tmp_path / "temp.json"),
            "FAKE_PHASE_TEMP_LOG": str(tmp_path / "phase-temp.txt"),
            "TMPDIR": str(inherited_tmp),
            "TMP": str(inherited_tmp),
            "TEMP": str(inherited_tmp),
            "TEMPDIR": str(inherited_tmp),
            "XDG_CACHE_HOME": str(inherited_cache),
            "PYTHONPYCACHEPREFIX": str(tmp_path / "inherited bytecode"),
        }
    )
    environment.pop("PYTHONDONTWRITEBYTECODE", None)
    result = subprocess.run(command, check=False, env=environment, capture_output=True, text=True)
    success = mode in ("success", "skip_sex")
    assert (result.returncode == 0) == success, result.stderr
    assert "[deepvariant_phasing_fused] scratch=" in result.stderr

    assert all(path.exists() == success for path in output_args.values())
    assert not (outputs / "tmp.g.vcf").exists()
    assert lease_log.read_text().splitlines() == (["status"] if mode == "deepvariant_failure" else ["status", "set"])
    if mode != "skip_sex":
        assert dv_env_log.read_text() == "|"
    metrics = json.loads((outputs / "metrics.json").read_text())
    assert metrics["success"] is success
    assert metrics["attempt"] == 3
    assert metrics["requested"]["low_memory_mb"] == 4000
    expected_phases = [
        "deepvariant_phasing_fused.deepvariant",
        "deepvariant_phasing_fused.phasing_merge",
    ]
    if mode == "deepvariant_failure":
        expected_phases.pop()
    assert [phase["label"] for phase in metrics["phases"]] == expected_phases
    assert metrics["scratch_removed"] is True
    assert list((scratch / "deepvariant_phasing_fused").iterdir()) == []
    assert list(inherited_tmp.iterdir()) == []
    assert not inherited_cache.exists()
    job_tmp = Path(metrics["scratch_job_directory"])
    assert job_tmp.is_relative_to(scratch)
    assert metrics["temporary_directory"] == str(job_tmp / "tmp")
    assert metrics["cache_directory"] == str(job_tmp / "cache")
    if mode != "skip_sex":
        observed = json.loads((tmp_path / "temp.json").read_text())
        for name in ("TMPDIR", "TMP", "TEMP", "TEMPDIR"):
            assert observed["env"][name] == str(job_tmp / "tmp")
        assert observed["env"]["XDG_CACHE_HOME"] == str(job_tmp / "cache")
        assert observed["env"]["PYTHONPYCACHEPREFIX"] == str(job_tmp / "cache" / "python")
        assert Path(observed["tempfile"]).is_relative_to(job_tmp / "tmp")
        assert Path(observed["bytecode"]).is_relative_to(job_tmp / "cache" / "python")
    if mode in ("success", "phasing_failure"):
        assert Path((tmp_path / "phase-temp.txt").read_text()).is_relative_to(job_tmp / "tmp")
