import gzip
import json
import os
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "scripts" / "run_fused_chrm_tail.py"


def _script(path, body):
    path.write_text("#!/usr/bin/env python3\n" + body)
    path.chmod(0o755)
    return path


def test_chrm_tail_keeps_intermediates_local_and_publishes_final_gvcfs(tmp_path):
    tools = tmp_path / "tools"
    tools.mkdir()
    gatk = _script(
        tools / "gatk",
        """import gzip
import sys
from pathlib import Path
args = sys.argv[1:]
tool = args[2] if args[:1] == ['--java-options'] else args[0]
output = Path(args[args.index('-O') + 1])
output.parent.mkdir(parents=True, exist_ok=True)
if tool == 'MergeMutectStats':
    output.write_text('stats\\n')
else:
    with gzip.open(output, 'wt') as handle:
        handle.write('##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n')
if tool == 'Mutect2':
    Path(str(output) + '.stats').write_text('stats\\n')
""",
    )
    tabix = _script(
        tools / "tabix",
        """import sys
from pathlib import Path
Path(sys.argv[-1] + '.tbi').write_bytes(b'index')
""",
    )
    bcftools = _script(
        tools / "bcftools",
        """import gzip
import sys
from pathlib import Path
args = sys.argv[1:]
output = Path(args[args.index('-o') + 1])
with gzip.open(output, 'wt') as handle:
    handle.write('##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n')
""",
    )
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    bams = []
    for name in ("chrm.bam", "chrm.bai", "shifted.bam", "shifted.bai", "numt.bam", "numt.bai", "shifted_numt.bam", "shifted_numt.bai"):
        path = inputs / name
        path.write_bytes(b"bam-or-index")
        bams.append(path)
    original = inputs / "original.fasta"
    shifted = inputs / "shifted.fasta"
    chain = inputs / "shift.chain"
    original.write_text(">chrM\nA\n")
    shifted.write_text(">chrM\nA\n")
    chain.write_text("chain\n")
    outputs = tmp_path / "outputs"
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    command = [
        sys.executable,
        str(RUNNER),
        "--bam-chrm", str(bams[0]), "--bai-chrm", str(bams[1]),
        "--bam-shifted-chrm", str(bams[2]), "--bai-shifted-chrm", str(bams[3]),
        "--bam-numts", str(bams[4]), "--bai-numts", str(bams[5]),
        "--bam-shifted-numts", str(bams[6]), "--bai-shifted-numts", str(bams[7]),
        "--original-reference", str(original), "--shifted-reference", str(shifted),
        "--chain", str(chain), "--sample", "S1",
        "--output-chrm-gvcf", str(outputs / "chrm.g.vcf.gz"),
        "--output-chrm-tbi", str(outputs / "chrm.g.vcf.gz.tbi"),
        "--output-numt-gvcf", str(outputs / "numt.g.vcf.gz"),
        "--output-numt-tbi", str(outputs / "numt.g.vcf.gz.tbi"),
        "--metrics", str(outputs / "metrics.json"),
        "--gatk", str(gatk), "--tabix", str(tabix), "--bcftools", str(bcftools),
        "--memory-mb", "5000", "--ssd-gb", "20",
        "--scratch-base", str(scratch), "--poll-interval", "0.05",
    ]
    subprocess.run(command, check=True, env=os.environ.copy())

    for name in ("chrm.g.vcf.gz", "chrm.g.vcf.gz.tbi", "numt.g.vcf.gz", "numt.g.vcf.gz.tbi"):
        assert (outputs / name).stat().st_size > 0
    assert gzip.open(outputs / "chrm.g.vcf.gz", "rt").read().startswith("##fileformat")
    metrics = json.loads((outputs / "metrics.json").read_text())
    assert metrics["success"] is True
    assert [phase["label"] for phase in metrics["phases"]] == [
        "chrm_tail_fused.stage",
        "chrm_tail_fused.mutect_filter",
        "chrm_tail_fused.bp_resolution",
    ]
    assert metrics["scratch_removed"] is True
    assert list((scratch / "chrm_tail_fused").iterdir()) == []
