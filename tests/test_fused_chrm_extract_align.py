import json
import os
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "scripts" / "run_fused_chrm_extract_align.py"


def _script(path, body):
    path.write_text("#!/usr/bin/env python3\n" + body)
    path.chmod(0o755)
    return path


def test_chrm_extract_align_keeps_fastqs_local_and_publishes_bams(tmp_path):
    tools = tmp_path / "tools"
    tools.mkdir()
    samtools = _script(
        tools / "samtools",
        """import gzip
import sys
from pathlib import Path
args = sys.argv[1:]
command = args[0]
if command == 'view':
    Path(args[args.index('-o') + 1]).write_bytes(b'view-bam')
elif command in ('sort', 'collate'):
    output = Path(args[args.index('-o') + 1])
    if not sys.stdin.isatty():
        data = sys.stdin.buffer.read()
    else:
        data = b''
    source = Path(args[-1])
    output.write_bytes(data or (source.read_bytes() if source.is_file() else b'sorted-bam'))
elif command == 'fastq':
    for flag, name in (('-1', b'r1'), ('-2', b'r2')):
        with gzip.open(args[args.index(flag) + 1], 'wb') as handle:
            handle.write(b'@' + name + b'\\nACGT\\n+\\n!!!!\\n')
elif command == 'index':
    Path(args[args.index('-o') + 1]).write_bytes(b'index')
else:
    raise SystemExit('unexpected samtools args: ' + repr(args))
""",
    )
    bwa = _script(tools / "bwa", "import sys\nsys.stdout.buffer.write(b'aligned-bam')\n")

    input_bam = tmp_path / "markdup.bam"
    input_bam.write_bytes(b"bam")
    numts = tmp_path / "numts.bed"
    numts.write_text("chr1\t1\t2\n")
    original = tmp_path / "original.fasta"
    shifted = tmp_path / "shifted.fasta"
    original.write_text(">chrM\nA\n")
    shifted.write_text(">chrM\nA\n")
    outputs = tmp_path / "outputs"
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    names = (
        "chrm.bam", "chrm.bai", "shifted_chrm.bam", "shifted_chrm.bai",
        "numts.bam", "numts.bai", "shifted_numts.bam", "shifted_numts.bai",
    )
    flags = (
        "--output-bam-chrm", "--output-bai-chrm",
        "--output-bam-shifted-chrm", "--output-bai-shifted-chrm",
        "--output-bam-numts", "--output-bai-numts",
        "--output-bam-shifted-numts", "--output-bai-shifted-numts",
    )
    command = [
        sys.executable, str(RUNNER), "--input-bam", str(input_bam),
        "--numts-bed", str(numts), "--original-reference", str(original),
        "--shifted-reference", str(shifted), "--sample", "S1",
    ]
    for flag, name in zip(flags, names):
        command.extend((flag, str(outputs / name)))
    command.extend(
        (
            "--metrics", str(outputs / "metrics.json"), "--samtools", str(samtools),
            "--bwa", str(bwa), "--threads", "2", "--memory-mb", "4000",
            "--ssd-gb", "20", "--scratch-base", str(scratch),
            "--poll-interval", "0.05",
        )
    )
    subprocess.run(command, check=True, env=os.environ.copy())

    for name in names:
        assert (outputs / name).stat().st_size > 0
    metrics = json.loads((outputs / "metrics.json").read_text())
    assert metrics["success"] is True
    assert [phase["label"] for phase in metrics["phases"]] == [
        "chrm_extract_align_fused.extract",
        "chrm_extract_align_fused.align",
    ]
    assert metrics["scratch_removed"] is True
    assert list((scratch / "chrm_extract_align_fused").iterdir()) == []
