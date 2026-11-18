import gzip
import json
import random
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
sys.path.insert(0, str(SCRIPTS))
from adapter_processing import prepare_adapters
from hash_utils import fnv1a64


def tools():
    found = {name: shutil.which(name) for name in ("AdapterRemoval", "pigz", "samtools")}
    if not all(found.values()):
        pytest.skip("requires AdapterRemoval, pigz and samtools")
    return found


def adapter_list(tmp_path):
    path = tmp_path / "known_adapters.txt"
    path.write_text(
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\tAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n"
    )
    return path


def sequences():
    rng = random.Random(93)
    fragment = "".join(rng.choices("ACGT", k=240))
    return fragment[:151], fragment[-151:].translate(str.maketrans("ACGT", "TGCA"))[::-1]


def test_extracted_single_pair_preserves_quality_and_checksums(tmp_path):
    binaries = tools()
    reads = sequences()
    quality = "I" * 151
    sam = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100000\n"
    for flag, sequence in zip((77, 141), reads):
        sam += f"pair\t{flag}\t*\t0\t0\t*\t*\t0\t0\t{sequence}\t{quality}\n"
    bam = tmp_path / "input.bam"
    subprocess.run(
        [binaries["samtools"], "view", "-b", "-o", str(bam), "-"],
        input=sam, text=True, capture_output=True, check=True, timeout=30,
    )
    outputs = {
        "raw-forward": "raw1.gz", "raw-reverse": "raw2.gz",
        "singletons": "singletons.gz", "forward": "cut1.gz", "reverse": "cut2.gz",
        "adapter-log": "adapter.log", "fastq-stats": "stats.tsv", "adapters": "adapters.out",
    }
    command = [
        sys.executable, str(SCRIPTS / "run_fused_external_adapter.py"),
        "--input-alignment", str(bam), "--sample", "sample", "--readgroup", "rg",
        "--adapter-list", str(adapter_list(tmp_path)),
        "--fastq-stats-script", str(SCRIPTS / "fastq_stats.py"),
        "--remove-duplicates-script", str(SCRIPTS / "remove_interleaved_duplicates.py"),
        "--pair-rescue-script", str(SCRIPTS / "fastq_pair_rescue.py"),
        "--error-file", str(tmp_path / "errors"),
        "--metrics", str(tmp_path / "metrics.json"),
        "--samtools", binaries["samtools"], "--pigz", binaries["pigz"],
        "--adapter-removal", binaries["AdapterRemoval"],
        "--initial-cores", "5", "--initial-memory-mb", "14250",
        "--lease-mode", "disabled", "--quota-mode", "disabled",
        "--ssd-gb", "1", "--scratch-base", str(tmp_path), "--poll-interval", "0.01",
    ]
    for flag, filename in outputs.items():
        command.extend(["--output-" + flag, str(tmp_path / filename)])
    result = subprocess.run(command, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    stats = dict(line.split("\t") for line in (tmp_path / "stats.tsv").read_text().splitlines())
    for mate in (1, 2):
        with gzip.open(tmp_path / f"raw{mate}.gz", "rt") as handle:
            raw = handle.read().splitlines()
        with gzip.open(tmp_path / f"cut{mate}.gz", "rt") as handle:
            prepared = handle.read().splitlines()
        assert len(raw) == len(prepared) == 4
        assert raw[1] == prepared[1] == reads[mate - 1]
        assert raw[3] == prepared[3] == quality
        assert fnv1a64(prepared[3]) == int(stats[f"compare_fastq{mate}_checksum_qual"])
    metrics = json.loads((tmp_path / "metrics.json").read_text())
    assert metrics["observed_fastq_encoding"]["samtools_ambiguous_phred33"] is True


@pytest.mark.parametrize(
    "from_samtools,qualities,expected",
    [
        (True, ("I", "I"), ("I", "I")),
        (True, ("@", "J"), ("@", "J")),
        (True, ("!", "I"), ("!", "I")),
        (False, ("I", "I"), ("*", "*")),
        (False, ("h", "h"), ("I", "I")),
        (False, ("!", "I"), ("!", "I")),
        (True, ("h", "h"), ("I", "I")),
        (True, ("J", "K"), ("+", ",")),
    ],
)
def test_guard_preserves_existing_decisions_outside_ambiguous_extracted_pair(
    tmp_path, from_samtools, qualities, expected
):
    binaries = tools()
    inputs = [tmp_path / "raw1.gz", tmp_path / "raw2.gz"]
    for mate, (path, quality, sequence) in enumerate(zip(inputs, qualities, sequences()), 1):
        with gzip.open(path, "wt") as handle:
            handle.write(f"@pair/{mate}\n{sequence}\n+\n{quality * 151}\n")
    originals = [path.read_bytes() for path in inputs]
    args = SimpleNamespace(
        sample="sample", attempt=1, error_file=str(tmp_path / "errors"),
        adapter_list=str(adapter_list(tmp_path)), remove_duplicated_reads=0,
        fastq_stats_script=str(SCRIPTS / "fastq_stats.py"),
        ssd_gb=1, adapter_memory_mb=768, poll_interval=0.01,
    )
    prepared, _ = prepare_adapters(
        args, *inputs, tmp_path, [], pigz=binaries["pigz"],
        adapter_removal=binaries["AdapterRemoval"], label="test",
        from_samtools_fastq=from_samtools,
    )
    for path, quality in zip(prepared[:2], expected):
        with gzip.open(path, "rt") as handle:
            records = handle.read().splitlines()
        assert len(records) == 4
        assert records[3] == quality * 151
    assert [path.read_bytes() for path in inputs] == originals
