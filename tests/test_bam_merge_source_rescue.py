import gzip
import shlex
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
from run_fused_alignment import merge_source_command, source_rescue_policy


@pytest.fixture(scope="module")
def merger(tmp_path_factory):
    if not shutil.which("gcc") or not shutil.which("pigz"):
        pytest.skip("requires gcc and pigz")
    binary = tmp_path_factory.mktemp("bam-merge") / "bam_merge"
    subprocess.run(["gcc", "-O2", "-std=c11", "-o", str(binary),
                    str(SCRIPTS / "bam_merge.c"), "-lz"], check=True)
    return binary


def raw_pair(tmp_path, *, corrupt=False):
    paths = [tmp_path / f"raw{i}.fq.gz" for i in (1, 2)]
    for mate, path in enumerate(paths, 1):
        with gzip.open(path, "wb") as handle:
            for name in ("aligned", "bad", "orphan", "discarded", "high"):
                if corrupt and mate == 2 and name == "orphan":
                    continue
                seq = b"ACGT" * (5 if name == "discarded" else 15)
                qual = b"!" + b"I" * (len(seq) - 1)
                if name == "high":
                    qual = qual[:-1] + b"~"
                if corrupt and mate == 1 and name == "bad":
                    qual = qual[:10] + b"\xc8" + qual[11:]
                handle.write(f"@{name}/{mate}\n".encode() + seq + b"\n+\n" + qual + b"\n")
    return paths


def merge_args(tmp_path, label):
    return ["-ua", str(tmp_path / f"{label}.bad1.gz"),
            "-ub", str(tmp_path / f"{label}.bad2.gz"),
            "-s", str(tmp_path / f"{label}.merge.tsv")]


@pytest.mark.parametrize("dedup", [False, True])
def test_rescued_adapter_reads_restore_to_matching_bam(tmp_path, merger, dedup):
    adapter, pigz, samtools = (shutil.which(tool) for tool in ("AdapterRemoval", "pigz", "samtools"))
    if not all((adapter, pigz, samtools)):
        pytest.skip("requires AdapterRemoval, pigz and samtools")
    raw1, raw2 = raw_pair(tmp_path, corrupt=True)
    if dedup:
        # A duplicate name with different sequence must not replace the first
        # accepted fragment when the untrimmed reads are restored.
        for mate, path in enumerate((raw1, raw2), 1):
            with gzip.open(path, "ab") as handle:
                handle.write(f"@aligned/{mate}\n".encode() + b"T" * 60
                             + b"\n+\n!" + b"I" * 59 + b"\n")
    original = [path.read_bytes() for path in (raw1, raw2)]
    adapters = tmp_path / "adapters.list"
    adapters.write_text("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\tAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n")
    args = SimpleNamespace(
        attempt=2, sample="sample", error_file=str(tmp_path / "errors"),
        pair_rescue_script=str(SCRIPTS / "fastq_pair_rescue.py"),
        remove_duplicates_script=str(SCRIPTS / "remove_interleaved_duplicates.py"),
        fastq_stats_script=str(SCRIPTS / "fastq_stats.py"),
        remove_duplicated_reads=int(dedup), adapter_list=str(adapters),
        ssd_gb=1, adapter_memory_mb=768, poll_interval=0.01,
        source_fastq1=str(raw1), source_fastq2=str(raw2),
    )
    prepared, _ = prepare_adapters(args, raw1, raw2, tmp_path, [],
                                   pigz=pigz, adapter_removal=adapter, label="test")
    stats = prepared[3]
    records = [gzip.open(path, "rt").read().splitlines() for path in prepared[:2]]
    sam = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100000\n@RG\tID:rg\tSM:sample\n"
    for i in range(0, len(records[0]), 4):
        for mate, lines in enumerate(records, 1):
            name = lines[i][1:].split("/")[0]
            sam += f"{name}\t{77 if mate == 1 else 141}\t*\t0\t0\t*\t*\t0\t0\t{lines[i+1]}\t{lines[i+3]}\tRG:Z:rg\n"
    assert "discarded\t" not in sam
    # Demonstrate the original failure: raw leftovers reintroduce invalid QUAL.
    legacy = subprocess.run([str(merger), "-a", str(raw1), "-b", str(raw2),
                             *merge_args(tmp_path, "legacy")], input=sam.encode(),
                            capture_output=True, timeout=30, check=True)
    bad = subprocess.run([samtools, "view", "-b", "-o", str(tmp_path / "bad.bam"), "-"],
                         input=legacy.stdout, capture_output=True, timeout=30)
    assert bad.returncode != 0
    assert b"invalid QUAL character" in bad.stderr
    command = shlex.split(merge_source_command(args, str(merger), source_rescue_policy(str(stats))))
    fixed = subprocess.run([*command, *merge_args(tmp_path, "fixed")], input=sam.encode(),
                           capture_output=True, timeout=30, check=True)
    lines = [line.split(b"\t") for line in fixed.stdout.splitlines() if not line.startswith(b"@")]
    assert {line[0] for line in lines} == {b"aligned", b"discarded", b"high"}
    assert len(lines) == 6
    assert all(len(line[9]) == len(line[10]) for line in lines)
    assert any(b"~" in line[10] for line in lines)
    bam = tmp_path / "restored.bam"
    subprocess.run([samtools, "view", "-b", "-o", str(bam), "-"],
                   input=fixed.stdout, capture_output=True, check=True)
    checked = tmp_path / "checked"
    result = subprocess.run([sys.executable, str(SCRIPTS / "bam_stats_compare_hts.py"),
                             "-i", str(bam), "--fastq-stats", str(stats),
                             "-s", str(tmp_path / "check.tsv"), "-c", str(checked)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert checked.read_text().strip() == "Match"
    assert [path.read_bytes() for path in (raw1, raw2)] == original


def test_clean_interleaved_and_separate_fastqs_restore_identically(tmp_path, merger):
    raw1, raw2 = raw_pair(tmp_path)
    sam = b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n"
    separate = subprocess.run([str(merger), "-a", str(raw1), "-b", str(raw2),
                               *merge_args(tmp_path, "separate")],
                              input=sam, capture_output=True, check=True, timeout=30)
    streamed = subprocess.run([sys.executable, str(SCRIPTS / "run_bam_merge_rescued.py"),
                               "--r1", str(raw1), "--r2", str(raw2), "--", str(merger),
                               *merge_args(tmp_path, "streamed")],
                              input=sam, capture_output=True, check=True, timeout=30)
    assert separate.stdout == streamed.stdout
    assert (tmp_path / "separate.merge.tsv").read_bytes() == (tmp_path / "streamed.merge.tsv").read_bytes()


@pytest.mark.parametrize("failure", ["missing", "gzip_crc"])
def test_producer_failure_is_not_hidden_by_successful_merger(tmp_path, failure):
    if not shutil.which("pigz"):
        pytest.skip("requires pigz")
    raw1, raw2 = raw_pair(tmp_path)
    if failure == "missing":
        raw1.unlink()
    else:
        data = bytearray(raw1.read_bytes())
        data[-8] ^= 1
        raw1.write_bytes(data)
    consumer = "import sys; f=open(sys.argv[-1], 'rb'); f.read()"
    result = subprocess.run([sys.executable, str(SCRIPTS / "run_bam_merge_rescued.py"),
                             "--r1", str(raw1), "--r2", str(raw2), "--",
                             sys.executable, "-c", consumer],
                            input=b"", capture_output=True, timeout=20)
    assert result.returncode != 0
    assert b"producer exited" in result.stderr


def test_early_merger_failure_does_not_wait_for_all_source_reads(tmp_path):
    if not shutil.which("pigz"):
        pytest.skip("requires pigz")
    raw1, raw2 = raw_pair(tmp_path)
    with gzip.open(raw1, "ab") as handle:
        handle.write(b"@extra/1\nACGT\n+\nIIII\n" * 200000)
    result = subprocess.run([sys.executable, str(SCRIPTS / "run_bam_merge_rescued.py"),
                             "--r1", str(raw1), "--r2", str(raw2), "--",
                             sys.executable, "-c", "raise SystemExit(17)"],
                            input=b"", capture_output=True, timeout=20)
    assert result.returncode == 17
