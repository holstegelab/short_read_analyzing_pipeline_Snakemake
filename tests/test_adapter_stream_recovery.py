import gzip
import io
import shutil
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import adapter_processing
import fastq_pair_rescue


@pytest.mark.parametrize("bad_quality", ["II\x00I", "II\ufffdI", "II\u00e9I", "II I"])
def test_rescue_rejects_quality_characters_that_lengths_alone_accept(bad_quality):
    with pytest.raises(ValueError, match="non-printable/non-ASCII"):
        fastq_pair_rescue._read_fastq4(
            io.StringIO(f"@r/1\nACGT\n+\n{bad_quality}\n"),
            fastq_pair_rescue._ReaderState(),
        )


def test_external_reader_stops_unread_producer_without_collecting_output(monkeypatch):
    proc = subprocess.Popen(
        [sys.executable, "-c",
         "import sys; sys.stdout.write('ready\\n'); sys.stdout.flush(); "
         "exec('while True: sys.stdout.write(\"x\" * 65536); sys.stdout.flush()')"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    reader = fastq_pair_rescue._ProcTextReader(proc)
    try:
        assert reader.readline() == "ready\n"

        def forbid_collecting(*args, **kwargs):
            pytest.fail("cleanup must not collect the unread FASTQ")

        monkeypatch.setattr(proc, "communicate", forbid_collecting)
        started = time.monotonic()
        reader.close()
        assert time.monotonic() - started < 3
        assert proc.poll() is not None
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()


def test_external_reader_reports_decompression_failure():
    proc = subprocess.Popen(
        [sys.executable, "-c", "import sys; sys.stderr.write('bad gzip'); sys.exit(2)"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    reader = fastq_pair_rescue._ProcTextReader(proc)
    assert reader.readline() == ""
    with pytest.raises(RuntimeError, match="bad gzip"):
        reader.close()


@pytest.mark.parametrize("attempt,previous_failure", [(1, False), (1, True), (2, False)])
def test_real_adapter_pipeline_preserves_high_qualities_and_rescues_invalid_pairs(
    tmp_path, capfd, attempt, previous_failure
):
    adapter = shutil.which("AdapterRemoval")
    pigz = shutil.which("pigz")
    if not adapter or not pigz:
        pytest.skip("requires AdapterRemoval and pigz")
    rescue = previous_failure or attempt >= 2
    bad = {1001: b"\xc8", 1002: b"\0", 1003: b"\xc3\xa9"} if rescue else {}
    inputs = [tmp_path / "R1.fastq.gz", tmp_path / "R2.fastq.gz"]
    for mate, path in enumerate(inputs, 1):
        with gzip.open(path, "wb") as handle:
            for i in range(1200):
                quality = b"!" + b"I" * 51
                if i == 1100:
                    quality = quality[:-1] + b"~"  # Phred 93, beyond the initial probe.
                if mate == 1 and i in bad:
                    quality = quality[:20] + bad[i] + quality[21:]
                handle.write(
                    f"@read{i}/{mate}\n".encode() + b"ACGT" * 13
                    + b"\n+\n" + quality + b"\n"
                )
    error_file = tmp_path / "errors"
    if previous_failure:
        error_file.write_text("sample\n")
    adapters = tmp_path / "known_adapters.txt"
    adapters.write_text("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\tAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n")
    args = SimpleNamespace(
        attempt=attempt, sample="sample", error_file=str(error_file),
        pair_rescue_script=str(REPO / "scripts/fastq_pair_rescue.py"),
        fastq_stats_script=str(REPO / "scripts/fastq_stats.py"),
        remove_duplicated_reads=0, adapter_list=str(adapters),
        ssd_gb=1, adapter_memory_mb=768, poll_interval=0.01,
    )
    prepared, _ = adapter_processing.prepare_adapters(
        args, *inputs, tmp_path, [], pigz=pigz, adapter_removal=adapter, label="test"
    )
    stats = dict(line.split("\t") for line in prepared[3].read_text().splitlines())
    assert stats["source_pair_rescue"] == str(int(rescue))
    assert stats["source_remove_duplicates"] == "0"
    for path in prepared[:2]:
        lines = gzip.open(path, "rb").read().splitlines()
        assert len(lines) == (1200 - len(bad)) * 4
        assert any(b"~" in quality for quality in lines[3::4])
        assert all(len(seq) == len(qual) for seq, qual in zip(lines[1::4], lines[3::4]))
        if rescue:
            assert not any(b"read1001/" in header for header in lines[::4])
    if rescue:
        assert "R1 total: 1197, corrupt: 3" in capfd.readouterr().err
