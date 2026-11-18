"""Opt-in real KMC regression tests; run on compute, not a login node.

Set KMC_TEST_BINARY and KMC_TEST_TOOLS to a post-deploy build. The sampling
regression streams about 1 GiB of synthetic FASTQ, with exact expected counts.
"""
from collections import Counter
import gzip
import os
from pathlib import Path
import subprocess

import pytest


KMC = os.environ.get("KMC_TEST_BINARY")
TOOLS = os.environ.get("KMC_TEST_TOOLS")
pytestmark = pytest.mark.skipif(
    not KMC or not TOOLS, reason="Opt-in real KMC tests require explicit binaries"
)
SEQUENCE = "ACGT" * 40
RECORD = ("@synthetic\n" + SEQUENCE + "\n+\n" + "I" * len(SEQUENCE) + "\n").encode()


def expected_counts(records):
    counts = Counter()
    for pos in range(len(SEQUENCE) - 31):
        kmer = SEQUENCE[pos : pos + 32]
        reverse = kmer.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        counts[min(kmer, reverse)] += records
    return dict(counts)


def run_kmc(tmp_path, inputs, readers=12):
    tmp_path.mkdir(parents=True, exist_ok=True)
    manifest = tmp_path / "inputs.lst"
    manifest.write_text("".join(str(p.resolve()) + "\n" for p in inputs))
    database = tmp_path / "db"
    result = subprocess.run(
        [KMC, "-fq", "-k32", "-ci1", "-cx4294967295", "-cs4294967295",
         f"-sf{readers}", "-sp12", "-sr1", "-m2", "@" + str(manifest),
         str(database), str(tmp_path)],
        capture_output=True, text=True, timeout=180,
    )
    (tmp_path / "stdout.log").write_text(result.stdout)
    (tmp_path / "stderr.log").write_text(result.stderr)
    return result, database


def dump_counts(tmp_path, database):
    dump = tmp_path / "counts.tsv"
    subprocess.run([TOOLS, "transform", str(database), "dump", str(dump)],
                   check=True, capture_output=True, timeout=60)
    return {fields[0]: int(fields[1]) for line in dump.read_text().splitlines()
            if (fields := line.split())}


def test_single_and_concatenated_gzip_members_have_exact_counts(tmp_path):
    for name, contents in {
        "single": gzip.compress(RECORD * 1001),
        "concatenated": gzip.compress(RECORD * 300) + gzip.compress(RECORD * 701),
    }.items():
        path = tmp_path / (name + ".fq.gz")
        path.write_bytes(contents)
        result, database = run_kmc(tmp_path / name, [path])
        assert result.returncode == 0, result.stderr
        assert dump_counts(tmp_path / name, database) == expected_counts(1001)


@pytest.mark.parametrize("damage", ["truncated_trailer", "bad_crc"])
def test_gzip_corruption_stays_fatal_during_sampling(tmp_path, damage):
    contents = bytearray(gzip.compress(RECORD * 1001))
    if damage == "truncated_trailer":
        contents = contents[:-6]
    else:
        contents[-8] ^= 0x80
    path = tmp_path / "corrupt.fq.gz"
    path.write_bytes(contents)
    result, _ = run_kmc(tmp_path / "run", [path])
    assert result.returncode != 0
    assert "gzip file" in result.stderr


@pytest.fixture
def large_member(tmp_path):
    # Stored DEFLATE forces multiple binary-reader packs; highly compressed
    # toy inputs can reach physical EOF before cancellation and miss the race.
    records = 262144
    path = tmp_path / "large.fq.gz"
    with gzip.open(path, "wb", compresslevel=0) as handle:
        block = RECORD * 8192
        for _ in range(records // 8192):
            handle.write(block)
    return path, records


def test_parallel_sampling_keeps_every_read_and_exact_kmer_counts(tmp_path, large_member):
    path, records = large_member
    inputs = []
    for index in range(12):
        linked = tmp_path / f"reads_{index}.fq.gz"
        os.link(path, linked)
        inputs.append(linked)
    # Both configurations run through Stage 2 to completion. The one-reader
    # run is a test oracle only, never the production workaround.
    for readers in (12, 1):
        output = tmp_path / f"run_{readers}"
        result, database = run_kmc(output, inputs, readers=readers)
        assert result.returncode == 0, result.stderr
        assert "Stage 2:" in result.stderr
        assert dump_counts(output, database) == expected_counts(records * 12)


def test_truncation_after_sampling_stays_fatal_in_full_counting(tmp_path, large_member):
    path, _ = large_member
    corrupt = tmp_path / "late_truncation.fq.gz"
    # Four members exceed the 256-MiB stats quota. Damage is only in the last
    # member: sampling should pass, but the full-counting reader must reject it.
    with corrupt.open("wb") as output:
        for _ in range(4):
            with path.open("rb") as source:
                while block := source.read(1024 * 1024):
                    output.write(block)
        output.truncate(output.tell() - 6)
    result, _ = run_kmc(tmp_path / "run", [corrupt])
    assert result.returncode != 0
    assert "Stage 1:" in result.stderr, result.stderr
    assert "Unexpected end of gzip file" in result.stderr
