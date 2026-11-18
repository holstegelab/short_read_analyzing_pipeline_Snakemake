"""End-to-end checks for opt-in rescue of RG-less Illumina uCRAMs."""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import read_samples
import readgroup_rescue


SAMTOOLS = shutil.which("samtools")
pytestmark = pytest.mark.skipif(SAMTOOLS is None, reason="samtools is required")


def run_samtools(*args):
    return subprocess.run(
        [SAMTOOLS, *map(str, args)], check=True, capture_output=True, text=True,
    ).stdout


def make_ucram(tmp_path):
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nACGTACGTACGT\n", encoding="utf-8")
    run_samtools("faidx", reference)
    names = [
        "LH00206:100:22JHY7LT4:5:1101:1825:1028",
        "LH00206:100:22JHY7LT4:6:1101:1826:1029",
    ]
    sam = tmp_path / "input.sam"
    rows = ["@HD\tVN:1.6\tSO:unsorted\tGO:query\n"]
    for name in names:
        for flag in (77, 141):
            barcode = "\tBC:Z:NTATTTGTTC+TNGNCGGTTA" if flag == 77 else ""
            rows.append(
                f"{name}\t{flag}\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII"
                f"{barcode}\n"
            )
    sam.write_text("".join(rows), encoding="utf-8")
    cram = tmp_path / "input.ucram"
    run_samtools(
        "view", "-O", "cram,version=3.1,no_ref=1", "-T", reference,
        "-o", cram, sam,
    )
    return cram, reference, names


def test_opt_in_discovers_and_splits_two_lanes(tmp_path):
    cram, reference, names = make_ucram(tmp_path)
    listing = tmp_path / "cohort.tsv"
    listing.write_text(
        f"study\tstudy_S1\tcram\tillumina_wgs\t\tF\t{cram}\t"
        f"{reference}\trescue_readgroups=true,cram_no_ref=true,"
        "rg_library=LIB1\n",
        encoding="utf-8",
    )
    sample = read_samples.read_samplefile(str(listing))[0]
    assert "readgroups" not in sample  # No expensive scan during DAG building.

    resolved, warnings = read_samples.get_readgroups(sample, str(tmp_path))
    assert warnings == []
    groups = resolved["readgroups"]
    assert len(groups) == 2
    assert {rg["info"]["ID"] for rg in groups} == {
        "study_S1_LH00206_100_22JHY7LT4_L5",
        "study_S1_LH00206_100_22JHY7LT4_L6",
    }
    assert {rg["source_read_count"] for rg in groups} == {2}

    output = tmp_path / "split"
    subprocess.run([
        "python", str(Path(readgroup_rescue.__file__)), "split",
        "--input", str(cram), "--output-dir", str(output),
        "--sample", "study_S1", "--groups-json", json.dumps(groups),
        "--reference", str(reference),
        "--output-fmt", "cram,version=3.1,no_ref=1",
        "--extension", "cram",
        "--threads", "0",
    ], check=True, capture_output=True, text=True)
    assert len(list(output.glob("*.cram"))) == 2
    for group in groups:
        rgid = group["info"]["ID"]
        output_cram = output / f"study_S1.{rgid}.cram"
        header = run_samtools("view", "-H", output_cram)
        assert f"@RG\tID:{rgid}\tSM:study_S1\tLB:LIB1\tPL:ILLUMINA" in header
        rows = run_samtools("view", "--reference", reference, output_cram).splitlines()
        assert len(rows) == 2
        assert {row.split("\t")[0] for row in rows} == {
            name for name in names if f"_L{name.split(':')[3]}" in rgid
        }
        assert all(f"RG:Z:{rgid}" in row for row in rows)
        assert "BC:Z:NTATTTGTTC+TNGNCGGTTA" in rows[0]
        assert "BC:Z:" not in rows[1]


def test_invalid_qname_fails_closed():
    with pytest.raises(ValueError, match="Cannot recover flowcell/lane"):
        readgroup_rescue.qname_group("unknown-name")


def test_changed_input_count_does_not_publish_partial_crams(tmp_path):
    cram, reference, _ = make_ucram(tmp_path)
    sample = {
        "sample": "study_S1", "sample_type": "illumina_wgs", "study": "study",
        "file1": [str(cram)], "file2": [], "file_type": "cram",
        "cram_refs": [str(reference)], "alt_name": set(),
        "rescue_readgroups": True,
    }
    groups = read_samples.get_readgroups(sample, str(tmp_path))[0]["readgroups"]
    groups[0]["source_read_count"] += 1
    output = tmp_path / "split"
    with pytest.raises(ValueError, match="Read count changed"):
        readgroup_rescue.split_rescued_cram(
            str(cram), str(output), "study_S1", groups, str(reference),
            "cram,version=3.1,no_ref=1", "cram", samtools=SAMTOOLS,
        )
    assert list(output.iterdir()) == []


def test_rgless_cram_requires_opt_in(tmp_path):
    cram, reference, _ = make_ucram(tmp_path)
    sample = {
        "sample": "study_S1", "sample_type": "illumina_wgs", "study": "study",
        "file1": [str(cram)], "file2": [], "file_type": "cram",
        "cram_refs": [str(reference)], "alt_name": set(),
    }
    with pytest.raises(ValueError, match="rescue_readgroups=true"):
        read_samples.get_readgroups(sample, str(tmp_path))


def test_pipeline_checkpoint_and_split_rule(tmp_path):
    pytest.importorskip("snakemake")
    cram, reference, _ = make_ucram(tmp_path)
    (tmp_path / "cohort.tsv").write_text(
        f"study\tstudy_S1\tcram\tillumina_wgs\tWGS\tM\t{cram.name}\t"
        f"{reference}\trescue_readgroups=true,cram_no_ref=true\n",
        encoding="utf-8",
    )
    (tmp_path / "source").mkdir()
    (tmp_path / "source" / "study_S1.started").touch()
    repo = Path(__file__).resolve().parents[1]
    environment = dict(os.environ, PYTHONPATH=str(repo))
    target = "readgroups/study_S1.sourcefile.input.checks_done"
    completed = subprocess.run(
        [sys.executable, "-m", "snakemake", "--snakefile", str(repo / "Snakefile"),
         "--cores", "3", "--nolock", target, "--config", "END_POINT=Align", "chrM=No"],
        cwd=tmp_path, env=environment, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90,
    )
    assert completed.returncode == 0, completed.stdout
    assert (tmp_path / target).is_file()
    assert len(list((tmp_path / "readgroups" / "study_S1.sourcefile.input").glob("*.cram"))) == 2

    dry_run = subprocess.run(
        [sys.executable, "-m", "snakemake", "--snakefile", str(repo / "Snakefile"),
         "--cores", "3", "--nolock", "--dry-run", "bams/study_S1.markdup.bam",
         "--config", "END_POINT=Align", "chrM=No"],
        cwd=tmp_path, env=environment, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=90,
    )
    assert dry_run.returncode == 0, dry_run.stdout
    assert "rule external_adapter_fused:" in dry_run.stdout
