"""Small real-Snakemake DAG/worker checks; never use a production run directory."""
import gzip
import os
import re
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


REPO = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(
    not Path("/gpfs/work3/0/qtholstg/hg38_res_v2/databases/Adapters_illumina.txt").is_file(),
    reason="site integration smoke tests require the currently configured reference bundle",
)


def command(root, *args, repo=REPO):
    env = dict(os.environ, PYTHONPATH=str(repo), SKIP_FASTQ_VALIDATION="1")
    env.pop("SHORT_READ_RESTART_MANIFEST", None)
    result = subprocess.run(
        [sys.executable, "-m", "snakemake", "--snakefile", str(repo / "Snakefile"),
         "--cores", "1", "--nolock", *args],
        cwd=root, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        timeout=90,
    )
    assert result.returncode == 0, result.stdout
    return result.stdout


def fastq_fixture(root, *, groups=1, sex="M", wgs=True):
    sample = "TEST_A"
    first, second = [], []
    for group in range(groups):
        for mate, files in ((1, first), (2, second)):
            name = f"lane{group}_R{mate}.fq.gz"
            with gzip.open(root / name, "wt") as handle:
                handle.write(f"@pair/{mate}\n{'ACGT' * 15}\n+\n{'!' + 'I' * 59}\n")
            files.append(name)
    fields = ["TEST", sample, "fastq_paired", "illumina_wgs" if wgs else "illumina_exome",
              "WGS" if wgs else "Agilent_V5", sex, ",".join(first), ",".join(second)]
    (root / "cohort.tsv").write_text("\t".join(fields) + "\n")
    return sample


@pytest.mark.parametrize("groups,sex,wgs", [(1, "M", True), (2, "F", False)])
def test_native_fastq_alignment_dag_has_only_fused_producers(tmp_path, groups, sex, wgs):
    sample = fastq_fixture(tmp_path, groups=groups, sex=sex, wgs=wgs)
    output = command(tmp_path, "--dry-run", f"bams/{sample}.markdup.bam",
                     "--config", "END_POINT=Align", "chrM=No")
    for name in ("adapter_removal", "align_reads_fused", "kmer_sex_fused", "markdup"):
        assert f"rule {name}:" in output
    assert ("rule merge_rgs:" in output) == (groups > 1)
    assert "rule external_adapter_fused:" not in output
    assert "rule align_reads:" not in output
    assert "rule kmer_reads:" not in output


@pytest.mark.skipif(not shutil.which("AdapterRemoval") or not shutil.which("pigz"),
                    reason="real AdapterRemoval and pigz required")
def test_native_adapter_worker_reparses_root_workflow_and_runs(tmp_path):
    sample = fastq_fixture(tmp_path)
    source = tmp_path / "source"
    source.mkdir()
    (source / f"{sample}.started").touch()
    target = f"fq/{sample}.{sample}_rg0.fastq.cut_1.fq.gz"
    command(tmp_path, target, "--allowed-rules", "adapter_removal",
            "--mode", "subprocess", "--force-use-threads", "--target-files-omit-workdir-adjustment",
            "--config", "END_POINT=Align", "chrM=No")
    assert gzip.open(tmp_path / target, "rt").read().startswith("@pair/1")
    assert (tmp_path / f"stats/{sample}.{sample}_rg0.fastq.stats.tsv").is_file()


@pytest.mark.parametrize("caller,chrm", [("Deepvariant", "Yes"), ("HaplotypeCaller", "No"), ("BOTH", "Yes")])
def test_root_workflow_lists_only_supported_rule_implementations(tmp_path, caller, chrm):
    fastq_fixture(tmp_path)
    output = command(tmp_path, "--list-rules", "--config", f"caller={caller}", f"chrM={chrm}")
    assert "align_reads_fused" in output
    assert "bam_qc_fused" in output
    assert "\nalign_reads\n" not in output
    assert "\nexternal_alignments_to_fastq\n" not in output
    if caller in {"Deepvariant", "BOTH"}:
        assert "deepvariant_phasing_fused" in output
        assert "DeepVariant_apptainer_all" in output


def alignment_fixture(root, kind, groups, erf=False):
    sample = "TEST_A"
    reference = root / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 200 + "\n")
    subprocess.run(["samtools", "faidx", str(reference)], check=True, capture_output=True)
    sam = root / "input.sam"
    sam.write_text("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:200\n" + "".join(
        f"@RG\tID:rg{i}\tSM:{sample}\tLB:lib\tPL:ILLUMINA\n" for i in range(groups)
    ) + "".join(
        f"read{i}\t{flag}\t*\t0\t0\t*\t*\t0\t0\t{'ACGT' * 15}\t{'!' + 'I' * 59}\tRG:Z:rg{i}\n"
        for i in range(groups) for flag in (77, 141)
    ))
    alignment = root / f"input.{kind}"
    subprocess.run(["samtools", "view", "-C" if kind == "cram" else "-b", "-T", str(reference),
                    "-o", str(alignment), str(sam)], check=True, capture_output=True)
    fields = ["TEST", sample, kind, "illumina_wgs", "WGS", "M", alignment.name,
              str(reference) if kind == "cram" else "", f"erf_correct={int(erf)}"]
    (root / "cohort.tsv").write_text("\t".join(fields) + "\n")
    return sample


@pytest.mark.skipif(not shutil.which("samtools"), reason="samtools required")
@pytest.mark.parametrize("kind,groups,erf", [("bam", 1, False), ("bam", 2, True), ("cram", 1, False), ("cram", 2, False)])
def test_alignment_input_dag_uses_external_fusion(tmp_path, kind, groups, erf):
    sample = alignment_fixture(tmp_path, kind, groups, erf)
    output = command(tmp_path, "--dry-run", f"bams/{sample}.markdup.bam", "--config", "END_POINT=Align", "chrM=No")
    assert "rule external_adapter_fused:" in output
    assert "rule adapter_removal:" not in output
    assert "rule split_alignments_by_readgroup:" in output
    assert "rule align_reads_fused:" in output


@pytest.mark.skipif(not os.environ.get("FUSED_CLEANUP_BASELINE"), reason="optional frozen pre-cleanup checkout")
@pytest.mark.parametrize("kind,groups", [("fastq", 1), ("fastq", 2), ("bam", 2), ("cram", 2)])
def test_complete_stage_dag_outputs_match_baseline(tmp_path, kind, groups):
    baseline = Path(os.environ["FUSED_CLEANUP_BASELINE"])
    if kind != "fastq" and not shutil.which("samtools"):
        pytest.skip("samtools required")
    sample = fastq_fixture(tmp_path, groups=groups) if kind == "fastq" else alignment_fixture(tmp_path, kind, groups)
    targets = [f"stats/{sample}.hs_metrics",
               f"deepvariant/gVCF/A0/{sample}.A0.wg.vcf.gz",
               f"chrM_analysis/variants/gvcf/{sample}.chrM_merged_BP_annotated.g.vcf.gz"]
    def rows(text):
        # Ignore dates/log locations; compare every DAG output and its producer.
        return sorted((fields[0], fields[2]) for line in text.splitlines()
                      if len(fields := line.split("\t")) >= 6 and fields[0] != "output_file")
    before = command(tmp_path, *targets, "--summary", repo=baseline)
    after = command(tmp_path, *targets, "--summary")
    assert rows(before)
    assert rows(after) == rows(before)


@pytest.mark.skipif(not os.environ.get("FUSED_CLEANUP_BASELINE") or not shutil.which("AdapterRemoval"),
                    reason="frozen baseline and real AdapterRemoval required")
@pytest.mark.parametrize("quality,dedup,attempt", [("!" + "I" * 59, 0, 1), ("h" * 60, 0, 1),
                                                 ("!" + "I" * 59, 1, 1), ("!" + "I" * 59, 0, 2)])
def test_native_adapter_outputs_match_original_inline_algorithm(tmp_path, quality, dedup, attempt):
    baseline = Path(os.environ["FUSED_CLEANUP_BASELINE"])
    observed = []
    for label, repo in (("before", baseline), ("after", REPO)):
        root = tmp_path / label
        root.mkdir()
        sample = fastq_fixture(root)
        # Real duplicate sequences and out-of-order pairs on the rescue route.
        for mate in (1, 2):
            records = [(f"pair{i}", "ACGT" * 15 if i < 2 else "TGCA" * 15) for i in range(3)]
            if mate == 2 and attempt == 2:
                records.reverse()
            with gzip.open(root / f"lane0_R{mate}.fq.gz", "wt") as handle:
                for name, sequence in records:
                    handle.write(f"@{name}/{mate}\n{sequence}\n+\n{quality}\n")
        listing = root / "cohort.tsv"
        listing.write_text(listing.read_text().rstrip() + f"\tremove_duplicated_reads={dedup}\n")
        (root / "source").mkdir()
        (root / "source" / f"{sample}.started").touch()
        prefix = f"{sample}.{sample}_rg0"
        command(root, f"fq/{prefix}.fastq.cut_1.fq.gz", "--allowed-rules", "adapter_removal",
                "--mode", "subprocess", "--force-use-threads", "--attempt", str(attempt),
                "--config", "END_POINT=Align", "chrM=No", repo=repo)
        def stable_settings(path):
            # Settings include executable/input/output filenames and elapsed time.
            lines = path.read_text().splitlines()
            return [line for line in lines if not re.search(r"/|time|date|command|started|finished", line, re.I)]
        observed.append((
            *(gzip.open(root / f"fq/{prefix}.fastq.cut_{mate}.fq.gz", "rt").read() for mate in (1, 2)),
            (root / f"stats/{prefix}.fastq.stats.tsv").read_text(),
            (root / f"stats/{prefix}.fastq.adapters").read_text(),
            stable_settings(root / f"stats/{prefix}.adapter_removal.log"),
            (root / "cohort.errors").read_text() if attempt >= 2 else "",
        ))
    assert observed[0] == observed[1]


@pytest.mark.skipif(not os.environ.get("FUSED_CLEANUP_BASELINE") or not shutil.which("AdapterRemoval")
                    or not shutil.which("samtools"), reason="frozen baseline and real tools required")
@pytest.mark.parametrize("kind", ["bam", "cram"])
def test_external_adapter_outputs_match_original_fused_runner(tmp_path, kind):
    baseline = Path(os.environ["FUSED_CLEANUP_BASELINE"])
    alignment_fixture(tmp_path, kind, 1)
    outputs = []
    for label, repo in (("before", baseline), ("after", REPO)):
        root = tmp_path / label
        root.mkdir()
        scratch = root / "scratch"
        scratch.mkdir()
        args = [sys.executable, str(repo / "scripts/run_fused_external_adapter.py"),
                "--input-alignment", str(tmp_path / f"input.{kind}"),
                "--cram-options", f"--reference {tmp_path / 'reference.fa'}" if kind == "cram" else "",
                "--sample", "TEST_A", "--readgroup", "rg0",
                "--adapter-list", "/gpfs/work3/0/qtholstg/hg38_res_v2/databases/Adapters_illumina.txt",
                "--fastq-stats-script", str(repo / "scripts/fastq_stats.py"),
                "--remove-duplicates-script", str(repo / "scripts/remove_interleaved_duplicates.py"),
                "--pair-rescue-script", str(repo / "scripts/fastq_pair_rescue.py"),
                "--error-file", str(root / "errors"), "--metrics", str(root / "metrics.json"),
                "--initial-cores", "5", "--initial-memory-mb", "14250",
                "--adapter-cores", "5", "--adapter-memory-mb", "768",
                "--lease-mode", "disabled", "--ssd-gb", "32",
                "--scratch-base", str(scratch), "--poll-interval", "0.05"]
        products = {"raw-forward": "raw1.fq.gz", "raw-reverse": "raw2.fq.gz",
                    "singletons": "singletons.fq.gz", "forward": "cut1.fq.gz",
                    "reverse": "cut2.fq.gz", "adapter-log": "adapter.log",
                    "fastq-stats": "stats.tsv", "adapters": "adapters.txt"}
        for option, name in products.items():
            args += ["--output-" + option, str(root / name)]
        result = subprocess.run(args, text=True, capture_output=True, timeout=60)
        assert result.returncode == 0, result.stdout + result.stderr
        outputs.append({name: gzip.open(root / name, "rb").read() if name.endswith(".gz")
                        else (root / name).read_bytes() for name in products.values() if name != "adapter.log"})
        assert list((scratch / "external_adapter_fused").iterdir()) == []
    assert outputs[0] == outputs[1]
