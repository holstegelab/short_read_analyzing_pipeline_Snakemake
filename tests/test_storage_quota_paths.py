"""Exercise the actual large-writer entry points using tiny mocked data."""
import json
import os
import shlex
import subprocess
import sys
import textwrap
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_fastq_adapter as native
import storage_quota


def test_shared_cli_passes_destinations_and_full_work_profile(monkeypatch):
    calls = []
    monkeypatch.setattr(storage_quota, "wait_for_storage_quota",
                        lambda *args, **kwargs: calls.append((args, kwargs)))
    assert storage_quota.main([
        "--destination", "/data/sample.bam", "/data/sample.bam.bai",
        "--cores", "8", "--memory-mb", "14000", "--quota-mode", "disabled",
        "--lease-command", "/pipeline/lease", "--phase", "restore_bam",
    ]) == 0
    assert calls == [(([Path("/data/sample.bam"), Path("/data/sample.bam.bai")],), dict(
        mode="disabled", resume_cores=8.0, resume_memory_mb=14000.0,
        lease_command="/pipeline/lease", resume_phase="restore_bam",
    ))]


def test_pre_upgrade_cli_uses_pipeline_local_lease_client(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(storage_quota, "wait_for_storage_quota",
                        lambda *args, **kwargs: calls.append(kwargs))
    storage_quota.main(["--destination", str(tmp_path / "out"), "--cores", "1", "--memory-mb", "512"])
    expected = str((REPO / "scripts/zslurm_lease_client.py").resolve())
    assert calls[0]["lease_command"] == expected
    # The live controller may already have expanded a native FASTQ command
    # without the new switches; its freshly started runner still gets them.
    fields = (
        "input-forward", "input-reverse", "sample", "readgroup", "adapter-list",
        "fastq-stats-script", "remove-duplicates-script", "pair-rescue-script",
        "error-file", "output-forward", "output-reverse", "output-adapter-log",
        "output-fastq-stats", "output-adapters", "metrics",
    )
    monkeypatch.setattr(sys, "argv", ["run_fastq_adapter.py"] + [
        token for field in fields for token in ("--" + field, "unused")
    ])
    args = native.parse_args()
    assert args.lease_command == expected
    assert args.quota_mode == "auto"
    assert args.initial_cores == 5.0


def test_native_fastq_waits_before_shared_storage_work(tmp_path, monkeypatch):
    input_file = tmp_path / "input"
    input_file.write_bytes(b"reads")
    output_dir = tmp_path / "fq"
    args = SimpleNamespace(**{name: str(input_file) for name in (
        "input_forward", "input_reverse", "adapter_list", "fastq_stats_script",
        "remove_duplicates_script", "pair_rescue_script",
    )}, **{name: str(output_dir / name) for name in (
        "output_forward", "output_reverse", "output_adapter_log", "output_fastq_stats",
        "output_adapters",
    )}, metrics=str(tmp_path / "metrics.json"), pigz="pigz", adapter_removal="AdapterRemoval",
        quota_mode="required", initial_cores=5.0, adapter_memory_mb=512.0,
        lease_command="/pipeline/lease", sample="S1", readgroup="RG1", attempt=1)
    events = []

    def wait(paths, **kwargs):
        assert not output_dir.exists()
        assert paths == [Path(args.output_forward), Path(args.output_reverse)]
        assert kwargs == dict(mode="required", resume_cores=5.0, resume_memory_mb=512.0,
                              lease_command="/pipeline/lease", resume_phase="native_fastq_prepare")
        events.append("wait")
        return {"enabled": True, "checks": 4}

    def prepare(_args, raw1, raw2, job_tmp, phases, **kwargs):
        assert events == ["wait"]
        assert job_tmp.parent == output_dir
        events.append("write")
        prepared = [job_tmp / str(i) for i in range(5)]
        for path in prepared:
            path.write_bytes(b"tiny result")
        return prepared, {"r1_base": 33}

    monkeypatch.setattr(native, "parse_args", lambda: args)
    monkeypatch.setattr(native, "executable", lambda value: value)
    monkeypatch.setattr(native, "wait_for_storage_quota", wait)
    monkeypatch.setattr(native, "prepare_adapters", prepare)
    assert native.main() == 0
    assert events == ["wait", "write"]
    assert Path(args.output_forward).read_bytes() == b"tiny result"
    assert json.loads(Path(args.metrics).read_text())["quota_wait"]["checks"] == 4
    assert not list(output_dir.glob(".adapter-*"))


class Wildcards(dict):
    def __getattr__(self, key):
        return self[key]


@pytest.mark.parametrize("groups,erf", [(1, False), (1, True), (2, False), (2, True)])
@pytest.mark.parametrize("kind", ["bam", "cram"])
def test_split_run_gates_real_writes_but_not_single_rg_hardlinks(groups, erf, kind):
    # Compile the real rule's Python run body; shell calls are recorded, never
    # executed. This tests its branches without touching a production dataset.
    text = (REPO / "Aligner.smk").read_text()
    block = text.split("rule split_alignments_by_readgroup:", 1)[1].split(
        "\ndef get_aligned_readgroup_folder", 1
    )[0]
    body = textwrap.dedent(block.split("    run:\n", 1)[1])
    readfile = f"/data/input.{kind}"
    sample = {"erf_correct": erf, "readgroups": [
        {"file": readfile, "file_type": kind, "reference_file": "/refs/genome.fa",
         "info": {"ID": f"rg{i}"}} for i in range(groups)
    ]}
    calls = []

    def shell(command, **kwargs):
        calls.append(command)
        return "/refs/genome.fa" if kwargs.get("read") else None

    namespace = dict(
        sampleinfo=lambda *_a, **_k: sample, SAMPLEINFO={"S1": sample},
        wildcards=Wildcards(sample="S1", filename="input"), os=os, sys=sys, shlex=shlex,
        output=SimpleNamespace(readgroups="/data/readgroups", done="/data/checks_done"),
        params=SimpleNamespace(fixer="/pipeline/fixer", reference_selector="/pipeline/select_ref",
                               quota_gate="/pipeline/storage_quota.py", quota_mode="auto",
                               lease_command="/pipeline/lease"),
        resources=SimpleNamespace(n=3, mem_mb=7500, use_threads=3), shell=shell,
        HG19_REFERENCE="hg19", HG19_B37_CHRY_REFERENCE="hg19_chry", HG38_CRAM_REFERENCE="hg38",
    )
    exec("def run_case():\n" + textwrap.indent(body, "    "), namespace)
    namespace["run_case"]()
    # Reference detection is read-only and may precede the check for CRAM.
    calls = [call for call in calls if "--primary-reference" not in call]
    if groups == 1 and not erf:
        assert len(calls) == 1
        assert f"ln {readfile}" in calls[0]
        assert "quota_gate" not in calls[0]
    else:
        assert "python {params.quota_gate:q}" in calls[0]
        assert "--cores {resources.n}" in calls[0]
        assert "--memory-mb {resources.mem_mb}" in calls[0]
        assert "--phase split_alignments" in calls[0]
        assert "mkdir -p" in calls[1]
        assert "samtools split" in calls[1] if groups > 1 else "/pipeline/fixer" in calls[1]


@pytest.mark.parametrize("rule,writer,phase,cores,memory", [
    ("download_encrypted_cram_from_tape", "{params.rclone:q}", "restore_download", 'n="0.2"', "mem_mb=2000"),
    ("decrypt_cram", "python {params.script}", "restore_decrypt", 'n="0.5"', "mem_mb=1000"),
    ("cram_to_markdup_bam", "samtools view", "restore_bam", "n=8", "mem_mb=14000"),
])
def test_restore_gates_precede_writers_and_do_not_overwrite_quota_logs(rule, writer, phase, cores, memory):
    text = (REPO / "Download_and_extract_bam.smk").read_text()
    block = text.split(f"rule {rule}:", 1)[1].split("\nrule ", 1)[0]
    shell = block.split("    shell:", 1)[1]
    assert shell.index(": > {log:q}") < shell.index("python {params.quota_gate:q}") < shell.index(writer)
    assert f"--phase {phase}" in shell
    assert "--cores {resources.n} --memory-mb {resources.mem_mb}" in shell
    assert cores in block and memory in block
    after_gate = shell.split(f"--phase {phase}", 1)[1]
    assert " > {log" not in after_gate
    assert " 2> {log" not in after_gate


@pytest.mark.parametrize("mode", ["auto", "disabled"])
def test_actual_restoration_rules_parse_and_format_in_snakemake(tmp_path, mode):
    # Supply only common's site values, then include the real restore rules.
    # No real data, credentials, downloads or conversions are used.
    snakefile = tmp_path / "Snakefile"
    snakefile.write_text(textwrap.dedent(f'''\
        import os
        import sys
        import types
        fake = types.ModuleType("common")
        values = dict(
            pj=os.path.join, sample_names=["TEST_A"], BAM="bams", CRAM="cram", LOG="logs",
            DCACHE_READ_CONFIG="test.conf", DCACHE_READ_PREFIX="test:/processed",
            DECRYPTION_PRIVATE_KEY="test.key", DECRYPTION_PASSPHRASE_FILE="",
            RCLONE="rclone", REF="genome.fa", CONDA_MAIN={str(REPO / 'envs/preprocess.yaml')!r},
            remote_base_for_sample=lambda sample: "cohort",
            zslurm_lease_command=lambda cfg: "/pipeline/lease",
            srcdir=lambda path: os.path.join({str(REPO)!r}, path),
        )
        for key, value in values.items():
            setattr(fake, key, value)
        sys.modules["common"] = fake
        include: {str(REPO / 'Download_and_extract_bam.smk')!r}
    '''))
    result = subprocess.run(
        [sys.executable, "-m", "snakemake", "--snakefile", str(snakefile),
         "--executor", "local", "--cores", "1", "--dry-run", "--printshellcmds",
         "bams/TEST_A.markdup.bam",
         "--config", f"storage_quota_mode={mode}"],
        cwd=tmp_path, capture_output=True, text=True, timeout=60,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert output.count("--quota-mode " + mode) == 3
    for phase in ("restore_download", "restore_decrypt", "restore_bam"):
        assert "--phase " + phase in output
    assert str((REPO / "scripts/storage_quota.py").resolve()) in output
    assert "--cores 8 --memory-mb 14000" in output
    assert not (tmp_path / "bams/TEST_A.markdup.bam").exists()
