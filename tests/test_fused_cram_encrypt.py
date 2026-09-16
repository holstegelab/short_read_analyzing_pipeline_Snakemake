import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "scripts" / "run_fused_cram_encrypt.py"
PREPROCESS_PYTHON = Path(
    "/home/hulsmanm/.snakemake/d6b82cba48e3c9dc9255474984326c51_/bin/python"
)
PREPROCESS_SAMTOOLS = PREPROCESS_PYTHON.with_name("samtools")


def _script(path: Path, body: str) -> Path:
    path.write_text("#!/usr/bin/env python3\n" + body, encoding="utf-8")
    path.chmod(0o755)
    return path


def _fixture_bam(root: Path) -> tuple[Path, Path, Path]:
    reference = root / "reference.fa"
    reference.write_text(">chr1\n" + "ACGT" * 2500 + "\n", encoding="utf-8")
    subprocess.run([PREPROCESS_SAMTOOLS, "faidx", reference], check=True)
    sequence = "ACGT" * 12 + "AC"
    quality = "I" * 50
    sam = root / "input.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:chr1\tLN:10000\n"
        "@RG\tID:rg0\tSM:S1\tLB:lib\tPL:ILLUMINA\n"
        f"primary\t0\tchr1\t101\t60\t50M\t*\t0\t0\t{sequence}\t{quality}\tRG:Z:rg0\n"
        f"secondary\t256\tchr1\t201\t20\t50M\t*\t0\t0\t{sequence}\t{quality}\tRG:Z:rg0\n"
        f"supplementary\t2048\tchr1\t301\t20\t50M\t*\t0\t0\t{sequence}\t{quality}\tRG:Z:rg0\n"
        f"unmapped\t4\t*\t0\t0\t*\t*\t0\t0\t{sequence}\t{quality}\tRG:Z:rg0\n",
        encoding="utf-8",
    )
    bam = root / "markdup.bam"
    subprocess.run([PREPROCESS_SAMTOOLS, "view", "-b", "-o", bam, sam], check=True)
    subprocess.run([PREPROCESS_SAMTOOLS, "index", bam], check=True)
    return bam, Path(str(bam) + ".bai"), reference


def _generate_keys(secret: Path, public: Path) -> None:
    subprocess.run(
        [
            PREPROCESS_PYTHON,
            "-c",
            "from crypt4gh.keys import c4gh; import sys; "
            "c4gh.generate(sys.argv[1], sys.argv[2], passphrase=None)",
            str(secret), str(public),
        ],
        check=True,
    )


def _decoded_records(cram: Path, reference: Path) -> list[str]:
    result = subprocess.run(
        [PREPROCESS_SAMTOOLS, "view", "-T", reference, cram],
        check=True, capture_output=True, text=True,
    )
    return result.stdout.splitlines()


@pytest.mark.skipif(
    not PREPROCESS_PYTHON.is_file() or not PREPROCESS_SAMTOOLS.is_file(),
    reason="preprocess environment with samtools and Crypt4GH required",
)
@pytest.mark.parametrize("failure", [None, "recipient", "crai_upload"])
def test_fused_cram_encrypt_upload_equivalence_and_failure_cleanup(tmp_path, failure):
    bam, bai, reference = _fixture_bam(tmp_path)
    baseline_cram = tmp_path / "baseline.cram"
    baseline_crai = tmp_path / "baseline.cram.crai"
    subprocess.run(
        [
            PREPROCESS_SAMTOOLS, "view", "--output-fmt", "cram,version=3.1,archive",
            "--reference", reference, "-@", "2", "--write-index", "-o",
            f"{baseline_cram}##idx##{baseline_crai}", bam,
        ],
        check=True,
    )
    secret = tmp_path / "secret.key"
    public = tmp_path / "public.key"
    _generate_keys(secret, public)
    if failure == "recipient":
        public.write_bytes(b"not a Crypt4GH public key\n")

    lease_log = tmp_path / "lease.log"
    lease = _script(
        tmp_path / "zslurm_lease",
        """import json
import os
import sys
args = sys.argv[1:]
command = args[1]
with open(os.environ['FAKE_LEASE_LOG'], 'a', encoding='utf-8') as handle:
    handle.write(command + '\\n')
if command == 'status':
    cores, memory = 1.85, 1800
else:
    cores = float(args[args.index('--cores') + 1])
    memory = float(args[args.index('--mem-mb') + 1])
print(json.dumps({'ok': True, 'held_cores': cores, 'held_mem_mb': memory}))
""",
    )
    upload_log = tmp_path / "upload.log"
    remote = tmp_path / "remote"
    remote.mkdir()
    upload = _script(
        tmp_path / "dcache_transfer.py",
        """import argparse
import os
import shutil
import zlib
from pathlib import Path
parser = argparse.ArgumentParser()
parser.add_argument('command')
parser.add_argument('--config')
parser.add_argument('--remote')
parser.add_argument('--source', type=Path)
parser.add_argument('--destination')
parser.add_argument('--checksum-output', type=Path)
parser.add_argument('--ada')
args = parser.parse_args()
with open(os.environ['FAKE_UPLOAD_LOG'], 'a', encoding='utf-8') as handle:
    handle.write(args.destination + '\\n')
if os.environ.get('FAKE_UPLOAD_FAIL_SUFFIX') and args.source.name.endswith(
    os.environ['FAKE_UPLOAD_FAIL_SUFFIX']
):
    raise SystemExit(19)
destination = Path(os.environ['FAKE_REMOTE_DIR']) / args.source.name
shutil.copyfile(args.source, destination)
checksum = 1
with args.source.open('rb') as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b''):
        checksum = zlib.adler32(chunk, checksum)
args.checksum_output.write_text(f'{checksum & 0xffffffff:08x}\\n')
""",
    )
    upload_config = tmp_path / "remote.conf"
    upload_config.write_text("[test]\ntype = webdav\n", encoding="utf-8")
    ada = _script(tmp_path / "ada", "raise SystemExit('not called by fake uploader')\n")
    outputs = tmp_path / "outputs"
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    command = [
        PREPROCESS_PYTHON, RUNNER, "--input-bam", bam, "--input-bai", bai,
        "--reference", reference, "--sample", "S1", "--private-key", secret,
        "--recipient-keys", public,
        "--upload-script", upload, "--upload-config", upload_config,
        "--upload-remote", "test", "--upload-directory", "/processed/cram",
        "--ada", ada,
        "--output-copied", outputs / "S1.mapped_hg38.cram.copied",
        "--output-checksum", outputs / "S1.mapped_hg38.cram.ADLER32",
        "--cram-log", outputs / "mCRAM.log", "--metrics", outputs / "metrics.json",
        "--samtools", PREPROCESS_SAMTOOLS, "--cram-threads", "2",
        "--initial-cores", "1.85", "--initial-memory-mb", "1800",
        "--encrypt-cores", "0.45", "--encrypt-memory-mb", "512",
        "--lease-command", lease, "--ssd-gb", "12",
        "--scratch-base", scratch, "--poll-interval", "0.05",
    ]
    environment = dict(
        os.environ,
        ZSLURM_LEASE_SOCKET="fake.socket",
        ZSLURM_LEASE_TOKEN="fake-token",
        ZSLURM_JOB_ID="123",
        FAKE_LEASE_LOG=str(lease_log),
        FAKE_UPLOAD_LOG=str(upload_log),
        FAKE_REMOTE_DIR=str(remote),
        FAKE_UPLOAD_FAIL_SUFFIX=(".crai" if failure == "crai_upload" else ""),
    )
    result = subprocess.run(list(map(str, command)), check=False, env=environment)
    metrics = json.loads((outputs / "metrics.json").read_text())
    assert metrics["success"] is (failure is None)
    assert metrics["scratch_removed"] is True
    assert list((scratch / "cram_encrypt_fused").iterdir()) == []
    assert lease_log.read_text().splitlines() == ["status", "set"]

    copied = outputs / "S1.mapped_hg38.cram.copied"
    checksum = outputs / "S1.mapped_hg38.cram.ADLER32"
    if failure is not None:
        assert result.returncode != 0
        assert not copied.exists()
        assert not checksum.exists()
        if failure == "recipient":
            assert not upload_log.exists()
        else:
            assert len(upload_log.read_text().splitlines()) == 2
            assert (remote / "S1.mapped_hg38.cram.c4gh").is_file()
        return

    assert result.returncode == 0
    assert copied.is_file() and checksum.is_file()
    assert copied.stat().st_size == 0
    assert len(upload_log.read_text().splitlines()) == 2
    encrypted = remote / "S1.mapped_hg38.cram.c4gh"
    result_crai = remote / "S1.mapped_hg38.cram.crai"
    assert encrypted.is_file() and result_crai.is_file()
    assert not list(outputs.glob("*.c4gh"))
    assert not list(outputs.glob("*.crai"))
    decrypted = tmp_path / "decrypted.cram"
    with encrypted.open("rb") as source, decrypted.open("wb") as destination:
        subprocess.run(
            [PREPROCESS_PYTHON, "-m", "crypt4gh", "decrypt", "--sk", secret],
            check=True, stdin=source, stdout=destination,
        )
    shutil.copyfile(result_crai, Path(str(decrypted) + ".crai"))
    subprocess.run([PREPROCESS_SAMTOOLS, "quickcheck", decrypted], check=True)
    baseline_records = _decoded_records(baseline_cram, reference)
    assert _decoded_records(decrypted, reference) == baseline_records
    baseline_region = subprocess.run(
        [PREPROCESS_SAMTOOLS, "view", "-T", reference, baseline_cram, "chr1:1-1000"],
        check=True, capture_output=True, text=True,
    ).stdout
    decrypted_region = subprocess.run(
        [PREPROCESS_SAMTOOLS, "view", "-T", reference, decrypted, "chr1:1-1000"],
        check=True, capture_output=True, text=True,
    ).stdout
    assert decrypted_region == baseline_region
    flags = [int(record.split("\t")[1]) for record in baseline_records]
    assert any(flag & 4 for flag in flags)
    assert any(flag & 256 for flag in flags)
    assert any(flag & 2048 for flag in flags)
    assert [phase["label"] for phase in metrics["phases"]] == [
        "cram_encrypt_fused.cram",
        "cram_encrypt_fused.encrypt",
        "cram_encrypt_fused.upload_cram",
        "cram_encrypt_fused.upload_crai",
    ]
