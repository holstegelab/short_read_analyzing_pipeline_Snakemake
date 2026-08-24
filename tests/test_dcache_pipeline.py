import sys
import subprocess
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

import read_samples
import dcache_transfer
from dcache_transfer import (
    DirectDCache,
    adler32_file,
    load_download_list,
    normalize_remote_path,
)
from update_dcache_listing import update_listing


def test_parse_dcache_uri_modern_and_legacy():
    assert read_samples.parse_dcache_uri("dcache:mine:/cohort/a.cram") == (
        "mine",
        "/cohort/a.cram",
    )
    assert read_samples.parse_dcache_uri("dcache://cohort/a.cram") == (
        "dcache",
        "/cohort/a.cram",
    )
    assert read_samples.parse_dcache_uri("/local/file") is None


def test_parse_dcache_uri_rejects_escape():
    with pytest.raises(ValueError, match="escapes"):
        read_samples.parse_dcache_uri("dcache:mine:/cohort/../secret")


def test_parse_s3_uri_and_reject_embedded_authentication():
    assert read_samples.parse_s3_uri(
        "s3://wanglab-dss-share/distribution/adsp/cram"
    ) == ("wanglab-dss-share", "/distribution/adsp/cram")
    assert read_samples.parse_s3_uri("/local/file") is None
    with pytest.raises(ValueError, match="credentials"):
        read_samples.parse_s3_uri("s3://user:secret@bucket/key")
    with pytest.raises(ValueError, match="escapes"):
        read_samples.parse_s3_uri("s3://bucket/root/../secret")


def test_external_cram_sidecars_are_parsed(tmp_path):
    listing = tmp_path / "cohort.tsv"
    source = tmp_path / "cohort.source"
    target = tmp_path / "cohort.target"
    (tmp_path / "mine_source.conf").write_text("[mine_source]\ntype = webdav\n", encoding="utf-8")
    (tmp_path / "mine_target.conf").write_text("[mine_target]\ntype = webdav\n", encoding="utf-8")
    source.write_text("dcache:mine_source:/input\n", encoding="utf-8")
    target.write_text("dcache:mine_target:/processed/cohort\n", encoding="utf-8")
    listing.write_text(
        "mine\tmine_sample1\tcram\tillumina_wgs\t\tF\t"
        "/run/sample.cram\t/reference/hg19.fa\t"
        "filesize=42.5,cram_no_ref=true\n",
        encoding="utf-8",
    )

    sample = read_samples.read_samplefile(str(listing))[0]

    assert sample["from_external"] == "dcache"
    assert sample["source_remote"] == "mine_source"
    assert sample["source_root"] == "/input"
    assert sample["target_remote"] == "mine_target"
    assert sample["target_root"] == "/processed/cohort"
    assert sample["file1"] == ["run/sample.cram"]
    assert sample["cram_refs"] == ["/reference/hg19.fa"]
    assert sample["filesize"] == 42.5
    assert sample["cram_no_ref"] is True
    assert read_samples.append_prefix(sample["prefix"], sample["file1"][0]) == (
        "dcache:mine_source:/input/run/sample.cram"
    )


def test_external_cram_no_ref_defaults_to_false(tmp_path):
    listing = tmp_path / "cohort.tsv"
    (tmp_path / "cohort.source").write_text(
        "dcache:mine_source:/input\n", encoding="utf-8"
    )
    (tmp_path / "mine_source.conf").write_text(
        "[mine_source]\ntype = webdav\n", encoding="utf-8"
    )
    listing.write_text(
        "mine\tmine_sample1\tcram\tillumina_wgs\t\tF\t"
        "/run/sample.cram\t/reference/hg19.fa\tfilesize=42.5\n",
        encoding="utf-8",
    )

    sample = read_samples.read_samplefile(str(listing))[0]

    assert sample["cram_no_ref"] is False


def test_s3_external_cram_sidecar_is_parsed_without_credentials(tmp_path):
    listing = tmp_path / "cohort.tsv"
    (tmp_path / "cohort.source").write_text(
        "s3://wanglab-dss-share/distribution/adsp/cram\n",
        encoding="utf-8",
    )
    listing.write_text(
        "study\tstudy_sample1\tcram\tillumina_wgs\t\tF\t"
        "/snd10000/sample.cram\t/reference/hg38.fa\tfilesize=2.5\n",
        encoding="utf-8",
    )

    sample = read_samples.read_samplefile(str(listing))[0]

    assert sample["from_external"] == "s3"
    assert sample["source_bucket"] == "wanglab-dss-share"
    assert sample["source_root"] == "/distribution/adsp/cram"
    assert sample["source_remote"] is None
    assert sample["source_config"] is None
    assert sample["file1"] == ["snd10000/sample.cram"]
    assert read_samples.append_prefix(sample["prefix"], sample["file1"][0]) == (
        "s3://wanglab-dss-share/distribution/adsp/cram/"
        "snd10000/sample.cram"
    )


def test_download_list_and_adler(tmp_path):
    file_list = tmp_path / "downloads.tsv"
    file_list.write_text("/remote/a.cram\tlocal/a.cram\n", encoding="utf-8")
    assert load_download_list(file_list) == [
        ("/remote/a.cram", Path("local/a.cram"))
    ]

    payload = tmp_path / "payload"
    payload.write_bytes(b"abc")
    assert adler32_file(payload) == "024d0127"


def test_remote_path_validation():
    assert normalize_remote_path("cohort/a.cram") == "/cohort/a.cram"
    with pytest.raises(ValueError):
        normalize_remote_path("../a.cram")


def test_update_listing_sets_reference_and_numeric_gib(tmp_path):
    listing = tmp_path / "cohort.tsv"
    sizes = tmp_path / "sizes.tsv"
    reference = tmp_path / "hg19.fa"
    reference.write_text(">chr1\nA\n", encoding="utf-8")
    listing.write_text(
        "mine\tmine_sample1\tcram\tillumina_wgs\t\tF\t"
        "/run/sample.cram\told.fa\t"
        "keep=yes,cram_no_ref=true,filesize=unknown\n",
        encoding="utf-8",
    )
    sizes.write_text(f"{42 * 1024**3}\trun/sample.cram\n", encoding="utf-8")

    updated, unused = update_listing(listing, sizes, reference)

    fields = listing.read_text(encoding="utf-8").rstrip("\n").split("\t")
    assert (updated, unused) == (1, 0)
    assert fields[7] == str(reference)
    assert fields[8] == "keep=yes,cram_no_ref=true,filesize=42.000000"


def test_ada_retries_http_429(monkeypatch, tmp_path):
    client = DirectDCache.__new__(DirectDCache)
    client.ada = "/managed/ada"
    client.config_path = tmp_path / "remote.conf"
    responses = [
        subprocess.CompletedProcess([], 1, b"", b"HTTP 429 Too Many Requests"),
        subprocess.CompletedProcess([], 0, b"ok", b""),
    ]

    monkeypatch.setattr(
        dcache_transfer, "_run", lambda *args, **kwargs: responses.pop(0)
    )
    monkeypatch.setattr(dcache_transfer.time, "sleep", lambda seconds: None)
    monkeypatch.setattr(
        dcache_transfer.random, "uniform", lambda start, end: 0.0
    )

    result = client._ada("--longlist", "/input/a.cram")

    assert result.returncode == 0
    assert responses == []


def test_stage_submits_pin_before_accepting_online_batch(tmp_path):
    client = DirectDCache.__new__(DirectDCache)
    client.remote = "mine"
    client.config_path = tmp_path / "mine.conf"
    client.config_path.touch()
    calls = []
    client._ada = lambda *args, **kwargs: calls.append(args)
    client._count_online = lambda path_list, expected: expected

    client.stage(
        ["/input/a.cram"],
        poll_seconds=0,
        timeout_seconds=1,
    )

    assert calls
    assert calls[0][0:2] == ("--stage", "--from-file")


def test_download_uses_literal_list_and_shared_slot(monkeypatch, tmp_path):
    client = DirectDCache.__new__(DirectDCache)
    client.remote = "mine"
    client.config_path = tmp_path / "mine.conf"
    client.config_path.touch()
    client.dcache_cp = "/managed/dcache_cp"
    client.copy_timeout = "600m"
    client.max_retries = 3
    client.retry_wait = 60
    commands = []
    monkeypatch.setattr(
        dcache_transfer,
        "_run",
        lambda command, **kwargs: commands.append(command),
    )

    client.download(
        [("/input/a.cram", tmp_path / "downloads" / "a.cram")],
        workers=1,
        transfer_slots=1,
        no_stage=True,
        no_destage=True,
        lifetime="7D",
        poll_seconds=1,
        timeout_seconds=1,
    )

    assert len(commands) == 1
    assert "--literal-file-list" in commands[0]
    assert (
        tmp_path / ".mine.download-slots" / "slot-0.lock"
    ).is_file()
