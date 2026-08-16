import json
import sys
from pathlib import Path

import pytest


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from start_sample_route import (  # noqa: E402
    StartSampleRouteError,
    expected_relative_files,
    make_partial_directory,
    promote_directory,
    quarantine_incomplete,
    safe_relative_path,
    sample_route,
    validate_materialized,
    write_completion_markers,
    write_manifest,
)
from start_sample_job import run_start_sample_job  # noqa: E402


def test_route_and_archive_conversion_are_explicit():
    sample = {"from_external": "archive", "file1": ["lane/R1.fq.bz2"], "file2": ["lane/R2.fq.bz2"]}
    assert sample_route(sample) == "archive"
    assert expected_relative_files(sample, "archive") == [
        Path("lane/R1.fq.gz"), Path("lane/R2.fq.gz")
    ]
    assert sample_route({"from_external": False}) == "active"
    with pytest.raises(StartSampleRouteError, match="unsafe"):
        safe_relative_path("../escape.fastq.gz")


def test_materialization_promote_validate_and_markers(tmp_path):
    sample = {"file1": ["R1.fastq.gz"], "file2": ["R2.fastq.gz"]}
    destination = tmp_path / "sample.data"
    partial = make_partial_directory(destination)
    (partial / "R1.fastq.gz").write_bytes(b"read1")
    (partial / "R2.fastq.gz").write_bytes(b"read2")
    records = validate_materialized(partial, sample, "dcache")
    write_manifest(partial, sample_name="S1", route="dcache", records=records)
    promote_directory(partial, destination)
    assert validate_materialized(destination, sample, "dcache") == records

    started = tmp_path / "source" / "S1.started"
    route_ready = tmp_path / "source" / "S1.route_ready"
    legacy = tmp_path / "source" / "S1.dcache_retrieved"
    write_completion_markers(
        started=started,
        route_ready=route_ready,
        legacy_marker=legacy,
        route="dcache",
        sample_name="S1",
        files=records,
    )
    assert json.loads(started.read_text())["route"] == "dcache"
    assert json.loads(route_ready.read_text())["files"] == records
    assert legacy.is_file()

    (destination / "R1.fastq.gz").write_bytes(b"changed")
    with pytest.raises(StartSampleRouteError, match="size changed"):
        validate_materialized(destination, sample, "dcache")
    quarantine = quarantine_incomplete(destination)
    assert quarantine is not None and quarantine.is_dir()
    assert not destination.exists()


def test_active_start_job_validates_without_owned_temp_directory(tmp_path):
    source = tmp_path / "input.cram"
    source.write_bytes(b"cram")
    started = tmp_path / "source" / "S1.started"
    route_ready = tmp_path / "source" / "S1.route_ready"

    run_start_sample_job(
        sample={
            "from_external": False,
            "prefix": str(tmp_path),
            "file1": [source.name],
            "file2": [],
        },
        sample_name="S1",
        expected_route="active",
        destination=None,
        started=started,
        route_ready=route_ready,
        append_prefix=lambda prefix, filename: str(Path(prefix) / filename),
        dcache_source_file=lambda _sample, _filename: (_ for _ in ()).throw(
            AssertionError("active route must not resolve dCache input")
        ),
        transfer_script=tmp_path / "unused.py",
        source_dir=tmp_path / "source",
        dcache_download_workers=1,
        dcache_download_lock_slots=1,
    )

    assert json.loads(started.read_text())["route"] == "active"
    assert json.loads(route_ready.read_text())["files"] == [
        {"bytes": 4, "path": str(source)}
    ]


def test_dcache_start_job_adopts_declared_materialization(tmp_path):
    sample = {
        "from_external": "dcache",
        "prefix": "dcache:remote:/cohort",
        "source_remote": "remote",
        "source_config": str(tmp_path / "rclone.conf"),
        "file1": ["sample.cram"],
        "file2": [],
    }
    destination = tmp_path / "source" / "S1.dcache_data"
    destination.mkdir(parents=True)
    (destination / "sample.cram").write_bytes(b"cram")
    records = validate_materialized(destination, sample, "dcache")
    write_manifest(
        destination,
        sample_name="S1",
        route="dcache",
        records=records,
    )

    started = tmp_path / "source" / "S1.started"
    route_ready = tmp_path / "source" / "S1.route_ready"
    run_start_sample_job(
        sample=sample,
        sample_name="S1",
        expected_route="dcache",
        destination=destination,
        started=started,
        route_ready=route_ready,
        append_prefix=lambda prefix, filename: f"{prefix}/{filename}",
        dcache_source_file=lambda _sample, _filename: (_ for _ in ()).throw(
            AssertionError("validated materialization must not download again")
        ),
        transfer_script=tmp_path / "unused.py",
        source_dir=tmp_path / "source",
        dcache_download_workers=2,
        dcache_download_lock_slots=4,
    )

    assert json.loads(route_ready.read_text())["route"] == "dcache"
    assert validate_materialized(destination, sample, "dcache") == records
