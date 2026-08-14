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
