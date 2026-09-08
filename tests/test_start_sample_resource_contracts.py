import re
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def _start_sample_cores_function():
    aligner = (REPO / "Aligner.smk").read_text()
    match = re.search(
        r"def _start_sample_cores\(wildcards\):\n(?P<body>.*?)(?=\n\ndef )",
        aligner,
        flags=re.DOTALL,
    )
    assert match is not None

    namespace = {}
    source = "def _start_sample_cores(wildcards):\n" + match.group("body")
    exec(source, namespace)
    return namespace["_start_sample_cores"], namespace


def test_start_sample_cores_preserve_fractional_zslurm_values():
    start_sample_cores, namespace = _start_sample_cores_function()
    namespace["_start_sample_route"] = lambda wildcards: wildcards["route"]

    expected = {
        "active": "0.1",
        "archive": "0.6",
        "dcache": "1.0",
        "s3": "0.5",
    }
    actual = {
        route: start_sample_cores({"route": route}) for route in expected
    }
    assert actual == expected
    assert all(isinstance(value, str) for value in actual.values())
