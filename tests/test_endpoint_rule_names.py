"""Every documented endpoint must parse with strict duplicate-rule checks."""

import os
from pathlib import Path
import subprocess
import sys

import pytest


REPO = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    "endpoint,caller,combine",
    [
        (endpoint, caller, combine)
        for endpoint in ("Genotype", "Combine", "VCF")
        for caller in ("Deepvariant", "HaplotypeCaller", "BOTH")
        for combine in ("GLnexus", "COMBINE_GVCF", "DBIMPORT")
    ]
    + [
        ("gVCF", caller, "GLnexus")
        for caller in ("Deepvariant", "HaplotypeCaller", "BOTH")
    ]
    + [("Align", "Deepvariant", "GLnexus"),
       ("PrepareRef", "Deepvariant", "GLnexus")],
)
def test_endpoint_rule_names_are_unambiguous(
    tmp_path, endpoint, caller, combine
):
    environment = dict(os.environ, PYTHONPATH=str(REPO))
    environment.pop("SHORT_READ_RESTART_MANIFEST", None)
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--snakefile",
            str(REPO / "Snakefile"),
            "--list-rules",
            "--cores",
            "1",
            "--config",
            f"END_POINT={endpoint}",
            f"caller={caller}",
            f"Combine_gVCF_method={combine}",
            "chrM=No",
        ],
        cwd=tmp_path,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=90,
    )
    assert result.returncode == 0, result.stdout
    assert "already used by another rule" not in result.stdout
