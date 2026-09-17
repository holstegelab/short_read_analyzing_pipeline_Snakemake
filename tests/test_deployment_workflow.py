import os
from pathlib import Path
import subprocess
import sys


REPO = Path(__file__).resolve().parents[1]
SNAKEFILE = REPO / "deployment" / "Snakefile"


def run_deployment(*arguments, cwd=REPO, expose_repo_on_pythonpath=True):
    environment = dict(os.environ)
    if expose_repo_on_pythonpath:
        environment["PYTHONPATH"] = str(REPO)
    else:
        environment.pop("PYTHONPATH", None)
    environment.pop("SHORT_READ_SITE_CONFIG", None)
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "--cores",
            "1",
            "--nolock",
            *arguments,
        ],
        cwd=cwd,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=90,
    )


def test_environment_dag_is_independent_of_samples():
    result = run_deployment("--dry-run", "--config", "END_POINT=gVCF",
                            "caller=Deepvariant")
    assert result.returncode == 0, result.stdout
    assert "8" in result.stdout
    for environment in (
        "preprocess",
        "align_fused",
        "qc_fused",
        "kmc",
        "kraken",
        "mosdepth",
        "capture_kit_finder",
        "vcf_handling",
    ):
        assert f"environment={environment}" in result.stdout


def test_environment_dag_loads_from_an_external_working_directory(tmp_path):
    result = run_deployment(
        "--dry-run",
        "--config",
        "END_POINT=gVCF",
        "caller=Deepvariant",
        cwd=tmp_path,
        expose_repo_on_pythonpath=False,
    )
    assert result.returncode == 0, result.stdout
    assert "environment=preprocess" in result.stdout


def test_optional_gcnv_environment_is_selected_explicitly():
    result = run_deployment(
        "--dry-run", "--config", "END_POINT=gVCF", "deployment_groups=gcnv"
    )
    assert result.returncode == 0, result.stdout
    assert "environment=gatk_gcnv" in result.stdout
    assert "environment=pca" in result.stdout


def test_unknown_environment_group_is_rejected():
    result = run_deployment(
        "--dry-run", "--config", "deployment_groups=unknown"
    )
    assert result.returncode != 0
    assert "unknown deployment_groups" in result.stdout
