import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


REPO = Path(__file__).resolve().parents[1]


def run_probe(site_file, cwd):
    code = (
        "import json, constants, os; "
        "print(json.dumps({'resources': constants.RESOURCES, "
        "'software': constants.SOFTWARE, 'tmp': constants.TMPDIR, "
        "'hg19': constants.HG19_REFERENCE, "
        "'site': constants.SITE_NAME, "
        "'site_file': os.environ.get('SHORT_READ_SITE_CONFIG')}))"
    )
    env = dict(os.environ, PYTHONPATH=str(REPO), SHORT_READ_SITE_CONFIG=str(site_file))
    result = subprocess.run(
        [sys.executable, "-c", code], cwd=cwd, env=env,
        text=True, capture_output=True, timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return json.loads(result.stdout)


def run_settings_probe(site_file, cwd):
    code = (
        "import json; from site_config import configure; "
        "settings = configure(); "
        "print(json.dumps(dict(settings.values)))"
    )
    env = dict(os.environ, PYTHONPATH=str(REPO), SHORT_READ_SITE_CONFIG=str(site_file))
    result = subprocess.run(
        [sys.executable, "-c", code], cwd=cwd, env=env,
        text=True, capture_output=True, timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return json.loads(result.stdout)


def test_site_config_reaches_constants_from_an_unrelated_workdir(tmp_path):
    root = tmp_path / "resources"
    shared = tmp_path / "shared tmp"
    site = tmp_path / "site.yaml"
    site.write_text(
        "schema_version: 1\n"
        "site: test-cluster\n"
        "paths:\n"
        f"  resource_root: {root}\n"
        f"  shared_tmp_root: {shared}\n"
        f"  hg19_b37_chry_reference: {root}/cram_refs/hg19_b37chrY.fa\n"
        "storage:\n"
        f"  dcache_processed_config: {tmp_path}/processed.conf\n"
        f"  dcache_read_config: {tmp_path}/read.conf\n"
        "encryption:\n"
        f"  sender_private_key: {tmp_path}/sender.key\n"
        "  recipient_public_keys:\n"
        f"    - {tmp_path}/recipient.pub\n"
        f"  decryption_private_key: {tmp_path}/recipient.key\n"
    )
    probe = run_probe(site, tmp_path)
    assert probe == {
        "resources": str(root),
        "software": str(root / "software"),
        "tmp": str(shared),
        "hg19": str(root / "cram_refs" / "hg19.fa"),
        "site": "test-cluster",
        "site_file": str(site.resolve()),
    }


def test_spider_example_contains_no_snellius_paths(tmp_path):
    probe = run_settings_probe(REPO / "config/sites/spider.example.yaml", tmp_path)
    assert probe["site_name"] == "spider"
    assert "/gpfs/" not in json.dumps(probe)


@pytest.mark.parametrize(
    "content,needle",
    [
        ("schema_version: 2\npaths:\n  resource_root: /x\n", "schema_version"),
        ("paths:\n  resource_root: /x\n  resorce_root: /typo\n", "resorce_root"),
        ("paths:\n  resource_root: $UNSET_FOR_SITE_TEST/x\n", "unresolved"),
    ],
)
def test_invalid_site_config_fails_before_constants_are_derived(tmp_path, content, needle):
    site = tmp_path / "bad.yaml"
    site.write_text(content)
    env = dict(os.environ, PYTHONPATH=str(REPO), SHORT_READ_SITE_CONFIG=str(site))
    env.pop("UNSET_FOR_SITE_TEST", None)
    result = subprocess.run(
        [sys.executable, "-c", "import constants"], cwd=tmp_path, env=env,
        text=True, capture_output=True, timeout=30,
    )
    assert result.returncode != 0
    assert needle in result.stderr
