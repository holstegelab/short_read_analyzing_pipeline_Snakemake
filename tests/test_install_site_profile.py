import subprocess
import sys
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts" / "install_site_profile.py"


def test_installed_profile_binds_site_and_pipeline(tmp_path):
    software = tmp_path / "software"
    site = tmp_path / "site.yaml"
    site.write_text(
        "\n".join(
            [
                "schema_version: 1",
                "site: test-spider",
                "paths:",
                f"  resource_root: {software / 'resources'}",
                f"  software_root: {software}",
                f"  conda_prefix: {software / 'conda'}",
                f"  apptainer_prefix: {software / 'apptainer'}",
                "encryption:",
                f"  sender_private_key: {tmp_path / 'sender'}",
                "  recipient_public_keys:",
                f"    - {tmp_path / 'recipient.pub'}",
                f"  decryption_private_key: {tmp_path / 'recipient'}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    destination = tmp_path / "profile"

    subprocess.run(
        [
            sys.executable,
            str(INSTALLER),
            "--site-config",
            str(site),
            "--profile-name",
            "test-profile",
            "--destination",
            str(destination),
        ],
        check=True,
    )

    profile = yaml.safe_load((destination / "config.yaml").read_text())
    assert profile["snakefile"] == str(ROOT / "Snakefile")
    assert profile["conda-prefix"] == str(software / "conda")
    assert profile["apptainer-prefix"] == str(software / "apptainer")
    assert profile["rerun-triggers"] == "mtime"
    assert profile["config"] == [f"site_config={site.resolve()}"]
