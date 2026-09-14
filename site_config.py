"""Load one site configuration before pipeline modules derive their paths.

The controller and every worker reparse the root Snakefile.  The selected
configuration path is therefore normalized and exported through
``SHORT_READ_SITE_CONFIG`` before ``common``/``constants`` are imported.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping, MutableMapping

import yaml


SITE_CONFIG_ENV = "SHORT_READ_SITE_CONFIG"


class SiteConfigError(ValueError):
    """Raised for an invalid or inconsistent site configuration."""


@dataclass(frozen=True)
class SiteSettings:
    source: Path | None
    values: Mapping[str, Any]

    def __getitem__(self, key: str) -> Any:
        return self.values[key]


_DEFAULTS: dict[str, Any] = {
    "site_name": "snellius-legacy-defaults",
    "resource_root": "/gpfs/work3/0/qtholstg/hg38_res_v2",
    "software_root": None,
    "shared_tmp_root": "tmp",
    "fallback_tmp_root": "/scratch-local",
    "conda_prefix": "/home/hulsmanm/.snakemake",
    "apptainer_prefix": "/home/hulsmanm/.apptainer",
    "deepvariant_native_prefix": None,
    "hg19_reference": None,
    "hg19_b37_chry_reference": "/gpfs/work3/0/qtholstg/marc/genome/hg19_b37chrY.fa",
    "hg38_cram_reference": None,
    "glnexus_preset": None,
    "check_empty_script": None,
    "dcache_processed_config": None,
    "dcache_read_config": None,
    "dcache_read_prefix": "agh_full_snellius:/tape/processed",
    "ada": None,
    "dcache_cp": None,
    "rclone": "rclone",
    "daget": "/opt/dacommands/bin/daget",
    "dals": "/opt/dacommands/bin/dals",
    "darelease": "/opt/dacommands/bin/darelease",
    "encryption_sender_private_key": None,
    "encryption_recipient_public_keys": None,
    "decryption_private_key": None,
    "decryption_passphrase_file": "",
}

_SECTIONS = {
    "paths": {key: key for key in {
        "resource_root", "software_root", "shared_tmp_root",
        "fallback_tmp_root", "conda_prefix", "apptainer_prefix",
        "deepvariant_native_prefix", "hg19_reference",
        "hg19_b37_chry_reference", "hg38_cram_reference",
        "glnexus_preset", "check_empty_script",
    }},
    "storage": {key: key for key in {
        "dcache_processed_config", "dcache_read_config", "dcache_read_prefix",
    }},
    "tools": {key: key for key in {
        "ada", "dcache_cp", "rclone", "daget", "dals", "darelease",
    }},
    "encryption": {
        "sender_private_key": "encryption_sender_private_key",
        "recipient_public_keys": "encryption_recipient_public_keys",
        "decryption_private_key": "decryption_private_key",
        "decryption_passphrase_file": "decryption_passphrase_file",
    },
}

_PATH_FIELDS = {
    "resource_root", "software_root", "shared_tmp_root", "fallback_tmp_root",
    "conda_prefix", "apptainer_prefix", "deepvariant_native_prefix",
    "hg19_reference", "hg19_b37_chry_reference", "hg38_cram_reference",
    "glnexus_preset", "check_empty_script", "dcache_processed_config",
    "dcache_read_config", "encryption_sender_private_key",
    "decryption_private_key", "decryption_passphrase_file",
}
_TOOL_FIELDS = {"ada", "dcache_cp", "rclone", "daget", "dals", "darelease"}

_CURRENT: SiteSettings | None = None


def _expand(value: str, *, key: str) -> str:
    expanded = os.path.expanduser(os.path.expandvars(value))
    if "$" in expanded:
        raise SiteConfigError(
            f"site setting {key!r} contains an unresolved environment variable: {value!r}"
        )
    return expanded


def _load(path: Path) -> Mapping[str, Any]:
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except OSError as error:
        raise SiteConfigError(f"cannot read site configuration {path}: {error}") from error
    except yaml.YAMLError as error:
        raise SiteConfigError(f"invalid YAML in site configuration {path}: {error}") from error
    if not isinstance(data, Mapping):
        raise SiteConfigError(f"site configuration {path} must contain a YAML mapping")
    return data


def _flatten(raw: Mapping[str, Any]) -> dict[str, Any]:
    allowed_top = {"schema_version", "site", *_SECTIONS}
    unknown_top = sorted(set(raw) - allowed_top)
    if unknown_top:
        raise SiteConfigError(f"unknown top-level site setting(s): {', '.join(unknown_top)}")
    if raw.get("schema_version", 1) != 1:
        raise SiteConfigError("only site configuration schema_version 1 is supported")

    result: dict[str, Any] = {}
    if "site" in raw:
        if not isinstance(raw["site"], str) or not raw["site"].strip():
            raise SiteConfigError("site must be a non-empty string")
        result["site_name"] = raw["site"].strip()
    for section, field_map in _SECTIONS.items():
        values = raw.get(section, {}) or {}
        if not isinstance(values, Mapping):
            raise SiteConfigError(f"{section} must be a mapping")
        unknown = sorted(set(values) - set(field_map))
        if unknown:
            raise SiteConfigError(
                f"unknown {section} setting(s): {', '.join(unknown)}"
            )
        result.update((field_map[key], value) for key, value in values.items())
    return result


def _derive(configured: Mapping[str, Any], *, explicit: bool) -> dict[str, Any]:
    values = dict(_DEFAULTS)
    values.update(configured)
    resource_root = values.get("resource_root")
    if not isinstance(resource_root, str) or not resource_root.strip():
        raise SiteConfigError("paths.resource_root must be a non-empty string")
    values["resource_root"] = _expand(resource_root, key="resource_root").rstrip(os.sep)

    if values.get("software_root") is None:
        values["software_root"] = os.path.join(values["resource_root"], "software")
    derived = {
        "hg19_reference": "cram_refs/hg19.fa",
        "hg38_cram_reference": "cram_refs/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        "check_empty_script": "scripts/check_empty.py",
        "dcache_processed_config": ".agh/agh_processed.conf",
        "dcache_read_config": ".agh/agh_processed.conf",
        "encryption_sender_private_key": ".c4gh/master_key_for_encryption",
        "decryption_private_key": ".c4gh/recipient1",
    }
    for key, relative in derived.items():
        if values.get(key) is None:
            values[key] = os.path.join(values["resource_root"], relative)
    if values.get("deepvariant_native_prefix") is None:
        values["deepvariant_native_prefix"] = os.path.join(
            values["software_root"], "deepvariant-1.9.0-native"
        )
    if values.get("glnexus_preset") is None:
        values["glnexus_preset"] = os.path.join(
            values["software_root"], "Glnexus_preset.yml"
        )
    if values.get("encryption_recipient_public_keys") is None:
        values["encryption_recipient_public_keys"] = [
            os.path.join(values["resource_root"], ".c4gh", "recipient1.pub"),
            os.path.join(values["resource_root"], ".c4gh", "recipient2.pub"),
        ]

    recipients = values["encryption_recipient_public_keys"]
    if not isinstance(recipients, list) or not recipients:
        raise SiteConfigError(
            "encryption.encryption_recipient_public_keys must be a non-empty list"
        )
    for key, value in list(values.items()):
        if key in _PATH_FIELDS:
            if value in (None, ""):
                continue
            if not isinstance(value, str):
                raise SiteConfigError(f"site setting {key!r} must be a string")
            values[key] = _expand(value, key=key)
        elif key in _TOOL_FIELDS:
            if value is None:
                continue
            if not isinstance(value, str) or not value:
                raise SiteConfigError(f"site setting {key!r} must be a non-empty string")
            values[key] = _expand(value, key=key)
    values["encryption_recipient_public_keys"] = [
        _expand(value, key="encryption_recipient_public_keys")
        if isinstance(value, str) and value else value
        for value in recipients
    ]
    if not all(isinstance(value, str) and value for value in values["encryption_recipient_public_keys"]):
        raise SiteConfigError("every encryption recipient public-key path must be non-empty")

    if explicit:
        absolute_fields = _PATH_FIELDS - {"shared_tmp_root", "fallback_tmp_root", "decryption_passphrase_file"}
        for key in sorted(absolute_fields):
            value = values.get(key)
            if value and not os.path.isabs(value):
                raise SiteConfigError(f"site setting {key!r} must be an absolute path")
        for value in values["encryption_recipient_public_keys"]:
            if not os.path.isabs(value):
                raise SiteConfigError("encryption recipient public-key paths must be absolute")
    return values


def _apply(workflow_config: MutableMapping[str, Any] | None, settings: SiteSettings) -> None:
    values = settings.values
    if workflow_config is not None:
        defaults = {
            "deepvariant_native_prefix": values["deepvariant_native_prefix"],
            "agh_processed": values["dcache_processed_config"],
            "dcache_read_token": values["dcache_read_config"],
            "token": values["dcache_read_config"],
            "prefix": values["dcache_read_prefix"],
            "path_to_private_key": values["encryption_sender_private_key"],
            "path_to_decryption_private_key": values["decryption_private_key"],
            "path_to_decryption_passphrase_file": values["decryption_passphrase_file"],
        }
        for key, value in defaults.items():
            workflow_config.setdefault(key, value)

    environment = {
        "SHORT_READ_RESOURCE_ROOT": values["resource_root"],
        "DEEPVARIANT_NATIVE_PREFIX": values["deepvariant_native_prefix"],
        "SHORT_READ_HG19_REFERENCE": values["hg19_reference"],
        "SHORT_READ_HG19_B37_CHRY_REFERENCE": values["hg19_b37_chry_reference"],
        "SHORT_READ_HG38_CRAM_REFERENCE": values["hg38_cram_reference"],
        "C4GH_DECRYPTION_PRIVATE_KEY": values["decryption_private_key"],
        "DAGET": values["daget"],
        "DALS": values["dals"],
        "DARELEASE": values["darelease"],
        "RCLONE": values["rclone"],
        "GATK_CNV_ROOT": os.path.join(values["software_root"], "gatk_4.4"),
    }
    if values.get("ada"):
        environment["ADA"] = values["ada"]
    if values.get("dcache_cp"):
        environment["DCACHE_CP"] = values["dcache_cp"]
    for key, value in environment.items():
        if value:
            os.environ[key] = str(value)


def configure(workflow_config: MutableMapping[str, Any] | None = None) -> SiteSettings:
    """Load, validate and apply the selected site configuration once."""
    global _CURRENT
    selected = workflow_config.get("site_config") if workflow_config is not None else None
    selected = selected or os.environ.get(SITE_CONFIG_ENV)
    source = Path(_expand(str(selected), key="site_config")).resolve() if selected else None

    if _CURRENT is not None:
        if source != _CURRENT.source:
            raise SiteConfigError(
                f"site configuration was already loaded from {_CURRENT.source}; cannot switch to {source}"
            )
        _apply(workflow_config, _CURRENT)
        return _CURRENT

    configured = _flatten(_load(source)) if source else {}
    values = _derive(configured, explicit=source is not None)
    _CURRENT = SiteSettings(source=source, values=MappingProxyType(values))
    if source is not None:
        os.environ[SITE_CONFIG_ENV] = str(source)
    _apply(workflow_config, _CURRENT)
    return _CURRENT


def current() -> SiteSettings:
    return configure()
