#!/usr/bin/env python3
"""Decrypt crypt4gh-encrypted CRAMs.

Usage:
    python scripts/decrypt_c4gh.py encrypted.c4gh [output.cram]

The script loads key locations from the pipeline configuration (same
locations used by Encrypt.smk). You can override them via CLI options.
"""

import argparse
import os
import sys

import yaml

pj = os.path.join
RESOURCES = "/gpfs/work3/0/qtholstg/hg38_res_v2/"

DEFAULT_CONFIG_PATH = "config.yaml"
DEFAULT_PRIVATE_KEY = pj(RESOURCES, ".c4gh/recipient1")
DEFAULT_PASSPHRASE = None


def load_config(path: str) -> dict:
    if not os.path.exists(path):
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def resolve_default(key: str, cfg: dict, fallback):
    return cfg.get(key, fallback) if isinstance(cfg, dict) else fallback


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="encrypted .c4gh file")
    parser.add_argument("output", nargs="?", help="output CRAM path")
    parser.add_argument(
        "--config", default=DEFAULT_CONFIG_PATH, help="path to config.yaml"
    )
    parser.add_argument(
        "--sk", dest="private_key", help="path to crypt4gh private key"
    )
    parser.add_argument(
        "--password-file", dest="password_file", help="path to passphrase file"
    )
    parser.add_argument(
        "--crypt4gh", default="python -m crypt4gh", help="crypt4gh executable"
    )
    parser.add_argument(
        "--replace", action="store_true", help="overwrite existing output"
    )
    return parser


def infer_output_path(input_path: str, provided: str | None) -> str:
    if provided:
        return provided
    if input_path.endswith(".c4gh"):
        return input_path[:-5]
    return f"{input_path}.decrypted"


def validate_file(path: str, description: str) -> None:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing {description}: {path}")


def decrypt(args) -> None:
    config = load_config(args.config) if args.config else {}
    key_path = args.private_key or resolve_default(
        "path_to_decryption_private_key", config, DEFAULT_PRIVATE_KEY
    )
    password_file = args.password_file or resolve_default(
        "path_to_decryption_passphrase_file", config, DEFAULT_PASSPHRASE
    )

    validate_file(args.input, "encrypted CRAM")
    validate_file(key_path, "private key")
    if password_file:
        validate_file(password_file, "passphrase file")

    output_path = infer_output_path(args.input, args.output)
    if os.path.exists(output_path) and not args.replace:
        raise FileExistsError(
            f"Output already exists: {output_path}. Use --replace to overwrite."
        )

    command = [args.crypt4gh, "decrypt", "--sk", key_path]
    if password_file:
        command.extend(["--password-file", password_file])
    command_str = " ".join(command)

    with open(args.input, "rb") as src, open(output_path, "wb") as dst:
        import subprocess

        process = subprocess.Popen(
            command_str,
            shell=True,
            stdin=src,
            stdout=dst,
            stderr=subprocess.PIPE,
        )
        _, stderr = process.communicate()
        if process.returncode != 0:
            sys.stderr.write(stderr.decode("utf-8", errors="replace"))
            raise RuntimeError(
                f"crypt4gh failed with exit code {process.returncode}".strip()
            )


def main():
    parser = build_parser()
    args = parser.parse_args()
    try:
        decrypt(args)
    except Exception as exc:  # noqa: BLE001
        import traceback

        traceback.print_exc()
        sys.stderr.write(f"Error: {exc}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
