#!/usr/bin/env python3
"""Build a deterministic CRAM decryption ZIP without credentials."""

import argparse
import hashlib
import os
import tempfile
import zipfile
from pathlib import Path


FIXED_TIMESTAMP = (2026, 1, 1, 0, 0, 0)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--readme", required=True)
    parser.add_argument("--decrypt-script", required=True)
    parser.add_argument("--bam-revert", required=True)
    parser.add_argument("--environment", required=True)
    parser.add_argument("--key-instructions", required=True)
    parser.add_argument("--short-readme", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_file(path):
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    return source.read_bytes()


def zip_info(name, executable=False):
    info = zipfile.ZipInfo(name, date_time=FIXED_TIMESTAMP)
    info.compress_type = zipfile.ZIP_DEFLATED
    info.create_system = 3
    mode = 0o755 if executable else 0o644
    info.external_attr = (mode & 0xFFFF) << 16
    return info


def main():
    args = parse_args()
    members = {
        "README.md": (read_file(args.readme), False),
        "decrypt_cram.py": (read_file(args.decrypt_script), True),
        "bam_revert.py": (read_file(args.bam_revert), True),
        "environment.yml": (read_file(args.environment), False),
        "KEY_PACKAGE_INSTRUCTIONS.md": (read_file(args.key_instructions), False),
        "KORTE_README_CRAM_DECRYPTIE.md": (read_file(args.short_readme), False),
    }
    manifest = "".join(
        f"{hashlib.sha256(content).hexdigest()}  {name}\n"
        for name, (content, _) in sorted(members.items())
    ).encode("ascii")
    members["MANIFEST.sha256"] = (manifest, False)

    destination = Path(args.output).absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = None
    try:
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
        )
        os.close(descriptor)
        with zipfile.ZipFile(temporary_name, mode="w") as archive:
            for name, (content, executable) in sorted(members.items()):
                archive.writestr(zip_info(name, executable), content)
        os.replace(temporary_name, destination)
        temporary_name = None
    finally:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass


if __name__ == "__main__":
    main()
