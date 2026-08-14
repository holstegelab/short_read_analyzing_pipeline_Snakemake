#!/usr/bin/env python3
"""Atomically decrypt one Crypt4GH-encrypted CRAM."""

import argparse
import os
import shlex
import stat
import subprocess
import sys
import tempfile
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="input .c4gh file")
    parser.add_argument("output", nargs="?", help="output CRAM (default: remove .c4gh)")
    parser.add_argument("--sk", required=True, help="recipient Crypt4GH private key")
    parser.add_argument(
        "--password-file",
        help="protected one-line passphrase file for unattended use",
    )
    parser.add_argument(
        "--crypt4gh",
        default="python -m crypt4gh",
        help="Crypt4GH command (default: python -m crypt4gh)",
    )
    parser.add_argument("--replace", action="store_true", help="replace existing output")
    return parser.parse_args()


def require_regular_file(path, label):
    candidate = Path(path).expanduser()
    if not candidate.is_file():
        raise FileNotFoundError(f"{label} is not a regular file: {candidate}")
    return candidate.resolve()


def check_private_permissions(path, label):
    mode = stat.S_IMODE(path.stat().st_mode)
    if mode & 0o077:
        raise PermissionError(
            f"{label} must not be accessible by group/others: {path} "
            f"(mode is {mode:04o}; use chmod 600)"
        )


def output_path(input_path, requested):
    if requested:
        return Path(requested).expanduser().absolute()
    name = str(input_path)
    if name.endswith(".c4gh"):
        name = name[:-5]
    else:
        name += ".decrypted"
    return Path(name).absolute()


def main():
    args = parse_args()
    encrypted = require_regular_file(args.input, "encrypted CRAM")
    private_key = require_regular_file(args.sk, "private key")
    check_private_permissions(private_key, "private key")

    destination = output_path(encrypted, args.output)
    if destination.exists() and not args.replace:
        raise FileExistsError(
            f"output already exists: {destination}; use --replace to replace it"
        )
    if not destination.parent.is_dir():
        raise FileNotFoundError(f"output directory does not exist: {destination.parent}")

    environment = os.environ.copy()
    if args.password_file:
        password_file = require_regular_file(args.password_file, "password file")
        check_private_permissions(password_file, "password file")
        passphrase = password_file.read_text(encoding="utf-8").rstrip("\r\n")
        if not passphrase:
            raise ValueError(f"password file is empty: {password_file}")
        environment["C4GH_PASSPHRASE"] = passphrase

    command = shlex.split(args.crypt4gh)
    if not command:
        raise ValueError("--crypt4gh must contain a command")
    command.extend(["decrypt", "--sk", str(private_key)])

    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{destination.name}.",
            suffix=".tmp",
            dir=destination.parent,
            delete=False,
        ) as temporary, encrypted.open("rb") as source:
            temporary_name = temporary.name
            subprocess.run(
                command,
                stdin=source,
                stdout=temporary,
                env=environment,
                check=True,
            )
            temporary.flush()
            os.fsync(temporary.fileno())
        os.replace(temporary_name, destination)
        temporary_name = None
    finally:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass

    print(destination)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"decrypt_cram: error: {error}", file=sys.stderr)
        sys.exit(1)
