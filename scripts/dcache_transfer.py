#!/usr/bin/env python3
"""Stage and transfer dCache files directly from Snellius.

Staging uses the managed ``ada`` executable. Downloads and uploads use the
lab's ``dcache_cp`` wrapper, which performs Adler-32 verification and installs
downloads atomically. Macaroon configs are passed by path and are never copied
or printed.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path, PurePosixPath


ONLINE_LOCALITIES = {"ONLINE", "ONLINE_AND_NEARLINE"}
REMOTE_RE = re.compile(r"^[A-Za-z0-9._-]+$")
ADLER_RE = re.compile(r"(?i)\badler32\s*=\s*([0-9a-f]{1,8})\b")


class TransferError(RuntimeError):
    pass


def log(message: str) -> None:
    print(f"[dcache] {message}", file=sys.stderr, flush=True)


def validate_remote(remote: str) -> str:
    if not REMOTE_RE.fullmatch(remote):
        raise ValueError(f"invalid dCache remote: {remote!r}")
    return remote


def normalize_remote_path(path: str) -> str:
    if "\n" in path or "\r" in path or "\0" in path:
        raise ValueError("remote path contains a control character")
    if any(part == ".." for part in PurePosixPath(path).parts):
        raise ValueError(f"remote path escapes its root: {path!r}")
    normalized = str(PurePosixPath("/" + path.lstrip("/")))
    if normalized == "/":
        raise ValueError("remote file path must not be the dCache root")
    return normalized


def adler32_file(path: Path) -> str:
    """Return local Adler-32; retained for focused tests and diagnostics."""
    import zlib

    checksum = 1
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            checksum = zlib.adler32(chunk, checksum)
    return f"{checksum & 0xFFFFFFFF:08x}"


def load_stage_list(path: Path) -> list[str]:
    result: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for lineno, line in enumerate(handle, 1):
            value = line.strip()
            if not value or value.startswith("#"):
                continue
            try:
                result.append(normalize_remote_path(value))
            except ValueError as exc:
                raise ValueError(f"{path}:{lineno}: {exc}") from exc
    if not result:
        raise ValueError(f"no remote paths in {path}")
    return result


def load_download_list(path: Path) -> list[tuple[str, Path]]:
    result: list[tuple[str, Path]] = []
    with path.open("r", encoding="utf-8") as handle:
        for lineno, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not row or row[0].strip().startswith("#"):
                continue
            if len(row) != 2:
                raise ValueError(f"{path}:{lineno}: expected two tab-separated columns")
            remote_path = normalize_remote_path(row[0].strip())
            local_path = Path(row[1].strip()).expanduser()
            result.append((remote_path, local_path))
    if not result:
        raise ValueError(f"no transfers in {path}")
    return result


def _find_tool(name: str, env_name: str) -> str:
    configured = os.environ.get(env_name)
    if configured:
        candidate = Path(configured).expanduser()
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())
        raise FileNotFoundError(f"{env_name} is not executable: {candidate}")

    if name == "ada":
        # Do not prefer /usr/hpc/bin/ada: it has produced incompatible
        # status/checksum behaviour on Snellius. dcache_cp manages this copy.
        candidates = [Path.home() / ".local/share/dcache_cp/ada"]
    else:
        candidates = [
            Path(sys.executable).with_name(name),
            Path.home()
            / "all_data/research/conda/marc/miniconda3/envs/snakemake/bin"
            / name,
        ]
    for candidate in candidates:
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())

    discovered = shutil.which(name)
    if discovered:
        return discovered
    raise FileNotFoundError(
        f"could not find {name}; set {env_name} or install dcache_cp in the "
        "Snakemake environment"
    )


def _run(argv: list[str], *, check: bool = True) -> subprocess.CompletedProcess:
    result = subprocess.run(
        argv,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if check and result.returncode != 0:
        detail = result.stderr.decode("utf-8", "replace").strip()
        raise TransferError(
            f"{Path(argv[0]).name} failed ({result.returncode}): {detail}"
        )
    return result


def _write_atomic(path: Path, content: str) -> None:
    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary_path = Path(temporary)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(content)
        os.replace(temporary_path, path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise


class DirectDCache:
    def __init__(
        self,
        config_path: Path,
        remote: str,
        *,
        dcache_cp: str | None = None,
        ada: str | None = None,
        copy_timeout: str = "600m",
        max_retries: int = 3,
        retry_wait: int = 60,
    ):
        self.remote = validate_remote(remote)
        self.config_path = config_path.expanduser().resolve()
        if not self.config_path.is_file():
            raise FileNotFoundError(
                f"dCache config does not exist: {self.config_path}"
            )
        self.dcache_cp = dcache_cp or _find_tool("dcache_cp", "DCACHE_CP")
        self.ada = ada or _find_tool("ada", "ADA")
        self.copy_timeout = str(copy_timeout)
        self.max_retries = int(max_retries)
        self.retry_wait = int(retry_wait)

    def _ada(self, *args: str, check: bool = True) -> subprocess.CompletedProcess:
        command = [
            self.ada,
            "--tokenfile",
            str(self.config_path),
            *map(str, args),
        ]
        for attempt in range(1, 9):
            result = _run(command, check=False)
            if result.returncode == 0:
                return result

            detail = result.stderr.decode("utf-8", "replace").strip()
            transient = any(
                marker in detail
                for marker in (
                    "429",
                    "Too Many Requests",
                    "502",
                    "503",
                    "504",
                )
            )
            if not transient or attempt == 8:
                if not check:
                    return result
                raise TransferError(
                    f"{Path(command[0]).name} failed "
                    f"({result.returncode}): {detail}"
                )

            wait_seconds = min(60.0, 2 ** attempt) + random.uniform(0.0, 5.0)
            log(
                f"transient ada API error; retry {attempt}/8 in "
                f"{wait_seconds:.1f}s"
            )
            time.sleep(wait_seconds)

        raise AssertionError("unreachable")

    def _count_online(self, path_list: Path, expected: int) -> int:
        result = self._ada("--longlist", "--from-file", str(path_list))
        listing = result.stdout.decode("utf-8", "replace")
        localities = [
            line.split()[-1]
            for line in listing.splitlines()
            if line.split()
        ]
        online = sum(value in ONLINE_LOCALITIES for value in localities)
        log(f"locality {online}/{expected} online")
        return online

    def stage(
        self,
        remote_paths: list[str],
        *,
        lifetime: str = "7D",
        poll_seconds: int = 60,
        timeout_seconds: int = 86400,
    ) -> None:
        paths = [normalize_remote_path(path) for path in remote_paths]
        if not paths:
            return

        fd, name = tempfile.mkstemp(prefix=".dcache-stage-", suffix=".txt")
        path_list = Path(name)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                handle.write("\n".join(paths) + "\n")

            # The managed ada supports bulk ``--from-file`` staging with its
            # 7D default, but (unlike single-file staging) does not accept a
            # separate --lifetime argument in this mode. Serialize only the
            # short API submission across compute nodes. Recall and polling
            # continue independently, so batches still stage in parallel.
            if str(lifetime).upper() != "7D":
                log(
                    f"warning: bulk staging uses ada's 7D lifetime; "
                    f"requested {lifetime}"
                )
            lock_path = self.config_path.parent / f".{self.remote}.stage.lock"
            with lock_path.open("a", encoding="utf-8") as lock_handle:
                fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
                try:
                    self._ada("--stage", "--from-file", str(path_list))
                finally:
                    fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)

            started = time.monotonic()
            while True:
                online = self._count_online(path_list, len(paths))
                if online == len(paths):
                    return
                if time.monotonic() - started >= int(timeout_seconds):
                    raise TransferError(
                        f"staging timed out after {timeout_seconds}s "
                        f"({online}/{len(paths)} online)"
                    )
                time.sleep(int(poll_seconds))
        finally:
            path_list.unlink(missing_ok=True)

    def unstage(self, remote_paths: list[str], *, check: bool = False) -> None:
        paths = [normalize_remote_path(path) for path in remote_paths]
        if not paths:
            return
        fd, name = tempfile.mkstemp(prefix=".dcache-unstage-", suffix=".txt")
        path_list = Path(name)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                handle.write("\n".join(paths) + "\n")
            self._ada(
                "--unstage",
                "--from-file",
                str(path_list),
                check=check,
            )
        finally:
            path_list.unlink(missing_ok=True)

    def remote_adler(
        self, remote_path: str, *, attempts: int = 10, wait_seconds: int = 2
    ) -> str:
        path = normalize_remote_path(remote_path)
        for attempt in range(1, attempts + 1):
            result = self._ada("--checksum", path)
            text = result.stdout.decode("utf-8", "replace")
            match = ADLER_RE.search(text)
            if match:
                return match.group(1).lower().zfill(8)
            if attempt < attempts:
                log(
                    f"checksum metadata not visible yet for {path}; "
                    f"retry {attempt}/{attempts}"
                )
                time.sleep(wait_seconds)
        raise TransferError(
            f"could not parse Adler-32 for {path} after {attempts} attempts"
        )

    def download(
        self,
        transfers: list[tuple[str, Path]],
        *,
        workers: int,
        transfer_slots: int,
        no_stage: bool,
        no_destage: bool,
        lifetime: str,
        poll_seconds: int,
        timeout_seconds: int,
    ) -> None:
        if workers < 1:
            raise ValueError("--workers must be >= 1")
        if transfer_slots < 1:
            raise ValueError("--transfer-slots must be >= 1")

        fd, name = tempfile.mkstemp(prefix=".dcache-download-", suffix=".tsv")
        dcache_file_list = Path(name)
        remote_paths = []
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                for remote_path, local_path in transfers:
                    normalized = normalize_remote_path(remote_path)
                    remote_paths.append(normalized)
                    local = local_path.expanduser().resolve()
                    local.parent.mkdir(parents=True, exist_ok=True)
                    handle.write(f"{self.remote}:{normalized}\t{local}\n")

            command = [
                self.dcache_cp,
                "--file-list",
                str(dcache_file_list),
                "--literal-file-list",
                "--config",
                str(self.config_path),
                "--remote",
                self.remote,
                "--workers",
                str(workers),
                "--copy-timeout",
                self.copy_timeout,
                "--max-retries",
                str(self.max_retries),
                "--retry-wait",
                str(self.retry_wait),
            ]
            if no_stage:
                command.append("--no-stage")
            else:
                command.extend(
                    [
                        "--stage-lifetime",
                        str(lifetime),
                        "--stage-poll",
                        str(poll_seconds),
                        "--stage-timeout",
                        str(timeout_seconds),
                    ]
                )
                if no_destage:
                    command.append("--no-destage")

            # Keep a cross-node hard limit in addition to the Snakemake
            # resource limit. The lock files live beside the shared macaroon
            # config, while advisory locks are released automatically if a
            # worker exits.
            lock_dir = (
                self.config_path.parent
                / f".{self.remote}.download-slots"
            )
            lock_dir.mkdir(mode=0o700, exist_ok=True)
            announced_wait = False
            while True:
                acquired = None
                for slot in range(int(transfer_slots)):
                    handle = (
                        lock_dir / f"slot-{slot}.lock"
                    ).open("a", encoding="utf-8")
                    try:
                        fcntl.flock(
                            handle.fileno(),
                            fcntl.LOCK_EX | fcntl.LOCK_NB,
                        )
                    except BlockingIOError:
                        handle.close()
                        continue
                    acquired = (slot, handle)
                    break

                if acquired is not None:
                    slot, handle = acquired
                    log(
                        f"using download slot {slot + 1}/"
                        f"{transfer_slots}"
                    )
                    try:
                        _run(command)
                    finally:
                        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
                        handle.close()
                    break

                if not announced_wait:
                    log(
                        f"waiting for one of {transfer_slots} shared "
                        "download slots"
                    )
                    announced_wait = True
                time.sleep(random.uniform(2.0, 5.0))

            # dcache_cp deliberately does not destage in --no-stage mode.
            # Here the pin belongs to the preceding stable batch stage, so
            # release it only after the verified download has succeeded.
            if no_stage and not no_destage:
                try:
                    self.unstage(remote_paths)
                except Exception as exc:
                    log(f"warning: failed to release stage pin: {exc}")
        finally:
            dcache_file_list.unlink(missing_ok=True)

    def upload(
        self,
        local_path: Path,
        remote_path: str,
        *,
        checksum_output: Path | None,
    ) -> None:
        local = local_path.expanduser().resolve()
        if not local.is_file():
            raise FileNotFoundError(f"upload source does not exist: {local}")
        destination = normalize_remote_path(remote_path)
        _run(
            [
                self.dcache_cp,
                str(local),
                f"{self.remote}:{destination}",
                "--config",
                str(self.config_path),
                "--remote",
                self.remote,
                "--workers",
                "1",
                "--copy-timeout",
                self.copy_timeout,
                "--max-retries",
                str(self.max_retries),
                "--retry-wait",
                str(self.retry_wait),
            ]
        )
        if checksum_output is not None:
            _write_atomic(checksum_output, self.remote_adler(destination) + "\n")


def add_transfer_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--remote", required=True)
    parser.add_argument("--dcache-cp")
    parser.add_argument("--ada")
    parser.add_argument("--copy-timeout", default="600m")
    parser.add_argument("--max-retries", type=int, default=3)
    parser.add_argument("--retry-wait", type=int, default=60)
    parser.add_argument("--dry-run", action="store_true")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    stage = subparsers.add_parser("stage", help="stage a one-column path list")
    add_transfer_options(stage)
    stage.add_argument("--file-list", type=Path, required=True)
    stage.add_argument("--lifetime", default="7D")
    stage.add_argument("--poll-seconds", type=int, default=60)
    stage.add_argument("--stage-timeout", type=int, default=86400)

    download = subparsers.add_parser(
        "download", help="download and verify a two-column file list"
    )
    add_transfer_options(download)
    download.add_argument("--file-list", type=Path, required=True)
    download.add_argument("--workers", type=int, default=2)
    download.add_argument(
        "--download-lock-slots",
        "--transfer-slots",
        dest="download_lock_slots",
        type=int,
        default=4,
        help=(
            "shared advisory download locks; separate from zslurm's "
            "dcache_download_slots scheduler resource"
        ),
    )
    download.add_argument("--no-stage", action="store_true")
    download.add_argument("--no-destage", action="store_true")
    download.add_argument("--lifetime", default="7D")
    download.add_argument("--poll-seconds", type=int, default=60)
    download.add_argument("--stage-timeout", type=int, default=86400)

    upload = subparsers.add_parser(
        "upload", help="upload and verify one local file"
    )
    add_transfer_options(upload)
    upload.add_argument("--source", type=Path, required=True)
    upload.add_argument("--destination", required=True)
    upload.add_argument("--checksum-output", type=Path)
    return parser


def make_client(args: argparse.Namespace) -> DirectDCache:
    return DirectDCache(
        args.config,
        args.remote,
        dcache_cp=args.dcache_cp,
        ada=args.ada,
        copy_timeout=args.copy_timeout,
        max_retries=args.max_retries,
        retry_wait=args.retry_wait,
    )


def main() -> int:
    args = build_parser().parse_args()

    if args.command == "stage":
        paths = load_stage_list(args.file_list)
        if args.dry_run:
            log(f"dry-run: would stage {len(paths)} file(s) directly")
            return 0
        make_client(args).stage(
            paths,
            lifetime=args.lifetime,
            poll_seconds=args.poll_seconds,
            timeout_seconds=args.stage_timeout,
        )
        return 0

    if args.command == "download":
        transfers = load_download_list(args.file_list)
        if args.dry_run:
            log(f"dry-run: would download {len(transfers)} file(s) directly")
            return 0
        make_client(args).download(
            transfers,
            workers=args.workers,
            transfer_slots=args.download_lock_slots,
            no_stage=args.no_stage,
            no_destage=args.no_destage,
            lifetime=args.lifetime,
            poll_seconds=args.poll_seconds,
            timeout_seconds=args.stage_timeout,
        )
        return 0

    if args.command == "upload":
        if args.dry_run:
            log(
                f"dry-run: would upload {args.source} to "
                f"{args.remote}:{args.destination}"
            )
            return 0
        make_client(args).upload(
            args.source,
            args.destination,
            checksum_output=args.checksum_output,
        )
        return 0

    raise AssertionError(f"unknown command: {args.command}")


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
    except (FileNotFoundError, TransferError, ValueError) as exc:
        log(f"error: {exc}")
        raise SystemExit(1)
