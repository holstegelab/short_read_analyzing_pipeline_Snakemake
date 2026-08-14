import argparse
import configparser
import concurrent.futures
import logging
import os
import posixpath
import shlex
import subprocess
import sys
import threading
import time
import zlib
from pathlib import Path

LOG = logging.getLogger("fs_to_dcache")
DEFAULT_COPY_TIMEOUT = "300m"
DEFAULT_RCLONE_CONFIG_CANDIDATES = [Path("~/config/rclone/rclone.conf"), Path("~/.config/rclone/rclone.conf")]  # first path is Snellius-specific


def setup_logging(verbose: bool):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s [%(levelname)s] %(message)s")


def format_bytes(value: int) -> str:
    units = ["B", "KiB", "MiB", "GiB", "TiB", "PiB"]
    size = float(value)
    for unit in units:
        if abs(size) < 1024 or unit == units[-1]:
            if unit == "B":
                return f"{int(size)}{unit}"
            return f"{size:.1f}{unit}"
        size /= 1024


def resolve_default_rclone_config() -> Path:
    env_value = os.environ.get("RCLONE_CONFIG")
    if env_value:
        return Path(env_value).expanduser()
    expanded = [candidate.expanduser() for candidate in DEFAULT_RCLONE_CONFIG_CANDIDATES]
    for candidate in expanded:
        if candidate.exists():
            return candidate
    return expanded[0]


def load_rclone_config(path: Path) -> configparser.ConfigParser:
    resolved = path.expanduser()
    if not resolved.exists():
        searched = ", ".join(str(candidate.expanduser()) for candidate in DEFAULT_RCLONE_CONFIG_CANDIDATES)
        raise FileNotFoundError(f"rclone config does not exist: {resolved}. Checked default locations: {searched}")
    parser = configparser.ConfigParser()
    with resolved.open("r", encoding="utf-8") as handle:
        parser.read_file(handle)
    if not parser.sections():
        raise ValueError(f"rclone config has no remotes: {resolved}")
    return parser


def resolve_remote_name(parser: configparser.ConfigParser, requested_remote: str | None) -> str:
    if requested_remote:
        if not parser.has_section(requested_remote):
            raise ValueError(f"remote {requested_remote!r} not found in config; available remotes: {', '.join(parser.sections())}")
        return requested_remote
    sections = parser.sections()
    if len(sections) == 1:
        return sections[0]
    raise ValueError("--remote is required when the config contains multiple remotes: " + ", ".join(sections))


def resolve_api_url(explicit_api: str | None, remote_config: configparser.SectionProxy) -> str | None:
    return explicit_api or os.environ.get("DCACHE_API") or os.environ.get("ADA_API") or remote_config.get("api", fallback=None)


def run_command(cmd: list[str]):
    LOG.debug("cmd: %s", " ".join(shlex.quote(str(x)) for x in cmd))
    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as exc:
        if exc.stdout:
            LOG.error("stdout: %s", exc.stdout.strip())
        if exc.stderr:
            LOG.error("stderr: %s", exc.stderr.strip())
        raise
    if result.stdout:
        LOG.debug("stdout: %s", result.stdout.strip())
    if result.stderr:
        LOG.debug("stderr: %s", result.stderr.strip())
    return result


def normalize_adler(value: str) -> str:
    s = str(value).strip().lower()
    if s.startswith("0x"):
        s = s[2:]
    s = "".join(ch for ch in s if ch in "0123456789abcdef")
    if not s:
        return s
    return s.zfill(8)[-8:]


def adler32_local(local_path: Path) -> str:
    adler = 1
    with local_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            adler = zlib.adler32(chunk, adler)
    adler &= 0xFFFFFFFF
    return f"{adler:08x}"


def plan_transfers(source: Path, destination: str, recursive: bool) -> list[dict]:
    source = source.expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(f"source does not exist: {source}")

    destination_is_dir = destination.endswith("/")
    cleaned_destination = destination.strip("/")
    if not cleaned_destination:
        raise ValueError("destination must not be empty")

    if source.is_file():
        resolved = source.resolve(strict=True)
        stat = resolved.stat()
        remote_path = posixpath.join(cleaned_destination, source.name) if destination_is_dir else cleaned_destination
        return [
            {
                "source": source,
                "resolved_source": resolved,
                "rel": source.name,
                "size": stat.st_size,
                "remote_path": remote_path,
            }
        ]

    if not source.is_dir():
        raise ValueError(f"source is not a regular file or directory: {source}")

    if not recursive:
        raise ValueError("source is a directory; use -R/--recursive to copy directories")

    remote_root = posixpath.join(cleaned_destination, source.name) if destination_is_dir else cleaned_destination
    discovered = []
    for path in sorted(source.rglob("*")):
        if not path.is_file():
            continue
        try:
            resolved = path.resolve(strict=True) if path.is_symlink() else path
        except FileNotFoundError as exc:
            raise FileNotFoundError(f"symlink target missing for {path}") from exc
        rel = path.relative_to(source)
        rel_posix = str(rel).replace(os.sep, "/")
        stat = resolved.stat()
        discovered.append(
            {
                "source": path,
                "resolved_source": resolved,
                "rel": rel_posix,
                "size": stat.st_size,
                "remote_path": posixpath.join(remote_root, rel_posix),
            }
        )
    return discovered


class Progress:
    def __init__(self, total_files: int, total_bytes: int):
        self.total_files = total_files
        self.total_bytes = total_bytes
        self.validated_files = 0
        self.validated_bytes = 0
        self.skipped_files = 0
        self.skipped_bytes = 0
        self.total_retries = 0
        self.failed = []
        self.lock = threading.Lock()
        self.start_time = time.monotonic()

    def success(self, rel_path: str, size: int, attempts: int = 1, skipped: bool = False):
        with self.lock:
            self.validated_files += 1
            self.validated_bytes += size
            if skipped:
                self.skipped_files += 1
                self.skipped_bytes += size
            if attempts > 1:
                self.total_retries += attempts - 1
            return self.validated_files, self.validated_bytes

    def failure(self, rel_path: str, exc: Exception):
        with self.lock:
            self.failed.append((rel_path, str(exc)))
            return len(self.failed)

    @property
    def done(self) -> int:
        return self.validated_files + len(self.failed)


class ProgressBar:
    """Thread-safe single-line progress bar written to stderr."""

    BAR_WIDTH = 30

    def __init__(self, progress: Progress, stream=None):
        self.progress = progress
        self.stream = stream or sys.stderr
        self._is_tty = hasattr(self.stream, "isatty") and self.stream.isatty()

    def update(self, last_file: str = ""):
        p = self.progress
        total = p.total_files
        done = p.done
        pct = done / total if total else 1.0
        filled = int(self.BAR_WIDTH * pct)
        bar = "\u2588" * filled + "\u2591" * (self.BAR_WIDTH - filled)
        elapsed = time.monotonic() - p.start_time
        elapsed_str = self._fmt_duration(elapsed)
        eta_str = self._fmt_duration(elapsed / pct - elapsed) if pct > 0 and pct < 1.0 else "--:--"
        errors_part = f" err:{len(p.failed)}" if p.failed else ""
        line = (
            f"\r  {bar}  {done}/{total} files  "
            f"{format_bytes(p.validated_bytes)}/{format_bytes(p.total_bytes)}  "
            f"{elapsed_str}<{eta_str}{errors_part}"
        )
        if last_file:
            max_name = 40
            display = last_file if len(last_file) <= max_name else "..." + last_file[-(max_name - 3):]
            line += f"  {display}"
        if self._is_tty:
            self.stream.write(f"\r\033[K{line}")
            self.stream.flush()

    def finish(self):
        if self._is_tty:
            self.stream.write("\r\033[K")
            self.stream.flush()

    @staticmethod
    def _fmt_duration(seconds: float) -> str:
        s = int(seconds)
        if s < 3600:
            return f"{s // 60:02d}:{s % 60:02d}"
        return f"{s // 3600}:{(s % 3600) // 60:02d}:{s % 60:02d}"


class DCacheCopier:
    def __init__(
        self,
        rclone_config: Path,
        remote: str,
        ada_cmd: str,
        api: str | None,
        max_retries: int,
        retry_wait: int,
        copy_timeout: str,
        skip_verified: bool = False,
    ):
        self.rclone_config = Path(rclone_config).expanduser()
        self.remote = remote
        self.ada_cmd = ada_cmd
        self.api = api
        self.max_retries = max_retries
        self.retry_wait = retry_wait
        self.copy_timeout = copy_timeout
        self.skip_verified = skip_verified
        self._seen_dirs = set()
        self._dirs_lock = threading.Lock()

        if not self.rclone_config.exists():
            raise FileNotFoundError(f"rclone config does not exist: {self.rclone_config}")

    def copy(self, file_entry: dict):
        local_path = Path(file_entry["resolved_source"])
        relative_path = file_entry["rel"].replace(os.sep, "/")
        remote_path = str(file_entry["remote_path"]).strip("/")
        if not remote_path:
            raise ValueError("remote path could not be derived")
        remote_dir = posixpath.dirname(remote_path)
        remote_name = posixpath.basename(remote_path)
        local_adler = adler32_local(local_path)
        remote_adler = ""

        if self.skip_verified:
            try:
                remote_adler = self._remote_adler(remote_path)
                if normalize_adler(local_adler) == normalize_adler(remote_adler):
                    LOG.info("skipping %s (remote checksum matches)", relative_path)
                    return {
                        "rel": relative_path,
                        "remote_path": remote_path,
                        "size": file_entry.get("size", 0),
                        "local_adler": local_adler,
                        "remote_adler": remote_adler,
                        "attempt": 0,
                        "skipped": True,
                    }
            except Exception:
                LOG.debug("remote checksum not available for %s; proceeding with upload", relative_path)

        for attempt in range(self.max_retries + 1):
            self._rclone_mkdir(remote_dir)
            self._rclone_copy(local_path, remote_dir, remote_name)
            remote_adler = self._remote_adler(remote_path)
            if normalize_adler(local_adler) == normalize_adler(remote_adler):
                return {
                    "rel": relative_path,
                    "remote_path": remote_path,
                    "size": file_entry.get("size", 0),
                    "local_adler": local_adler,
                    "remote_adler": remote_adler,
                    "attempt": attempt + 1,
                    "skipped": False,
                }
            LOG.warning(
                "checksum mismatch for %s: local=%s remote=%s (attempt %d/%d)",
                relative_path,
                local_adler,
                remote_adler,
                attempt + 1,
                self.max_retries + 1,
            )
            self._rclone_delete(remote_dir, remote_name)
            if attempt < self.max_retries:
                time.sleep(self.retry_wait)

        raise RuntimeError(f"checksum mismatch for {relative_path}: local={local_adler} remote={remote_adler}")

    def _rclone_mkdir(self, remote_dir: str):
        if not remote_dir:
            return
        with self._dirs_lock:
            if remote_dir in self._seen_dirs:
                return
            self._seen_dirs.add(remote_dir)
        run_command(["rclone", "--config", str(self.rclone_config), "mkdir", f"{self.remote}:{remote_dir}"])

    def _rclone_copy(self, local_path: Path, remote_dir: str, remote_name: str):
        remote_path = f"{remote_dir}/{remote_name}" if remote_dir else remote_name
        run_command(
            [
                "rclone",
                "--config",
                str(self.rclone_config),
                "-v",
                "--timeout",
                self.copy_timeout,
                "copyto",
                str(local_path),
                f"{self.remote}:{remote_path}",
            ]
        )

    def _rclone_delete(self, remote_dir: str, remote_name: str):
        remote_path = f"{remote_dir}/{remote_name}" if remote_dir else remote_name
        run_command(["rclone", "--config", str(self.rclone_config), "-v", "deletefile", f"{self.remote}:{remote_path}"])

    def _remote_adler(self, remote_path: str) -> str:
        cmd = [self.ada_cmd, "--tokenfile", str(self.rclone_config)]
        if self.api:
            cmd.extend(["--api", self.api])
        cmd.extend(["--checksum", remote_path])
        result = run_command(cmd)
        stdout = result.stdout.strip().split()
        for token in stdout:
            if "=" in token:
                key, value = token.split("=", 1)
                if key.lower().startswith("adler"):
                    return value.strip()
        for line in result.stdout.splitlines():
            if "adler32" not in line.lower():
                continue
            parts = line.replace(",", " ").split()
            for part in parts:
                if part.lower().startswith("adler32="):
                    return part.split("=", 1)[1]
        raise RuntimeError(f"unable to parse remote checksum for {remote_path}; ada output: {result.stdout.strip()!r}")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Copy local files to dCache and validate every upload with Adler-32.",
        epilog=(
            "File copy: a destination ending in '/' is treated as a remote directory.\n"
            "File copy: otherwise the destination is treated as the exact remote file path.\n"
            "Directory copy requires -R/--recursive. If the destination ends in '/', the source directory name is preserved below it.\n"
            "Otherwise the destination is treated as the exact target directory path."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("source", type=Path, help="Local file or directory to copy")
    parser.add_argument("destination", help="Remote dCache path")
    parser.add_argument("-R", "--recursive", action="store_true", help="Copy directories recursively")
    parser.add_argument(
        "--config",
        "--rclone-config",
        dest="config",
        type=Path,
        help="rclone config file. Defaults to RCLONE_CONFIG, then ~/config/rclone/rclone.conf, then ~/.config/rclone/rclone.conf",
    )
    parser.add_argument("--remote", default=os.environ.get("RCLONE_REMOTE"), help="rclone remote name. If omitted, the only config section is used")
    parser.add_argument("--ada", default=os.environ.get("ADA", "ada"), help="ada executable used to fetch remote checksums")
    parser.add_argument("--api", help="Optional dCache API URL for ada. If omitted, ada configuration or env vars are used")
    parser.add_argument("--dry-run", action="store_true", help="Show planned transfers without copying")
    parser.add_argument("--no-skip-verified", dest="skip_verified", action="store_false", default=True,
                        help="Re-upload files even if remote checksum already matches (default: skip verified)")
    parser.add_argument("--workers", type=int, default=4, help="Number of concurrent upload threads (default: 4)")
    parser.add_argument("--max-retries", type=int, default=3, help="Max retries on checksum mismatch (default: 3)")
    parser.add_argument("--retry-wait", type=int, default=60, help="Seconds to wait between retries (default: 60)")
    parser.add_argument("--copy-timeout", default=DEFAULT_COPY_TIMEOUT, help="rclone --timeout value (default: %(default)s)")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()
    setup_logging(args.verbose)

    rclone_config = args.config.expanduser() if args.config else resolve_default_rclone_config()
    config = load_rclone_config(rclone_config)
    remote = resolve_remote_name(config, args.remote)
    api = resolve_api_url(args.api, config[remote])
    if args.workers < 1:
        raise ValueError("--workers must be >= 1")
    if args.max_retries < 0:
        raise ValueError("--max-retries must be >= 0")
    if args.retry_wait < 0:
        raise ValueError("--retry-wait must be >= 0")

    files = plan_transfers(args.source, args.destination, args.recursive)
    if not files:
        LOG.info("no files to process")
        return 0

    total_bytes = sum(entry["size"] for entry in files)

    if args.dry_run:
        LOG.info("dry run: %d file(s), %s", len(files), format_bytes(total_bytes))
        for entry in files:
            LOG.info("  %s -> %s (%s)", entry["rel"], entry["remote_path"], format_bytes(entry["size"]))
        return 0

    progress = Progress(total_files=len(files), total_bytes=total_bytes)
    copier = DCacheCopier(
        rclone_config=rclone_config,
        remote=remote,
        ada_cmd=args.ada,
        api=api,
        max_retries=args.max_retries,
        retry_wait=args.retry_wait,
        copy_timeout=args.copy_timeout,
        skip_verified=args.skip_verified,
    )

    LOG.info("config  : %s", rclone_config)
    LOG.info("remote  : %s", remote)
    if api:
        LOG.info("api     : %s", api)
    LOG.info("source  : %s", args.source)
    LOG.info("dest    : %s", args.destination)
    LOG.info("files   : %d (%s)", len(files), format_bytes(total_bytes))
    LOG.info("workers : %d  retries: %d  skip-verified: %s",
             args.workers, args.max_retries, "yes" if args.skip_verified else "no")
    LOG.info("")

    bar = ProgressBar(progress)

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_to_entry = {executor.submit(copier.copy, entry): entry for entry in files}
        for future in concurrent.futures.as_completed(future_to_entry):
            entry = future_to_entry[future]
            try:
                result = future.result()
            except Exception as exc:
                failures = progress.failure(entry["rel"], exc)
                bar.finish()
                LOG.error("FAIL %s: %s", entry["rel"], exc)
                bar.update(entry["rel"])
                continue
            skipped = result.get("skipped", False)
            validated_files, validated_bytes = progress.success(
                entry["rel"], entry["size"],
                attempts=result.get("attempt", 1),
                skipped=skipped,
            )
            if skipped:
                LOG.debug("skip %s (verified)", result["rel"])
            else:
                bar.finish()
                LOG.info("ok   %s -> %s (%s)", result["rel"], result["remote_path"], format_bytes(entry["size"]))
            bar.update(result["rel"])

    bar.finish()

    elapsed = time.monotonic() - progress.start_time
    uploaded_bytes = progress.validated_bytes - progress.skipped_bytes
    uploaded_files = progress.validated_files - progress.skipped_files

    LOG.info("")
    LOG.info("====== transfer summary ======")
    LOG.info("  planned   : %d files, %s", progress.total_files, format_bytes(progress.total_bytes))
    LOG.info("  uploaded  : %d files, %s", uploaded_files, format_bytes(uploaded_bytes))
    if progress.skipped_files:
        LOG.info("  skipped   : %d files, %s (already verified)", progress.skipped_files, format_bytes(progress.skipped_bytes))
    if progress.total_retries:
        LOG.warning("  retries   : %d", progress.total_retries)
    if progress.failed:
        LOG.error("  FAILED    : %d", len(progress.failed))
        for rel_path, message in progress.failed:
            LOG.error("    %s | %s", rel_path, message)
    LOG.info("  elapsed   : %s", ProgressBar._fmt_duration(elapsed))
    if elapsed > 0 and uploaded_bytes > 0:
        LOG.info("  speed     : %s/s", format_bytes(int(uploaded_bytes / elapsed)))
    status = "COMPLETED" if not progress.failed else "COMPLETED WITH ERRORS"
    LOG.info("  status    : %s", status)
    LOG.info("==============================")

    return 1 if progress.failed else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        LOG.warning("interrupted")
        raise SystemExit(130)
    except (FileNotFoundError, ValueError) as exc:
        LOG.error("%s", exc)
        raise SystemExit(1)
