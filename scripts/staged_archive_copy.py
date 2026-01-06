import argparse
import json
import logging
import os
import queue
import shlex
import signal
import subprocess
import sys
import threading
import time
import traceback
import zlib
from pathlib import Path
import posixpath

LOG = logging.getLogger("staged_copy")

DCACHE_API = "https://dcacheview.grid.surfsara.nl:22880/api/v1"
DEFAULT_CONFIG_PATH = Path.home() / "macaroons" / "staged_archive_copy.conf"
DEFAULT_MACAROON_PATH = Path.home() / "macaroons" / "agh_full_snellius.conf"
DEFAULT_ADA_PATH = Path.home() / "projects" / "short_read_analyzing_pipeline_Snakemake" / "ada"
DEFAULT_RCLONE_REMOTE = "agh_full_snellius"

_PROGRESS_LOCK = threading.Lock()
_LAST_PROGRESS_LEN = 0
_LAST_PROGRESS_LOG_TS = 0.0
_LAST_PROGRESS_BYTES = 0
_LAST_PROGRESS_TIME = None
_LAST_SPEED = None
_PROGRESS_LOG_INTERVAL = 30.0
_PROGRESS_SPEED_WINDOW = 15.0
_IS_TTY = sys.stdout.isatty()


def _format_bytes(value: int) -> str:
    units = ["B", "KiB", "MiB", "GiB", "TiB", "PiB"]
    size = float(value)
    for unit in units:
        if abs(size) < 1024 or unit == units[-1]:
            if unit == "B":
                return f"{int(size)}{unit}"
            return f"{size:.1f}{unit}"
        size /= 1024
    return f"{size:.1f}{units[-1]}"


def _format_rate(value: float) -> str:
    if value is None or value <= 0:
        return "0B/s"
    units = ["B/s", "KiB/s", "MiB/s", "GiB/s", "TiB/s"]
    rate = float(value)
    for unit in units:
        if rate < 1024 or unit == units[-1]:
            if unit == "B/s":
                return f"{int(rate)}{unit}"
            return f"{rate:.1f}{unit}"
        rate /= 1024
    return f"{rate:.1f}{units[-1]}"


def _log_progress(summary: dict, context: str = "", force: bool = False, line_text: str | None = None):
    global _LAST_PROGRESS_LOG_TS
    total = summary.get("total_files", 0)
    if total <= 0:
        return
    now = time.monotonic()
    if not force and now - _LAST_PROGRESS_LOG_TS < _PROGRESS_LOG_INTERVAL:
        return
    counts = summary.get("counts", {})
    validated = counts.get(STATUS_VALIDATED, summary.get("validated_files", 0))
    pending = counts.get(STATUS_PENDING, 0)
    staged = counts.get(STATUS_STAGED, 0)
    copying = counts.get(STATUS_COPYING, 0)
    failed = counts.get(STATUS_FAILED, 0)
    total_bytes = summary.get("total_bytes", 0)
    validated_bytes = summary.get("validated_bytes", 0)
    remaining = summary.get("remaining_files", max(total - validated, 0))
    remaining_bytes = summary.get("remaining_bytes", max(total_bytes - validated_bytes, 0))
    percent = (validated / total) * 100 if total else 0.0
    byte_percent = (validated_bytes / total_bytes) * 100 if total_bytes else 0.0
    meta = summary.get("meta") or {}
    source = meta.get("source", "-")
    dest = meta.get("dest", "-")
    stage_active = summary.get("staging_bytes_active", 0)
    stage_cap = summary.get("staging_cap_bytes", 0)
    stage_workers = (summary.get("stage_workers_active", 0), summary.get("stage_workers_total", 0))
    copy_workers = (summary.get("copy_workers_active", 0), summary.get("copy_workers_total", 0))
    speed_label = summary.get("_speed_label", "0B/s")
    stage_queue = summary.get("pending_batches", 0)
    copy_queue = summary.get("copy_queue", 0)
    label = f"[{context}] " if context else ""
    if line_text is None:
        line_text = (
            f"{source} -> {dest} | stage {_format_bytes(stage_active)}/{_format_bytes(stage_cap)} | "
            f"workers s={stage_workers[0]}/{stage_workers[1]} c={copy_workers[0]}/{copy_workers[1]} | "
            f"files {validated}/{total} ({percent:.1f}%) todo={remaining} | "
            f"bytes {_format_bytes(validated_bytes)}/{_format_bytes(total_bytes)} ({byte_percent:.1f}%) | "
            f"remain {_format_bytes(remaining_bytes)} | speed {speed_label} | "
            f"queues stage={stage_queue} copy={copy_queue} | failed={failed}"
        )
    LOG.info("%s%s", label, line_text)
    _LAST_PROGRESS_LOG_TS = now


def update_progress_line(summary: dict, context: str = ""):
    global _LAST_PROGRESS_LEN, _LAST_PROGRESS_BYTES, _LAST_PROGRESS_TIME, _LAST_SPEED
    total = summary.get("total_files", 0)
    if total <= 0:
        return ""
    validated = summary.get("validated_files", 0)
    counts = summary.get("counts", {})
    pending = counts.get(STATUS_PENDING, 0)
    staged = counts.get(STATUS_STAGED, 0)
    copying = counts.get(STATUS_COPYING, 0)
    failed = counts.get(STATUS_FAILED, 0)
    total_bytes = summary.get("total_bytes", 0)
    validated_bytes = summary.get("validated_bytes", 0)
    remaining = summary.get("remaining_files", max(total - validated, 0))
    remaining_bytes = summary.get("remaining_bytes", max(total_bytes - validated_bytes, 0))
    percent = (validated / total) * 100 if total else 0.0
    byte_percent = (validated_bytes / total_bytes) * 100 if total_bytes else 0.0
    now = time.monotonic()
    current_bytes = validated_bytes
    speed = None
    delta_time = None
    if _LAST_PROGRESS_TIME is not None and now > _LAST_PROGRESS_TIME:
        delta_bytes = current_bytes - _LAST_PROGRESS_BYTES
        delta_time = now - _LAST_PROGRESS_TIME
        if delta_time > 0:
            speed = max(delta_bytes / delta_time, 0.0)
    if speed is not None:
        if _LAST_SPEED is None:
            smoothed = speed
        else:
            window = max(_PROGRESS_SPEED_WINDOW, 1.0)
            alpha = min(delta_time / window, 1.0) if delta_time is not None else 1.0
            smoothed = (_LAST_SPEED * (1.0 - alpha)) + (speed * alpha)
        _LAST_SPEED = smoothed
    speed_value = _LAST_SPEED if _LAST_SPEED is not None else 0.0
    _LAST_PROGRESS_TIME = now
    _LAST_PROGRESS_BYTES = current_bytes
    speed_label = _format_rate(speed_value)
    summary["_speed_label"] = speed_label
    meta = summary.get("meta") or {}
    source = meta.get("source", "-")
    dest = meta.get("dest", "-")
    stage_active = summary.get("staging_bytes_active", 0)
    stage_cap = summary.get("staging_cap_bytes", 0)
    stage_workers = (summary.get("stage_workers_active", 0), summary.get("stage_workers_total", 0))
    copy_workers = (summary.get("copy_workers_active", 0), summary.get("copy_workers_total", 0))
    stage_queue = summary.get("pending_batches", 0)
    copy_queue = summary.get("copy_queue", 0)
    label = f"[{context}] " if context else ""
    stage_str = f"{_format_bytes(stage_active)}/{_format_bytes(stage_cap)}"
    line_body = (
        f"{source} -> {dest} | stage {stage_str} | "
        f"workers s={stage_workers[0]}/{stage_workers[1]} c={copy_workers[0]}/{copy_workers[1]} | "
        f"files {validated}/{total} ({percent:.1f}%) todo={remaining} | "
        f"bytes {_format_bytes(validated_bytes)}/{_format_bytes(total_bytes)} ({byte_percent:.1f}%) | "
        f"remain {_format_bytes(remaining_bytes)} | speed {speed_label} | "
        f"queues stage={stage_queue} copy={copy_queue} | failed={failed}"
    )
    line = label + line_body
    with _PROGRESS_LOCK:
        if _IS_TTY:
            padding = max(_LAST_PROGRESS_LEN - len(line), 0)
            sys.stdout.write("\r" + line + (" " * padding))
            sys.stdout.flush()
            _LAST_PROGRESS_LEN = len(line)
        else:
            _LAST_PROGRESS_LEN = 0
    _log_progress(summary, context, line_text=line_body)
    return line_body


def finalize_progress_line(summary: dict, context: str = ""):
    global _LAST_PROGRESS_LEN
    if summary.get("total_files", 0) > 0:
        line_body = update_progress_line(summary, context)
    else:
        line_body = ""
        with _PROGRESS_LOCK:
            if _IS_TTY and _LAST_PROGRESS_LEN:
                sys.stdout.write("\n")
                sys.stdout.flush()
                _LAST_PROGRESS_LEN = 0
            else:
                _LAST_PROGRESS_LEN = 0
    _log_progress(summary, context, force=True, line_text=line_body)


STATUS_PENDING = "pending"
STATUS_STAGED = "staged"
STATUS_COPYING = "copying"
STATUS_VALIDATED = "validated"
STATUS_FAILED = "failed"

def load_config_file(path: Path):
    cfg = {}
    resolved = path.expanduser()
    if not resolved.exists():
        return cfg
    with resolved.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                continue
            key, value = [part.strip() for part in line.split("=", 1)]
            if not key:
                continue
            if key == "daget_flags":
                tokens = [v for v in value.split() if v]
                cfg[key] = tokens
            else:
                cfg[key] = value
    return cfg

class Ledger:
    def __init__(self, path: Path):
        self.path = path
        self.lock = threading.Lock()
        self.state = {"files": {}, "meta": {}}
        if self.path.exists():
            with self.path.open("r", encoding="utf-8") as handle:
                self.state = json.load(handle)
        else:
            self.path.parent.mkdir(parents=True, exist_ok=True)
        self._dirty = False
        self._save_interval = 5.0
        self._stop = False
        self._flusher = threading.Thread(target=self._auto_flush, daemon=True)
        self._flusher.start()

    def _auto_flush(self):
        while not self._stop:
            time.sleep(self._save_interval)
            self.flush()

    def flush(self):
        with self.lock:
            if not self._dirty:
                return
            tmp = self.path.with_suffix(".tmp")
            with tmp.open("w", encoding="utf-8") as handle:
                json.dump(self.state, handle, indent=2, sort_keys=True)
            tmp.replace(self.path)
            self._dirty = False

    def stop(self):
        self._stop = True
        self._flusher.join(timeout=2.0)
        self.flush()

    def info(self, rel_path: str):
        with self.lock:
            entry = self.state["files"].setdefault(rel_path, {"status": STATUS_PENDING, "size": 0, "retries": 0})
            return json.loads(json.dumps(entry))

    def update(self, rel_path: str, **fields):
        with self.lock:
            entry = self.state["files"].setdefault(rel_path, {"status": STATUS_PENDING, "size": 0, "retries": 0})
            entry.update(fields)
            entry.setdefault("history", []).append({"time": time.time(), **fields})
            self._dirty = True

    def set_meta(self, key: str, value):
        with self.lock:
            self.state["meta"][key] = value
            self._dirty = True

    def meta(self):
        with self.lock:
            return self.state.get("meta", {}).copy()

    def iter_files(self):
        with self.lock:
            for rel_path, entry in self.state["files"].items():
                yield rel_path, entry.copy()

    def summary(self):
        with self.lock:
            counts = {}
            total_bytes = 0
            validated_bytes = 0
            for entry in self.state["files"].values():
                status = entry.get("status", STATUS_PENDING)
                counts[status] = counts.get(status, 0) + 1
                size = entry.get("size", 0) or 0
                total_bytes += size
                if status == STATUS_VALIDATED:
                    validated_bytes += size
            total_files = len(self.state["files"])
            remaining_files = max(total_files - counts.get(STATUS_VALIDATED, 0), 0)
            remaining_bytes = max(total_bytes - validated_bytes, 0)
            return {
                "total_files": total_files,
                "total_bytes": total_bytes,
                "validated_files": counts.get(STATUS_VALIDATED, 0),
                "validated_bytes": validated_bytes,
                "counts": counts,
                "remaining_files": remaining_files,
                "remaining_bytes": remaining_bytes,
                "meta": self.state.get("meta", {}).copy(),
            }

    def log_summary(self, context: str = ""):
        summary = self.summary()
        total = summary["total_files"]
        label = f"[{context}] " if context else ""
        if total == 0:
            LOG.info("%sprogress: no files tracked", label)
            return summary
        counts = summary["counts"]
        validated = counts.get(STATUS_VALIDATED, 0)
        pending = counts.get(STATUS_PENDING, 0)
        staged = counts.get(STATUS_STAGED, 0)
        copying = counts.get(STATUS_COPYING, 0)
        failed = counts.get(STATUS_FAILED, 0)
        remaining = max(total - validated, 0)
        percent = (validated / total) * 100 if total else 0.0
        total_bytes = summary["total_bytes"]
        validated_bytes = summary["validated_bytes"]
        byte_percent = (validated_bytes / total_bytes) * 100 if total_bytes else 0.0
        LOG.info(
            "%sprogress: validated=%d/%d (%.1f%%), pending=%d, staged=%d, copying=%d, failed=%d, remaining=%d | bytes=%d/%d (%.1f%%)",
            label,
            validated,
            total,
            percent,
            pending,
            staged,
            copying,
            failed,
            remaining,
            validated_bytes,
            total_bytes,
            byte_percent,
        )
        if remaining == 0 and total > 0:
            LOG.info("%sprogress: all files validated", label)
        return summary

class Batch:
    def __init__(self, files):
        self.files = files
        self.size = sum(f["size"] for f in files)

class Scheduler:
    def __init__(self, source_root: Path, dest_root: str, ledger: Ledger, batch_size_bytes: int, staging_cap_bytes: int, daget_flags: list[str], use_daget: bool):
        self.source_root = source_root
        self.dest_root = dest_root
        self.ledger = ledger
        self.batch_size_bytes = batch_size_bytes
        self.staging_cap_bytes = staging_cap_bytes
        self.daget_flags = list(daget_flags or [])
        self.use_daget = use_daget
        self.pending_batches = queue.Queue()
        self.copy_queue = queue.Queue()
        self.active_staged_bytes = 0
        self.active_lock = threading.Lock()
        self.worker_lock = threading.Lock()
        self.active_stage_workers = 0
        self.active_copy_workers = 0
        self.stage_workers_total = 0
        self.copy_workers_total = 0

    def enqueue(self, batch: Batch):
        self.pending_batches.put(batch)

    def start(self, stage_workers, copy_workers, copier):
        stage_threads = []
        self.stage_workers_total = max(stage_workers, 1)
        for idx in range(self.stage_workers_total):
            t = threading.Thread(target=self._stage_loop, name=f"stage-{idx}", daemon=True)
            t.start()
            stage_threads.append(t)

        copy_threads = []
        self.copy_workers_total = max(copy_workers, 1)
        for idx in range(self.copy_workers_total):
            t = threading.Thread(target=self._copy_loop, name=f"copy-{idx}", args=(copier,), daemon=True)
            t.start()
            copy_threads.append(t)

        # Wait for all scheduled work (including retries) to be processed.
        # We loop because copy failures can re-enqueue into pending_batches.
        try:
            while True:
                self.pending_batches.join()
                self.copy_queue.join()
                if self.pending_batches.empty() and self.copy_queue.empty():
                    break
        finally:
            finalize_progress_line(self._progress_summary(), "copy")

        # Now it is safe to stop workers
        for _ in stage_threads:
            self.pending_batches.put(None)
        for t in stage_threads:
            t.join()

        for _ in copy_threads:
            self.copy_queue.put(None)
        for t in copy_threads:
            t.join()

    def _stage_loop(self):
        while True:
            batch = self.pending_batches.get()
            if batch is None:
                self.pending_batches.task_done()
                return
            with self.worker_lock:
                self.active_stage_workers += 1
            while True:
                with self.active_lock:
                    if self.active_staged_bytes + batch.size <= self.staging_cap_bytes:
                        self.active_staged_bytes += batch.size
                        break
                time.sleep(2)
            try:
                if self.use_daget:
                    paths = [str(batch_file.get("resolved_source") or batch_file["source"]) for batch_file in batch.files]
                    cmd = ["daget", "-av", *self.daget_flags, *paths]
                    run_command(cmd)
                now = time.time()
                for file_entry in batch.files:
                    self.ledger.update(file_entry["rel"], status=STATUS_STAGED, stage_time=now)
                self.copy_queue.put(batch)
                update_progress_line(self._progress_summary(), "copy")
            except Exception:
                LOG.error("staging failed", exc_info=True)
                for file_entry in batch.files:
                    self.ledger.update(file_entry["rel"], status=STATUS_FAILED, last_error="staging")
                with self.active_lock:
                    self.active_staged_bytes -= batch.size
            finally:
                with self.worker_lock:
                    self.active_stage_workers = max(self.active_stage_workers - 1, 0)
                self.pending_batches.task_done()

    def _copy_loop(self, copier):
        while True:
            batch = self.copy_queue.get()
            if batch is None:
                self.copy_queue.task_done()
                return
            with self.worker_lock:
                self.active_copy_workers += 1
            failed = False
            retry_files = []
            succeeded_files = []
            succeeded_size = 0
            for file_entry in batch.files:
                status = self.ledger.info(file_entry["rel"]) ["status"]
                if status == STATUS_VALIDATED:
                    continue
                try:
                    self.ledger.update(file_entry["rel"], status=STATUS_COPYING)
                    update_progress_line(self._progress_summary(), "copy")
                    copier.copy(file_entry)
                    self.ledger.update(file_entry["rel"], status=STATUS_VALIDATED, validate_time=time.time())
                    succeeded_files.append(file_entry)
                    succeeded_size += file_entry.get("size", 0)
                except Exception:
                    LOG.error("copy failed for %s", file_entry["rel"], exc_info=True)
                    info = self.ledger.info(file_entry["rel"])
                    self.ledger.update(file_entry["rel"], status=STATUS_FAILED, retries=info.get("retries", 0) + 1, last_error="copy")
                    update_progress_line(self._progress_summary(), "copy")
                    failed = True
                    retry_files.append(file_entry)
            try:
                if succeeded_files and self.use_daget:
                    release_paths = [str(file_entry.get("resolved_source") or file_entry["source"]) for file_entry in succeeded_files]
                    run_command(["darelease", *release_paths])
            except Exception:
                LOG.warning("darelease failed", exc_info=True)
            with self.active_lock:
                self.active_staged_bytes -= succeeded_size
            thread_name = threading.current_thread().name
            if succeeded_files or retry_files:
                if LOG.isEnabledFor(logging.DEBUG):
                    if succeeded_files:
                        LOG.debug("%s: validated %d file(s) (%d bytes)", thread_name, len(succeeded_files), succeeded_size)
                    if retry_files:
                        LOG.debug("%s: %d file(s) queued for retry", thread_name, len(retry_files))
                update_progress_line(self._progress_summary(), "copy")
            if failed:
                LOG.warning("batch failed, requeueing after delay")
                time.sleep(60)
                if retry_files:
                    retry_batch = Batch(retry_files)
                    # Retry directly; files remain staged because we did not release them
                    self.copy_queue.put(retry_batch)
            self.copy_queue.task_done()
            with self.worker_lock:
                self.active_copy_workers = max(self.active_copy_workers - 1, 0)

    def _progress_summary(self):
        summary = self.ledger.summary()
        summary.update(
            {
                "staging_bytes_active": self.active_staged_bytes,
                "staging_cap_bytes": self.staging_cap_bytes,
                "stage_workers_active": self.active_stage_workers,
                "stage_workers_total": self.stage_workers_total,
                "copy_workers_active": self.active_copy_workers,
                "copy_workers_total": self.copy_workers_total,
                "pending_batches": self.pending_batches.qsize(),
                "copy_queue": self.copy_queue.qsize(),
            }
        )
        return summary

def run_command(cmd):
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

class Copier:
    def __init__(self, dest_root: str, config_path: Path, remote: str, ada_script: Path, max_retries: int, sleep_seconds: int):
        cleaned_root = dest_root.strip("/")
        if not cleaned_root:
            raise ValueError("dest_root must not be empty")
        self.dest_root = cleaned_root
        self.config_path = Path(config_path).expanduser()
        self.remote = remote
        self.ada_script = Path(ada_script).expanduser()
        self.max_retries = max_retries
        self.sleep_seconds = sleep_seconds

    def copy(self, file_entry):
        local_path = Path(file_entry.get("resolved_source") or file_entry["source"])
        relative_path = file_entry["rel"].replace(os.sep, "/")
        remote_path = posixpath.join(self.dest_root, relative_path) if relative_path else self.dest_root
        remote_path = remote_path.strip("/")
        if not remote_path:
            raise ValueError("remote path could not be derived")
        remote_dir = posixpath.dirname(remote_path)
        if not remote_dir:
            remote_dir = self.dest_root
        remote_name = posixpath.basename(remote_path)
        local_adler = self._adler_local(local_path)
        retries = 0
        while retries <= self.max_retries:
            self._rclone_mkdir(remote_dir)
            self._rclone_copy(local_path, remote_dir, remote_name)
            remote_adler = self._remote_adler(remote_path)
            if self._normalize_adler(local_adler) == self._normalize_adler(remote_adler):
                return
            self._rclone_delete(remote_dir, remote_name)
            retries += 1
            if retries > self.max_retries:
                break
            time.sleep(self.sleep_seconds)
        raise RuntimeError(f"checksum mismatch for {relative_path}: local={local_adler} remote={remote_adler}")

    def _adler_local(self, local_path: Path):
        adler = 1
        with local_path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
                adler = zlib.adler32(chunk, adler)
        adler &= 0xffffffff
        return f"{adler:08x}"

    def _normalize_adler(self, value: str):
        s = str(value).strip().lower()
        if s.startswith("0x"):
            s = s[2:]
        s = "".join(ch for ch in s if ch in "0123456789abcdef")
        if not s:
            return s
        return s.zfill(8)[-8:]

    def _rclone_mkdir(self, remote_dir: str):
        if remote_dir:
            run_command(["rclone", "--config", str(self.config_path), "mkdir", f"{self.remote}:{remote_dir}"])

    def _rclone_copy(self, local_path: Path, remote_dir: str, remote_name: str):
        run_command(["rclone", "--config", str(self.config_path), "-v", "--timeout", "600m", "copyto", str(local_path), f"{self.remote}:{remote_dir}/{remote_name}"])

    def _rclone_delete(self, remote_dir: str, remote_name: str):
        run_command(["rclone", "--config", str(self.config_path), "-v", "deletefile", f"{self.remote}:{remote_dir}/{remote_name}"])

    def _remote_adler(self, remote_path: str):
        cmd = [str(self.ada_script), "--tokenfile", str(self.config_path), "--api", DCACHE_API, "--checksum", remote_path]
        result = run_command(cmd)
        stdout = result.stdout.strip().split()
        for token in stdout:
            if "=" in token:
                key, value = token.split("=", 1)
                if key.lower().startswith("adler"):
                    return value.strip()
        for line in result.stdout.splitlines():
            if "adler32" in line:
                parts = line.replace(",", " ").split()
                for part in parts:
                    if part.lower().startswith("adler32="):
                        return part.split("=", 1)[1]
        raise RuntimeError(f"unable to parse remote checksum for {remote_path}")

def discover_files(source_root: Path, ledger: Ledger):
    discovered = []
    for path in sorted(source_root.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(source_root)
        rel_str = str(rel)
        stat = path.stat()
        resolved = path
        if path.is_symlink():
            try:
                resolved = path.resolve(strict=True)
            except FileNotFoundError:
                LOG.error("symlink target missing for %s", rel_str)
                ledger.update(rel_str, status=STATUS_FAILED, last_error="missing symlink target")
                continue
        entry = ledger.info(rel_str)
        if entry.get("status") == STATUS_VALIDATED:
            continue
        record = {"source": path, "resolved_source": resolved, "rel": rel_str, "size": stat.st_size}
        ledger.update(rel_str, size=stat.st_size)
        discovered.append(record)
    return discovered

def make_batches(files, batch_size_bytes):
    batches = []
    current = []
    current_size = 0
    for entry in files:
        size = entry["size"]
        if size > batch_size_bytes:
            batches.append(Batch([entry]))
            continue
        if current_size + size > batch_size_bytes and current:
            batches.append(Batch(current))
            current = []
            current_size = 0
        current.append(entry)
        current_size += size
    if current:
        batches.append(Batch(current))
    return batches

def setup_logging(verbose: bool):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")

class GracefulExit(Exception):
    pass

def install_signal_handlers():
    def handler(signum, frame):
        raise GracefulExit()
    signal.signal(signal.SIGINT, handler)
    signal.signal(signal.SIGTERM, handler)

def parse_size(value: str):
    units = {"k": 1024, "m": 1024**2, "g": 1024**3, "t": 1024**4}
    if value is None:
        raise ValueError("size value is None")
    value = value.strip().lower()
    if not value:
        raise ValueError("size value is empty")
    if value[-1] in units:
        return int(float(value[:-1]) * units[value[-1]])
    return int(value)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("dest_root")
    parser.add_argument("state", type=Path)
    parser.add_argument("--config", type=Path, default=DEFAULT_MACAROON_PATH)
    parser.add_argument("--ada", type=Path, default=DEFAULT_ADA_PATH)
    parser.add_argument("--batch-size", default="1t")
    parser.add_argument("--staging-cap", default="2t")
    parser.add_argument("--stage-workers", type=int, default=1)
    parser.add_argument("--copy-workers", type=int, default=4)
    parser.add_argument("--max-retries", type=int, default=3)
    parser.add_argument("--retry-wait", type=int, default=60)
    parser.add_argument("--daget-flag", action="append")
    parser.add_argument("--config-file", type=Path, default=DEFAULT_CONFIG_PATH)
    parser.add_argument("--rclone-remote")
    parser.add_argument("--source-mode", choices=["tape", "fs"], default="tape")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    defaults = load_config_file(args.config_file)
    rclone_config = Path(args.config).expanduser() if args.config else (Path(defaults.get("rclone_config")).expanduser() if defaults.get("rclone_config") else None)
    ada_script = Path(args.ada).expanduser() if args.ada else (Path(defaults.get("ada_script")).expanduser() if defaults.get("ada_script") else None)
    daget_flags = args.daget_flag if args.daget_flag is not None else defaults.get("daget_flags", [])
    rclone_remote = args.rclone_remote or defaults.get("rclone_remote") or DEFAULT_RCLONE_REMOTE
    if not rclone_config or not rclone_config.exists():
        raise ValueError("rclone config path must be provided via CLI or config file")
    if not ada_script or not ada_script.exists():
        raise ValueError("ADA script path must be provided via CLI or config file")

    setup_logging(args.verbose)
    install_signal_handlers()

    ledger = Ledger(args.state)
    ledger.set_meta("source", str(args.source))
    ledger.set_meta("dest", args.dest_root)
    ledger.set_meta("source_mode", args.source_mode)

    files = discover_files(args.source, ledger)
    if not files:
        LOG.info("no files to process")
        ledger.stop()
        return
    batch_bytes = parse_size(args.batch_size)
    batches = make_batches(files, batch_bytes)
    staging_cap = parse_size(args.staging_cap)
    if staging_cap < batch_bytes:
        raise ValueError("staging cap must be >= batch size")

    copier = Copier(args.dest_root, rclone_config, rclone_remote, ada_script, args.max_retries, args.retry_wait)
    use_daget = args.source_mode == "tape"
    scheduler = Scheduler(args.source, args.dest_root, ledger, batch_bytes, staging_cap, daget_flags, use_daget)
    for batch in batches:
        scheduler.enqueue(batch)

    LOG.info("starting transfer: %d files across %d batches", len(files), len(batches))
    try:
        scheduler.start(args.stage_workers, args.copy_workers, copier)
    except GracefulExit:
        LOG.warning("received termination signal, flushing state")
    except Exception:
        LOG.error("fatal error: %s", traceback.format_exc())
        raise
    finally:
        ledger.stop()

if __name__ == "__main__":
    main()
