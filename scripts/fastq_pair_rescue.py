#!/usr/bin/env python3
"""
Interleave two FASTQ files robustly, rescuing pairs when possible.

- Reads R1 and R2 in lockstep and outputs an interleaved FASTQ stream.
- Skips corrupt records (malformed headers, missing '+', missing qualities, length mismatch).
- If read names mismatch between sides, drops the pair and continues.
- Stops output when either side reaches EOF (to avoid emitting unpaired trailing reads).

Exit code is 0 on success; errors opening files will raise normally.
Summary statistics are printed to stderr unless --quiet is set.
"""
from __future__ import annotations

import argparse
import bz2
import gzip
import sys
from dataclasses import dataclass
from typing import Optional, TextIO, Tuple
from collections import OrderedDict
import shutil
import subprocess


@dataclass
class _ReaderState:
    pending_header: Optional[str] = None


class _ProcTextReader:
    def __init__(self, proc: subprocess.Popen):
        self._proc = proc
        self._stream = proc.stdout  # type: ignore[assignment]
        self._closed = False

    def readline(self) -> str:
        return self._stream.readline()  # type: ignore[no-any-return]

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        try:
            # Drain and wait, then report errors but do not raise
            _out, err = self._proc.communicate()
            rc = self._proc.returncode
            if rc not in (0, None):
                if err:
                    print(f"[fastq_pair_rescue] decompressor exited with code {rc}:\n{err}", file=sys.stderr)
                else:
                    print(f"[fastq_pair_rescue] decompressor exited with code {rc}", file=sys.stderr)
        except Exception:
            # Best-effort cleanup; do not raise
            pass


def _should_use_external(path: str, decompress: str) -> bool:
    if decompress == 'external':
        return True
    if decompress == 'internal':
        return False
    # auto
    if path.endswith('.gz') and shutil.which('pigz'):
        return True
    if path.endswith('.bz2') and (shutil.which('pbzip2') or shutil.which('bzip2')):
        return True
    return False


def _describe_input_source(path: str, decompress: str, threads: int) -> str:
    if path == '-':
        return 'stdin'
    use_ext = _should_use_external(path, decompress)
    if use_ext and path.endswith('.gz'):
        return f"pigz -cd {'-p '+str(threads) if threads and threads>0 else ''}"
    if use_ext and path.endswith('.bz2'):
        if shutil.which('pbzip2'):
            return f"pbzip2 -cd {'-p '+str(threads) if threads and threads>0 else ''}"
        return "bzip2 -cd"
    # internal
    if path.endswith('.gz'):
        return 'python-gzip'
    if path.endswith('.bz2'):
        return 'python-bz2'
    return 'plain-text'


def _open_text_file(path: str, mode: str, *, decompress: str = 'auto', threads: int = 0) -> Tuple[TextIO, bool]:
    if 'b' in mode:
        raise ValueError("Binary mode not supported; use text mode like 'rt' or 'wt'.")
    if path == '-':
        if 'r' in mode:
            return sys.stdin, False
        if 'w' in mode or 'a' in mode:
            return sys.stdout, False
        raise ValueError(f"Unsupported mode '{mode}' for standard streams")
    # Readers
    if 'r' in mode:
        use_ext = _should_use_external(path, decompress)
        if use_ext and path.endswith('.gz'):
            tflag = ['-p', str(threads)] if threads and threads > 0 else []
            proc = subprocess.Popen(['pigz', '-cd', *tflag, '--', path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding='utf-8', errors='replace')
            return _ProcTextReader(proc), True
        if use_ext and path.endswith('.bz2'):
            if shutil.which('pbzip2'):
                tflag = ['-p', str(threads)] if threads and threads > 0 else []
                proc = subprocess.Popen(['pbzip2', '-cd', *tflag, '--', path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding='utf-8', errors='replace')
            else:
                proc = subprocess.Popen(['bzip2', '-cd', '--', path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding='utf-8', errors='replace')
            return _ProcTextReader(proc), True
        # Fallback to internal readers
        if path.endswith('.gz'):
            return gzip.open(path, mode, encoding='utf-8', errors='replace'), True
        if path.endswith('.bz2'):
            return bz2.open(path, mode, encoding='utf-8', errors='replace'), True
        return open(path, mode, encoding='utf-8', errors='replace'), True
    # Writers
    if path.endswith('.gz'):
        return gzip.open(path, mode, encoding='utf-8', errors='replace'), True
    if path.endswith('.bz2'):
        return bz2.open(path, mode, encoding='utf-8', errors='replace'), True
    return open(path, mode, encoding='utf-8', errors='replace'), True


def _read_fastq4(handle: TextIO, state: _ReaderState) -> Optional[Tuple[str, str, str, str]]:
    # Find next header line starting with '@'
    header: Optional[str] = None
    if state.pending_header is not None:
        header = state.pending_header
        state.pending_header = None
    else:
        while True:
            line = handle.readline()
            if not line:
                return None
            line = line.rstrip('\r\n')
            if not line:
                continue
            if line.startswith('@'):
                header = line
                break
            # Skip stray lines quietly; caller logs summary counts.

    # Read sequence, plus, qualities
    seq_raw = handle.readline()
    plus_raw = handle.readline()
    qual_raw = handle.readline()

    if not seq_raw or not plus_raw or not qual_raw:
        # Try to resync: search for next header start and buffer it for next call
        _resync_to_header(handle, state, first_line=plus_raw if plus_raw else qual_raw)
        raise ValueError("Incomplete FASTQ record (missing one or more lines)")

    seq = seq_raw.rstrip('\r\n')
    plus = plus_raw.rstrip('\r\n')
    qual = qual_raw.rstrip('\r\n')

    if not plus.startswith('+'):
        # The plus line might actually be the next header; attempt to resync
        if plus.startswith('@'):
            state.pending_header = plus
        else:
            _resync_to_header(handle, state)
        raise ValueError("Malformed FASTQ record: '+' line missing")

    if len(qual) != len(seq):
        # Qual length mismatch; drop
        raise ValueError("Malformed FASTQ record: sequence/quality length mismatch")

    return header, seq, plus, qual


def _resync_to_header(handle: TextIO, state: _ReaderState, *, first_line: Optional[str] = None) -> None:
    # Try to find the next header and store in state.pending_header
    if first_line and first_line.startswith('@'):
        state.pending_header = first_line.rstrip('\r\n')
        return
    while True:
        line = handle.readline()
        if not line:
            return
        if line.startswith('@'):
            state.pending_header = line.rstrip('\r\n')
            return


def _norm(header: str) -> str:
    core = header[1:].strip()
    token = core.split()[0]
    if token.endswith('/1') or token.endswith('/2'):
        token = token[:-2]
    return token


@dataclass
class Stats:
    r1_total: int = 0
    r2_total: int = 0
    r1_corrupt: int = 0
    r2_corrupt: int = 0
    mismatched: int = 0
    pairs_written: int = 0


def interleave_rescue(r1_path: str, r2_path: str, out_path: str, quiet: bool = False, buffer_size: int = 2048, decompress: str = 'auto', threads: int = 0) -> Stats:
    if not quiet:
        print(f"[fastq_pair_rescue] Mode: decompress={decompress} threads={threads} buffer_size={buffer_size}", file=sys.stderr)
        print(f"[fastq_pair_rescue] R1: {r1_path} via {_describe_input_source(r1_path, decompress, threads)}", file=sys.stderr)
        print(f"[fastq_pair_rescue] R2: {r2_path} via {_describe_input_source(r2_path, decompress, threads)}", file=sys.stderr)
    in1, close1 = _open_text_file(r1_path, 'rt', decompress=decompress, threads=threads)
    in2, close2 = _open_text_file(r2_path, 'rt', decompress=decompress, threads=threads)
    out, close_out = _open_text_file(out_path, 'wt')

    st = Stats()
    s1 = _ReaderState()
    s2 = _ReaderState()
    buf1: "OrderedDict[str, Tuple[str, str, str, str]]" = OrderedDict()
    buf2: "OrderedDict[str, Tuple[str, str, str, str]]" = OrderedDict()
    eof1 = False
    eof2 = False

    try:
        while True:
            # Try to read next valid from R1 if we can grow buffer
            if not eof1 and len(buf1) < buffer_size:
                while True:
                    try:
                        rec1 = _read_fastq4(in1, s1)
                    except ValueError as e:
                        st.r1_corrupt += 1
                        if not quiet:
                            print(f"[fastq_pair_rescue] R1 corrupt: {e}", file=sys.stderr)
                        continue
                    if rec1 is None:
                        eof1 = True
                        break
                    st.r1_total += 1
                    n1 = _norm(rec1[0])
                    if n1 not in buf1:
                        buf1[n1] = rec1
                    break

            # Try to read next valid from R2 if we can grow buffer
            if not eof2 and len(buf2) < buffer_size:
                while True:
                    try:
                        rec2 = _read_fastq4(in2, s2)
                    except ValueError as e:
                        st.r2_corrupt += 1
                        if not quiet:
                            print(f"[fastq_pair_rescue] R2 corrupt: {e}", file=sys.stderr)
                        continue
                    if rec2 is None:
                        eof2 = True
                        break
                    st.r2_total += 1
                    n2 = _norm(rec2[0])
                    if n2 not in buf2:
                        buf2[n2] = rec2
                    break

            # If both EOF and no buffers, stop
            if eof1 and eof2 and not buf1 and not buf2:
                break

            # Find a matching name across buffers
            match_name: Optional[str] = None
            for name in buf1.keys():
                if name in buf2:
                    match_name = name
                    break

            if match_name is not None:
                r1 = buf1.pop(match_name)
                r2 = buf2.pop(match_name)
                # Write interleaved pair
                h1, s1q, p1, q1 = r1
                h2, s2q, p2, q2 = r2
                out.write(f"{h1}\n{s1q}\n{p1}\n{q1}\n")
                out.write(f"{h2}\n{s2q}\n{p2}\n{q2}\n")
                st.pairs_written += 1
                continue

            # No match yet. If both buffers are at capacity or one side hit EOF, evict one oldest to progress
            if (len(buf1) >= buffer_size or eof1) and buf1 and ((len(buf2) >= buffer_size or eof2) and buf2):
                # Drop oldest from the larger buffer (or R1 if equal)
                if len(buf1) >= len(buf2):
                    name, _rec = buf1.popitem(last=False)
                    st.mismatched += 1
                    if not quiet:
                        print(f"[fastq_pair_rescue] Dropping unmatched R1: {name}", file=sys.stderr)
                else:
                    name, _rec = buf2.popitem(last=False)
                    st.mismatched += 1
                    if not quiet:
                        print(f"[fastq_pair_rescue] Dropping unmatched R2: {name}", file=sys.stderr)
                continue

            # Otherwise, continue filling buffers
            if (not eof1 and len(buf1) < buffer_size) or (not eof2 and len(buf2) < buffer_size):
                continue

            # If we reach here, we cannot progress further
            break

    finally:
        if close1:
            in1.close()
        if close2:
            in2.close()
        if close_out:
            out.flush()
            out.close()

    if not quiet:
        print(
            (
                f"R1 total: {st.r1_total}, corrupt: {st.r1_corrupt}. "
                f"R2 total: {st.r2_total}, corrupt: {st.r2_corrupt}. "
                f"Mismatched dropped: {st.mismatched}. "
                f"Pairs written: {st.pairs_written}."
            ),
            file=sys.stderr,
        )
    return st


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Rescue and interleave paired FASTQ files.")
    p.add_argument('--r1', required=True, help='Path to R1 FASTQ (.gz/.bz2 supported)')
    p.add_argument('--r2', required=True, help='Path to R2 FASTQ (.gz/.bz2 supported)')
    p.add_argument('-o', '--output', default='-', help="Output path (use '-' for stdout)")
    p.add_argument('--quiet', action='store_true', help='Suppress progress messages on stderr')
    p.add_argument('--buffer-size', type=int, default=2048, help='Lookahead buffer size per side for resynchronization')
    p.add_argument('--decompress', choices=['auto','external','internal'], default='auto', help='Decompression mode for inputs')
    p.add_argument('--threads', type=int, default=0, help='Threads for external decompression tools (pigz/pbzip2)')
    return p.parse_args()


def main() -> None:
    args = parse_args()
    interleave_rescue(args.r1, args.r2, args.output, quiet=args.quiet, buffer_size=args.buffer_size, decompress=args.decompress, threads=args.threads)


if __name__ == '__main__':
    main()
