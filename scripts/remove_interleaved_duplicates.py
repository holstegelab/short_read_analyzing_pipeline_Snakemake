#!/usr/bin/env python3
"""Remove duplicate fragments from an interleaved FASTQ stream."""

from __future__ import annotations

import argparse
import bz2
import gzip
import sys
from dataclasses import dataclass
from typing import Optional, TextIO, Tuple


@dataclass
class _ReaderState:
    """Keeps track of FASTQ layout (2-line vs 4-line) and buffered headers."""

    layout: Optional[int] = None
    pending_header: Optional[str] = None


def _open_text_file(path: str, mode: str) -> Tuple[TextIO, bool]:
    """Open a text file (optionally gz/bz2 compressed). Returns handle and close flag."""

    if 'b' in mode:
        raise ValueError("Binary mode not supported; use text mode like 'rt' or 'wt'.")

    if path == '-':
        if 'r' in mode:
            return sys.stdin, False
        if 'w' in mode or 'a' in mode:
            return sys.stdout, False
        raise ValueError(f"Unsupported mode '{mode}' for standard streams")

    if path.endswith('.gz'):
        return gzip.open(path, mode, encoding='utf-8'), True
    if path.endswith('.bz2'):
        return bz2.open(path, mode, encoding='utf-8'), True
    return open(path, mode, encoding='utf-8'), True


def _read_record(handle: TextIO, state: _ReaderState) -> Optional[Tuple[str, str, Optional[str], Optional[str]]]:
    """Read the next FASTQ record, respecting layout detection and buffered headers."""

    header: Optional[str] = None

    while True:
        if state.pending_header is not None:
            header = state.pending_header
            state.pending_header = None
        else:
            header_line = handle.readline()
            if not header_line:
                return None
            header = header_line.rstrip('\r\n')

        if header == '':
            # Skip blank lines quietly.
            continue
        break

    if not header.startswith('@'):
        raise ValueError(f"Malformed FASTQ: header lines must start with '@' (got: {header!r})")

    seq_line = handle.readline()
    if not seq_line:
        raise ValueError(f"Incomplete FASTQ record: missing sequence for {header!r}")
    sequence = seq_line.rstrip('\r\n')

    plus_line: Optional[str] = None
    qual_line: Optional[str] = None

    if state.layout == 4:
        plus_raw = handle.readline()
        if not plus_raw:
            raise ValueError(f"Incomplete FASTQ record: missing '+' line for {header!r}")
        if not plus_raw.startswith('+'):
            raise ValueError(f"Malformed FASTQ: expected '+' line after sequence for {header!r}")
        plus_line = plus_raw.rstrip('\r\n')

        qual_raw = handle.readline()
        if not qual_raw:
            raise ValueError(f"Incomplete FASTQ record: missing qualities for {header!r}")
        qual_line = qual_raw.rstrip('\r\n')

    elif state.layout == 2:
        # Two-line FASTQ (no quality lines) – nothing more to read.
        pass

    else:
        # Determine layout lazily based on the next line encountered.
        third_raw = handle.readline()
        if not third_raw:
            raise ValueError(f"Incomplete FASTQ record: unable to determine layout for {header!r}")
        third = third_raw.rstrip('\r\n')

        if third.startswith('+'):
            state.layout = 4
            plus_line = third
            qual_raw = handle.readline()
            if not qual_raw:
                raise ValueError(f"Incomplete FASTQ record: missing qualities for {header!r}")
            qual_line = qual_raw.rstrip('\r\n')
        else:
            state.layout = 2
            state.pending_header = third

    return header, sequence, plus_line, qual_line


def _normalise_read_name(header: str) -> str:
    """Normalise a FASTQ header to a fragment identifier (strip adornments)."""

    core = header[1:].strip()
    if not core:
        raise ValueError("Encountered empty FASTQ header")
    token = core.split()[0]
    if token.endswith('/1') or token.endswith('/2'):
        token = token[:-2]
    return token


def remove_duplicates(input_path: str, output_path: str) -> Tuple[int, int]:
    """Process an interleaved FASTQ stream and drop duplicate fragments by read name."""

    infile, close_in = _open_text_file(input_path, 'rt')
    outfile, close_out = _open_text_file(output_path, 'wt')

    state = _ReaderState()
    seen = set()
    total_pairs = 0
    duplicate_pairs = 0

    try:
        while True:
            first = _read_record(infile, state)
            if first is None:
                break

            second = _read_record(infile, state)
            if second is None:
                raise ValueError(
                    "Interleaved FASTQ terminated unexpectedly: missing mate for fragment"
                )

            h1, s1, p1, q1 = first
            h2, s2, p2, q2 = second

            name1 = _normalise_read_name(h1)
            name2 = _normalise_read_name(h2)
            if name1 != name2:
                raise ValueError(
                    f"Mismatched read names in fragment: {name1!r} vs {name2!r}. "
                    "Input must be interleaved pairs in order."
                )

            total_pairs += 1

            if name1 in seen:
                duplicate_pairs += 1
                continue
            seen.add(name1)

            outfile.write(f"{h1}\n{s1}\n")
            if state.layout == 4:
                assert p1 is not None and q1 is not None
                outfile.write(f"{p1}\n{q1}\n")

            outfile.write(f"{h2}\n{s2}\n")
            if state.layout == 4:
                assert p2 is not None and q2 is not None
                outfile.write(f"{p2}\n{q2}\n")

    finally:
        if close_in:
            infile.close()
        if close_out:
            outfile.close()

    return total_pairs, duplicate_pairs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Remove duplicate fragments (by read name) from an interleaved FASTQ file."
    )
    parser.add_argument(
        '-i',
        '--input',
        default='-',
        help="Input FASTQ path (supports .gz/.bz2). Use '-' for stdin (default).",
    )
    parser.add_argument(
        '-o',
        '--output',
        default='-',
        help="Output FASTQ path (supports .gz/.bz2). Use '-' for stdout (default).",
    )
    parser.add_argument(
        '--quiet',
        action='store_true',
        help='Suppress progress summary on stderr.',
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    total_pairs, duplicate_pairs = remove_duplicates(args.input, args.output)

    if not args.quiet:
        unique_pairs = total_pairs - duplicate_pairs
        msg = (
            f"Fragments processed: {total_pairs}. "
            f"Unique fragments written: {unique_pairs}. "
            f"Duplicates removed: {duplicate_pairs}."
        )
        print(msg, file=sys.stderr)


if __name__ == '__main__':
    main()
