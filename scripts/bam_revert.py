#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import subprocess
import sys
import tempfile
import time
from collections import OrderedDict
from contextlib import contextmanager
from pathlib import Path

RESTORE_TAGS = {"RG", "YB", "YQ", "ZB", "ZQ"}
_DNA_COMPLEMENT = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn",
    "TGCAYRKMVHDBNtgcayrkmvhdbn",
)


class BamRecord:
    """Minimal SAM record implementation required by this standalone tool."""

    def __init__(self, row):
        if len(row) < 11:
            raise ValueError("SAM record has fewer than 11 fields")
        self.qname = row[0][:-2] if row[0].endswith(("/1", "/2")) else row[0]
        self.flag = int(row[1])
        self.seq = row[9]
        self.qual = row[10]
        self.tags = {}
        for tag_field in row[11:]:
            fields = tag_field.split(":", 2)
            if len(fields) != 3:
                continue
            tag, value_type, value = fields
            if value_type == "i":
                value = int(value)
            elif value_type == "f":
                value = float(value)
            self.tags[tag] = value

    def getTagValue(self, tag, default=None):
        return self.tags.get(tag, default)

    def is_primary(self):
        return not bool(self.flag & (0x100 | 0x800))

    def is_unmapped(self):
        return bool(self.flag & 0x4)

    def restored_payload(self):
        """Restore hard-clipped bases/qualities and the original orientation."""
        has_restore_tags = any(tag in self.tags for tag in ("YB", "YQ", "ZB", "ZQ"))
        if self.is_unmapped() and not has_restore_tags:
            return self.seq, self.qual

        sequence = (
            str(self.getTagValue("YB", ""))
            + self.seq
            + str(self.getTagValue("ZB", ""))
        )
        quality = (
            str(self.getTagValue("YQ", ""))
            + self.qual
            + str(self.getTagValue("ZQ", ""))
        )
        if self.flag & 0x10:
            sequence = sequence.translate(_DNA_COMPLEMENT)[::-1]
            quality = quality[::-1]
        return sequence, quality


class BamRevertError(RuntimeError):
    pass


def parse_record(row):
    if len(row) < 11:
        raise BamRevertError(f"Malformed SAM record with {len(row)} fields")

    # BamRecord only needs these tags for FASTQ reconstruction. Ignoring all
    # other tags also avoids rejecting valid B/H auxiliary tag types that are
    # irrelevant to reversion.
    restore_fields = []
    for field in row[11:]:
        parts = field.split(":", 2)
        if len(parts) == 3 and parts[0] in RESTORE_TAGS:
            restore_fields.append(field)
    return BamRecord(row[:11] + restore_fields)


def primary_pairs(querygroup, stats):
    if not querygroup:
        return []

    stats["alignment_counter"] += len(querygroup)
    primaries_by_read_group = OrderedDict()
    for record in querygroup:
        if not record.is_primary():
            continue
        read_group = str(record.getTagValue("RG", ""))
        primaries_by_read_group.setdefault(read_group, []).append(record)

    if not primaries_by_read_group:
        raise BamRevertError(f"Read {querygroup[0].qname!r} has no primary alignments")

    pairs = []
    for read_group, primaries in primaries_by_read_group.items():
        reads1 = [record for record in primaries if record.flag & 0x40 and not record.flag & 0x80]
        reads2 = [record for record in primaries if record.flag & 0x80 and not record.flag & 0x40]
        label = f"read {querygroup[0].qname!r}"
        if read_group:
            label += f" in read group {read_group!r}"
        if len(reads1) != 1 or len(reads2) != 1:
            raise BamRevertError(
                f"{label} must have exactly one primary R1 and R2; "
                f"found R1={len(reads1)}, R2={len(reads2)}"
            )
        pairs.append((reads1[0], reads2[0], read_group))
        stats["fragment_counter"] += 1
    return pairs


def to_fastq_record(record, read_number, read_group=None):
    sequence, quality = record.restored_payload()
    if sequence == "*" or quality == "*":
        raise BamRevertError(f"Read {record.qname!r}/{read_number} has missing sequence or quality")
    if len(sequence) != len(quality):
        raise BamRevertError(
            f"Read {record.qname!r}/{read_number} has sequence length {len(sequence)} "
            f"but quality length {len(quality)}"
        )
    name = record.qname
    if read_group is not None:
        name = f"{name}:{read_group}"
    return f"@{name}/{read_number}\n{sequence}\n+\n{quality}\n"


class SamtoolsSamReader:
    def __init__(self, args):
        self.args = args
        self.processes = []
        self.stream = None

    def __enter__(self):
        input_stream = None
        if self.args.input == "-":
            input_stream = getattr(sys.stdin, "buffer", sys.stdin)

        try:
            if self.args.name_collated:
                command = [self.args.samtools, "view", "-h", "--threads", str(self.args.threads)]
                if self.args.reference:
                    command.extend(["--reference", self.args.reference])
                command.append(self.args.input)
                view_process = subprocess.Popen(
                    command,
                    stdin=input_stream,
                    stdout=subprocess.PIPE,
                    text=True,
                )
                self.processes.append(("samtools view", view_process))
            else:
                collate_command = [
                    self.args.samtools,
                    "collate",
                    "-@",
                    str(self.args.threads),
                    "-Ou",
                ]
                if self.args.reference:
                    collate_command.extend(["--reference", self.args.reference])
                collate_command.append(self.args.input)
                collate_process = subprocess.Popen(
                    collate_command,
                    stdin=input_stream,
                    stdout=subprocess.PIPE,
                )
                self.processes.append(("samtools collate", collate_process))

                view_process = subprocess.Popen(
                    [self.args.samtools, "view", "-h", "-"],
                    stdin=collate_process.stdout,
                    stdout=subprocess.PIPE,
                    text=True,
                )
                collate_process.stdout.close()
                self.processes.append(("samtools view", view_process))
        except Exception:
            self._terminate()
            raise

        self.stream = view_process.stdout
        return self.stream

    def _terminate(self):
        if self.stream is not None:
            self.stream.close()
        for _, process in reversed(self.processes):
            if process.poll() is None:
                process.terminate()
        for _, process in reversed(self.processes):
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is not None:
            self._terminate()
            return False

        if self.stream is not None:
            self.stream.close()
        failures = []
        for name, process in reversed(self.processes):
            return_code = process.wait()
            if return_code != 0:
                failures.append(f"{name} exited with status {return_code}")
        if failures:
            raise BamRevertError("; ".join(failures))
        return False


@contextmanager
def atomic_output_paths(targets):
    absolute_targets = [Path(target).expanduser().absolute() for target in targets]
    if len(set(absolute_targets)) != len(absolute_targets):
        raise BamRevertError("FASTQ and statistics outputs must be different paths")

    temporary_paths = []
    try:
        for target in absolute_targets:
            if not target.parent.is_dir():
                raise BamRevertError(f"Output directory does not exist: {target.parent}")
            descriptor, temporary_name = tempfile.mkstemp(
                prefix=f".{target.name}.", suffix=".tmp", dir=target.parent
            )
            os.close(descriptor)
            temporary_paths.append(Path(temporary_name))
        yield temporary_paths
        for temporary_path, target in zip(temporary_paths, absolute_targets):
            os.replace(temporary_path, target)
    finally:
        for temporary_path in temporary_paths:
            try:
                temporary_path.unlink()
            except FileNotFoundError:
                pass


def revert_stream(stream, out_file1, out_file2, stats, include_read_group=False):
    reader = csv.reader(stream, delimiter="\t", quoting=csv.QUOTE_NONE)
    current_name = None
    querygroup = []
    rows_seen = 0
    last_report = time.monotonic()

    def write_group(group):
        for read1, read2, read_group in primary_pairs(group, stats):
            output_read_group = read_group if include_read_group and read_group else None
            out_file1.write(to_fastq_record(read1, 1, output_read_group))
            out_file2.write(to_fastq_record(read2, 2, output_read_group))

    for row in reader:
        if not row:
            continue
        if row[0].startswith("@"):
            continue
        rows_seen += 1
        if rows_seen % 100000 == 0:
            now = time.monotonic()
            elapsed = max(now - last_report, 1e-9)
            print(
                f"{rows_seen} alignments read at {int(100000 / elapsed)} alignments/sec; "
                f"{stats['fragment_counter']} fragments written",
                file=sys.stderr,
            )
            last_report = now

        record = parse_record(row)
        raw_name = row[0]
        if current_name is None or raw_name == current_name:
            querygroup.append(record)
        else:
            write_group(querygroup)
            querygroup = [record]
        current_name = raw_name

    if querygroup:
        write_group(querygroup)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Restore paired FASTQ payloads from a SAM, BAM, or CRAM file."
    )
    parser.add_argument("-i", "--input", default="-", help="input SAM/BAM/CRAM (default: stdin)")
    parser.add_argument("--f1", required=True, help="output R1 FASTQ.gz")
    parser.add_argument("--f2", required=True, help="output R2 FASTQ.gz")
    parser.add_argument("-s", "--stats", required=True, help="output statistics TSV")
    parser.add_argument("-T", "--reference", help="reference FASTA used to decode CRAM input")
    parser.add_argument("--threads", type=int, default=2, help="samtools threads (default: 2)")
    parser.add_argument("--samtools", default="samtools", help="samtools executable")
    parser.add_argument(
        "--name-collated",
        action="store_true",
        help="skip automatic samtools collate when input is already grouped by query name",
    )
    parser.add_argument(
        "--include-read-group-in-name",
        action="store_true",
        help="append the RG value to reconstructed FASTQ names",
    )
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    return args


def run(args):
    output_paths = [Path(path).expanduser().absolute() for path in (args.f1, args.f2, args.stats)]
    protected_paths = []
    if args.input != "-":
        protected_paths.append(Path(args.input).expanduser().absolute())
    if args.reference:
        protected_paths.append(Path(args.reference).expanduser().absolute())
    overlap = set(output_paths).intersection(protected_paths)
    if overlap:
        raise BamRevertError(
            "Output paths must not overwrite the input or reference: "
            + ", ".join(str(path) for path in sorted(overlap))
        )

    stats = {"alignment_counter": 0, "fragment_counter": 0}
    with atomic_output_paths([args.f1, args.f2, args.stats]) as temporary_paths:
        temporary_f1, temporary_f2, temporary_stats = temporary_paths
        with gzip.open(temporary_f1, "wt", encoding="utf-8", newline="\n") as out_file1, gzip.open(
            temporary_f2, "wt", encoding="utf-8", newline="\n"
        ) as out_file2:
            with SamtoolsSamReader(args) as stream:
                revert_stream(
                    stream,
                    out_file1,
                    out_file2,
                    stats,
                    include_read_group=args.include_read_group_in_name,
                )

        with open(temporary_stats, "wt", encoding="utf-8", newline="\n") as stats_file:
            for key in ("alignment_counter", "fragment_counter"):
                stats_file.write(f"{key}\t {stats[key]}\n")
    print(
        f"Done: restored {stats['fragment_counter']} read pairs from "
        f"{stats['alignment_counter']} alignments",
        file=sys.stderr,
    )


def main(argv=None):
    args = parse_args(argv)
    try:
        run(args)
    except (BamRevertError, OSError, subprocess.SubprocessError) as error:
        print(f"bam_revert: error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
