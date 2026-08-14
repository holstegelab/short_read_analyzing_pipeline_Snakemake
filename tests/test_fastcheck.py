import gzip
import importlib.util
import os
import shutil
import subprocess
import sys
import sysconfig
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPTS = REPOSITORY / "scripts"
sys.path.insert(0, str(SCRIPTS))
from fastcheck_loader import _local_extension_is_stale


FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3


def fnv1a64(value):
    result = FNV_OFFSET
    for byte in value.encode("ascii"):
        result ^= byte
        result = (result * FNV_PRIME) & 0xFFFFFFFFFFFFFFFF
    return result


def normalized_quality(sequence, quality):
    return "".join(
        "!" if base in "Nn" else ("!" if score == "#" else score)
        for base, score in zip(sequence, quality)
    )


class FastcheckPairedFastqTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if shutil.which("pigz") is None:
            raise unittest.SkipTest("pigz is required by fastcheck.fastq_stats")

        cls.build_dir = tempfile.TemporaryDirectory(prefix="fastcheck-test-")
        build_path = Path(cls.build_dir.name)
        shutil.copy2(SCRIPTS / "fastcheck.c", build_path / "fastcheck.c")
        shutil.copy2(SCRIPTS / "setup.py", build_path / "setup.py")
        subprocess.run(
            [sys.executable, "setup.py", "build_ext", "--inplace"],
            cwd=build_path,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        extension_path = next(build_path.glob("fastcheck*.so"))
        spec = importlib.util.spec_from_file_location("fastcheck", extension_path)
        cls.fastcheck = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.fastcheck)

    @classmethod
    def tearDownClass(cls):
        cls.build_dir.cleanup()

    @staticmethod
    def records(count):
        records = []
        quality_alphabet = "!$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI"
        for index in range(count):
            name = f"read{index:06d}"
            sequence1 = "ACGT"
            sequence2 = "TGCA"
            quality1 = "I" * 4
            marker = quality_alphabet[index % len(quality_alphabet)]
            quality2 = marker + "HGF"
            records.append((name, sequence1, quality1, sequence2, quality2))
        return records

    @staticmethod
    def write_fastqs(directory, records):
        fastq1 = directory / "reads_R1.fastq.gz"
        fastq2 = directory / "reads_R2.fastq.gz"
        with gzip.open(fastq1, "wt", encoding="ascii", newline="\n") as out1, gzip.open(
            fastq2, "wt", encoding="ascii", newline="\n"
        ) as out2:
            for name, sequence1, quality1, sequence2, quality2 in records:
                out1.write(f"@{name}/1\n{sequence1}\n+\n{quality1}\n")
                out2.write(f"@{name}/2\n{sequence2}\n+\n{quality2}\n")
        return fastq1, fastq2

    @staticmethod
    def expected(records):
        seq1_hash = qual1_hash = seq2_hash = qual2_hash = 0
        bases1 = bases2 = 0
        for _, sequence1, quality1, sequence2, quality2 in records:
            seq1_hash ^= fnv1a64(sequence1)
            qual1_hash ^= fnv1a64(normalized_quality(sequence1, quality1))
            seq2_hash ^= fnv1a64(sequence2)
            qual2_hash ^= fnv1a64(normalized_quality(sequence2, quality2))
            bases1 += len(sequence1)
            bases2 += len(sequence2)
        return seq1_hash, qual1_hash, seq2_hash, qual2_hash, len(records), bases1, bases2

    def run_fastcheck(self, records):
        with tempfile.TemporaryDirectory(prefix="fastcheck-input-") as directory:
            fastq1, fastq2 = self.write_fastqs(Path(directory), records)
            return self.fastcheck.fastq_stats(str(fastq1), str(fastq2))

    def test_quality_buffer_boundaries_include_every_r2_record(self):
        for count in (1999, 2000, 2001):
            with self.subTest(count=count):
                records = self.records(count)
                self.assertEqual(self.run_fastcheck(records), self.expected(records))

    def test_checksums_do_not_depend_on_which_record_is_at_buffer_boundary(self):
        records = self.records(2001)
        reordered = list(reversed(records))
        self.assertEqual(self.run_fastcheck(records), self.run_fastcheck(reordered))


class FastcheckLoaderTests(unittest.TestCase):
    def test_local_extension_is_rebuilt_when_source_is_newer(self):
        extension_suffix = sysconfig.get_config_var("EXT_SUFFIX")
        if not extension_suffix:
            self.skipTest("Python does not expose an extension suffix")
        with tempfile.TemporaryDirectory(prefix="fastcheck-loader-test-") as directory:
            directory = Path(directory)
            source = directory / "fastcheck.c"
            extension = directory / f"fastcheck{extension_suffix}"
            source.write_text("/* source */\n", encoding="ascii")
            extension.write_bytes(b"extension")
            os.utime(extension, (1, 1))
            os.utime(source, (2, 2))
            self.assertTrue(_local_extension_is_stale(str(directory)))

            os.utime(extension, (3, 3))
            self.assertFalse(_local_extension_is_stale(str(directory)))


if __name__ == "__main__":
    unittest.main()
