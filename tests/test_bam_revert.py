import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
BAM_REVERT = REPOSITORY / "scripts" / "bam_revert.py"


REFERENCE_SEQUENCE = "A" * 200
COORDINATE_SORTED_SAM = """\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:200
@RG\tID:rg1\tSM:test
pairA\t99\tchr1\t1\t60\t4M2H\t=\t100\t103\tACGT\tIIII\tRG:Z:rg1\tZB:Z:TT\tZQ:Z:##\tXA:B:i,1,2
pairB\t99\tchr1\t10\t60\t4M\t=\t90\t84\tCCCC\tJJJJ\tRG:Z:rg1
pairB\t147\tchr1\t90\t60\t4M\t=\t10\t-84\tGGGG\tHHHH\tRG:Z:rg1
pairA\t147\tchr1\t100\t60\t2H4M\t=\t1\t-103\tGGTT\tMLKJ\tRG:Z:rg1\tYB:Z:CC\tYQ:Z:ON\tXB:H:ABCD
"""


def read_fastq(path):
    records = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        while True:
            header = handle.readline().rstrip("\n")
            if not header:
                break
            sequence = handle.readline().rstrip("\n")
            plus = handle.readline().rstrip("\n")
            quality = handle.readline().rstrip("\n")
            if plus != "+":
                raise AssertionError(f"Malformed FASTQ plus line: {plus!r}")
            records.append((header, sequence, quality))
    return sorted(records)


class BamRevertTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.samtools = shutil.which("samtools")
        if cls.samtools is None:
            raise unittest.SkipTest("samtools is required")

    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory(prefix="bam-revert-test-")
        self.directory = Path(self.temporary_directory.name)
        self.reference = self.directory / "reference.fa"
        self.reference.write_text(f">chr1\n{REFERENCE_SEQUENCE}\n", encoding="ascii")
        subprocess.run([self.samtools, "faidx", self.reference], check=True)
        self.sam = self.directory / "coordinate_sorted.sam"
        self.sam.write_text(COORDINATE_SORTED_SAM, encoding="ascii")

    def tearDown(self):
        self.temporary_directory.cleanup()

    def outputs(self, prefix):
        return (
            self.directory / f"{prefix}_R1.fastq.gz",
            self.directory / f"{prefix}_R2.fastq.gz",
            self.directory / f"{prefix}.stats.tsv",
        )

    def run_revert(
        self, input_path, prefix, reference=None, check=True, include_read_group=False
    ):
        fastq1, fastq2, stats = self.outputs(prefix)
        command = [
            sys.executable,
            BAM_REVERT,
            "--input",
            str(input_path),
            "--f1",
            str(fastq1),
            "--f2",
            str(fastq2),
            "--stats",
            str(stats),
            "--threads",
            "1",
            "--samtools",
            self.samtools,
        ]
        if reference is not None:
            command.extend(["--reference", str(reference)])
        if include_read_group:
            command.append("--include-read-group-in-name")
        result = subprocess.run(command, check=check, text=True, capture_output=True)
        return result, fastq1, fastq2, stats

    def assert_restored_payload(self, fastq1, fastq2, stats):
        self.assertEqual(
            read_fastq(fastq1),
            sorted(
                [
                    ("@pairA/1", "ACGTTT", "IIII##"),
                    ("@pairB/1", "CCCC", "JJJJ"),
                ]
            ),
        )
        self.assertEqual(
            read_fastq(fastq2),
            sorted(
                [
                    ("@pairA/2", "AACCGG", "JKLMNO"),
                    ("@pairB/2", "CCCC", "HHHH"),
                ]
            ),
        )
        self.assertEqual(
            stats.read_text(encoding="utf-8"),
            "alignment_counter\t 4\nfragment_counter\t 2\n",
        )

    def test_coordinate_sorted_sam_is_collated_automatically(self):
        _, fastq1, fastq2, stats = self.run_revert(self.sam, "sam")
        self.assert_restored_payload(fastq1, fastq2, stats)

    def test_script_runs_without_repository_python_modules(self):
        standalone = self.directory / "standalone" / "bam_revert.py"
        standalone.parent.mkdir()
        shutil.copy2(BAM_REVERT, standalone)
        original = globals()["BAM_REVERT"]
        try:
            globals()["BAM_REVERT"] = standalone
            _, fastq1, fastq2, stats = self.run_revert(self.sam, "standalone")
        finally:
            globals()["BAM_REVERT"] = original
        self.assert_restored_payload(fastq1, fastq2, stats)

    def test_cram_is_decoded_with_explicit_reference(self):
        cram = self.directory / "reads.cram"
        subprocess.run(
            [self.samtools, "view", "-C", "--reference", self.reference, "-o", cram, self.sam],
            check=True,
        )
        _, fastq1, fastq2, stats = self.run_revert(cram, "cram", reference=self.reference)
        self.assert_restored_payload(fastq1, fastq2, stats)

    def test_samtools_failure_does_not_replace_existing_outputs(self):
        fastq1, fastq2, stats = self.outputs("failure")
        sentinels = {
            fastq1: b"existing R1\n",
            fastq2: b"existing R2\n",
            stats: b"existing stats\n",
        }
        for path, content in sentinels.items():
            path.write_bytes(content)

        result, _, _, _ = self.run_revert(self.directory / "missing.cram", "failure", check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("bam_revert: error:", result.stderr)
        for path, content in sentinels.items():
            self.assertEqual(path.read_bytes(), content)
        self.assertEqual(list(self.directory.glob(".*.tmp")), [])

    def test_duplicate_query_names_are_separated_by_read_group(self):
        multi_group_sam = self.directory / "multiple_read_groups.sam"
        multi_group_sam.write_text(
            """\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:200
@RG\tID:rg1\tSM:test
@RG\tID:rg2\tSM:test
duplicate\t99\tchr1\t1\t60\t4M\t=\t100\t103\tAAAA\tIIII\tRG:Z:rg1
duplicate\t99\tchr1\t2\t60\t4M\t=\t101\t103\tCCCC\tJJJJ\tRG:Z:rg2
duplicate\t147\tchr1\t100\t60\t4M\t=\t1\t-103\tTTTT\tHHHH\tRG:Z:rg1
duplicate\t147\tchr1\t101\t60\t4M\t=\t2\t-103\tGGGG\tGGGG\tRG:Z:rg2
""",
            encoding="ascii",
        )
        _, fastq1, fastq2, stats = self.run_revert(
            multi_group_sam,
            "multiple_read_groups",
            include_read_group=True,
        )
        self.assertEqual(
            read_fastq(fastq1),
            sorted(
                [
                    ("@duplicate:rg1/1", "AAAA", "IIII"),
                    ("@duplicate:rg2/1", "CCCC", "JJJJ"),
                ]
            ),
        )
        self.assertEqual(
            read_fastq(fastq2),
            sorted(
                [
                    ("@duplicate:rg1/2", "AAAA", "HHHH"),
                    ("@duplicate:rg2/2", "CCCC", "GGGG"),
                ]
            ),
        )
        self.assertEqual(
            stats.read_text(encoding="utf-8"),
            "alignment_counter\t 4\nfragment_counter\t 2\n",
        )


if __name__ == "__main__":
    unittest.main()
