"""Regression tests for paired FASTQ handling and SAM conversion."""

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "R_scripts"))
from extract_read_tags import parse_sam  # noqa: E402
from sam_to_tables import convert  # noqa: E402
from select_read_pairs import select  # noqa: E402


def fastq(name, sequence):
    return f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n"


class HelperTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_gzipped_pairs_stay_synchronized_and_adapter_is_required(self):
        sam = self.root / "r1.sam"
        sam.write_text("@HD\tVN:1.6\nread_a/1\t0\ttranscript\t1\t42\t100M\t*\t0\t0\tAC\tII\n"
                       "read_b/1\t0\ttranscript\t1\t42\t100M\t*\t0\t0\tAC\tII\n")
        r1 = self.root / "r1.fq.gz"
        r2 = self.root / "r2.fq.gz"
        with gzip.open(r1, "wt") as handle:
            handle.write(fastq("read_a/1", "ACGT"))
            handle.write(fastq("read_b/1", "ACGT"))
        with gzip.open(r2, "wt") as handle:
            handle.write(fastq("read_a/2", "A" * 15 + "GTCAG" + "T" * 80))
            handle.write(fastq("read_b/2", "A" * 15 + "AAAAA" + "T" * 80))
        out1, out2 = self.root / "selected1.fq", self.root / "selected2.fq"
        select(sam, str(r1), str(r2), out1, out2)
        self.assertEqual(out1.read_text().count("@"), 1)
        self.assertTrue(out2.read_text().startswith("@read_a/2\n"))

    def test_mismatched_pair_names_fail(self):
        sam = self.root / "r1.sam"
        sam.write_text("read_a\t0\ttranscript\t1\t42\t4M\t*\t0\t0\tACGT\tIIII\n")
        r1, r2 = self.root / "r1.fq", self.root / "r2.fq"
        r1.write_text(fastq("read_a", "ACGT"))
        r2.write_text(fastq("wrong_name", "A" * 15 + "GTCAG"))
        with self.assertRaisesRegex(ValueError, "out of sync"):
            select(sam, str(r1), str(r2), self.root / "o1", self.root / "o2")

    def test_sam_tables_normalize_mate_suffix_around_umi(self):
        r1, r2 = self.root / "r1.sam", self.root / "r2.sam"
        r1.write_text("read_a/1\t0\tENST1|ENSG1|demo\t1\t42\t4M\t*\t0\t0\tACGT\tIIII\n")
        r2.write_text("read_a/2:AAAAAAAAAAAAAAA\t4\t*\t0\t0\t*\t*\t0\t0\tTTTT\tIIII\n")
        o1, o2 = self.root / "o1.tsv", self.root / "o2.tsv"
        convert(r1, r2, o1, o2)
        self.assertTrue(o1.read_text().startswith("read_a\tENSG1\tENST1\t"))
        self.assertTrue(o2.read_text().startswith("read_a:AAAAAAAAAAAAAAA\t4\t"))

    def test_genome_hit_labels_ignore_unmapped_records(self):
        sam, out = self.root / "genome.sam", self.root / "tags.tsv"
        sam.write_text("single\t0\tchr1\t1\n"
                       "multi\t0\tchr1\t1\n"
                       "multi\t256\tchr2\t1\n"
                       "unmapped\t4\t*\t0\n")
        parse_sam(sam, out)
        self.assertEqual(out.read_text(), "single\tunique\nmulti\tmulti\n")


if __name__ == "__main__":
    unittest.main()
