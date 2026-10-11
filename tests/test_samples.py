"""Checks for multi-sample FASTQ pairing and manifest validation."""

import tempfile
import unittest
from pathlib import Path

from scripts.make_samplesheet import discover_pairs, write_sheet
from scripts.sample_sheet import load_samples


class SampleSheetTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.r1 = self.root / "R1"
        self.r2 = self.root / "R2"
        self.r1.mkdir()
        self.r2.mkdir()

    def test_two_folders_pair_by_exact_shared_id(self):
        for sample in ("patient7", "patient8"):
            (self.r1 / f"{sample}_R1_001.fastq.gz").touch()
            (self.r2 / f"{sample}_R2_001.fastq.gz").touch()
        pairs = discover_pairs(self.r1, self.r2)
        output = self.root / "samples.tsv"
        write_sheet(pairs, output)
        config = {"samples": str(output), "sample": "legacy", "R1": "old1", "R2": "old2"}
        loaded = load_samples(config)
        self.assertEqual(list(loaded), ["patient7", "patient8"])
        self.assertTrue(loaded["patient7"]["R1"].endswith("patient7_R1_001.fastq.gz"))
        with self.assertRaises(FileExistsError):
            write_sheet(pairs, output)

    def test_user_filename_pattern_in_one_folder(self):
        (self.r1 / "may2026_PB25_2_R1.fastq.gz").touch()
        (self.r1 / "may2026_PB25_2_R2.fastq.gz").touch()
        pairs = discover_pairs(self.r1, self.r1)
        self.assertEqual(list(pairs), ["may2026_PB25_2"])
        self.assertEqual(pairs["may2026_PB25_2"][0].name, "may2026_PB25_2_R1.fastq.gz")
        self.assertEqual(pairs["may2026_PB25_2"][1].name, "may2026_PB25_2_R2.fastq.gz")

    def test_partial_overlap_does_not_pair_wrong_samples(self):
        (self.r1 / "patient7_R1.fq").touch()
        (self.r2 / "patient70_R2.fq").touch()
        with self.assertRaisesRegex(ValueError, "Unpaired samples"):
            discover_pairs(self.r1, self.r2)

    def test_ambiguous_duplicate_id_is_rejected(self):
        (self.r1 / "patient7_R1.fq").touch()
        (self.r1 / "patient7_R1_001.fq").touch()
        (self.r2 / "patient7_R2.fq").touch()
        with self.assertRaisesRegex(ValueError, "Multiple R1 FASTQs"):
            discover_pairs(self.r1, self.r2)

    def test_custom_regex_extracts_shared_part(self):
        (self.r1 / "patient7_read1_runA.fastq.gz").touch()
        (self.r2 / "patient7_read2_runA.fastq.gz").touch()
        pairs = discover_pairs(
            self.r1, self.r2,
            r"^(?P<sample>.+)_read1_runA\.fastq\.gz$",
            r"^(?P<sample>.+)_read2_runA\.fastq\.gz$",
        )
        self.assertEqual(list(pairs), ["patient7"])

    def test_duplicate_manifest_row_is_rejected(self):
        sheet = self.root / "samples.tsv"
        sheet.write_text("sample\tR1\tR2\np1\ta.fq\tb.fq\np1\tc.fq\td.fq\n")
        with self.assertRaisesRegex(ValueError, "duplicate sample"):
            load_samples({"samples": str(sheet)})

    def test_legacy_single_sample_still_works(self):
        self.assertEqual(
            load_samples({"sample": "p1", "R1": "a.fq", "R2": "b.fq"}),
            {"p1": {"R1": "a.fq", "R2": "b.fq"}},
        )


if __name__ == "__main__":
    unittest.main()
