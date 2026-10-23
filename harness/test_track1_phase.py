"""Verify the pair phase reader and the pair level rates on hand-written fixtures.

Every record of this test is written by hand. The contig chrT does not exist, so no fixture can be
mistaken for a real sample. The test reads no alignment, no release file and no patient data.
"""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "pipelines" / "phase_pair" / "bin"))

import phase_pairs  # noqa: E402
import track1_phase_control as control  # noqa: E402

HEADER = ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")


def vcf_line(pos: int, ref: str, alt: str, genotype: str, phase_set: str) -> str:
    """Build one VCF data line with a phase set.

    Args:
        pos: The one-based position.
        ref: The REF allele.
        alt: The ALT allele.
        genotype: The GT value.
        phase_set: The PS value.

    Returns:
        The line with its newline.
    """
    return (f"chrT\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT:PS\t{genotype}:{phase_set}\n")


def write_vcf(path: Path, lines: list[str]) -> None:
    """Write a small VCF.

    Args:
        path: The output path.
        lines: The data lines.
    """
    path.write_text("##fileformat=VCFv4.2\n" + HEADER + "".join(lines))


def write_submission(path: Path, pairs: list[tuple[int, int]]) -> None:
    """Write a submission file with one row per pair.

    Args:
        path: The output path.
        pairs: The position pairs.
    """
    fields = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2",
              "alt_2", "epcr", "finding_type", "notes"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for rank, (first, second) in enumerate(pairs, start=1):
            writer.writerow({"proband_id": "PROBAND01", "chrom_1": "chrT", "pos_1": first,
                             "ref_1": "A", "alt_1": "T", "chrom_2": "chrT", "pos_2": second,
                             "ref_2": "A", "alt_2": "T", "epcr": f"{1.0 - 0.1 * rank:.2f}",
                             "finding_type": "primary", "notes": "fixture"})


class HaplotypeTest(unittest.TestCase):
    """The haplotype reader accepts a phased heterozygote alone."""

    def test_phased_heterozygote(self) -> None:
        """A phased heterozygote names the haplotype of the alternate allele."""
        self.assertEqual(phase_pairs.haplotype_of_alt("1|0"), 0)
        self.assertEqual(phase_pairs.haplotype_of_alt("0|1"), 1)

    def test_unphased_or_homozygous(self) -> None:
        """An unphased or a homozygous genotype returns None."""
        for genotype in ("0/1", "1/1", "1|1", "0|0", "."):
            self.assertIsNone(phase_pairs.haplotype_of_alt(genotype))

    def test_the_control_agrees_with_the_reader(self) -> None:
        """The validation script and the pipeline reader use the same rule."""
        for genotype in ("1|0", "0|1", "0/1", "1|1"):
            self.assertEqual(phase_pairs.haplotype_of_alt(genotype),
                             control.haplotype_of_alt(genotype))


class ReaderTest(unittest.TestCase):
    """The reader states the phase of every submitted pair."""

    def setUp(self) -> None:
        """Build one phased VCF and one submission file."""
        self.dir = Path(tempfile.mkdtemp())
        write_vcf(self.dir / "phased.vcf", [
            vcf_line(100, "A", "T", "1|0", "100"),
            vcf_line(250, "A", "T", "0|1", "100"),
            vcf_line(400, "A", "T", "1|0", "100"),
            vcf_line(9000, "A", "T", "0/1", "."),
        ])
        write_submission(self.dir / "sub.csv", [(100, 250), (100, 400), (100, 9000), (250, 8000)])
        self.results = [phase_pairs.evaluate_pair(pair,
                                                  phase_pairs.read_phased(self.dir / "phased.vcf"))
                        for pair in phase_pairs.submitted_pairs(self.dir / "sub.csv")]

    def test_opposite_haplotypes_are_trans(self) -> None:
        """Two alternate alleles on opposite haplotypes give trans."""
        self.assertEqual(self.results[0]["phase_state"], "PHASED")
        self.assertEqual(self.results[0]["orientation"], "trans")
        self.assertTrue(self.results[0]["explains_a_recessive_disease"])

    def test_one_haplotype_is_cis(self) -> None:
        """Two alternate alleles on one haplotype give cis and explain no recessive disease."""
        self.assertEqual(self.results[1]["orientation"], "cis")
        self.assertFalse(self.results[1]["explains_a_recessive_disease"])

    def test_an_unphased_variant_gives_no_orientation(self) -> None:
        """A variant with no phase set leaves the pair unphased."""
        self.assertEqual(self.results[2]["phase_state"], "NOT PHASED")
        self.assertIsNone(self.results[2]["orientation"])

    def test_a_variant_outside_the_input_is_named(self) -> None:
        """A submitted variant that the phasing input does not hold is reported as such."""
        self.assertEqual(self.results[3]["phase_state"], "NOT IN THE INPUT")
        self.assertFalse(self.results[3]["second_in_the_phasing_input"])

    def test_the_band_follows_the_distance(self) -> None:
        """The distance band of each pair follows the position difference."""
        self.assertEqual(self.results[0]["band"], "0_to_300")
        self.assertEqual(self.results[2]["band"], "1001_to_10000")


class ClassAGateTest(unittest.TestCase):
    """The class A gate file carries no result of the sample."""

    def test_holds_only_the_fixed_keys(self) -> None:
        """The record holds the fixed key list and nothing else."""
        work = Path(tempfile.mkdtemp())
        write_vcf(work / "phased.vcf", [vcf_line(100, "A", "T", "1|0", "100"),
                                        vcf_line(250, "A", "T", "0|1", "100")])
        write_submission(work / "sub.csv", [(100, 250)])
        out = work / "out"
        out.mkdir()
        results = [phase_pairs.evaluate_pair(pair, phase_pairs.read_phased(work / "phased.vcf"))
                   for pair in phase_pairs.submitted_pairs(work / "sub.csv")]
        table, note = phase_pairs.write_rows(results, out)
        record = phase_pairs.build_gates("fixture_v1", work / "phased.vcf", work / "sub.csv",
                                         table, note, {"whatshap": "2.8"})
        self.assertEqual(tuple(record), phase_pairs.CLASS_A_KEYS)

    def test_holds_no_count_and_no_coordinate(self) -> None:
        """No value of the record is a count, a coordinate, a gene or a genotype."""
        work = Path(tempfile.mkdtemp())
        write_vcf(work / "phased.vcf", [vcf_line(100, "A", "T", "1|0", "100"),
                                        vcf_line(250, "A", "T", "0|1", "100")])
        write_submission(work / "sub.csv", [(100, 250)])
        out = work / "out"
        out.mkdir()
        results = [phase_pairs.evaluate_pair(pair, phase_pairs.read_phased(work / "phased.vcf"))
                   for pair in phase_pairs.submitted_pairs(work / "sub.csv")]
        table, note = phase_pairs.write_rows(results, out)
        record = phase_pairs.build_gates("fixture_v1", work / "phased.vcf", work / "sub.csv",
                                         table, note, {"whatshap": "2.8"})
        text = json.dumps(record)
        for forbidden in ("chrT", "trans", "cis", "PHASED", "\"pairs", "orientation"):
            self.assertNotIn(forbidden, text)


class PairRateTest(unittest.TestCase):
    """The pair level rates count the right pairs."""

    def test_a_pair_needs_a_truth_genotype_on_both_sides(self) -> None:
        """A call with no truth heterozygote produces no pair."""
        calls = [{"chrom": "chrT", "pos": 100, "ref": "A", "alt": "T", "interval": "GENE",
                  "ps": "100", "haplotype": 0},
                 {"chrom": "chrT", "pos": 250, "ref": "A", "alt": "T", "interval": "GENE",
                  "ps": "100", "haplotype": 1}]
        truth = {("chrT", 100, "A", "T"): 0}
        self.assertEqual(control.build_pairs(calls, truth), [])

    def test_the_orientation_is_invariant_to_a_block_flip(self) -> None:
        """A flip of both haplotypes keeps the same orientation."""
        calls = [{"chrom": "chrT", "pos": 100, "ref": "A", "alt": "T", "interval": "GENE",
                  "ps": "100", "haplotype": 0},
                 {"chrom": "chrT", "pos": 250, "ref": "A", "alt": "T", "interval": "GENE",
                  "ps": "100", "haplotype": 1}]
        flipped = [call | {"haplotype": 1 - call["haplotype"]} for call in calls]
        truth = {("chrT", 100, "A", "T"): 0, ("chrT", 250, "A", "T"): 1}
        first = control.build_pairs(calls, truth)[0]
        second = control.build_pairs(flipped, truth)[0]
        self.assertEqual(first["called_orientation"], "trans")
        self.assertEqual(second["called_orientation"], "trans")
        self.assertTrue(first["agrees"])
        self.assertTrue(second["agrees"])

    def test_two_blocks_leave_the_pair_unphased(self) -> None:
        """Two phase sets leave the pair without a called orientation."""
        calls = [{"chrom": "chrT", "pos": 100, "ref": "A", "alt": "T", "interval": "GENE",
                  "ps": "100", "haplotype": 0},
                 {"chrom": "chrT", "pos": 250, "ref": "A", "alt": "T", "interval": "GENE",
                  "ps": "200", "haplotype": 1}]
        truth = {("chrT", 100, "A", "T"): 0, ("chrT", 250, "A", "T"): 1}
        pair = control.build_pairs(calls, truth)[0]
        self.assertFalse(pair["phased"])
        self.assertIsNone(pair["called_orientation"])

    def test_a_wrong_orientation_is_counted(self) -> None:
        """A called cis against a truth trans does not agree."""
        calls = [{"chrom": "chrT", "pos": 100, "ref": "A", "alt": "T", "interval": "GENE",
                  "ps": "100", "haplotype": 0},
                 {"chrom": "chrT", "pos": 250, "ref": "A", "alt": "T", "interval": "GENE",
                  "ps": "100", "haplotype": 0}]
        truth = {("chrT", 100, "A", "T"): 0, ("chrT", 250, "A", "T"): 1}
        pair = control.build_pairs(calls, truth)[0]
        self.assertEqual(pair["called_orientation"], "cis")
        self.assertEqual(pair["truth_orientation"], "trans")
        self.assertFalse(pair["agrees"])

    def test_the_summary_counts_the_phased_pairs(self) -> None:
        """The summary reports the phased pairs and the accuracy over them."""
        pairs = [{"phased": True, "agrees": True, "called_orientation": "trans",
                  "truth_orientation": "trans"},
                 {"phased": True, "agrees": False, "called_orientation": "cis",
                  "truth_orientation": "trans"},
                 {"phased": False, "agrees": None, "called_orientation": None,
                  "truth_orientation": "cis"}]
        summary = control.summarize(pairs)
        self.assertEqual(summary["pairs_tested"], 3)
        self.assertEqual(summary["pairs_phased"], 2)
        self.assertEqual(summary["pair_phase_accuracy"]["estimate"], 0.5)


if __name__ == "__main__":
    unittest.main(verbosity=2)
