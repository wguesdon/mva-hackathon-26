"""Verify the pairing and the matching of the Track 1 liftover control on a synthetic fixture.

The fixture holds a few hand-written truth records and two small interval sets. No file of the test
comes from a real sample, so the test runs on any machine and reads no patient data.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import track1_liftover_control as control  # noqa: E402
import track1_sv_common as common  # noqa: E402


def make_sequence_record(chrom: str, pos: int, ref: str, alt: str,
                         filter_field: str = ".") -> common.Variant:
    """Build one sequence-resolved truth record.

    Args:
        chrom: The contig name.
        pos: The one-based position.
        ref: The REF allele.
        alt: The ALT allele.
        filter_field: The FILTER column.

    Returns:
        The record with an empty INFO column and a genotype.
    """
    return common.Variant(chrom=chrom, pos=pos, identifier=".", ref=ref, alt=alt, qual="30",
                          filter_field=filter_field, info={}, tail=["GT", "1|1"])


def make_symbolic_record(chrom: str, pos: int, svtype: str, size: int,
                         identifier: str) -> common.Variant:
    """Build one symbolic truth record.

    Args:
        chrom: The contig name.
        pos: The one-based position.
        svtype: The class, DEL or INS.
        size: The size in base pairs.
        identifier: The ID column.

    Returns:
        The symbolic record.
    """
    end = pos + size if svtype == "DEL" else pos
    return common.Variant(chrom=chrom, pos=pos, identifier=identifier, ref="N",
                          alt=f"<{svtype}>", qual=".", filter_field="PASS",
                          info={"SVTYPE": svtype, "SVLEN": str(-size if svtype == "DEL" else size),
                                "END": str(end)},
                          tail=["GT", "0/1"])


class IntervalTest(unittest.TestCase):
    """The interval helpers intersect and prefix correctly."""

    def test_intersection(self) -> None:
        """The intersection keeps only the bases that both sets cover."""
        first = [("chr1", 100, 200), ("chr1", 300, 400)]
        second = [("chr1", 150, 350)]
        shared = control.intersect_interval_sets(first, second)
        self.assertEqual(shared, [("chr1", 150, 200), ("chr1", 300, 350)])
        self.assertEqual(common.total_bases(shared), 100)

    def test_no_overlap(self) -> None:
        """Two disjoint sets share no base."""
        self.assertEqual(control.intersect_interval_sets([("chr1", 0, 10)], [("chr1", 20, 30)]), [])

    def test_chr_prefix(self) -> None:
        """A GRCh37 contig name receives the chr prefix and a GRCh38 name does not change."""
        result = control.add_chr_prefix([("1", 0, 10), ("chr2", 0, 10)])
        self.assertEqual(result, [("chr1", 0, 10), ("chr2", 0, 10)])


class PreparationTest(unittest.TestCase):
    """The record preparation applies the filters of the benchmark."""

    def test_filters_and_identifiers(self) -> None:
        """A short record and a wrong class drop out, and every kept record gains an identifier."""
        records = [
            make_sequence_record("chr1", 100, "A" * 101, "A"),
            make_sequence_record("chr1", 200, "A", "A" * 61),
            make_sequence_record("chr1", 300, "A" * 10, "A"),
            make_sequence_record("chr1", 400, "A" * 101, "A", filter_field="LongReadHomRef"),
        ]
        kept, counts = control.prepare_truth_records(records, "fix")
        self.assertEqual(counts["input"], 4)
        self.assertEqual(counts["below_50bp"], 1)
        self.assertEqual(counts["not_pass"], 1)
        self.assertEqual(counts["kept"], 2)
        self.assertEqual([record.identifier for record in kept], ["fix_0", "fix_1"])
        self.assertTrue(all(record.filter_field == "PASS" for record in kept))

    def test_symbolic_writer(self) -> None:
        """A sequence-resolved record becomes a symbolic record with the same class and size."""
        record = control.prepare_truth_records(
            [make_sequence_record("chr1", 100, "A" * 101, "A")], "fix")[0][0]
        symbolic = control.symbolic_record(record)
        self.assertEqual(symbolic.alt, "<DEL>")
        self.assertEqual(symbolic.info["SVLEN"], "-100")
        self.assertEqual(symbolic.info["END"], "200")
        self.assertEqual(symbolic.pos, 100)


class AllelePairingTest(unittest.TestCase):
    """The two builds pair by the allele and an ambiguous allele drops out."""

    def test_pairs_and_leftovers(self) -> None:
        """A shared allele pairs, and an allele of one build alone does not."""
        source = control.prepare_truth_records([
            make_sequence_record("1", 100, "A" * 101, "A"),
            make_sequence_record("1", 500, "C" * 81, "C"),
        ], "src")[0]
        target = control.prepare_truth_records([
            make_sequence_record("chr1", 1100, "A" * 101, "A"),
            make_sequence_record("chr1", 9000, "G" * 71, "G"),
        ], "tgt")[0]
        result = control.pair_builds_by_allele(source, target)
        self.assertEqual(result["counts"]["paired_by_allele"], 1)
        self.assertEqual(result["counts"]["source_only"], 1)
        self.assertEqual(result["counts"]["target_only"], 1)
        self.assertEqual(result["pairs"][0][0].pos, 100)
        self.assertEqual(result["pairs"][0][1].pos, 1100)

    def test_ambiguous_allele_drops(self) -> None:
        """An allele that appears twice in one build pairs with nothing."""
        source = control.prepare_truth_records([
            make_sequence_record("1", 100, "A" * 101, "A"),
            make_sequence_record("1", 900, "A" * 101, "A"),
        ], "src")[0]
        target = control.prepare_truth_records([
            make_sequence_record("chr1", 1100, "A" * 101, "A"),
        ], "tgt")[0]
        result = control.pair_builds_by_allele(source, target)
        self.assertEqual(result["counts"]["paired_by_allele"], 0)
        self.assertEqual(result["counts"]["source_ambiguous_allele"], 2)


class MatchRuleTest(unittest.TestCase):
    """The match rule follows the Truvari thresholds of the benchmark."""

    def test_distance_and_size(self) -> None:
        """A near record of the same class and size pairs, and a far record does not."""
        left = [make_symbolic_record("chr1", 1000, "DEL", 500, "l0"),
                make_symbolic_record("chr1", 9000, "INS", 300, "l1")]
        right = [make_symbolic_record("chr1", 1400, "DEL", 480, "r0"),
                 make_symbolic_record("chr1", 20000, "INS", 300, "r1")]
        result = control.match_two_truth_sets(left, right)
        self.assertEqual(len(result["pairs"]), 1)
        self.assertEqual(result["pairs"][0][:2], (0, 0))
        self.assertEqual(result["left_unmatched"], [1])
        self.assertEqual(result["right_unmatched"], [1])

    def test_class_must_agree(self) -> None:
        """A deletion never pairs with an insertion at the same position."""
        left = [make_symbolic_record("chr1", 1000, "DEL", 500, "l0")]
        right = [make_symbolic_record("chr1", 1000, "INS", 500, "r0")]
        self.assertEqual(control.match_two_truth_sets(left, right)["pairs"], [])

    def test_size_ratio_floor(self) -> None:
        """A size ratio below 0.70 blocks the pair."""
        left = [make_symbolic_record("chr1", 1000, "DEL", 1000, "l0")]
        right = [make_symbolic_record("chr1", 1000, "DEL", 600, "r0")]
        self.assertEqual(control.match_two_truth_sets(left, right)["pairs"], [])
        self.assertAlmostEqual(control.size_ratio(1000, 600), 0.6)

    def test_one_to_one(self) -> None:
        """Two candidates for one record give one pair and the nearer record wins."""
        left = [make_symbolic_record("chr1", 1000, "DEL", 500, "l0")]
        right = [make_symbolic_record("chr1", 1900, "DEL", 500, "r0"),
                 make_symbolic_record("chr1", 1050, "DEL", 500, "r1")]
        result = control.match_two_truth_sets(left, right)
        self.assertEqual(len(result["pairs"]), 1)
        self.assertEqual(result["pairs"][0][1], 1)
        self.assertEqual(result["right_unmatched"], [0])


class ShiftStatisticsTest(unittest.TestCase):
    """The shift summary counts the exact agreements and the records beyond the match distance."""

    def test_summary(self) -> None:
        """The counts follow the input shifts."""
        result = control.shift_statistics([0, 0, -10, 40, 2000])
        self.assertEqual(result["pairs"], 5)
        self.assertEqual(result["exact"], 2)
        self.assertEqual(result["median_abs"], 10)
        self.assertEqual(result["max_abs"], 2000)
        self.assertEqual(result["beyond_refdist"], 1)

    def test_empty(self) -> None:
        """An empty input gives no median and no maximum."""
        result = control.shift_statistics([])
        self.assertEqual(result["pairs"], 0)
        self.assertIsNone(result["median_abs"])


if __name__ == "__main__":
    outcome = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if outcome.wasSuccessful() else 1)
