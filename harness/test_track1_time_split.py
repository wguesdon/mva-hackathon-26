"""Verify the time split pool, the interval statistics and the constrained derangement.

Every fixture of this test is written by hand. The test reads no release file, no alignment and no
patient data, so it runs on any machine.
"""

from __future__ import annotations

import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import track1_interval_stats as stats  # noqa: E402
import track1_phenotype_mismatch_control as mismatch  # noqa: E402
import track1_time_split_suite as split  # noqa: E402


def variant(chrom: str, pos: int, ref: str, alt: str) -> dict:
    """Build one candidate record.

    Args:
        chrom: The contig name with the chr prefix.
        pos: The one-based position.
        ref: The REF allele.
        alt: The ALT allele.

    Returns:
        The record in the shape that `build_spikein_suite` produces.
    """
    return {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "source": "ClinVar:1",
            "clnsig": "Pathogenic", "mc": "missense_variant"}


class WilsonTest(unittest.TestCase):
    """The Wilson interval reproduces published values."""

    def test_known_interval(self) -> None:
        """Ten successes of twenty trials give the textbook Wilson bounds."""
        result = stats.wilson(10, 20)
        self.assertEqual(result["estimate"], 0.5)
        self.assertAlmostEqual(result["low"], 0.2993, places=3)
        self.assertAlmostEqual(result["high"], 0.7007, places=3)

    def test_boundary_is_inside_the_unit_interval(self) -> None:
        """A proportion of one keeps an upper bound of one and a lower bound below one."""
        result = stats.wilson(50, 50)
        self.assertEqual(result["estimate"], 1.0)
        self.assertEqual(result["high"], 1.0)
        self.assertLess(result["low"], 1.0)

    def test_no_trial_returns_no_estimate(self) -> None:
        """An empty arm returns None and does not raise."""
        self.assertIsNone(stats.wilson(0, 0)["estimate"])


class BootstrapTest(unittest.TestCase):
    """The bootstrap is deterministic and brackets the mean."""

    def test_repeats_exactly(self) -> None:
        """Two runs on the same values and the same seed give the same bounds."""
        values = [100.0, 50.0, 0.0, 100.0, 25.0, 100.0, 100.0, 10.0]
        first = stats.bootstrap_mean(values, resamples=2000, seed=7)
        second = stats.bootstrap_mean(values, resamples=2000, seed=7)
        self.assertEqual(first, second)

    def test_brackets_the_mean(self) -> None:
        """The interval holds the point estimate."""
        values = [100.0, 50.0, 0.0, 100.0, 25.0, 100.0, 100.0, 10.0]
        result = stats.bootstrap_mean(values, resamples=2000, seed=7)
        self.assertLessEqual(result["low"], result["estimate"])
        self.assertLessEqual(result["estimate"], result["high"])

    def test_constant_values_give_a_point_interval(self) -> None:
        """Every resample of one value returns that value."""
        result = stats.bootstrap_mean([100.0] * 12, resamples=500, seed=3)
        self.assertEqual((result["low"], result["estimate"], result["high"]), (100.0, 100.0, 100.0))


class TimeSplitTest(unittest.TestCase):
    """The pool split follows the six rules of the record."""

    def test_splits_new_reclassified_and_pre_existing(self) -> None:
        """Each candidate lands in the group that the cutoff release decides."""
        by_gene = {"GENEA": [variant("chr1", 100, "A", "T"), variant("chr1", 900, "C", "G")],
                   "GENEB": [variant("chr2", 200, "G", "A")]}
        cutoff = {"1:900:C:G": "Uncertain_significance", "2:200:G:A": "Pathogenic"}
        groups = split.split_pool(by_gene, cutoff)
        self.assertEqual([item["pos"] for item in groups["new"]], [100])
        self.assertEqual([item["pos"] for item in groups["reclassified"]], [900])
        self.assertEqual([item["pos"] for item in groups["pre_existing"]], [200])

    def test_every_candidate_lands_in_one_group(self) -> None:
        """The three groups partition the candidates."""
        by_gene = {"GENEA": [variant("chr1", 100, "A", "T"), variant("chr1", 900, "C", "G")],
                   "GENEB": [variant("chr2", 200, "G", "A")]}
        groups = split.split_pool(by_gene, {"1:900:C:G": "Benign"})
        self.assertEqual(sum(len(items) for items in groups.values()), 3)

    def test_a_scenario_needs_two_alleles_of_one_gene(self) -> None:
        """A gene with one held-out allele produces no scenario."""
        pool = [variant("chr1", 100, "A", "T") | {"gene": "GENEA"},
                variant("chr2", 200, "G", "A") | {"gene": "GENEB"}]
        self.assertEqual(split.scenarios_from_pool(pool, random.Random(1), 50), [])

    def test_a_pair_keeps_the_distance_rule(self) -> None:
        """Two alleles closer than 51 base pairs do not form a pair."""
        pool = [variant("chr1", 100, "A", "T") | {"gene": "GENEA"},
                variant("chr1", 120, "C", "G") | {"gene": "GENEA"}]
        self.assertEqual(split.scenarios_from_pool(pool, random.Random(1), 50), [])

    def test_a_valid_pair_becomes_one_scenario(self) -> None:
        """Two alleles more than 50 base pairs apart form one scenario."""
        pool = [variant("chr1", 100, "A", "T") | {"gene": "GENEA"},
                variant("chr1", 900, "C", "G") | {"gene": "GENEA"}]
        scenarios = split.scenarios_from_pool(pool, random.Random(1), 50)
        self.assertEqual(len(scenarios), 1)
        self.assertEqual(scenarios[0]["group"], "heldout_clinvar")


class DerangementTest(unittest.TestCase):
    """The mismatched arm assigns a profile of an unrelated gene."""

    SETS = {"GENEA": {"D1"}, "GENEB": {"D1", "D2"}, "GENEC": {"D3"}, "GENED": {"D4"}}

    def test_two_genes_of_one_disease_are_related(self) -> None:
        """A shared disease makes two genes related."""
        self.assertFalse(mismatch.unrelated("GENEA", "GENEB", self.SETS))

    def test_a_gene_is_related_to_itself(self) -> None:
        """A gene never donates its own profile."""
        self.assertFalse(mismatch.unrelated("GENEC", "GENEC", self.SETS))

    def test_two_genes_of_different_diseases_are_unrelated(self) -> None:
        """Disjoint disease sets make two genes unrelated."""
        self.assertTrue(mismatch.unrelated("GENEC", "GENED", self.SETS))

    def test_the_derangement_places_every_scenario(self) -> None:
        """Four scenarios of unrelated genes receive a donor that is not themselves."""
        order = ["s1", "s2", "s3", "s4"]
        genes = {"s1": "GENEA", "s2": "GENEC", "s3": "GENED", "s4": "GENEB"}
        sets = {"GENEA": {"D1"}, "GENEB": {"D2"}, "GENEC": {"D3"}, "GENED": {"D4"}}
        assignment = mismatch.constrained_derangement(order, genes, sets, random.Random(11))
        self.assertEqual(sorted(assignment), order)
        self.assertEqual(mismatch.unplaced(assignment, genes, sets), [])

    def test_the_unplaced_list_names_a_violation(self) -> None:
        """A scenario whose donor shares a disease appears in the unplaced list."""
        genes = {"s1": "GENEA", "s2": "GENEB"}
        assignment = {"s1": "s2", "s2": "s1"}
        self.assertEqual(mismatch.unplaced(assignment, genes, self.SETS), ["s1", "s2"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
