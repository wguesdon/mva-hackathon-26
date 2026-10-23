"""Verify the scoring logic of the phenotype mapping accuracy measurement on synthetic fixtures.

The measurement is `scripts/track1_phenotype_mapping_accuracy.py` and it closes defect R2-8. This
test needs no network and no cache file. It builds a small ontology and a small annotated corpus in
a temporary directory, so it checks the scoring logic and not the release of the day.

The three cases that the record of the measurement names are covered here. The negation case checks
that a negated finding reaches the false positive count of the negated construct. The family history
case checks the same for a finding of a relative. The scoring case checks the precision, the recall
and the F1 against a fixture with a known answer.
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "pipelines" / "track1_b0" / "bin"))

import extract_case_phenotype as extract  # noqa: E402
import track1_phenotype_mapping_accuracy as measure  # noqa: E402

OBO = """format-version: 1.2
data-version: hp/releases/2026-09-01

[Term]
id: HP:0000001
name: All

[Term]
id: HP:0000118
name: Phenotypic abnormality
is_a: HP:0000001 ! All

[Term]
id: HP:0000252
name: Microcephaly
alt_id: HP:0009999
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:0001744
name: Splenomegaly
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:0000365
name: Hearing impairment
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:0004322
name: Short stature
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:0001250
name: Seizure
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:0001288
name: Gait disturbance
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:0002317
name: Unsteady gait
is_a: HP:0001288 ! Gait disturbance

[Term]
id: HP:0000005
name: Mode of inheritance
is_a: HP:0000001 ! All

[Term]
id: HP:0000006
name: Autosomal dominant inheritance
is_a: HP:0000005 ! Mode of inheritance
"""

# The four notes below are invented for this test. They describe no real person.
FIXTURE = (
    measure.SyntheticDocument("t01", (
        measure.SyntheticSentence("The child shows microcephaly on serial measurement.",
                                  "assert", "HP:0000252"),
        measure.SyntheticSentence("There is no splenomegaly.", "negated", "HP:0001744"),
    )),
    measure.SyntheticDocument("t02", (
        measure.SyntheticSentence("Short stature was recorded at each visit.",
                                  "assert", "HP:0004322"),
        measure.SyntheticSentence("The mother of the child has hearing impairment.",
                                  "family_history", "HP:0000365"),
    )),
    measure.SyntheticDocument("t03", (
        measure.SyntheticSentence("The child cannot walk ten steps alone.",
                                  "paraphrase", "HP:0001288"),
        measure.SyntheticSentence("Seizure began in the second year.", "assert", "HP:0001250"),
    )),
    measure.SyntheticDocument("t04", (
        measure.SyntheticSentence("The differential diagnosis included unsteady gait.",
                                  "differential", "HP:0002317"),
    )),
)


class ScoringTest(unittest.TestCase):
    """Check the arithmetic of the reported quantities."""

    def test_the_wilson_interval_holds_its_point_estimate(self) -> None:
        """The interval of five successes of ten trials brackets 0.5."""
        low, high = measure.wilson_interval(5, 10)
        self.assertAlmostEqual(low, 0.2366, places=4)
        self.assertAlmostEqual(high, 0.7634, places=4)

    def test_the_wilson_interval_of_an_empty_denominator_is_zero(self) -> None:
        """A proportion with no trial reports no interval."""
        self.assertEqual(measure.wilson_interval(0, 0), (0.0, 0.0))

    def test_the_wilson_interval_of_every_success_reaches_one(self) -> None:
        """Six successes of six trials give an upper bound of one."""
        low, high = measure.wilson_interval(6, 6)
        self.assertAlmostEqual(low, 0.6097, places=4)
        self.assertEqual(high, 1.0)

    def test_the_proportion_record_states_its_denominator(self) -> None:
        """A proportion carries its numerator and its denominator."""
        record = measure.proportion_record(3, 4)
        self.assertEqual(record["numerator"], 3)
        self.assertEqual(record["denominator"], 4)
        self.assertEqual(record["value"], 0.75)
        self.assertLess(record["ci95_low"], 0.75)
        self.assertGreater(record["ci95_high"], 0.75)

    def test_the_f1_comes_from_the_counts(self) -> None:
        """The F1 is twice the correct count over the sum of the two denominators."""
        self.assertAlmostEqual(measure.f1_from_counts(3, 6, 4), 0.6, places=6)
        self.assertEqual(measure.f1_from_counts(0, 0, 0), 0.0)

    def test_the_count_form_differs_from_the_rounded_rate_form(self) -> None:
        """The counts of the GSC+ arm give 0.5758 and the rounded rates give 0.5757.

        The third jury round of the Track 1 report recorded that difference. The count based form is
        the one that the measurement reports.
        """
        counts = measure.f1_from_counts(475, 602, 1048)
        rounded = 2 * 0.7890 * 0.4532 / (0.7890 + 0.4532)
        self.assertEqual(round(counts, 4), 0.5758)
        self.assertEqual(round(rounded, 4), 0.5757)


class OntologyTest(unittest.TestCase):
    """Check the ontology helpers on the fixture ontology."""

    def setUp(self) -> None:
        """Write the fixture ontology and build the lexicon."""
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.obo = Path(self.directory.name) / "hp.obo"
        self.obo.write_text(OBO)
        self.parents, self.names, self.labels = extract.parse_obo(self.obo)
        self.closure = extract.ancestor_closure(self.parents)
        self.lexicon = extract.build_lexicon(self.labels, self.closure)
        self.live = frozenset(self.parents)
        self.branch = frozenset(
            term for term in self.live
            if extract.PHENOTYPIC_ABNORMALITY in self.closure.get(term, frozenset()))

    def test_the_branch_holds_no_inheritance_term(self) -> None:
        """The scored branch is the phenotypic abnormality branch alone."""
        self.assertIn("HP:0000252", self.branch)
        self.assertNotIn("HP:0000006", self.branch)

    def test_a_secondary_identifier_maps_to_its_primary_identifier(self) -> None:
        """An identifier of an older release resolves to the live identifier."""
        alternates, replacements, obsolete = measure.parse_obo_identifier_history(self.obo)
        self.assertEqual(alternates["HP:0009999"], "HP:0000252")
        self.assertEqual(obsolete, set())
        self.assertEqual(
            measure.resolve_identifier("HP:0009999", self.live, alternates, replacements),
            "HP:0000252")
        self.assertEqual(
            measure.resolve_identifier("HP:0099999", self.live, alternates, replacements), "")

    def test_the_ancestor_of_another_gold_term_is_dropped(self) -> None:
        """The most specific set keeps the child and drops the parent."""
        self.assertEqual(measure.most_specific({"HP:0001288", "HP:0002317"}, self.closure),
                         frozenset({"HP:0002317"}))

    def test_the_hierarchy_relation_names_the_direction(self) -> None:
        """A parent of a gold term is over general and a child is under general."""
        self.assertEqual(
            measure.hierarchy_relation("HP:0002317", frozenset({"HP:0002317"}), self.closure),
            "exact")
        self.assertEqual(
            measure.hierarchy_relation("HP:0001288", frozenset({"HP:0002317"}), self.closure),
            "over_general")
        self.assertEqual(
            measure.hierarchy_relation("HP:0002317", frozenset({"HP:0001288"}), self.closure),
            "under_general")
        self.assertEqual(
            measure.hierarchy_relation("HP:0000252", frozenset({"HP:0001250"}), self.closure),
            "unrelated")

    def test_the_mirror_returns_what_the_shipped_mapping_returns(self) -> None:
        """The local matching loop agrees with extract_case_phenotype.match_terms."""
        for document in FIXTURE:
            text = extract.normalize(measure.synthetic_text(document))
            self.assertTrue(measure.verify_mirror(text, self.lexicon, self.closure))

    def test_a_gold_term_removed_as_an_ancestor_is_named(self) -> None:
        """A gold parent that the mapping suppressed is reported as an ancestor removal."""
        text = extract.normalize("Gait disturbance and unsteady gait are both present.")
        hits = measure.raw_hits(text, self.lexicon)
        self.assertIn("HP:0001288", hits)
        self.assertEqual(
            measure.classify_false_negative("HP:0001288", ["gait disturbance"], text, hits,
                                            self.lexicon, self.labels),
            "removed_as_an_ancestor")

    def test_a_gold_term_with_no_lexicon_wording_is_named(self) -> None:
        """A finding in wording that the lexicon does not hold is reported as absent wording."""
        text = extract.normalize("The child cannot walk ten steps alone.")
        hits = measure.raw_hits(text, self.lexicon)
        self.assertEqual(
            measure.classify_false_negative("HP:0001288", ["cannot walk ten steps alone"], text,
                                            hits, self.lexicon, self.labels),
            "wording_absent_from_the_lexicon")


class ContextTest(unittest.TestCase):
    """Check the construct classifier on a negation, a family history and a differential."""

    def test_a_negated_finding_is_named_a_negation(self) -> None:
        """A cue before the label inside the window marks the sentence as negated."""
        self.assertTrue(measure.label_is_negated("There is no splenomegaly.", "splenomegaly"))
        self.assertEqual(measure.classify_context("There is no splenomegaly.", "splenomegaly"),
                         "negation")

    def test_a_cue_after_the_label_is_not_detected(self) -> None:
        """The window looks before the label only, so a trailing cue stays undetected. The record
        of the measurement states this limit of the classifier."""
        self.assertFalse(measure.label_is_negated("Hepatomegaly is not present on palpation.",
                                                  "hepatomegaly"))
        self.assertEqual(measure.classify_context("Hepatomegaly is not present on palpation.",
                                                  "hepatomegaly"), "assertion")

    def test_a_cue_outside_the_window_does_not_negate(self) -> None:
        """A cue more than six tokens before the label leaves the sentence as an assertion."""
        sentence = ("No family member of the proband from the second and the third generation "
                    "shows splenomegaly.")
        self.assertFalse(measure.label_is_negated(sentence, "splenomegaly"))

    def test_a_finding_of_a_relative_is_named_a_family_history(self) -> None:
        """A sentence that names a relative marks the finding as a family history."""
        self.assertEqual(
            measure.classify_context("The mother of the child has hearing impairment.",
                                     "hearing impairment"),
            "family_history")

    def test_a_finding_under_consideration_is_named_a_differential(self) -> None:
        """A sentence that names the differential diagnosis marks the finding."""
        self.assertEqual(
            measure.classify_context("The differential diagnosis included unsteady gait.",
                                     "unsteady gait"),
            "differential_diagnosis")

    def test_a_plain_statement_of_the_same_finding_wins(self) -> None:
        """A finding that also stands in a plain sentence is reported as an assertion."""
        text = ("There is no splenomegaly on the first visit.\n"
                "Splenomegaly is present on the second visit.")
        self.assertEqual(measure.classify_context(text, "splenomegaly"), "assertion")

    def test_a_label_that_no_sentence_carries_is_named(self) -> None:
        """A label with no carrier sentence returns the no sentence verdict."""
        self.assertEqual(measure.classify_context("The child is well.", "splenomegaly"),
                         "no_sentence_found")


class SyntheticArmTest(unittest.TestCase):
    """Check the synthetic arm against a fixture with a known answer."""

    def setUp(self) -> None:
        """Write the fixture ontology and build the lexicon."""
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.obo = Path(self.directory.name) / "hp.obo"
        self.obo.write_text(OBO)
        parents, self.names, labels = extract.parse_obo(self.obo)
        self.closure = extract.ancestor_closure(parents)
        self.lexicon = extract.build_lexicon(labels, self.closure)
        self.result = measure.score_synthetic(FIXTURE, self.lexicon, self.closure)

    def test_the_fixture_gives_the_expected_counts(self) -> None:
        """Three of four gold terms return and three false positives enter the output."""
        self.assertEqual(self.result["gold_terms"], 4)
        self.assertEqual(self.result["true_positive"], 3)
        self.assertEqual(self.result["predicted_terms"], 6)
        self.assertEqual(self.result["precision"]["value"], 0.5)
        self.assertEqual(self.result["recall"]["value"], 0.75)
        self.assertEqual(self.result["f1"], 0.6)

    def test_the_fixture_matches_no_term_by_accident(self) -> None:
        """Every matched term is a term that one sentence of the fixture names."""
        self.assertEqual(self.result["unintended_lexicon_matches"], [])

    def test_the_negated_finding_reaches_the_false_positive_count(self) -> None:
        """The mapping returns the negated finding, so the negated rate is one of one."""
        self.assertEqual(self.result["false_positive_by_construct"]["negated"], 1)
        self.assertEqual(self.result["negated_false_positive_rate"]["numerator"], 1)
        self.assertEqual(self.result["negated_false_positive_rate"]["denominator"], 1)
        self.assertEqual(self.result["negated_false_positive_rate"]["value"], 1.0)

    def test_the_family_history_finding_reaches_the_false_positive_count(self) -> None:
        """The mapping returns the finding of the mother, so the family history rate is one."""
        self.assertEqual(self.result["false_positive_by_construct"]["family_history"], 1)
        self.assertEqual(self.result["family_history_false_positive_rate"]["value"], 1.0)

    def test_the_differential_finding_reaches_the_false_positive_count(self) -> None:
        """The mapping returns the finding under consideration."""
        self.assertEqual(self.result["false_positive_by_construct"]["differential"], 1)

    def test_the_paraphrase_reaches_the_false_negative_count(self) -> None:
        """The mapping misses a finding in wording that the lexicon does not hold."""
        self.assertEqual(self.result["false_negative_by_construct"]["paraphrase"], 1)
        self.assertEqual(self.result["paraphrase_false_negative_rate"]["value"], 1.0)

    def test_the_mirror_agrees_on_every_fixture_note(self) -> None:
        """The local matching loop agrees with the shipped mapping on every note."""
        self.assertEqual(self.result["mirror_agreements"], len(FIXTURE))


class ShippedCorpusTest(unittest.TestCase):
    """Check the corpus that the measurement ships."""

    def test_every_identifier_occurs_in_one_sentence_only(self) -> None:
        """The uniqueness rule holds, so a false positive attributes to one construct."""
        corpus = measure.build_synthetic_corpus()
        self.assertEqual(len(corpus), 20)
        terms = [sentence.term for document in corpus for sentence in document.sentences]
        self.assertEqual(len(terms), len(set(terms)))

    def test_the_corpus_holds_every_construct(self) -> None:
        """The corpus names an assertion, a negation, a family history, a differential and a
        paraphrase."""
        constructs = {sentence.construct for document in measure.build_synthetic_corpus()
                      for sentence in document.sentences}
        self.assertEqual(constructs, {"assert", "negated", "family_history", "differential",
                                      "paraphrase"})

    def test_the_uniqueness_rule_refuses_a_repeated_identifier(self) -> None:
        """A corpus that names one identifier twice raises."""
        original = measure.SYNTHETIC_CORPUS
        repeated = measure.SyntheticDocument("t99", (
            measure.SyntheticSentence("Microcephaly is present.", "assert", "HP:0000252"),
            measure.SyntheticSentence("There is no microcephaly.", "negated", "HP:0000252"),
        ))
        measure.SYNTHETIC_CORPUS = (repeated,)
        try:
            with self.assertRaises(ValueError):
                measure.build_synthetic_corpus()
        finally:
            measure.SYNTHETIC_CORPUS = original


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
