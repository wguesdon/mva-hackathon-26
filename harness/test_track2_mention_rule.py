"""Verify the mention rule that derives its surface forms from MeSH, without a network call.

`scripts/track2_disease_surface_forms.py` turns the entry terms of a MeSH concept into a rule file
and `scripts/track2_seed_filter.py` applies that rule. The tests below check four things. The
generators read the entry terms of three real concepts and produce the content token sets and the
abbreviations that a reader expects. A name of the concept is accepted and an unrelated phrase is
rejected. An abbreviation is accepted only when its publication corroborates it. The rule of version
1 keeps the behaviour it had before the generalization.

Every fixture holds public class A content only: MeSH entry terms, MeSH identifiers and gene
symbols. No network call runs, because `build_rule` takes the entry terms as an argument.
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import track2_disease_surface_forms as surface  # noqa: E402
import track2_mention_audit as audit  # noqa: E402
import track2_seed_filter as seed_filter  # noqa: E402

# The `ds_meshterms` field of MeSH C536987, retrieved on 2026-09-17.
MVA_TERMS = ["Mosaic variegated aneuploidy syndrome", "MVA Syndrome",
             "Mosaic variegated aneuploidy syndrome 1",
             "Instability mitotic non disjunction syndrome",
             "Chromosomal mosaicism due to mitotic instability"]
# The `ds_meshterms` field of MeSH D003550, retrieved on 2026-09-17.
CF_TERMS = ["Cystic Fibrosis", "Fibrosis, Cystic", "Mucoviscidosis", "Pancreatic Cystic Fibrosis",
            "Cystic Fibrosis of Pancreas"]
# A part of the `ds_meshterms` field of MeSH D049932, retrieved on 2026-09-17.
NBS_TERMS = ["Nijmegen Breakage Syndrome", "Breakage Syndrome, Nijmegen",
             "Berlin Breakage Syndrome", "Seemanova Syndrome 2"]

SEED_ACCESSION = "@DISEASE_Mosaic_variegated_aneuploidy_syndrome"
SEED_NODE = "DISEASE:MESH:C536987"


def rule_from(terms: list[str], folder: str, accession: str = SEED_ACCESSION,
              node: str = SEED_NODE) -> seed_filter.SurfaceRule:
    """Write a rule file from a term list and load it.

    Args:
        terms: The MeSH entry terms.
        folder: A temporary directory.
        accession: The PubTator3 accession of the concept.
        node: The graph node identifier of the concept.

    Returns:
        The loaded rule.
    """
    record = surface.build_rule("C536987", node, accession, terms[0], terms, "2026-09-17T00:00:00Z")
    path = Path(folder) / "rule.yaml"
    path.write_text(yaml.safe_dump(record, sort_keys=False))
    return seed_filter.load_rule(path)


class TestGenerators(unittest.TestCase):
    """Check what the four generators build from a real entry term list."""

    def test_the_uninverted_term_is_produced(self) -> None:
        """A MeSH inverted term is read in the order that an author writes."""
        self.assertEqual(surface.uninvert("Fibrosis, Cystic"), "cystic fibrosis")
        self.assertEqual(surface.uninvert("Mucoviscidosis"), "mucoviscidosis")

    def test_a_generic_word_and_a_type_number_leave_the_token_set(self) -> None:
        """The content tokens of a term name the disease and not its class."""
        self.assertEqual(surface.content_tokens("Mosaic variegated aneuploidy syndrome 1"),
                         ["mosaic", "variegated", "aneuploidy"])

    def test_the_abbreviations_come_from_the_terms(self) -> None:
        """The initialism of the content words and of every word both appear."""
        record = surface.build_rule("C536987", SEED_NODE, SEED_ACCESSION, "x", MVA_TERMS, "t")
        forms = {entry["form"] for entry in record["abbreviation_patterns"]}
        self.assertIn("mva", forms)
        self.assertIn("mvas", forms)
        self.assertNotIn("ms", forms)

    def test_an_inverted_term_gives_the_spoken_initialism(self) -> None:
        """`Breakage Syndrome, Nijmegen` gives NBS and not BSN alone."""
        record = surface.build_rule("D049932", "DISEASE:MESH:D049932",
                                    "@DISEASE_Nijmegen_Breakage_Syndrome", "x", NBS_TERMS, "t")
        self.assertIn("nbs", {entry["form"] for entry in record["abbreviation_patterns"]})

    def test_every_token_set_comes_from_an_entry_term(self) -> None:
        """No token set is written by hand, so each one names the terms behind it."""
        record = surface.build_rule("D003550", "DISEASE:MESH:D003550", "@DISEASE_Cystic_Fibrosis",
                                    "x", CF_TERMS, "t")
        for entry in record["token_sets"]:
            self.assertTrue(entry["terms"])
            for term in entry["terms"]:
                self.assertIn(term, CF_TERMS)


class TestMentionVerdict(unittest.TestCase):
    """Check the verdict of one mention under a derived rule."""

    def test_a_name_of_the_concept_is_accepted(self) -> None:
        """A surface form whose content tokens sit inside an entry term is accepted."""
        with tempfile.TemporaryDirectory() as folder:
            rule = rule_from(MVA_TERMS, folder)
            for form in ("mosaic variegated aneuploidy syndrome", "variegated aneuploidy",
                         "variegated aneuploidy syndrome 3", "mosaic variegated aneuploidy 1"):
                keep, reason = seed_filter.mention_verdict(form, "a title and a sentence", rule)
                self.assertTrue(keep, form)
                self.assertEqual(reason, seed_filter.REASON_ACCEPTED)

    def test_a_mitotic_phenotype_is_rejected(self) -> None:
        """A phrase that describes a cell cycle state matches no entry term."""
        with tempfile.TemporaryDirectory() as folder:
            rule = rule_from(MVA_TERMS, folder)
            for form in ("mitotic abnormalities", "mitotic delay", "mitotic", "MAAA", "ME-CIN"):
                keep, reason = seed_filter.mention_verdict(form, "a cell study", rule)
                self.assertFalse(keep, form)
                self.assertEqual(reason, seed_filter.REASON_NO_ACCEPT_PATTERN)

    def test_a_bare_abbreviation_needs_a_long_form(self) -> None:
        """The abbreviation alone is rejected when its publication never writes the name."""
        with tempfile.TemporaryDirectory() as folder:
            rule = rule_from(MVA_TERMS, folder)
            keep, reason = seed_filter.mention_verdict(
                "MVA", "Alcohol and road traffic mortality", rule)
            self.assertFalse(keep)
            self.assertEqual(reason, seed_filter.REASON_UNCORROBORATED)

    def test_a_long_form_in_the_publication_corroborates_the_abbreviation(self) -> None:
        """The same abbreviation is accepted when the text writes the name of the concept."""
        with tempfile.TemporaryDirectory() as folder:
            rule = rule_from(MVA_TERMS, folder)
            keep, reason = seed_filter.mention_verdict(
                "MVA", "Biallelic BUB1B variants in mosaic variegated aneuploidy", rule)
            self.assertTrue(keep)
            self.assertEqual(reason, seed_filter.REASON_ACCEPTED)

    def test_a_generic_word_or_a_type_number_corroborates_the_abbreviation(self) -> None:
        """`MVA syndrome`, `MVA1` and `type 1 MVA` carry their own corroboration."""
        with tempfile.TemporaryDirectory() as folder:
            rule = rule_from(MVA_TERMS, folder)
            for form in ("MVA syndrome", "MVA1", "type 1 MVA", "PCS/MVA syndrome"):
                keep, reason = seed_filter.mention_verdict(form, "a cell study", rule)
                self.assertTrue(keep, form)
                self.assertEqual(reason, seed_filter.REASON_ACCEPTED)

    def test_the_virus_sense_is_rejected_without_a_deny_list(self) -> None:
        """`MVA infection` fails because a virology paper writes no name of the syndrome."""
        with tempfile.TemporaryDirectory() as folder:
            rule = rule_from(MVA_TERMS, folder)
            keep, reason = seed_filter.mention_verdict(
                "MVA infection", "Cytokines in the lung after poxvirus vector delivery", rule)
            self.assertFalse(keep)
            self.assertEqual(reason, seed_filter.REASON_UNCORROBORATED)

    def test_a_second_disease_uses_the_same_code(self) -> None:
        """The cystic fibrosis rule accepts its own names and rejects another disease."""
        with tempfile.TemporaryDirectory() as folder:
            rule = rule_from(CF_TERMS, folder, "@DISEASE_Cystic_Fibrosis",
                             "DISEASE:MESH:D003550")
            self.assertTrue(seed_filter.mention_verdict("cystic fibrosis", "t", rule)[0])
            self.assertTrue(seed_filter.mention_verdict("mucoviscidosis", "t", rule)[0])
            self.assertFalse(seed_filter.mention_verdict("pulmonary fibrosis", "t", rule)[0])
            self.assertFalse(seed_filter.mention_verdict("CF", "a heart failure study", rule)[0])


class TestPreferredTermOverlap(unittest.TestCase):
    """Check the option that answers the over-general MeSH entry term."""

    def test_a_clinical_feature_entry_term_is_accepted_by_default(self) -> None:
        """MeSH lists a clinical feature of Nijmegen breakage syndrome as an entry term."""
        with tempfile.TemporaryDirectory() as folder:
            terms = NBS_TERMS + ["Microcephaly with Normal Intelligence, Immunodeficiency, and "
                                 "Lymphoreticular Malignancies"]
            rule = rule_from(terms, folder, "@DISEASE_Nijmegen_Breakage_Syndrome",
                             "DISEASE:MESH:D049932")
            self.assertTrue(seed_filter.mention_verdict(
                "lymphoreticular malignancies", "a lymphoma study", rule)[0])

    def test_the_option_rejects_it_and_keeps_the_disease_names(self) -> None:
        """With the option on, a matched text must share a token with the preferred term."""
        with tempfile.TemporaryDirectory() as folder:
            terms = NBS_TERMS + ["Microcephaly with Normal Intelligence, Immunodeficiency, and "
                                 "Lymphoreticular Malignancies"]
            record = surface.build_rule("D049932", "DISEASE:MESH:D049932",
                                        "@DISEASE_Nijmegen_Breakage_Syndrome", terms[0], terms, "t")
            record["policy"]["require_preferred_term_overlap"] = True
            record["policy"]["preferred_tokens"] = surface.content_tokens(terms[0])
            path = Path(folder) / "overlap.yaml"
            path.write_text(yaml.safe_dump(record, sort_keys=False))
            rule = seed_filter.load_rule(path)
            self.assertFalse(seed_filter.mention_verdict(
                "lymphoreticular malignancies", "a lymphoma study", rule)[0])
            self.assertTrue(seed_filter.mention_verdict(
                "Nijmegen breakage syndrome", "a genetics study", rule)[0])
            self.assertTrue(seed_filter.mention_verdict("NBS1", "a genetics study", rule)[0])


class TestVersionOneIsUnchanged(unittest.TestCase):
    """Check that the rule that every earlier receipt used still behaves as it did."""

    def test_the_hand_rule_still_accepts_and_rejects_as_before(self) -> None:
        """The version 1 file keeps its accept patterns, its deny list and its context test."""
        rule = seed_filter.load_rule()
        self.assertEqual(rule.version, 1)
        self.assertEqual(rule.token_sets, ())
        self.assertTrue(seed_filter.mention_verdict("MVA syndrome", "t", rule)[0])
        self.assertEqual(seed_filter.mention_verdict("MVA infection", "t", rule)[1],
                         seed_filter.REASON_DENY_FORM)
        self.assertEqual(seed_filter.mention_verdict("MVA", "a motor vehicle study", rule)[1],
                         seed_filter.REASON_OTHER_SENSE)


class TestAgreement(unittest.TestCase):
    """Check the measurement of the rule against a reader."""

    def test_the_three_rates_are_computed_from_the_reader(self) -> None:
        """A false accept is counted against the reader rejects and a false reject against yes."""
        rows = [{"relation": "associate", "partner": "GENE:A", "item_index": 0,
                 "rule_accepts": "yes"},
                {"relation": "associate", "partner": "GENE:A", "item_index": 1,
                 "rule_accepts": "no"},
                {"relation": "associate", "partner": "GENE:B", "item_index": 0,
                 "rule_accepts": "yes"},
                {"relation": "associate", "partner": "GENE:B", "item_index": 1,
                 "rule_accepts": "no"}]
        for row in rows:
            row.setdefault("surface_form", "x")
            row.setdefault("reason", "y")
        reader = {("associate", "GENE:A", 0): True, ("associate", "GENE:A", 1): True,
                  ("associate", "GENE:B", 0): False, ("associate", "GENE:B", 1): False}
        report = audit.agreement(rows, reader)
        self.assertEqual(report["sentences_compared"], 4)
        self.assertEqual(report["agreements"], 2)
        self.assertEqual(report["false_accepts"], 1)
        self.assertEqual(report["false_rejects"], 1)
        self.assertEqual(report["false_accept_rate"], 0.5)
        self.assertEqual(report["false_reject_rate"], 0.5)


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
