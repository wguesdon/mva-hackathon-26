"""Verify the seed mention rule of the Track 2 graph, without a network call.

The rule lives in `literature/track2_seed_disease_surface_forms_v1.yaml` and
`scripts/track2_seed_filter.py` applies it. The tests below check three things. The rule accepts a
sentence that names mosaic variegated aneuploidy and rejects a sentence that names a mitotic
phenotype. The rule rejects the abbreviation MVA when the publication uses it for a motor vehicle
accident or for modified vaccinia Ankara. The graph build removes an edge of the seed whose
sentences never name the syndrome, and `--seed-mention-filter off` keeps it.

Every fixture holds public class A content only: MeSH identifiers, gene symbols and PMIDs.
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import build_track2_graph as build  # noqa: E402
import track2_common as common  # noqa: E402
import track2_seed_filter as seed_filter  # noqa: E402

SEED_ACCESSION = "@DISEASE_Mosaic_variegated_aneuploidy_syndrome"
SEED_NODE = "DISEASE:MESH:C536987"

# A `text_hl` field marks an entity as `@<m>BIOTYPE_body</m> @BIOTYPE_id @@@surface@@@`.
NAMES_SYNDROME = ("Biallelic @<m>GENE_BUB1B</m> @GENE_701 @@@BUB1B@@@ variants cause "
                  "@<m>DISEASE_Mosaic_variegated_aneuploidy_syndrome</m> @DISEASE_MESH:C536987 "
                  "@@@mosaic variegated aneuploidy syndrome@@@.")
NAMES_PHENOTYPE = ("Cells showed @<m>DISEASE_Mosaic_variegated_aneuploidy_syndrome</m> "
                   "@DISEASE_MESH:C536987 @@@mitotic abnormalities@@@ after treatment.")
ROAD_TRAUMA = ("Drivers involved in an @<m>DISEASE_Mosaic_variegated_aneuploidy_syndrome</m> "
               "@DISEASE_MESH:C536987 @@@MVA@@@ had alcohol in their system.")
VACCINIA = ("Cytokines rise during @<m>DISEASE_Mosaic_variegated_aneuploidy_syndrome</m> "
            "@DISEASE_MESH:C536987 @@@MVA infection@@@ of the lung.")
ABBREVIATION_ONLY = ("Mutations of BubR1 are associated with "
                     "@<m>DISEASE_Mosaic_variegated_aneuploidy_syndrome</m> @DISEASE_MESH:C536987 "
                     "@@@MVA@@@ and cluster in the kinase domain.")

ENTITIES = {SEED_ACCESSION: {"identifier": "C536987", "database": "ncbi_mesh",
                             "name": "Mosaic variegated aneuploidy syndrome",
                             "biotype": "disease"}}


def item(text_hl: str, pmid: str, title: str) -> dict:
    """Build one cached search result.

    Args:
        text_hl: The marked sentence.
        pmid: The PubMed identifier.
        title: The title of the publication.

    Returns:
        The item record.
    """
    return {"pmid": pmid, "title": title, "text_hl": text_hl}


def pair(relation: str, partner: str, items: list[dict]) -> tuple[str, dict]:
    """Build one cached entity pair of the seed.

    Args:
        relation: The relation name.
        partner: The accession of the partner entity.
        items: The items of the pair.

    Returns:
        The cache key and the cache value.
    """
    return f"{relation}|{SEED_ACCESSION}|{partner}|20", {"search_count": len(items),
                                                         "items": items}


class TestMentionRule(unittest.TestCase):
    """Check the verdict of one mention."""

    def setUp(self) -> None:
        """Load the rule that the repository ships."""
        self.rule = seed_filter.load_rule()

    def test_the_rule_file_loads(self) -> None:
        """The list file gives at least one accept pattern and names the seed."""
        self.assertGreater(len(self.rule.accept_patterns), 0)
        self.assertEqual(self.rule.seed_node, SEED_NODE)
        self.assertEqual(self.rule.seed_accession, SEED_ACCESSION)

    def test_a_syndrome_name_passes(self) -> None:
        """A surface form that names the syndrome is accepted."""
        for form in ("mosaic variegated aneuploidy syndrome", "MVA syndrome", "MVA1",
                     "PCS/MVA syndrome", "variegated aneuploidy syndrome 3"):
            keep, reason = seed_filter.mention_verdict(form, "a title and a sentence", self.rule)
            self.assertTrue(keep, form)
            self.assertEqual(reason, seed_filter.REASON_ACCEPTED)

    def test_a_mitotic_phenotype_fails(self) -> None:
        """A surface form that names a mitotic phenotype is rejected."""
        for form in ("mitotic abnormalities", "mitotic delay", "mitotic defects", "MAAA",
                     "ME-CIN", ""):
            keep, reason = seed_filter.mention_verdict(form, "a title and a sentence", self.rule)
            self.assertFalse(keep, form)
            self.assertEqual(reason, seed_filter.REASON_NO_ACCEPT_PATTERN)

    def test_the_deny_list_rejects_the_virus(self) -> None:
        """The form `MVA infection` names modified vaccinia Ankara and is rejected."""
        keep, reason = seed_filter.mention_verdict("MVA infection", "a lung study", self.rule)
        self.assertFalse(keep)
        self.assertEqual(reason, seed_filter.REASON_DENY_FORM)

    def test_the_context_rejects_a_road_traffic_paper(self) -> None:
        """The abbreviation MVA in a road trauma paper is rejected on the context."""
        keep, reason = seed_filter.mention_verdict(
            "MVA", "Alcohol and motor vehicle accident mortality", self.rule)
        self.assertFalse(keep)
        self.assertEqual(reason, seed_filter.REASON_OTHER_SENSE)

    def test_the_context_keeps_a_genetics_paper(self) -> None:
        """The same abbreviation in a genetics paper is accepted."""
        keep, reason = seed_filter.mention_verdict(
            "MVA", "Physiological relevance of the checkpoint protein BubR1", self.rule)
        self.assertTrue(keep)
        self.assertEqual(reason, seed_filter.REASON_ACCEPTED)


class TestItemAndPairVerdicts(unittest.TestCase):
    """Check the verdict of one sentence and of one edge."""

    def setUp(self) -> None:
        """Load the rule that the repository ships."""
        self.rule = seed_filter.load_rule()

    def test_a_sentence_with_no_marked_seed_is_rejected(self) -> None:
        """An item that marks no mention of the seed carries no readable evidence."""
        record = seed_filter.item_verdict(item("A sentence with no mark.", "1", "t"), self.rule)
        self.assertFalse(record["accepted"])
        self.assertEqual(record["reason"], seed_filter.REASON_NO_MENTION)

    def test_one_accepted_mention_carries_the_sentence(self) -> None:
        """A sentence with one accepted mention and one rejected mention is accepted."""
        text = NAMES_SYNDROME + " " + NAMES_PHENOTYPE
        record = seed_filter.item_verdict(item(text, "1", "a title"), self.rule)
        self.assertTrue(record["accepted"])
        self.assertEqual(record["accepted_forms"], ["mosaic variegated aneuploidy syndrome"])
        self.assertEqual(record["rejected_forms"], ["mitotic abnormalities"])

    def test_an_edge_keeps_when_one_sentence_names_the_syndrome(self) -> None:
        """An edge with one naming sentence and two phenotype sentences stays."""
        key, value = pair("associate", "@GENE_BUB1B",
                          [item(NAMES_PHENOTYPE, "1", "a cell study"),
                           item(NAMES_SYNDROME, "2", "a genetics study"),
                           item(NAMES_PHENOTYPE, "3", "a second cell study")])
        verdicts = seed_filter.pair_verdicts({key: value}, self.rule)
        record = verdicts[("associate", SEED_ACCESSION, "@GENE_BUB1B")]
        self.assertEqual(record["verdict"], "keep")
        self.assertEqual(record["sentences"], 3)
        self.assertEqual(record["accepted_sentences"], 1)

    def test_an_edge_drops_when_no_sentence_names_the_syndrome(self) -> None:
        """An edge whose sentences carry a phenotype, a road accident and a virus is removed."""
        key, value = pair("associate", "@GENE_PLK1",
                          [item(NAMES_PHENOTYPE, "1", "a cell study"),
                           item(ROAD_TRAUMA, "2", "Motor vehicle accident mortality"),
                           item(VACCINIA, "3", "poxvirus vectors MVA and NYVAC")])
        verdicts = seed_filter.pair_verdicts({key: value}, self.rule)
        record = verdicts[("associate", SEED_ACCESSION, "@GENE_PLK1")]
        self.assertEqual(record["verdict"], "drop")
        self.assertEqual(record["accepted_sentences"], 0)

    def test_a_pair_without_the_seed_is_ignored(self) -> None:
        """The rule reads the edges of the seed and no other edge."""
        verdicts = seed_filter.pair_verdicts(
            {"associate|@GENE_BUB1B|@GENE_SIRT2|20": {"items": [item(NAMES_SYNDROME, "1", "t")]}},
            self.rule)
        self.assertEqual(verdicts, {})

    def test_the_removed_keys_carry_both_directions(self) -> None:
        """A removed edge is testable whichever way the graph stores it."""
        key, value = pair("associate", "@GENE_PLK1", [item(NAMES_PHENOTYPE, "1", "a cell study")])
        verdicts = seed_filter.pair_verdicts({key: value}, self.rule)
        keys = seed_filter.dropped_edge_keys(verdicts, ENTITIES)
        self.assertIn((SEED_NODE, "GENE:PLK1", "associate"), keys)
        self.assertIn(("GENE:PLK1", SEED_NODE, "associate"), keys)

    def test_an_unlisted_form_is_reported(self) -> None:
        """A surface form that the list file does not record is named in the summary."""
        key, value = pair("associate", "@GENE_PLK1", [item(NAMES_PHENOTYPE, "1", "a cell study")])
        verdicts = seed_filter.pair_verdicts({key: value}, self.rule)
        report = seed_filter.summary(verdicts, self.rule, set())
        self.assertIn("mitotic abnormalities", report["forms_not_in_the_list_file"])
        self.assertEqual(seed_filter.summary(verdicts, self.rule,
                                             {"mitotic abnormalities"})[
                             "forms_not_in_the_list_file"], [])


class TestBuildIntegration(unittest.TestCase):
    """Check that the graph build removes the rejected edges and keeps the option to skip."""

    def setUp(self) -> None:
        """Build the edge rows and the cached evidence of a small seed neighbourhood."""
        self.edges = [
            {"source_id": SEED_NODE, "target_id": "GENE:BUB1B", "relation": "associate",
             "evidence_count": 27, "specificity_npmi": "", "pmids": "", "source": "pubtator3"},
            {"source_id": SEED_NODE, "target_id": "GENE:BUB1B", "relation": "associate",
             "evidence_count": 27, "specificity_npmi": "", "pmids": "",
             "source": "pubtator3_consequence"},
            {"source_id": SEED_NODE, "target_id": "GENE:PLK1", "relation": "associate",
             "evidence_count": 11, "specificity_npmi": "", "pmids": "",
             "source": "pubtator3_consequence"},
            {"source_id": SEED_NODE, "target_id": "CHEMBL:CHEMBL1", "relation": "indication",
             "evidence_count": 1, "specificity_npmi": "", "pmids": "", "source": "chembl"},
            {"source_id": "GENE:BUB1B", "target_id": "GENE:SIRT2", "relation": "associate",
             "evidence_count": 5, "specificity_npmi": "", "pmids": "", "source": "pubtator3"},
        ]
        good = pair("associate", "@GENE_BUB1B", [item(NAMES_SYNDROME, "1", "a genetics study"),
                                                 item(ABBREVIATION_ONLY, "2", "BubR1 regulation")])
        bad = pair("associate", "@GENE_PLK1", [item(NAMES_PHENOTYPE, "3", "a cell study")])
        self.cache = {"fetched_at": "2026-09-15T00:00:00+00:00", "endpoint": "test",
                      "pairs": dict([good, bad])}

    def test_the_strict_filter_removes_only_the_rejected_edge(self) -> None:
        """Two rows of the PLK1 edge go and the BUB1B edge and the curated edge stay."""
        with tempfile.TemporaryDirectory() as folder:
            evidence_path = Path(folder) / "evidence.json"
            common.write_json(evidence_path, self.cache)
            kept, record = build.apply_seed_filter(
                self.edges, SEED_NODE, ENTITIES, seed_filter.DEFAULT_RULE, evidence_path,
                Path(folder) / "verdicts.tsv")
        self.assertEqual(record["seed_literature_edges_before"], 2)
        self.assertEqual(record["seed_literature_edges_after"], 1)
        self.assertEqual(record["seed_literature_edges_removed"], 1)
        self.assertEqual(record["edge_rows_removed"], 1)
        self.assertEqual(len(kept), 4)
        self.assertNotIn("GENE:PLK1", {row["target_id"] for row in kept})
        self.assertIn("CHEMBL:CHEMBL1", {row["target_id"] for row in kept})
        self.assertEqual(record["mode"], build.SEED_FILTER_STRICT)

    def test_a_missing_evidence_cache_stops_the_build(self) -> None:
        """A build with no cached sentence fails rather than keep an unreadable edge."""
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaises(FileNotFoundError):
                build.apply_seed_filter(self.edges, SEED_NODE, ENTITIES,
                                        seed_filter.DEFAULT_RULE, Path(folder) / "missing.json",
                                        Path(folder) / "verdicts.tsv")

    def test_the_seed_edge_count_reads_distinct_keys(self) -> None:
        """Two rows of one edge from two literature layers count once."""
        self.assertEqual(build.seed_edge_count(self.edges, SEED_NODE), 2)


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
