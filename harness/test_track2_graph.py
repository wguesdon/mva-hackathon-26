"""Verify the Track 2 graph parsing and assembly on small fixtures, without a network call.

The fixtures below copy the exact shape of a real PubTator3 and Open Targets response. They hold
public class A content only: gene symbols, MeSH identifiers, ChEMBL identifiers and PMIDs. The
tests check the neighbour derivation, the two tier selection, the markup removal, the identifier
mapping, the ChEMBL indication join, the specificity weight, the consequence layer and the output
tables. They never call an endpoint.
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import build_track2_graph as build  # noqa: E402
import fetch_chembl_indications as indications  # noqa: E402
import fetch_consequence_relations as consequences  # noqa: E402
import fetch_open_targets_drugs as drugs  # noqa: E402
import fetch_pubtator_entities as entities  # noqa: E402
import fetch_pubtator_evidence as evidence  # noqa: E402
import fetch_pubtator_relations as relations  # noqa: E402
import track2_common as common  # noqa: E402

# One relation record has a type, a source accession, a target accession and a publication count.
BUB1B_RELATIONS = [
    {"type": "associate", "source": "@DISEASE_Mosaic_variegated_aneuploidy_syndrome",
     "target": "@GENE_BUB1B", "publications": 27},
    {"type": "associate", "source": "@GENE_BUB1B", "target": "@GENE_CDC20", "publications": 41},
    {"type": "associate", "source": "@DISEASE_Aneuploidy", "target": "@GENE_BUB1B",
     "publications": 48},
    {"type": "associate", "source": "@GENE_AURKB", "target": "@GENE_BUB1B", "publications": 13},
    {"type": "associate", "source": "@CHEMICAL_Nocodazole", "target": "@GENE_BUB1B",
     "publications": 13},
    {"type": "negative_correlate", "source": "@GENE_RARE1", "target": "@GENE_BUB1B",
     "publications": 1},
    {"type": "associate", "source": "@DISEASE_Breast_Neoplasms", "target": "@GENE_BUB1B",
     "publications": 31},
]
TRIP13_RELATIONS = [
    {"type": "associate", "source": "@DISEASE_Mosaic_variegated_aneuploidy_syndrome",
     "target": "@GENE_TRIP13", "publications": 3},
    {"type": "associate", "source": "@GENE_TRIP13", "target": "@GENE_CDC20", "publications": 5},
    {"type": "inhibit", "source": "@GENE_TRIP13", "target": "@GENE_AURKB", "publications": 2},
    {"type": "associate", "source": "@DISEASE_Aneuploidy", "target": "@GENE_TRIP13",
     "publications": 9},
]
# A neighbour gene brings in a hub disease and a strong compound, plus a weak compound to drop.
AURKB_RELATIONS = [
    {"type": "associate", "source": "@GENE_AURKB", "target": "@GENE_CDC20", "publications": 30},
    {"type": "inhibit", "source": "@CHEMICAL_Barasertib", "target": "@GENE_AURKB",
     "publications": 25},
    {"type": "associate", "source": "@CHEMICAL_Weak_compound", "target": "@GENE_AURKB",
     "publications": 2},
    {"type": "associate", "source": "@DISEASE_Unrelated_hub_disease", "target": "@GENE_AURKB",
     "publications": 900},
    {"type": "associate", "source": "@DISEASE_Breast_Neoplasms", "target": "@GENE_AURKB",
     "publications": 40},
]
CDC20_RELATIONS = [
    {"type": "associate", "source": "@GENE_CDC20", "target": "@GENE_BUB1B", "publications": 41},
]

RELATION_CACHE = {"BUB1B": BUB1B_RELATIONS, "TRIP13": TRIP13_RELATIONS,
                  "AURKB": AURKB_RELATIONS, "CDC20": CDC20_RELATIONS}
SEEDS = ["BUB1B", "TRIP13"]

ENTITY_CACHE = {
    "@DISEASE_Mosaic_variegated_aneuploidy_syndrome": {
        "accession": "@DISEASE_Mosaic_variegated_aneuploidy_syndrome", "identifier": "C536987",
        "database": "ncbi_mesh", "name": "Mosaic variegated aneuploidy syndrome"},
    "@CHEMICAL_Nocodazole": {"accession": "@CHEMICAL_Nocodazole", "identifier": "D015739",
                             "database": "ncbi_mesh", "name": "Nocodazole"},
    "@CHEMICAL_Barasertib": {"accession": "@CHEMICAL_Barasertib", "identifier": "",
                             "database": "", "name": "Barasertib"},
    "@DISEASE_Breast_Neoplasms": {"accession": "@DISEASE_Breast_Neoplasms",
                                  "identifier": "D001943", "database": "ncbi_mesh",
                                  "name": "Breast Neoplasms", "biotype": "disease"},
}

# A search result marks each entity as `@<m>ACCESSION</m> @TYPE_IDENTIFIER @@@surface form@@@`.
TEXT_HL = ("Furthermore, mutations of @<m>GENE_BUB1B</m> @GENE_12236 @@@BubR1@@@ are associated "
           "with @<m>DISEASE_Mosaic_variegated_aneuploidy_syndrome</m> "
           "@DISEASE_MESH:C536987 @@@MVA@@@ in a @SPECIES_9606 @@@patient@@@ cohort.")
SEARCH_PAGE = {"count": 27, "total_pages": 3, "page_size": 10,
               "results": [{"pmid": 35804254, "title": "A reappraisal", "date": "2022-11-01T00:00:00Z",
                            "text_hl": TEXT_HL},
                           {"pmid": 21389084, "title": "Another paper", "date": "2011-03-01T00:00:00Z",
                            "text_hl": "@<m>GENE_BUB1B</m> @GENE_701 @@@BUB1B@@@ is mutated."},
                           {"title": "A record with no PMID", "text_hl": ""}]}

EVIDENCE_CACHE = {
    "associate|@DISEASE_Mosaic_variegated_aneuploidy_syndrome|@GENE_BUB1B": {
        "evidence": [{"pmid": "35804254", "sentence": "A sentence."},
                     {"pmid": "21389084", "sentence": "Another sentence."},
                     {"pmid": "35804254", "sentence": "A repeat of the first PMID."}]},
}

OPEN_TARGETS_PAYLOAD = {"data": {"targets": [
    {"id": "ENSG00000178999", "approvedSymbol": "AURKB",
     "pathways": [{"pathwayId": "R-HSA-2500257", "pathway": "Resolution of Sister Chromatid Cohesion",
                   "topLevelTerm": "Cell Cycle"}],
     "drugAndClinicalCandidates": {"count": 1, "rows": [
         {"maxClinicalStage": "PHASE_2_3",
          "drug": {"id": "CHEMBL415049", "name": "BARASERTIB", "drugType": "Small molecule",
                   "maximumClinicalStage": "PHASE_2_3",
                   "mechanismsOfAction": {"rows": [
                       {"actionType": "INHIBITOR", "targetName": "Aurora kinase B",
                        "mechanismOfAction": "Aurora-B inhibitor"}]}}}]}},
    {"id": "ENSG00000156970", "approvedSymbol": "BUB1B", "pathways": [],
     "drugAndClinicalCandidates": {"count": 0, "rows": []}}]}}

# A ChEMBL indication record carries a MeSH identifier, an EFO identifier and a clinical phase.
# `D001943` stands for a disease that the literature layer already holds. `D999999` does not.
CHEMBL_INDICATIONS = [
    {"molecule_chembl_id": "CHEMBL415049", "mesh_id": "D001943", "mesh_heading": "Breast Neoplasms",
     "efo_id": "EFO:0000305", "efo_term": "breast carcinoma", "max_phase_for_ind": 2.0,
     "indication_refs": [{"ref_type": "ClinicalTrials", "ref_id": "NCT00001,NCT00002"}]},
    {"molecule_chembl_id": "CHEMBL415049", "mesh_id": "D001943", "mesh_heading": "Breast Neoplasms",
     "efo_id": "EFO:0003869", "efo_term": "breast neoplasm", "max_phase_for_ind": 4.0,
     "indication_refs": [{"ref_type": "ClinicalTrials", "ref_id": "NCT00002,NCT00003"}]},
    {"molecule_chembl_id": "CHEMBL415049", "mesh_id": "D999999", "mesh_heading": "Not in the graph",
     "efo_id": "EFO:0009999", "efo_term": "absent", "max_phase_for_ind": 1.0,
     "indication_refs": []},
]

# A total publication count per entity, and a corpus size, drive the specificity weight.
# AKT1 is the hub case: it is mentioned in far more papers than BUB1B.
ENTITY_COUNTS = {
    "@GENE_BUB1B": 4301, "@GENE_TRIP13": 900, "@GENE_CDC20": 3000, "@GENE_AURKB": 2500,
    "@GENE_AKT1": 355552, "@GENE_RARE1": 40,
    "@DISEASE_Mosaic_variegated_aneuploidy_syndrome": 3734, "@DISEASE_Aneuploidy": 38482,
    "@DISEASE_Breast_Neoplasms": 400000, "@DISEASE_Unrelated_hub_disease": 500000,
    "@CHEMICAL_Nocodazole": 5000, "@CHEMICAL_Barasertib": 300,
    "@CHEMICAL_Weak_compound": 100,
}
CORPUS_SIZE = 41148566

# A consequence entity returns relation types that a gene query never returns, such as `treat`.
ANEUPLOIDY_RELATIONS = [
    {"type": "cause", "source": "@CHEMICAL_Colchicine", "target": "@DISEASE_Aneuploidy",
     "publications": 38},
    {"type": "treat", "source": "@CHEMICAL_Barasertib", "target": "@DISEASE_Aneuploidy",
     "publications": 24},
    {"type": "associate", "source": "@DISEASE_Aneuploidy", "target": "@GENE_BUB1B",
     "publications": 48},
    {"type": "cause", "source": "@CHEMICAL_Tiny_trace_compound", "target": "@DISEASE_Aneuploidy",
     "publications": 1},
    {"type": "associate", "source": "@DISEASE_Aneuploidy", "target": "@SPECIES_9606",
     "publications": 12},
]
CONSEQUENCES = [{"accession": "@DISEASE_Aneuploidy", "label": "Aneuploidy", "npmi": 0.19,
                 "seed_count": 3, "publications": 57}]

MAP_PAYLOAD = {"data": {"mapIds": {"mappings": [
    {"term": "BUB1B", "hits": [{"id": "ENSG00000156970", "name": "BUB1B", "entity": "target"},
                               {"id": "ENSG00000259288", "name": "BUB1B-PAK6", "entity": "target"}]},
    {"term": "TRP53", "hits": []}]}}}


def build_targets_cache() -> dict:
    """Parse the Open Targets fixture into the cache shape that the build script reads.

    Returns:
        The target records, keyed by gene symbol.
    """
    return {record["symbol"]: record for record in drugs.parse_targets(OPEN_TARGETS_PAYLOAD)}


class AccessionTest(unittest.TestCase):
    """Check the accession parsing and the node identifier mapping."""

    def test_split_and_kind(self) -> None:
        """An accession splits into a biotype and a body, and maps to a node kind."""
        self.assertEqual(common.split_accession("@GENE_BUB1B"), ("GENE", "BUB1B"))
        self.assertEqual(common.accession_kind("@DISEASE_Aneuploidy"), "disease")
        self.assertEqual(common.accession_kind("@CHEMICAL_Nocodazole"), "chemical")
        self.assertEqual(common.accession_kind("@NOT_AN_ACCESSION_TYPE"), "")
        self.assertEqual(common.split_accession("no prefix here"), ("", ""))

    def test_a_gene_keeps_its_symbol(self) -> None:
        """A gene node identifier uses the symbol, so it stays readable and stable."""
        self.assertEqual(common.node_id_for_accession("@GENE_BUB1B"), "GENE:BUB1B")

    def test_a_resolved_entity_uses_its_database_identifier(self) -> None:
        """A MeSH identifier enters the node identifier with its namespace."""
        self.assertEqual(
            common.node_id_for_accession("@CHEMICAL_Nocodazole", "D015739", "ncbi_mesh"),
            "CHEMICAL:MESH:D015739")

    def test_an_identifier_that_already_has_a_namespace_is_not_doubled(self) -> None:
        """An identifier such as `MESH:C536987` keeps its own namespace."""
        self.assertEqual(
            common.node_id_for_accession("@DISEASE_X", "MESH:C536987", "ncbi_mesh"),
            "DISEASE:MESH:C536987")

    def test_an_unresolved_entity_falls_back_to_the_accession(self) -> None:
        """An unresolved entity never receives an invented identifier."""
        self.assertEqual(common.node_id_for_accession("@CHEMICAL_Barasertib"),
                         "CHEMICAL:PUBTATOR:Barasertib")

    def test_an_unsupported_biotype_has_no_node(self) -> None:
        """A species accession is not a graph node kind, so it returns no identifier."""
        self.assertEqual(common.node_id_for_accession("@SPECIES_9606"), "")


class HighlightTest(unittest.TestCase):
    """Check the removal of the PubTator3 annotation markup."""

    def test_the_sentence_keeps_only_the_surface_form(self) -> None:
        """The plain sentence holds the words of the paper and no markup."""
        sentence = common.plain_sentence(TEXT_HL)
        self.assertEqual(sentence, "Furthermore, mutations of BubR1 are associated with MVA "
                                   "in a patient cohort.")
        self.assertNotIn("@", sentence)
        self.assertNotIn("<m>", sentence)

    def test_an_empty_highlight_gives_an_empty_sentence(self) -> None:
        """A result with no highlight gives an empty sentence and never raises."""
        self.assertEqual(common.plain_sentence(""), "")

    def test_the_identifiers_come_out_of_the_highlight(self) -> None:
        """Each marked entity yields its database identifier."""
        found = common.highlight_identifiers(TEXT_HL)
        self.assertEqual(found["@DISEASE_Mosaic_variegated_aneuploidy_syndrome"], "MESH:C536987")
        self.assertEqual(found["@GENE_BUB1B"], "12236")


class NeighbourDerivationTest(unittest.TestCase):
    """Check that the neighbour list comes from the data and not from a fixed list."""

    def test_the_counts_add_up_over_the_seed_genes(self) -> None:
        """A gene that both seeds link to sums both publication counts."""
        totals = relations.gene_neighbour_counts(RELATION_CACHE, SEEDS)
        self.assertEqual(totals["@GENE_CDC20"], 46)
        self.assertEqual(totals["@GENE_AURKB"], 15)
        self.assertNotIn("@GENE_BUB1B", totals)
        self.assertNotIn("@GENE_TRIP13", totals)

    def test_the_threshold_drops_a_weak_neighbour(self) -> None:
        """A gene below the threshold never enters the neighbour list."""
        chosen = relations.derive_neighbours(RELATION_CACHE, SEEDS, 5, 40)
        self.assertEqual(chosen, ["CDC20", "AURKB"])
        self.assertNotIn("RARE1", chosen)

    def test_the_cap_limits_the_neighbour_count(self) -> None:
        """The cap keeps the best supported neighbour only."""
        self.assertEqual(relations.derive_neighbours(RELATION_CACHE, SEEDS, 5, 1), ["CDC20"])

    def test_a_disease_is_never_a_gene_neighbour(self) -> None:
        """Only a gene accession enters the neighbour list."""
        chosen = relations.derive_neighbours(RELATION_CACHE, SEEDS, 1, 40)
        self.assertNotIn("Mosaic variegated aneuploidy syndrome", chosen)


class SelectionTest(unittest.TestCase):
    """Check the two tier rule that keeps the graph inside the MVA neighbourhood."""

    def setUp(self) -> None:
        """Run the selection once for every test in this class."""
        self.neighbours = ["CDC20", "AURKB"]
        self.nodes, self.edges = common.select_graph(RELATION_CACHE, SEEDS, self.neighbours, 20)

    def test_every_seed_relation_is_kept(self) -> None:
        """The first tier keeps each relation of a seed gene, however weak."""
        keys = {(e["source"], e["target"], e["type"]) for e in self.edges}
        self.assertIn(("@DISEASE_Mosaic_variegated_aneuploidy_syndrome", "@GENE_BUB1B",
                       "associate"), keys)
        self.assertIn(("@GENE_RARE1", "@GENE_BUB1B", "negative_correlate"), keys)

    def test_a_strong_neighbour_compound_enters_the_graph(self) -> None:
        """A compound above the threshold becomes a node and keeps its edge."""
        self.assertIn("@CHEMICAL_Barasertib", self.nodes)
        keys = {(e["source"], e["target"], e["type"]) for e in self.edges}
        self.assertIn(("@CHEMICAL_Barasertib", "@GENE_AURKB", "inhibit"), keys)

    def test_a_weak_neighbour_compound_is_dropped(self) -> None:
        """A compound below the threshold never enters the graph."""
        self.assertNotIn("@CHEMICAL_Weak_compound", self.nodes)

    def test_a_hub_disease_of_a_neighbour_does_not_enter(self) -> None:
        """A second tier disease is dropped, so a hub gene cannot pull in the whole of PubMed."""
        self.assertNotIn("@DISEASE_Unrelated_hub_disease", self.nodes)

    def test_a_second_tier_edge_between_two_graph_nodes_is_kept(self) -> None:
        """A relation that links two genes already in the graph closes the subgraph."""
        keys = {(e["source"], e["target"], e["type"]) for e in self.edges}
        self.assertIn(("@GENE_AURKB", "@GENE_CDC20", "associate"), keys)

    def test_a_repeated_relation_appears_once(self) -> None:
        """The same triple from two gene queries yields one record."""
        keys = [(e["source"], e["target"], e["type"]) for e in self.edges]
        self.assertEqual(len(keys), len(set(keys)))


class SearchParsingTest(unittest.TestCase):
    """Check the parsing of a search result page."""

    def test_a_result_without_a_pmid_is_dropped(self) -> None:
        """Only a result that carries a PMID becomes evidence."""
        parsed = evidence.parse_results(SEARCH_PAGE)
        self.assertEqual([item["pmid"] for item in parsed], ["35804254", "21389084"])

    def test_the_evidence_carries_a_sentence_and_a_year(self) -> None:
        """A parsed result holds the plain sentence and the publication year."""
        first = evidence.parse_results(SEARCH_PAGE)[0]
        self.assertEqual(first["year"], "2022")
        self.assertTrue(first["sentence"].startswith("Furthermore, mutations of BubR1"))

    def test_the_pair_query_has_the_endpoint_shape(self) -> None:
        """The query text matches the form that the search endpoint expects."""
        self.assertEqual(common.pair_query("@GENE_BUB1B", "@DISEASE_X", "associate"),
                         "relations:associate|@GENE_BUB1B|@DISEASE_X")

    def test_the_ranking_puts_the_best_supported_pair_first(self) -> None:
        """Without a priority set, a pair with more publications is fetched first."""
        ranked = evidence.rank_records(BUB1B_RELATIONS)
        self.assertEqual(ranked[0]["publications"], 48)
        self.assertEqual(ranked[-1]["publications"], 1)

    def test_only_a_disease_that_every_seed_gene_shares_gets_priority(self) -> None:
        """MVA links to both seed genes, so it ranks first. Breast Neoplasms links to one only."""
        chosen = evidence.priority_accessions(RELATION_CACHE, SEEDS, len(SEEDS))
        self.assertIn("@DISEASE_Mosaic_variegated_aneuploidy_syndrome", chosen)
        self.assertNotIn("@DISEASE_Breast_Neoplasms", chosen)
        self.assertIn("@GENE_BUB1B", chosen)
        self.assertIn("@GENE_TRIP13", chosen)

    def test_the_rare_disease_edge_outranks_a_hub_edge(self) -> None:
        """A small MVA literature must not lose its provenance to a large cancer literature."""
        hub = {"type": "associate", "source": "@DISEASE_Unrelated_hub_disease",
               "target": "@GENE_AURKB", "publications": 900}
        mva = {"type": "associate", "source": "@DISEASE_Mosaic_variegated_aneuploidy_syndrome",
               "target": "@GENE_BUB1B", "publications": 27}
        chosen = evidence.priority_accessions(RELATION_CACHE, SEEDS, len(SEEDS))
        ranked = evidence.rank_records([hub, mva], chosen)
        self.assertEqual(ranked[0]["source"],
                         "@DISEASE_Mosaic_variegated_aneuploidy_syndrome")

    def test_a_looser_threshold_admits_a_shared_cancer(self) -> None:
        """One seed link is enough when the caller lowers the threshold."""
        chosen = evidence.priority_accessions(RELATION_CACHE, SEEDS, 1)
        self.assertIn("@DISEASE_Breast_Neoplasms", chosen)


class EntityResolutionTest(unittest.TestCase):
    """Check that only an exact accession match becomes an identifier."""

    def test_only_the_exact_accession_is_accepted(self) -> None:
        """A near match such as `@GENE_BUB1BB` never replaces `@GENE_BUB1B`."""
        hits = [{"_id": "@GENE_BUB1BB", "db_id": "567006"},
                {"_id": "@GENE_BUB1B", "db_id": "701"}]
        self.assertEqual(entities.pick_hit("@GENE_BUB1B", hits)["db_id"], "701")

    def test_no_match_gives_no_identifier(self) -> None:
        """An accession with no exact hit returns an empty mapping."""
        self.assertEqual(entities.pick_hit("@GENE_BUB1B", [{"_id": "@GENE_OTHER"}]), {})
        self.assertEqual(entities.pick_hit("@GENE_BUB1B", []), {})

    def test_only_a_disease_or_a_chemical_needs_a_call(self) -> None:
        """A gene and a species are not sent to the autocomplete endpoint."""
        wanted = entities.accessions_to_resolve(
            {"@GENE_BUB1B", "@DISEASE_Aneuploidy", "@CHEMICAL_Nocodazole", "@SPECIES_9606"})
        self.assertEqual(wanted, ["@CHEMICAL_Nocodazole", "@DISEASE_Aneuploidy"])


class OpenTargetsParsingTest(unittest.TestCase):
    """Check the parsing of the Open Targets GraphQL responses."""

    def test_a_read_through_locus_never_replaces_the_gene(self) -> None:
        """Only a hit whose name equals the term becomes the Ensembl identifier."""
        mapping = drugs.parse_mappings(MAP_PAYLOAD)
        self.assertEqual(mapping, {"BUB1B": "ENSG00000156970"})
        self.assertNotIn("TRP53", mapping)

    def test_the_compound_row_keeps_the_phase_and_the_action(self) -> None:
        """A candidate row carries the ChEMBL identifier, the action type and the phase."""
        targets = build_targets_cache()
        drug = targets["AURKB"]["drugs"][0]
        self.assertEqual(drug["drug_id"], "CHEMBL415049")
        self.assertEqual(drug["action_type"], "inhibitor")
        self.assertEqual(drug["max_clinical_stage_on_target"], "PHASE_2_3")

    def test_a_gene_with_no_compound_gives_an_empty_list(self) -> None:
        """BUB1B has no compound, and the parser reports that without raising."""
        targets = build_targets_cache()
        self.assertEqual(targets["BUB1B"]["drugs"], [])
        self.assertEqual(targets["BUB1B"]["drug_count"], 0)

    def test_batching_covers_every_gene(self) -> None:
        """The batches hold each symbol exactly once."""
        batches = drugs.batched(["A", "B", "C", "D", "E"], 2)
        self.assertEqual(batches, [["A", "B"], ["C", "D"], ["E"]])


class IndicationParsingTest(unittest.TestCase):
    """Check the ChEMBL indication join and the merge of repeated records."""

    def test_an_indication_outside_the_graph_is_dropped(self) -> None:
        """A disease that the literature layer never saw teaches the control nothing."""
        kept = indications.keep_indications(CHEMBL_INDICATIONS, {"D001943"})
        self.assertEqual([entry["mesh_id"] for entry in kept], ["D001943"])

    def test_two_records_of_one_pair_merge_and_keep_the_highest_phase(self) -> None:
        """The same molecule and disease under two EFO terms give one edge at the highest phase."""
        kept = indications.keep_indications(CHEMBL_INDICATIONS, {"D001943"})
        self.assertEqual(len(kept), 1)
        self.assertEqual(kept[0]["max_phase_for_ind"], 4.0)

    def test_the_reference_count_counts_the_distinct_trials(self) -> None:
        """A reference field holds a comma separated trial list, and each identifier counts once."""
        self.assertEqual(indications.reference_count(CHEMBL_INDICATIONS[0]), 2)
        self.assertEqual(indications.reference_count(CHEMBL_INDICATIONS[2]), 1)

    def test_the_graph_disease_identifiers_come_from_the_entity_cache(self) -> None:
        """Only a disease resolved against MeSH can join a ChEMBL indication."""
        self.assertEqual(indications.graph_disease_mesh(ENTITY_CACHE), {"C536987", "D001943"})

    def test_the_compound_list_comes_from_the_open_targets_cache(self) -> None:
        """The indication fetch asks only about compounds that the graph already holds."""
        self.assertEqual(indications.graph_drug_ids(build_targets_cache()), ["CHEMBL415049"])


class SpecificityTest(unittest.TestCase):
    """Check the normalized pointwise mutual information weight."""

    def test_a_specific_pair_scores_above_a_hub_pair(self) -> None:
        """The MVA pair is rarer and tighter than a pair that only shares a large literature."""
        specific = common.npmi(27, 3734, 4301, CORPUS_SIZE)
        hub = common.npmi(27, 400000, 355552, CORPUS_SIZE)
        self.assertGreater(specific, hub)
        self.assertGreater(specific, 0)
        self.assertLess(hub, 0)

    def test_the_value_stays_inside_the_defined_range(self) -> None:
        """The measure is normalized, so it never leaves the range -1 to 1."""
        for pair, a, b in [(27, 3734, 4301), (1, 40, 4301), (48, 38482, 4301),
                           (2, 500000, 355552)]:
            value = common.npmi(pair, a, b, CORPUS_SIZE)
            self.assertIsNotNone(value)
            self.assertGreaterEqual(value, -1.0)
            self.assertLessEqual(value, 1.0)

    def test_a_perfect_association_approaches_one(self) -> None:
        """When every mention of both entities is shared, the weight reaches its maximum."""
        self.assertAlmostEqual(common.npmi(500, 500, 500, CORPUS_SIZE), 1.0, places=6)

    def test_an_impossible_or_missing_count_gives_no_value(self) -> None:
        """A missing count never becomes a zero weight, because zero has a meaning here."""
        self.assertIsNone(common.npmi(0, 10, 10, CORPUS_SIZE))
        self.assertIsNone(common.npmi(5, 0, 10, CORPUS_SIZE))
        self.assertIsNone(common.npmi(5, 10, 0, CORPUS_SIZE))
        self.assertIsNone(common.npmi(5, 10, 10, 0))
        self.assertIsNone(common.npmi(20, 10, 10, CORPUS_SIZE))

    def test_the_counts_of_one_pair_are_added_over_the_relation_types(self) -> None:
        """A pair with two relation types contributes both counts to one co-mention total."""
        records = [{"source": "@GENE_A", "target": "@GENE_B", "publications": 13},
                   {"source": "@GENE_B", "target": "@GENE_A", "publications": 3}]
        totals = common.pair_total(records)
        self.assertEqual(totals[("@GENE_A", "@GENE_B")], 16)

    def test_both_edges_of_a_pair_carry_the_same_weight(self) -> None:
        """Specificity describes the entity pair, so the relation type does not change it."""
        records = [{"source": "@GENE_BUB1B", "target": "@GENE_CDC20", "type": "associate",
                    "publications": 41},
                   {"source": "@GENE_BUB1B", "target": "@GENE_CDC20", "type": "inhibit",
                    "publications": 4}]
        totals = common.pair_total(records)
        first = build.specificity_for(records[0], totals, ENTITY_COUNTS, CORPUS_SIZE)
        second = build.specificity_for(records[1], totals, ENTITY_COUNTS, CORPUS_SIZE)
        self.assertEqual(first, second)
        self.assertTrue(first)

    def test_an_entity_with_no_count_gives_an_empty_weight(self) -> None:
        """An unseen entity leaves the column empty instead of taking a made up value."""
        record = {"source": "@GENE_UNCOUNTED", "target": "@GENE_BUB1B", "type": "associate",
                  "publications": 5}
        totals = common.pair_total([record])
        self.assertEqual(build.specificity_for(record, totals, ENTITY_COUNTS, CORPUS_SIZE), "")


class ConsequenceSelectionTest(unittest.TestCase):
    """Check that the consequence entities are derived from the data, not named."""

    def setUp(self) -> None:
        """Score the candidates once for every test in this class."""
        self.links = consequences.seed_disease_links(RELATION_CACHE, SEEDS)
        self.scored = consequences.score_consequences(self.links, ENTITY_COUNTS, SEEDS,
                                                      CORPUS_SIZE, 2)

    def test_only_a_disease_linked_to_enough_seeds_is_scored(self) -> None:
        """A disease that one seed gene mentions is not a shared consequence."""
        labels = {record["label"] for record in self.scored}
        self.assertIn("Mosaic variegated aneuploidy syndrome", labels)
        self.assertIn("Aneuploidy", labels)

    def test_the_mechanism_term_outranks_the_hub_term(self) -> None:
        """A term that appears beside every gene scores below a term that names the mechanism."""
        by_label = {record["label"]: record["npmi"] for record in self.scored}
        self.assertGreater(by_label["Mosaic variegated aneuploidy syndrome"],
                           by_label["Aneuploidy"])
        if "Breast Neoplasms" in by_label:
            self.assertGreater(by_label["Aneuploidy"], by_label["Breast Neoplasms"])

    def test_the_publication_floor_drops_a_coincidence(self) -> None:
        """A rare term that shares two papers scores well but must not enter the graph."""
        scored = [{"accession": "@DISEASE_Fluke", "label": "Fluke", "npmi": 0.5,
                   "publications": 2},
                  {"accession": "@DISEASE_Real", "label": "Real", "npmi": 0.2,
                   "publications": 40}]
        kept = consequences.select_consequences(scored, 0.07, 5, 12)
        self.assertEqual([record["label"] for record in kept], ["Real"])

    def test_the_specificity_threshold_drops_a_hub(self) -> None:
        """A term below the specificity threshold never becomes a consequence."""
        scored = [{"accession": "@DISEASE_Hub", "label": "Hub", "npmi": 0.01,
                   "publications": 900}]
        self.assertEqual(consequences.select_consequences(scored, 0.07, 5, 12), [])

    def test_the_limit_caps_the_set(self) -> None:
        """The cap keeps the most specific consequence only."""
        scored = [{"accession": "@DISEASE_A", "label": "A", "npmi": 0.4, "publications": 40},
                  {"accession": "@DISEASE_B", "label": "B", "npmi": 0.3, "publications": 40}]
        self.assertEqual([r["label"] for r in consequences.select_consequences(scored, 0.07, 5, 1)],
                         ["A"])


class ConsequenceLayerTest(unittest.TestCase):
    """Check the edges that the consequence layer adds to the graph."""

    def setUp(self) -> None:
        """Run the consequence selection once for every test in this class."""
        self.base, _ = common.select_graph(RELATION_CACHE, SEEDS, ["CDC20", "AURKB"], 20)
        self.added, self.records = common.select_consequence_layer(
            CONSEQUENCES, {"@DISEASE_Aneuploidy": ANEUPLOIDY_RELATIONS}, self.base, 20)

    def test_a_compound_reaches_the_graph_through_the_consequence(self) -> None:
        """Colchicine has no tie to a seed gene and enters only through aneuploidy."""
        self.assertIn("@CHEMICAL_Colchicine", self.added)
        self.assertNotIn("@CHEMICAL_Colchicine", self.base)

    def test_the_treat_relation_survives(self) -> None:
        """A disease query returns `cause` and `treat`, which a gene query never returns."""
        kinds = {record["type"] for record in self.records}
        self.assertIn("treat", kinds)
        self.assertIn("cause", kinds)

    def test_a_trace_compound_is_dropped(self) -> None:
        """A compound below the publication threshold does not enter the graph."""
        self.assertNotIn("@CHEMICAL_Tiny_trace_compound", self.added)

    def test_an_entity_already_in_the_graph_keeps_its_edge(self) -> None:
        """The link from aneuploidy back to BUB1B is kept because both nodes exist."""
        pairs = {(r["source"], r["target"]) for r in self.records}
        self.assertIn(("@DISEASE_Aneuploidy", "@GENE_BUB1B"), pairs)

    def test_a_species_never_enters(self) -> None:
        """A species is not a graph node kind, so its relation is dropped."""
        self.assertNotIn("@SPECIES_9606", self.added)

    def test_no_consequence_cache_gives_an_empty_layer(self) -> None:
        """A build that runs before the consequence fetch produces no consequence edge."""
        with tempfile.TemporaryDirectory() as folder:
            self.assertEqual(common.consequence_nodes(Path(folder), 20), set())
            self.assertEqual(common.consequence_records(Path(folder), 20), [])
            self.assertEqual(common.consequence_accessions(Path(folder)), set())


class GraphAssemblyTest(unittest.TestCase):
    """Check the node table, the edge table and the manifest counts."""

    def setUp(self) -> None:
        """Assemble the graph once from the fixtures."""
        self.neighbours = ["CDC20", "AURKB"]
        accessions, selected = common.select_graph(RELATION_CACHE, SEEDS, self.neighbours, 20)
        self.nodes = build.build_literature_nodes(accessions, ENTITY_CACHE)
        self.totals = common.pair_total(selected)
        self.edges = build.build_literature_edges(selected, self.nodes, ENTITY_CACHE,
                                                  EVIDENCE_CACHE, self.totals, ENTITY_COUNTS,
                                                  CORPUS_SIZE)
        targets = build_targets_cache()
        drug_edges, self.phases = build.build_drug_layer(targets, self.nodes)
        pathway_edges = build.build_pathway_layer(targets, self.nodes)
        kept = indications.keep_indications(CHEMBL_INDICATIONS, {"D001943"})
        self.indication_edges, self.indication_phases = build.build_indication_layer(
            kept, self.nodes)
        self.edges.extend(drug_edges)
        self.edges.extend(pathway_edges)
        self.edges.extend(self.indication_edges)

    def test_the_node_kinds_are_inside_the_contract(self) -> None:
        """Every node kind is one of gene, disease, chemical or pathway."""
        kinds = {row["kind"] for row in self.nodes.values()}
        self.assertTrue(kinds <= {"gene", "disease", "chemical", "pathway"}, kinds)

    def test_a_species_never_becomes_a_node(self) -> None:
        """No node identifier starts with a species prefix."""
        self.assertFalse([n for n in self.nodes if n.startswith("SPECIES")])

    def test_the_disease_node_carries_its_mesh_identifier(self) -> None:
        """The MVA disease node uses the identifier that the entity endpoint returned."""
        node = self.nodes["DISEASE:MESH:C536987"]
        self.assertEqual(node["kind"], "disease")
        self.assertEqual(node["source"], "pubtator3")
        self.assertEqual(node["label"], "Mosaic variegated aneuploidy syndrome")

    def test_the_pmids_reach_the_edge_and_repeat_once(self) -> None:
        """An edge carries the PMIDs of its pair, with no duplicate."""
        edge = [e for e in self.edges
                if e["source_id"] == "DISEASE:MESH:C536987" and e["target_id"] == "GENE:BUB1B"][0]
        self.assertEqual(edge["pmids"], "35804254,21389084")
        self.assertEqual(edge["evidence_count"], 27)
        self.assertEqual(edge["source"], "pubtator3")

    def test_an_edge_with_no_evidence_record_has_an_empty_pmid_field(self) -> None:
        """A pair that was never fetched keeps an empty PMID list instead of an invented one."""
        edge = [e for e in self.edges
                if e["source_id"] == "GENE:BUB1B" and e["target_id"] == "GENE:CDC20"][0]
        self.assertEqual(edge["pmids"], "")

    def test_the_compound_edge_points_at_the_gene(self) -> None:
        """An Open Targets compound edge runs from the compound to its target gene."""
        edge = [e for e in self.edges if e["source_id"] == "CHEMBL:CHEMBL415049"][0]
        self.assertEqual(edge["target_id"], "GENE:AURKB")
        self.assertEqual(edge["relation"], "inhibitor")
        self.assertEqual(edge["source"], "open_targets")
        self.assertEqual(self.nodes["CHEMBL:CHEMBL415049"]["kind"], "chemical")

    def test_a_compound_of_a_gene_outside_the_graph_is_dropped(self) -> None:
        """A compound only enters through a gene that the graph already holds."""
        targets = build_targets_cache()
        targets["NOT_IN_GRAPH"] = {"symbol": "NOT_IN_GRAPH", "drugs": [
            {"drug_id": "CHEMBL999", "drug_name": "Outsider", "action_type": "inhibitor"}],
            "pathways": []}
        extra_edges, _ = build.build_drug_layer(targets, dict(self.nodes))
        self.assertFalse([e for e in extra_edges if e["source_id"] == "CHEMBL:CHEMBL999"])

    def test_the_pathway_edge_and_node_exist(self) -> None:
        """A Reactome pathway becomes a node and a gene links to it."""
        self.assertEqual(self.nodes["PATHWAY:R-HSA-2500257"]["kind"], "pathway")
        edge = [e for e in self.edges if e["target_id"] == "PATHWAY:R-HSA-2500257"][0]
        self.assertEqual(edge["source_id"], "GENE:AURKB")
        self.assertEqual(edge["relation"], "participates_in")

    def test_the_clinical_phase_reaches_the_side_table(self) -> None:
        """The edge schema has no phase column, so the phase goes to drug_evidence.tsv."""
        row = [p for p in self.phases if p["drug_id"] == "CHEMBL:CHEMBL415049"][0]
        self.assertEqual(row["max_clinical_stage_on_target"], "PHASE_2_3")
        self.assertEqual(row["target_gene_id"], "GENE:AURKB")

    def test_the_indication_edge_runs_from_the_disease_to_the_compound(self) -> None:
        """The control looks for an edge with the relation `indication` on a disease node."""
        edge = [e for e in self.edges if e["relation"] == "indication"][0]
        self.assertEqual(edge["source_id"], "DISEASE:MESH:D001943")
        self.assertEqual(edge["target_id"], "CHEMBL:CHEMBL415049")
        self.assertEqual(edge["source"], "chembl")
        self.assertEqual(edge["evidence_count"], 3)
        self.assertEqual(edge["pmids"], "")

    def test_the_indication_phase_reaches_the_side_table(self) -> None:
        """The clinical phase of an indication goes to indication_evidence.tsv."""
        row = self.indication_phases[0]
        self.assertEqual(row["disease_id"], "DISEASE:MESH:D001943")
        self.assertEqual(row["max_phase_for_ind"], 4.0)

    def test_the_case_keeps_a_gene_route_after_the_direct_edge_is_removed(self) -> None:
        """Breast Neoplasms reaches BARASERTIB through AURKB, so the held-out case is eligible."""
        routes = build.gene_route_cases(self.edges, self.indication_edges)
        self.assertEqual(routes["cases_with_a_gene_route"], 1)
        self.assertEqual(routes["cases_without_a_gene_route"], 0)
        self.assertEqual(routes["diseases_with_a_recoverable_case"], 1)

    def test_a_case_with_no_shared_gene_is_not_eligible(self) -> None:
        """A compound that shares no gene with the disease is unreachable once the edge goes."""
        orphan = [{"source_id": "DISEASE:MESH:C536987", "target_id": "CHEMBL:CHEMBL415049",
                   "relation": "indication", "evidence_count": 1, "pmids": "", "source": "chembl"}]
        routes = build.gene_route_cases(self.edges + orphan, orphan)
        self.assertEqual(routes["cases_with_a_gene_route"], 0)
        self.assertEqual(routes["diseases_with_a_recoverable_case"], 0)

    def test_an_indication_on_a_missing_node_never_becomes_an_edge(self) -> None:
        """A compound that the graph does not hold produces no indication edge."""
        absent = [{"molecule_chembl_id": "CHEMBL000", "mesh_id": "D001943",
                   "max_phase_for_ind": 4.0, "reference_count": 3}]
        made, _ = build.build_indication_layer(absent, dict(self.nodes))
        self.assertEqual(made, [])

    def test_every_edge_endpoint_is_a_node(self) -> None:
        """The graph holds no edge that points at a missing node."""
        for edge in self.edges:
            self.assertIn(edge["source_id"], self.nodes, edge)
            self.assertIn(edge["target_id"], self.nodes, edge)


class TableWritingTest(unittest.TestCase):
    """Check that the written tables keep the contract shape."""

    def test_the_header_and_the_column_count_hold(self) -> None:
        """Each row has exactly as many columns as the header."""
        rows = [{"node_id": "GENE:BUB1B", "kind": "gene", "label": "BUB1B", "source": "pubtator3"}]
        with tempfile.TemporaryDirectory() as folder:
            path = common.write_tsv(Path(folder) / "nodes.tsv", build.NODE_HEADER, rows)
            lines = path.read_text().strip().splitlines()
        self.assertEqual(lines[0].split("\t"), build.NODE_HEADER)
        self.assertEqual(len(lines[1].split("\t")), len(build.NODE_HEADER))

    def test_a_tab_inside_a_label_never_adds_a_column(self) -> None:
        """A label that holds a tab or a newline stays inside one field."""
        rows = [{"node_id": "DISEASE:MESH:X", "kind": "disease",
                 "label": "a\tlabel\nwith breaks", "source": "pubtator3"}]
        with tempfile.TemporaryDirectory() as folder:
            path = common.write_tsv(Path(folder) / "nodes.tsv", build.NODE_HEADER, rows)
            lines = path.read_text().strip().splitlines()
        self.assertEqual(len(lines), 2)
        self.assertEqual(len(lines[1].split("\t")), len(build.NODE_HEADER))

    def test_the_counts_group_the_rows(self) -> None:
        """The manifest counter groups the rows of one column."""
        rows = [{"kind": "gene"}, {"kind": "gene"}, {"kind": "disease"}]
        self.assertEqual(build.count_by(rows, "kind"), {"gene": 2, "disease": 1})

    def test_the_digest_changes_with_the_content(self) -> None:
        """The SHA-256 of a table follows its content."""
        with tempfile.TemporaryDirectory() as folder:
            first = Path(folder) / "a.tsv"
            first.write_text("one")
            digest = common.sha256_of(first)
            first.write_text("two")
            self.assertNotEqual(digest, common.sha256_of(first))
            self.assertEqual(common.sha256_of(Path(folder) / "missing.tsv"), "")


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
