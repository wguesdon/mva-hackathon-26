"""Verify the mention rule that runs on every disease node of a graph, without a network call.

`scripts/track2_graph_mention_rule.py` derives one rule per disease node from the MeSH entry terms
of its concept, decides every cached sentence of every literature edge of that node, and removes the
edges that the rule rejects. ADR 0010 fixes its scope and its treatment of an edge that carries no
cached sentence.

The tests below check five things. The node table and the edge table are read into one group per
disease node. The rule of one node decides the edges of that node and no other node. A disease node
with no rule keeps every edge. The summary of a node gives the share of its edges that no sentence
supports. The drop of an edge whose only rule-matched mentions are uncorroborated abbreviations is
counted separately, because that class carries the false reject cost.

Every fixture holds public class A content only: MeSH entry terms, MeSH identifiers, gene symbols
and PMIDs. No network call runs, because `build_rule` takes the entry terms as an argument.
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import track2_common as common  # noqa: E402
import track2_disease_surface_forms as surface  # noqa: E402
import track2_graph_mention_rule as graph_rule  # noqa: E402
import track2_seed_filter as seed_filter  # noqa: E402

SEED_NODE = "DISEASE:MESH:C536987"
SEED_ACCESSION = "@DISEASE_Mosaic_variegated_aneuploidy_syndrome"
OTHER_NODE = "DISEASE:MESH:D049932"
OTHER_ACCESSION = "@DISEASE_Nijmegen_Breakage_Syndrome"
NO_RULE_NODE = "DISEASE:PUBTATOR:176430"
NO_RULE_ACCESSION = "@DISEASE_176430"

# The `ds_meshterms` field of MeSH C536987, retrieved on 2026-09-17.
MVA_TERMS = ["Mosaic variegated aneuploidy syndrome", "MVA Syndrome",
             "Mosaic variegated aneuploidy syndrome 1",
             "Instability mitotic non disjunction syndrome",
             "Chromosomal mosaicism due to mitotic instability"]
# A part of the `ds_meshterms` field of MeSH D049932, retrieved on 2026-09-17.
NBS_TERMS = ["Nijmegen Breakage Syndrome", "Breakage Syndrome, Nijmegen",
             "Berlin Breakage Syndrome", "Seemanova Syndrome 2"]

ENTITIES = {
    SEED_ACCESSION: {"identifier": "C536987", "database": "ncbi_mesh",
                     "name": "Mosaic variegated aneuploidy syndrome", "biotype": "disease"},
    OTHER_ACCESSION: {"identifier": "D049932", "database": "ncbi_mesh",
                      "name": "Nijmegen Breakage Syndrome", "biotype": "disease"},
    NO_RULE_ACCESSION: {"identifier": "", "database": "", "name": "176430", "biotype": "disease"},
}

NAMES_MVA = ("Biallelic @<m>GENE_BUB1B</m> @GENE_701 @@@BUB1B@@@ variants cause "
             "@<m>DISEASE_Mosaic_variegated_aneuploidy_syndrome</m> @DISEASE_MESH:C536987 "
             "@@@mosaic variegated aneuploidy syndrome@@@.")
NAMES_PHENOTYPE = ("Cells showed @<m>DISEASE_Mosaic_variegated_aneuploidy_syndrome</m> "
                   "@DISEASE_MESH:C536987 @@@mitotic abnormalities@@@ after treatment.")
BARE_ABBREVIATION = ("Checkpoint work in @<m>DISEASE_Mosaic_variegated_aneuploidy_syndrome</m> "
                     "@DISEASE_MESH:C536987 @@@MVA@@@ continues.")
NAMES_NBS = ("Loss of @<m>GENE_NBN</m> @GENE_4683 @@@NBN@@@ causes "
             "@<m>DISEASE_Nijmegen_Breakage_Syndrome</m> @DISEASE_MESH:D049932 "
             "@@@Nijmegen breakage syndrome@@@.")

NODE_ROWS = [
    {"node_id": SEED_NODE, "kind": "disease", "label": "Mosaic variegated aneuploidy syndrome",
     "source": "pubtator3"},
    {"node_id": OTHER_NODE, "kind": "disease", "label": "Nijmegen Breakage Syndrome",
     "source": "pubtator3"},
    {"node_id": NO_RULE_NODE, "kind": "disease", "label": "176430", "source": "pubtator3"},
    {"node_id": "GENE:BUB1B", "kind": "gene", "label": "BUB1B", "source": "pubtator3"},
]

EDGE_ROWS = [
    {"source_id": SEED_NODE, "target_id": "GENE:BUB1B", "relation": "associate",
     "evidence_count": 27, "specificity_npmi": "", "pmids": "", "source": "pubtator3"},
    {"source_id": SEED_NODE, "target_id": "GENE:BUB1B", "relation": "associate",
     "evidence_count": 27, "specificity_npmi": "", "pmids": "",
     "source": "pubtator3_consequence"},
    {"source_id": SEED_NODE, "target_id": "GENE:PLK1", "relation": "associate",
     "evidence_count": 11, "specificity_npmi": "", "pmids": "", "source": "pubtator3"},
    {"source_id": SEED_NODE, "target_id": "GENE:CDC20", "relation": "associate",
     "evidence_count": 4, "specificity_npmi": "", "pmids": "", "source": "pubtator3"},
    {"source_id": OTHER_NODE, "target_id": "GENE:NBN", "relation": "associate",
     "evidence_count": 60, "specificity_npmi": "", "pmids": "", "source": "pubtator3"},
    {"source_id": OTHER_NODE, "target_id": "GENE:TP53", "relation": "associate",
     "evidence_count": 9, "specificity_npmi": "", "pmids": "", "source": "pubtator3"},
    {"source_id": NO_RULE_NODE, "target_id": "GENE:BUB1B", "relation": "associate",
     "evidence_count": 3, "specificity_npmi": "", "pmids": "", "source": "pubtator3"},
    {"source_id": SEED_NODE, "target_id": "CHEMBL:CHEMBL1", "relation": "indication",
     "evidence_count": 1, "specificity_npmi": "", "pmids": "", "source": "chembl"},
]


def item(text_hl: str, pmid: str, title: str) -> dict:
    """Build one cached search result.

    Args:
        text_hl: The marked passage.
        pmid: The PubMed identifier.
        title: The title of the publication.

    Returns:
        The item record.
    """
    return {"pmid": pmid, "title": title, "text_hl": text_hl}


def pair(relation: str, accession: str, partner: str, items: list[dict]) -> tuple[str, dict]:
    """Build one cached entity pair.

    Args:
        relation: The relation name.
        accession: The accession of the disease node.
        partner: The accession of the partner entity.
        items: The items of the pair.

    Returns:
        The cache key and the cache value.
    """
    key = f"{relation}|{accession}|{partner}|{graph_rule.MAX_PAGES_ALL}"
    return key, {"search_count": len(items), "items": items}


def write_rule(folder: Path, node: str, accession: str, mesh_id: str,
               terms: list[str]) -> Path:
    """Write the rule file of one disease node from a term list.

    Args:
        folder: The rule directory.
        node: The graph node identifier.
        accession: The PubTator3 accession.
        mesh_id: The MeSH identifier.
        terms: The MeSH entry terms.

    Returns:
        The rule file path.
    """
    record = surface.build_rule(mesh_id, node, accession, terms[0], terms, "2026-09-18T00:00:00Z")
    path = graph_rule.rule_path_for(folder, node)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(record, sort_keys=False))
    return path


def write_graph(folder: Path) -> Path:
    """Write the node table and the edge table of the fixture graph.

    Args:
        folder: The graph directory.

    Returns:
        The graph directory.
    """
    common.write_tsv(folder / "nodes.tsv", ["node_id", "kind", "label", "source"], NODE_ROWS)
    common.write_tsv(folder / "edges.tsv",
                     ["source_id", "target_id", "relation", "evidence_count", "specificity_npmi",
                      "pmids", "source"], EDGE_ROWS)
    return folder


def write_caches(folder: Path) -> Path:
    """Write the sentence cache of each disease node of the fixture.

    Args:
        folder: The cache directory.

    Returns:
        The cache directory.
    """
    folder.mkdir(parents=True, exist_ok=True)
    seed_pairs = dict([
        pair("associate", SEED_ACCESSION, "@GENE_BUB1B",
             [item(NAMES_MVA, "1", "a genetics study"),
              item(NAMES_PHENOTYPE, "2", "a cell study")]),
        pair("associate", SEED_ACCESSION, "@GENE_PLK1",
             [item(NAMES_PHENOTYPE, "3", "a cell study")]),
        pair("associate", SEED_ACCESSION, "@GENE_CDC20",
             [item(BARE_ABBREVIATION, "4", "Checkpoint silencing")]),
    ])
    other_pairs = dict([
        pair("associate", OTHER_ACCESSION, "@GENE_NBN",
             [item(NAMES_NBS, "5", "a genetics study")]),
    ])
    common.write_json(graph_rule.cache_path_for(folder, SEED_NODE),
                      {"fetched_at": "2026-09-18T00:00:00+00:00", "endpoint": "test",
                       "node": SEED_NODE, "unresolved_edges": [], "refused_pairs": [],
                       "pairs": seed_pairs})
    common.write_json(graph_rule.cache_path_for(folder, OTHER_NODE),
                      {"fetched_at": "2026-09-18T00:00:00+00:00", "endpoint": "test",
                       "node": OTHER_NODE, "unresolved_edges": [], "refused_pairs": [],
                       "pairs": other_pairs})
    return folder


class TestIdentifiers(unittest.TestCase):
    """Check the two mappings that name a rule file and a cache file."""

    def test_the_mesh_identifier_comes_from_the_node(self) -> None:
        """A MeSH node gives its identifier and a PubTator node gives none."""
        self.assertEqual(graph_rule.mesh_id_of(SEED_NODE), "C536987")
        self.assertEqual(graph_rule.mesh_id_of(OTHER_NODE), "D049932")
        self.assertEqual(graph_rule.mesh_id_of(NO_RULE_NODE), "")

    def test_the_file_stem_holds_no_colon(self) -> None:
        """A node identifier becomes a file name with no separator character."""
        self.assertEqual(graph_rule.safe_stem(SEED_NODE), "disease_mesh_c536987")


class TestGraphReading(unittest.TestCase):
    """Check that the tables are read into one group per disease node."""

    def test_only_the_disease_rows_are_read(self) -> None:
        """The gene row of the node table is not a disease node."""
        with tempfile.TemporaryDirectory() as folder:
            graph = write_graph(Path(folder))
            nodes = graph_rule.disease_nodes(graph / "nodes.tsv")
        self.assertEqual(set(nodes), {SEED_NODE, OTHER_NODE, NO_RULE_NODE})

    def test_a_curated_edge_and_a_second_layer_row_count_once(self) -> None:
        """The `indication` row comes from ChEMBL, and the two literature rows of the BUB1B edge
        are one distinct source, target and relation key."""
        with tempfile.TemporaryDirectory() as folder:
            graph = write_graph(Path(folder))
            nodes = graph_rule.disease_nodes(graph / "nodes.tsv")
            by_node = graph_rule.literature_edges_by_node(graph / "edges.tsv", nodes)
        self.assertEqual(len(by_node[SEED_NODE]), 3)
        self.assertEqual(len(by_node[OTHER_NODE]), 2)
        self.assertEqual(len(by_node[NO_RULE_NODE]), 1)
        self.assertNotIn("indication", {row["relation"] for row in by_node[SEED_NODE]})


class TestGraphFilter(unittest.TestCase):
    """Check what the rule of each node removes from the whole edge set."""

    def setUp(self) -> None:
        """Write the fixture graph, the two rule files and the two sentence caches."""
        self.folder = tempfile.TemporaryDirectory()
        base = Path(self.folder.name)
        self.graph = write_graph(base / "graph")
        self.rules = base / "rules"
        write_rule(self.rules, SEED_NODE, SEED_ACCESSION, "C536987", MVA_TERMS)
        write_rule(self.rules, OTHER_NODE, OTHER_ACCESSION, "D049932", NBS_TERMS)
        self.cache = write_caches(base / "cache")

    def tearDown(self) -> None:
        """Remove the fixture directory."""
        self.folder.cleanup()

    def test_each_node_is_decided_by_its_own_rule(self) -> None:
        """The seed loses two edges and the second disease node keeps its named edge."""
        kept, record = graph_rule.apply_graph_filter(list(EDGE_ROWS), self.rules, self.cache,
                                                     ENTITIES, self.graph)
        self.assertEqual(record["mode"], "graph")
        self.assertEqual(record["disease_nodes"], 3)
        self.assertEqual(record["disease_nodes_with_a_rule"], 2)
        self.assertEqual(record["literature_edges_decided"], 4)
        self.assertEqual(record["literature_edges_kept"], 2)
        self.assertEqual(record["literature_edges_dropped"], 2)
        targets = {(row["source_id"], row["target_id"]) for row in kept}
        self.assertIn((SEED_NODE, "GENE:BUB1B"), targets)
        self.assertIn((OTHER_NODE, "GENE:NBN"), targets)
        self.assertNotIn((SEED_NODE, "GENE:PLK1"), targets)
        self.assertNotIn((SEED_NODE, "GENE:CDC20"), targets)

    def test_an_edge_with_no_cached_sentence_stays(self) -> None:
        """The TP53 edge of the second node has no cached pair, so the build keeps it."""
        kept, _ = graph_rule.apply_graph_filter(list(EDGE_ROWS), self.rules, self.cache,
                                                ENTITIES, self.graph)
        self.assertIn((OTHER_NODE, "GENE:TP53"),
                      {(row["source_id"], row["target_id"]) for row in kept})

    def test_a_node_with_no_rule_keeps_every_edge(self) -> None:
        """A disease node whose identifier carries no MeSH concept is not filtered."""
        kept, _ = graph_rule.apply_graph_filter(list(EDGE_ROWS), self.rules, self.cache,
                                                ENTITIES, self.graph)
        self.assertIn((NO_RULE_NODE, "GENE:BUB1B"),
                      {(row["source_id"], row["target_id"]) for row in kept})

    def test_the_curated_layer_is_untouched(self) -> None:
        """The ChEMBL indication row of the seed survives the rule."""
        kept, record = graph_rule.apply_graph_filter(list(EDGE_ROWS), self.rules, self.cache,
                                                     ENTITIES, self.graph)
        self.assertIn("chembl", {row["source"] for row in kept})
        self.assertEqual(record["edge_rows_removed"], len(EDGE_ROWS) - len(kept))

    def test_a_missing_rule_directory_stops_the_build(self) -> None:
        """A build with no rule file fails rather than keep every unread edge in silence."""
        with tempfile.TemporaryDirectory() as empty:
            with self.assertRaises(FileNotFoundError):
                graph_rule.apply_graph_filter(list(EDGE_ROWS), Path(empty), self.cache,
                                              ENTITIES, self.graph)


class TestNodeSummary(unittest.TestCase):
    """Check the per-node measurement and the distribution that the report reads."""

    def setUp(self) -> None:
        """Write the fixture and decide the seed node."""
        self.folder = tempfile.TemporaryDirectory()
        base = Path(self.folder.name)
        self.graph = write_graph(base / "graph")
        self.rules = base / "rules"
        write_rule(self.rules, SEED_NODE, SEED_ACCESSION, "C536987", MVA_TERMS)
        self.other_rules = self.rules
        write_rule(self.other_rules, OTHER_NODE, OTHER_ACCESSION, "D049932", NBS_TERMS)
        self.cache = write_caches(base / "cache")
        self.rule = seed_filter.load_rule(graph_rule.rule_path_for(self.rules, SEED_NODE))
        cached = common.read_json(graph_rule.cache_path_for(self.cache, SEED_NODE))
        self.verdicts = graph_rule.node_pair_verdicts(cached, self.rule)

    def tearDown(self) -> None:
        """Remove the fixture directory."""
        self.folder.cleanup()

    def test_the_share_before_the_rule_counts_the_decided_edges(self) -> None:
        """The seed holds three distinct literature edges and two of them lose every sentence."""
        nodes = graph_rule.disease_nodes(self.graph / "nodes.tsv")
        by_node = graph_rule.literature_edges_by_node(self.graph / "edges.tsv", nodes)
        row = graph_rule.node_summary(SEED_NODE, nodes[SEED_NODE], by_node[SEED_NODE],
                                      self.rule, self.verdicts)
        self.assertEqual(row["literature_edges"], 3)
        self.assertEqual(row["edges_with_sentence"], 3)
        self.assertEqual(row["edges_without_sentence"], 0)
        self.assertEqual(row["edges_kept"], 1)
        self.assertAlmostEqual(row["share_unsupported_before"], 0.6667, places=4)
        self.assertEqual(row["share_unsupported_after"], 0.0)

    def test_an_edge_with_no_cached_sentence_leaves_the_denominator(self) -> None:
        """The second disease node carries two edges and the rule read one of them."""
        nodes = graph_rule.disease_nodes(self.graph / "nodes.tsv")
        by_node = graph_rule.literature_edges_by_node(self.graph / "edges.tsv", nodes)
        rule = seed_filter.load_rule(graph_rule.rule_path_for(self.other_rules, OTHER_NODE))
        cached = common.read_json(graph_rule.cache_path_for(self.cache, OTHER_NODE))
        row = graph_rule.node_summary(OTHER_NODE, nodes[OTHER_NODE], by_node[OTHER_NODE], rule,
                                      graph_rule.node_pair_verdicts(cached, rule))
        self.assertEqual(row["literature_edges"], 2)
        self.assertEqual(row["edges_with_sentence"], 1)
        self.assertEqual(row["edges_without_sentence"], 1)
        self.assertEqual(row["share_unsupported_before"], 0.0)

    def test_the_abbreviation_drop_is_counted_on_its_own(self) -> None:
        """The CDC20 edge falls only on an abbreviation with no long form and the PLK1 edge does
        not."""
        records = {key[2]: value for key, value in self.verdicts.items()}
        self.assertTrue(graph_rule.abbreviation_only(records["@GENE_CDC20"]))
        self.assertFalse(graph_rule.abbreviation_only(records["@GENE_PLK1"]))
        self.assertFalse(graph_rule.abbreviation_only(records["@GENE_BUB1B"]))

    def test_the_distribution_carries_the_quartiles(self) -> None:
        """A list of rates gives the median, the two quartiles and the two extremes."""
        block = graph_rule.quantile_block([0.0, 0.25, 0.5, 0.75, 1.0])
        self.assertEqual(block["n"], 5)
        self.assertEqual(block["min"], 0.0)
        self.assertEqual(block["max"], 1.0)
        self.assertEqual(block["median"], 0.5)
        self.assertEqual(block["q1"], 0.25)
        self.assertEqual(block["q3"], 0.75)
        self.assertEqual(graph_rule.quantile_block([]), {"n": 0})


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
