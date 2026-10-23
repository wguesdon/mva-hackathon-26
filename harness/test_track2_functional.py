"""Verify the optional STRING functional gene layer of the Track 2 graph, without a network call.

The fixtures hold public class A content only: gene symbols, STRING protein identifiers and one
channel row copied from the STRING v12.0 human file. The tests check the channel combination, the
symbol mapping, the neighbour derivation, the edge table, the three gene layers of the ranker, and
that the default layer builds the same graph as before the option existed.
"""

from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import build_track2_functional_layer as layer  # noqa: E402
import track2_functional_compare as compare  # noqa: E402
import track2_rank as ranker  # noqa: E402

LINKS_HEADER = ("protein1 protein2 neighborhood fusion cooccurence coexpression experimental "
                "database textmining combined_score\n")
LINKS_ROWS = [
    "9606.P_BUB1B 9606.P_CDC20 0 0 0 300 950 900 990 999",
    "9606.P_CDC20 9606.P_BUB1B 0 0 0 300 950 900 990 999",
    "9606.P_BUB1B 9606.P_MAD2 0 0 0 0 800 0 900 850",
    "9606.P_MAD2 9606.P_BUB1B 0 0 0 0 800 0 900 850",
    "9606.P_BUB1B 9606.P_TEXT 0 0 0 0 0 0 990 990",
    "9606.P_TEXT 9606.P_BUB1B 0 0 0 0 0 0 990 990",
    "9606.P_TRIP13 9606.P_MAD2 0 0 0 0 0 900 0 900",
    "9606.P_MAD2 9606.P_TRIP13 0 0 0 0 0 900 0 900",
]

NODES = """node_id\tkind\tlabel\tsource
DISEASE:D1\tdisease\tseed disease\tpubtator3
GENE:BUB1B\tgene\tBUB1B\tpubtator3
GENE:CDC20\tgene\tCDC20\tpubtator3
GENE:AKT1\tgene\tAKT1\tpubtator3
CHEMBL:C1\tchemical\tdrug one\topen_targets
"""
EDGES = """source_id\ttarget_id\trelation\tevidence_count\tspecificity_npmi\tpmids\tsource
DISEASE:D1\tGENE:BUB1B\tassociate\t20\t0.3\tP1,P2\tpubtator3
GENE:BUB1B\tGENE:AKT1\tassociate\t40\t0.01\tP3\tpubtator3
GENE:BUB1B\tGENE:CDC20\tinteract\t9\t0.2\tP4\tpubtator3_consequence
GENE:AKT1\tCHEMBL:C1\tinhibitor\t1\t\t\topen_targets
"""
FUNCTIONAL_NODES = """node_id\tkind\tlabel\tsource
GENE:MAD1L1\tgene\tMAD1L1\tstring_v12_0
"""
FUNCTIONAL_EDGES = """source_id\ttarget_id\trelation\tevidence_count\tspecificity_npmi\tpmids\tsource\tstring_expdb_score
GENE:BUB1B\tGENE:CDC20\tfunctional_interaction\t1\t\t\tstring_v12_0_expdb\t995
GENE:BUB1B\tGENE:MAD1L1\tfunctional_interaction\t1\t\t\tstring_v12_0_expdb\t960
"""


def write_tables(folder: Path) -> tuple[Path, Path, Path]:
    """Write the default graph fixture and the functional layer fixture.

    Args:
        folder: The output folder.

    Returns:
        The nodes path, the edges path and the functional layer directory.
    """
    functional = folder / "functional"
    functional.mkdir(parents=True, exist_ok=True)
    (folder / "nodes.tsv").write_text(NODES)
    (folder / "edges.tsv").write_text(EDGES)
    (functional / "nodes.tsv").write_text(FUNCTIONAL_NODES)
    (functional / "edges.tsv").write_text(FUNCTIONAL_EDGES)
    return folder / "nodes.tsv", folder / "edges.tsv", functional


class ChannelCombinationTest(unittest.TestCase):
    """The STRING rule combines channel scores with the prior correction."""

    def test_no_evidence_gives_zero(self) -> None:
        """A pair with no channel evidence scores zero."""
        self.assertEqual(layer.combine_channels([0, 0]), 0)

    def test_certain_channel_gives_one_thousand(self) -> None:
        """A channel score of 1000 makes the combined score 1000."""
        self.assertEqual(layer.combine_channels([1000, 0]), 1000)

    def test_real_string_row_is_reproduced(self) -> None:
        """All seven channels of a STRING v12.0 row reproduce its combined score within rounding."""
        self.assertLessEqual(abs(layer.combine_channels([0, 0, 0, 45, 134, 0, 81]) - 173), 2)

    def test_two_equal_channels(self) -> None:
        """Two channels of 500 give 739 under the rule."""
        self.assertEqual(layer.combine_channels([500, 500]), 739)


class MappingTest(unittest.TestCase):
    """A symbol maps through the preferred name first and through a unique alias second."""

    def test_mapping_tiers(self) -> None:
        """A preferred name, a previous symbol, an ambiguous alias and an unknown symbol."""
        names = {"9606.P1": "SEPTIN9", "9606.P2": "BUB1B", "9606.P3": "DUP", "9606.P4": "DUP"}
        aliases = {"current_symbol": {}, "previous_symbol": {"SEPT9": {"9606.P1"},
                                                             "AMBIG": {"9606.P2", "9606.P3"}}}
        mapped, method, unmapped = layer.map_symbols(["BUB1B", "SEPT9", "AMBIG", "DUP", "X"],
                                                     names, aliases)
        self.assertEqual(mapped, {"BUB1B": "9606.P2", "SEPT9": "9606.P1"})
        self.assertEqual(method, {"BUB1B": "preferred_name", "SEPT9": "previous_symbol"})
        self.assertEqual(unmapped, ["AMBIG", "DUP", "X"])


class LinksTest(unittest.TestCase):
    """The layer keeps one unordered pair and ignores the text mining channel."""

    def setUp(self) -> None:
        """Write a small gzip links file."""
        self.folder = tempfile.TemporaryDirectory()
        self.path = Path(self.folder.name) / "links.txt.gz"
        with gzip.open(self.path, "wt") as handle:
            handle.write(LINKS_HEADER + "\n".join(LINKS_ROWS) + "\n")

    def tearDown(self) -> None:
        """Remove the temporary directory."""
        self.folder.cleanup()

    def test_text_mining_pair_is_dropped(self) -> None:
        """A pair that only text mining supports never reaches the threshold."""
        links = layer.read_functional_links(self.path)
        self.assertEqual(set(links), {("9606.P_BUB1B", "9606.P_CDC20"),
                                      ("9606.P_BUB1B", "9606.P_MAD2"),
                                      ("9606.P_MAD2", "9606.P_TRIP13")})

    def test_neighbours_sum_over_seeds(self) -> None:
        """MAD2 links to two seeds, so it ranks above CDC20 with one seed."""
        links = layer.read_functional_links(self.path)
        names = {"9606.P_MAD2": "MAD2L1", "9606.P_CDC20": "CDC20"}
        rows = layer.derive_functional_neighbours(
            links, {"BUB1B": "9606.P_BUB1B", "TRIP13": "9606.P_TRIP13"}, names)
        self.assertEqual([row["symbol"] for row in rows], ["MAD2L1", "CDC20"])
        self.assertEqual(rows[0]["seeds"], ["BUB1B", "TRIP13"])

    def test_edges_stay_inside_the_universe(self) -> None:
        """An edge to a protein outside the universe is not written."""
        links = layer.read_functional_links(self.path)
        rows = layer.functional_edges(links, {"9606.P_BUB1B": "BUB1B", "9606.P_CDC20": "CDC20"})
        self.assertEqual([(row["source_id"], row["target_id"]) for row in rows],
                         [("GENE:BUB1B", "GENE:CDC20")])
        self.assertEqual(rows[0]["relation"], "functional_interaction")


class GeneLayerTest(unittest.TestCase):
    """The ranker loads one of three gene layers, and the default layer is unchanged."""

    def setUp(self) -> None:
        """Write the fixture tables."""
        self.folder = tempfile.TemporaryDirectory()
        self.nodes, self.edges, self.functional = write_tables(Path(self.folder.name))

    def tearDown(self) -> None:
        """Remove the temporary directory."""
        self.folder.cleanup()

    def edge_set(self, graph: ranker.Graph) -> set[tuple[str, str, str]]:
        """Return the endpoints and the relation of every edge.

        Args:
            graph: The graph.

        Returns:
            The edge triples.
        """
        return {(edge.source_id, edge.target_id, edge.relation) for edge in graph.edges}

    def test_default_layer_equals_the_plain_loader(self) -> None:
        """The default layer reads no functional file and builds the same graph."""
        plain = ranker.load_graph(self.nodes, self.edges)
        layered = ranker.load_layered_graph(self.nodes, self.edges,
                                            functional_dir=Path(self.folder.name) / "absent")
        self.assertEqual(self.edge_set(plain), self.edge_set(layered))
        self.assertEqual(plain.stats, layered.stats)

    def test_functional_layer_replaces_literature_gene_edges(self) -> None:
        """Both co-mention gene edges leave, and the disease and drug edges stay."""
        graph = ranker.load_layered_graph(self.nodes, self.edges, "functional", self.functional)
        self.assertEqual(self.edge_set(graph), {
            ("DISEASE:D1", "GENE:BUB1B", "associate"),
            ("GENE:AKT1", "CHEMBL:C1", "inhibitor"),
            ("GENE:BUB1B", "GENE:CDC20", "functional_interaction"),
            ("GENE:BUB1B", "GENE:MAD1L1", "functional_interaction")})
        self.assertEqual(graph.kind("GENE:MAD1L1"), "gene")

    def test_union_layer_keeps_both(self) -> None:
        """The union holds the four default edges and the two functional edges."""
        graph = ranker.load_layered_graph(self.nodes, self.edges, "union", self.functional)
        self.assertEqual(len(graph.edges), 6)

    def test_functional_relation_has_its_own_class(self) -> None:
        """A functional edge is neither a co-mention nor a pharmacological claim."""
        self.assertEqual(ranker.relation_class("functional_interaction"), ranker.FUNCTIONAL)

    def test_unknown_layer_is_refused(self) -> None:
        """An unknown layer name raises an error."""
        with self.assertRaises(ValueError):
            ranker.apply_gene_layer({}, [], "citation", {}, [])

    def test_drug_route_changes_with_the_layer(self) -> None:
        """The drug on AKT1 is reachable through co-mention and not through the functional layer."""
        options = ranker.RankOptions()
        default = ranker.load_layered_graph(self.nodes, self.edges)
        functional = ranker.load_layered_graph(self.nodes, self.edges, "functional",
                                               self.functional)
        evidence = ranker.DrugEvidence(stage={"CHEMBL:C1": "PHASE_3"})
        self.assertEqual([c.node_id for c in ranker.rank_candidates(default, ["DISEASE:D1"],
                                                                    options, evidence)],
                         ["CHEMBL:C1"])
        self.assertEqual(ranker.rank_candidates(functional, ["DISEASE:D1"], options, evidence), [])


class CompareHelpersTest(unittest.TestCase):
    """The comparison helpers compute a rank correlation and a recovery at another top k."""

    def test_spearman_of_identical_and_reversed_orders(self) -> None:
        """Identical ranks give 1 and reversed ranks give -1."""
        self.assertEqual(compare.spearman([1, 2, 3, 4], [1, 2, 3, 4]), 1.0)
        self.assertEqual(compare.spearman([1, 2, 3, 4], [4, 3, 2, 1]), -1.0)
        self.assertIsNone(compare.spearman([1, 2], [2, 1]))

    def test_average_ranks_share_a_tie(self) -> None:
        """Two tied values share the average of their two ranks."""
        self.assertEqual(compare.average_ranks([10, 20, 20, 30]), [1.0, 2.5, 2.5, 4.0])

    def test_rows_at_k(self) -> None:
        """A best rank of 3 is a recovery at top 3 and not at top 1, and rank 0 is never one."""
        rows = [{"best_rank": 3, "recovered": True}, {"best_rank": 0, "recovered": False}]
        self.assertEqual([row["recovered"] for row in compare.rows_at_k(rows, 1)], [False, False])
        self.assertEqual([row["recovered"] for row in compare.rows_at_k(rows, 3)], [True, False])


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
