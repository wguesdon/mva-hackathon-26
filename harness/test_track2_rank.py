"""Verify the Track 2 path ranker and its two controls on synthetic class A input.

Every graph of this test is written here and its answer is known by hand. No patient data and no
real identifier enters the test. The hand graph carries one relation of every class, so the test
separates a co-mention from a pharmacological action and a restoring action from a reinforcing one.

    seed_disease -- strong_gene    associate      evidence 20, PMIDs P1 P2 P3
    seed_disease -- weak_gene      associate      evidence 2,  PMID  P4
    seed_disease -- partner_gene   associate      evidence 6,  PMID  P13
    seed_disease -- known_drug     indication     evidence 5,  PMID  P10
    strong_gene  -- good_drug      agonist        evidence 10, PMIDs P5 P6
    strong_gene  -- blocker_drug   inhibitor      evidence 10, PMIDs P5 P6
    strong_gene  -- reagent        associate      evidence 50, 30 PMIDs
    strong_gene  -- shared_pathway participates_in
    partner_gene -- shared_pathway participates_in
    partner_gene -- pathway_blocker inhibitor     evidence 12, PMIDs P14 P15
    weak_gene    -- far_gene       interact       evidence 4,  PMID  P7
    weak_gene    -- known_drug     agonist        evidence 8,  PMIDs P11 P12
    far_gene     -- far_drug       inhibitor      evidence 9,  PMIDs P8 P9
    other_disease -- strong_gene   associate      evidence 30, PMIDs P16 to P20
    other_disease -- other_drug    indication     evidence 100

The defect of the seed is a loss of function of strong_gene. The good drug therefore restores and
the blocker drug reinforces, and both carry the same evidence score.
"""

from __future__ import annotations

import json
import math
import random
import statistics
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import track2_controls as controls  # noqa: E402
import track2_null_correction as correction  # noqa: E402
import track2_make_graph_fixture as fixture  # noqa: E402
import track2_rank as ranker  # noqa: E402

SEED = "DISEASE:SYNSEED"
OTHER_DISEASE = "DISEASE:SYNOTHERDIS"
STRONG = "GENE:SYNSTRONG"
WEAK = "GENE:SYNWEAK"
FAR_GENE = "GENE:SYNFARGENE"
PARTNER = "GENE:SYNPART"
GOOD = "CHEMBL:SYNGOOD"
BLOCKER = "CHEMBL:SYNBLOCK"
PATHWAY_BLOCKER = "CHEMBL:SYNPBLOCK"
REAGENT = "CHEMICAL:SYNREAGENT"
FAR = "CHEMBL:SYNFAR"
KNOWN = "CHEMBL:SYNKNOWN"
OTHER_DRUG = "CHEMBL:SYNOTHER"

NODES = """node_id\tkind\tlabel\tsource
DISEASE:SYNSEED\tdisease\tseed_disease\ttest
DISEASE:SYNOTHERDIS\tdisease\tother_disease\ttest
GENE:SYNSTRONG\tgene\tstrong_gene\ttest
GENE:SYNWEAK\tgene\tweak_gene\ttest
GENE:SYNFARGENE\tgene\tfar_gene\ttest
GENE:SYNPART\tgene\tpartner_gene\ttest
PATHWAY:SYNPATH\tpathway\tshared_pathway\ttest
CHEMBL:SYNGOOD\tchemical\tgood_drug\ttest
CHEMBL:SYNBLOCK\tchemical\tblocker_drug\ttest
CHEMBL:SYNPBLOCK\tchemical\tpathway_blocker\ttest
CHEMICAL:SYNREAGENT\tchemical\treagent\ttest
CHEMBL:SYNFAR\tchemical\tfar_drug\ttest
CHEMBL:SYNKNOWN\tchemical\tknown_drug\ttest
CHEMBL:SYNOTHER\tchemical\tother_drug\ttest
"""

REAGENT_PMIDS = ",".join(f"R{number:02d}" for number in range(30))
OTHER_PMIDS = "P16,P17,P18,P19,P20"

EDGES = f"""source_id\ttarget_id\trelation\tevidence_count\tpmids\tsource
DISEASE:SYNSEED\tGENE:SYNSTRONG\tassociate\t20\tP1,P2,P3\tpubtator3
DISEASE:SYNSEED\tGENE:SYNWEAK\tassociate\t2\tP4\tpubtator3
DISEASE:SYNSEED\tGENE:SYNPART\tassociate\t6\tP13\tpubtator3
DISEASE:SYNSEED\tCHEMBL:SYNKNOWN\tindication\t5\tP10\tchembl
GENE:SYNSTRONG\tCHEMBL:SYNGOOD\tagonist\t10\tP5,P6\topen_targets
GENE:SYNSTRONG\tCHEMBL:SYNBLOCK\tinhibitor\t10\tP5,P6\topen_targets
GENE:SYNSTRONG\tCHEMICAL:SYNREAGENT\tassociate\t50\t{REAGENT_PMIDS}\tpubtator3
GENE:SYNSTRONG\tPATHWAY:SYNPATH\tparticipates_in\t1\t\topen_targets_reactome
GENE:SYNPART\tPATHWAY:SYNPATH\tparticipates_in\t1\t\topen_targets_reactome
GENE:SYNPART\tCHEMBL:SYNPBLOCK\tinhibitor\t12\tP14,P15\topen_targets
GENE:SYNWEAK\tGENE:SYNFARGENE\tinteract\t4\tP7\tpubtator3
GENE:SYNWEAK\tCHEMBL:SYNKNOWN\tagonist\t8\tP11,P12\topen_targets
GENE:SYNFARGENE\tCHEMBL:SYNFAR\tinhibitor\t9\tP8,P9\topen_targets
DISEASE:SYNOTHERDIS\tGENE:SYNSTRONG\tassociate\t30\t{OTHER_PMIDS}\tpubtator3
DISEASE:SYNOTHERDIS\tCHEMBL:SYNOTHER\tindication\t100\t\tchembl
"""

DRUGS = """drug_id\tdrug_name\tdrug_type\ttarget_gene_id\taction_type\tmax_clinical_stage_on_target\tmax_clinical_stage_overall\tmechanism_of_action
CHEMBL:SYNGOOD\tgood_drug\tSmall molecule\tGENE:SYNSTRONG\tagonist\tPHASE_2\tPHASE_2\tstrong_gene agonist
CHEMBL:SYNBLOCK\tblocker_drug\tSmall molecule\tGENE:SYNSTRONG\tinhibitor\tPHASE_3\tPHASE_3\tstrong_gene inhibitor
CHEMBL:SYNPBLOCK\tpathway_blocker\tSmall molecule\tGENE:SYNPART\tinhibitor\tPHASE_1\tPHASE_1\tpartner_gene inhibitor
CHEMBL:SYNFAR\tfar_drug\tSmall molecule\tGENE:SYNFARGENE\tinhibitor\tAPPROVAL\tAPPROVAL\tfar_gene inhibitor
CHEMBL:SYNKNOWN\tknown_drug\tSmall molecule\tGENE:SYNWEAK\tagonist\tAPPROVAL\tAPPROVAL\tweak_gene agonist
CHEMBL:SYNOTHER\tother_drug\tSmall molecule\tGENE:SYNWEAK\tinhibitor\tAPPROVAL\tAPPROVAL\tweak_gene inhibitor
"""

LOOSE = ranker.RankOptions(terminal_classes=ranker.RELATION_CLASSES, require_drug_evidence=False,
                           blocked_intermediate_kinds=())
DEFECT = ranker.RankOptions(defect="loss_of_function", defect_genes=(STRONG,))


def write_graph(folder: Path, nodes: str = NODES, edges: str = EDGES,
                drugs: str = DRUGS) -> tuple[Path, Path, Path]:
    """Write one graph to a temporary directory.

    Args:
        folder: The temporary directory.
        nodes: The content of nodes.tsv.
        edges: The content of edges.tsv.
        drugs: The content of drug_evidence.tsv.

    Returns:
        The nodes path, the edges path and the drug evidence path.
    """
    nodes_path = folder / "nodes.tsv"
    edges_path = folder / "edges.tsv"
    drugs_path = folder / "drug_evidence.tsv"
    nodes_path.write_text(nodes, encoding="utf-8")
    edges_path.write_text(edges, encoding="utf-8")
    drugs_path.write_text(drugs, encoding="utf-8")
    return nodes_path, edges_path, drugs_path


class RelationClassTest(unittest.TestCase):
    """Check the declared class of every relation the graph uses."""

    def test_a_pubtator_relation_is_a_co_mention(self) -> None:
        """PubTator3 states that two entities share a publication and nothing more."""
        for name in ("associate", "positive_correlate", "negative_correlate", "stimulate",
                     "interact"):
            self.assertEqual(ranker.relation_class(name), ranker.CO_MENTION)

    def test_inhibit_and_inhibitor_differ(self) -> None:
        """The two names differ by one letter and by class, and the class decides the ranking."""
        self.assertEqual(ranker.relation_class("inhibit"), ranker.CO_MENTION)
        self.assertEqual(ranker.relation_class("inhibitor"), ranker.PHARMACOLOGICAL)

    def test_an_unknown_relation_takes_the_conservative_class(self) -> None:
        """A relation that the table does not name must not nominate a drug."""
        self.assertEqual(ranker.relation_class("brand_new_relation"), ranker.CO_MENTION)
        self.assertEqual(ranker.relation_class("brand_new_relation", ("brand_new_relation",)),
                         ranker.PHARMACOLOGICAL)

    def test_an_unknown_relation_is_reported(self) -> None:
        """A new action type of a source database must not enter the co-mention class in silence."""
        folder = tempfile.TemporaryDirectory()
        edges = EDGES.replace("\tagonist\t", "\tsuper_new_agonist\t")
        nodes_path, edges_path, _ = write_graph(Path(folder.name), NODES, edges)
        graph = ranker.load_graph(nodes_path, edges_path)
        self.assertIn("super_new_agonist", ranker.unclassified_relations(graph))
        self.assertNotIn("inhibitor", ranker.unclassified_relations(graph))
        folder.cleanup()

    def test_every_declared_action_direction_is_known(self) -> None:
        """An action type with a direction must map to increase or decrease and nothing else."""
        self.assertEqual(ranker.action_direction_of("inhibitor"), ranker.DECREASE)
        self.assertEqual(ranker.action_direction_of("agonist"), ranker.INCREASE)
        self.assertEqual(ranker.action_direction_of("modulator"), ranker.UNCLEAR)
        self.assertEqual(ranker.action_direction_of(""), ranker.UNCLEAR)


class RankerTest(unittest.TestCase):
    """Check the paths, the score, the constraints and the output of the ranker."""

    def test_the_reagent_wins_without_the_constraints(self) -> None:
        """A co-mention with 30 PMIDs beats every drug, which is the failure the constraints fix."""
        order = [item.node_id for item in ranker.rank_candidates(self.graph, [SEED], LOOSE)]
        self.assertEqual(order[0], REAGENT)
        self.assertIn(OTHER_DRUG, order)

    def test_the_terminal_constraint_alone_removes_the_reagent(self) -> None:
        """A co-mention may support a step inside a path but may not nominate a drug."""
        options = ranker.RankOptions(require_drug_evidence=False, blocked_intermediate_kinds=())
        reached = {item.node_id for item in ranker.rank_candidates(self.graph, [SEED], options)}
        self.assertNotIn(REAGENT, reached)
        self.assertIn(GOOD, reached)

    def test_the_eligibility_constraint_alone_removes_the_reagent(self) -> None:
        """A compound that was never a drug is not a repurposing candidate."""
        options = ranker.RankOptions(terminal_classes=ranker.RELATION_CLASSES,
                                     blocked_intermediate_kinds=())
        reached = {item.node_id for item in
                   ranker.rank_candidates(self.graph, [SEED], options, self.evidence)}
        self.assertNotIn(REAGENT, reached)
        self.assertIn(GOOD, reached)

    def test_a_path_may_not_hop_through_a_second_disease(self) -> None:
        """The drugs of a co-mentioned disease are not evidence about the seed disease."""
        default = {item.node_id for item in
                   ranker.rank_candidates(self.graph, [SEED], evidence=self.evidence)}
        self.assertNotIn(OTHER_DRUG, default)
        allowed = ranker.RankOptions(blocked_intermediate_kinds=())
        opened = {item.node_id for item in
                  ranker.rank_candidates(self.graph, [SEED], allowed, self.evidence)}
        self.assertIn(OTHER_DRUG, opened)

    def test_the_candidate_set_and_the_path_length(self) -> None:
        """Five drugs remain, at the distance the hand graph gives them."""
        candidates = self.by_id()
        self.assertEqual(set(candidates), {GOOD, BLOCKER, PATHWAY_BLOCKER, FAR, KNOWN})
        self.assertEqual(candidates[KNOWN].path_length, 1)
        self.assertEqual(candidates[GOOD].path_length, 2)
        self.assertEqual(candidates[FAR].path_length, 3)

    def test_the_hand_computed_scores_match(self) -> None:
        """Each score is the formula of the module docstring on the hand graph."""
        expected = {KNOWN: math.log1p(5) + math.log1p(1),
                    GOOD: 0.5 * (math.log1p(10) + math.log1p(5)),
                    BLOCKER: 0.5 * (math.log1p(10) + math.log1p(5)),
                    PATHWAY_BLOCKER: 0.5 * (math.log1p(6) + math.log1p(3)),
                    FAR: 0.25 * (math.log1p(2) + math.log1p(4))}
        for node_id, candidate in self.by_id().items():
            self.assertAlmostEqual(candidate.score, expected[node_id], places=10)

    def test_the_best_path_carries_the_pmids_of_every_step(self) -> None:
        """The path of the good drug joins both edges and both PMID sets."""
        candidate = self.by_id()[GOOD]
        self.assertEqual(candidate.best_path.node_ids, (SEED, STRONG, GOOD))
        self.assertEqual(candidate.pmids, ("P1", "P2", "P3", "P5", "P6"))
        self.assertEqual(candidate.min_evidence, 10)
        self.assertEqual(ranker.path_string(self.graph, candidate.best_path),
                         "seed_disease -> strong_gene -> good_drug")
        self.assertEqual(ranker.relation_string(self.graph, candidate.best_path),
                         "associate -> agonist")
        self.assertEqual(candidate.terminal_class, ranker.PHARMACOLOGICAL)
        self.assertEqual(candidate.clinical_stage, "PHASE_2")
        self.assertEqual(candidate.mechanism, "strong_gene agonist")
        self.assertEqual(candidate.eligibility, "drug_record")

    def test_a_longer_pharmacological_route_beats_a_closer_co_mention(self) -> None:
        """The far drug keeps its three edge route although a co-mention sits closer."""
        candidate = self.by_id()[FAR]
        self.assertEqual(candidate.best_path.node_ids, (SEED, WEAK, FAR_GENE, FAR))
        self.assertEqual(candidate.terminal_relation, "inhibitor")

    def test_a_seed_set_records_the_seed_of_each_path(self) -> None:
        """Ranking from a consequence is a different question, so each row names its own seed."""
        candidates = {item.node_id: item for item in
                      ranker.rank_candidates(self.graph, [SEED, FAR_GENE], evidence=self.evidence)}
        self.assertEqual(candidates[FAR].seed_id, FAR_GENE)
        self.assertEqual(candidates[FAR].path_length, 1)
        self.assertEqual(candidates[GOOD].seed_id, SEED)

    def test_the_output_table_holds_the_path_the_reason_and_the_pmids(self) -> None:
        """The written file carries the contract columns and the explanation of each row."""
        folder = Path(self.folder.name)
        candidates = ranker.rank_candidates(self.graph, [SEED], DEFECT, self.evidence)
        path = ranker.write_ranking(self.graph, candidates, folder / "ranking",
                                    stamp="2026_01_01_000000")
        self.assertEqual(path.name, "candidates_2026_01_01_000000.tsv")
        lines = path.read_text(encoding="utf-8").strip().splitlines()
        self.assertEqual(lines[0].split("\t"), list(ranker.OUTPUT_COLUMNS))
        self.assertEqual(len(lines), 6)
        row = dict(zip(ranker.OUTPUT_COLUMNS, lines[1].split("\t")))
        self.assertEqual(row["chembl_id"], GOOD)
        self.assertEqual(row["path"], "seed_disease -> strong_gene -> good_drug")
        self.assertEqual(row["pmids"], "P1,P2,P3,P5,P6")
        self.assertEqual(row["relations"], "associate -> agonist")
        self.assertEqual(row["direction_verdict"], "restores_target")
        self.assertEqual(row["seed"], SEED)

    def test_a_missing_column_stops_the_load(self) -> None:
        """A table that breaks the input contract raises instead of ranking half a graph."""
        broken = Path(self.folder.name) / "broken.tsv"
        broken.write_text("source_id\ttarget_id\trelation\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            ranker.read_edges(broken)

    def test_a_chemical_prefix_of_either_builder_is_a_chemical(self) -> None:
        """An edge endpoint that the node table misses keeps its kind from the identifier prefix."""
        self.assertEqual(ranker.kind_from_id("CHEMICAL:MESH:D015739"), "chemical")
        self.assertEqual(ranker.kind_from_id("CHEMBL:CHEMBL123"), "chemical")
        self.assertEqual(ranker.kind_from_id("GENE:BUB1B"), "gene")
        self.assertEqual(ranker.kind_from_id("OTHER:1"), "unknown")

    def test_an_unknown_seed_stops_the_run(self) -> None:
        """A seed that is not a node of the graph raises."""
        with self.assertRaises(ValueError):
            ranker.resolve_seeds(self.graph, ["GENE:ABSENT"])

    def test_an_unknown_terminal_class_stops_the_run(self) -> None:
        """A typed class name that does not exist raises instead of ranking with no constraint."""
        with self.assertRaises(ValueError):
            ranker.parse_terminal_classes("pharmacologial")
        self.assertEqual(ranker.parse_terminal_classes("all"), ranker.RELATION_CLASSES)

    def by_id(self) -> dict[str, ranker.Candidate]:
        """Return the constrained candidates of the hand graph by identifier.

        Returns:
            Every candidate, keyed by node identifier.
        """
        return {item.node_id: item for item in
                ranker.rank_candidates(self.graph, [SEED], evidence=self.evidence)}

    def setUp(self) -> None:
        """Write the hand graph and load it."""
        self.folder = tempfile.TemporaryDirectory()
        nodes_path, edges_path, drugs_path = write_graph(Path(self.folder.name))
        self.graph = ranker.load_graph(nodes_path, edges_path)
        self.evidence = ranker.read_drug_evidence(drugs_path)

    def tearDown(self) -> None:
        """Remove the temporary directory."""
        self.folder.cleanup()


class DirectionTest(unittest.TestCase):
    """Check the direction of effect against a declared loss of function defect."""

    def test_an_agonist_of_the_defect_gene_restores(self) -> None:
        """The needed direction of a loss of function is an increase of function."""
        candidate = self.by_id()[GOOD]
        self.assertEqual(candidate.direction_verdict, "restores_target")
        self.assertEqual(candidate.action_direction, ranker.INCREASE)
        self.assertIn("strong_gene", candidate.direction_reason)

    def test_an_inhibitor_of_the_defect_gene_reinforces(self) -> None:
        """A drug that pushes the same way as the defect is a harm and not a candidate."""
        candidate = self.by_id()[BLOCKER]
        self.assertEqual(candidate.direction_verdict, "reinforces_target")
        self.assertEqual(candidate.action_direction, ranker.DECREASE)

    def test_an_inhibitor_of_a_pathway_partner_reinforces_the_pathway(self) -> None:
        """The shared pathway comes from the graph, so no compound and no gene is written here."""
        candidate = self.by_id()[PATHWAY_BLOCKER]
        self.assertEqual(candidate.direction_verdict, "reinforces_pathway")
        self.assertIn("shared_pathway", candidate.direction_reason)

    def test_a_drug_on_an_unrelated_gene_does_not_compose(self) -> None:
        """A co-mention step carries no sign, so the direction does not travel along the path."""
        candidate = self.by_id()[FAR]
        self.assertEqual(candidate.direction_verdict, "unclear_composition")

    def test_an_indication_edge_names_no_target_function(self) -> None:
        """A direct clinical edge ends on the disease, so it states nothing about a gene."""
        candidate = self.by_id()[KNOWN]
        self.assertEqual(candidate.direction_verdict, "unclear_target")

    def test_the_reinforcing_drug_sinks_below_every_other_candidate(self) -> None:
        """The blocker keeps the score of the good drug and still ranks last."""
        order = [item.node_id for item in
                 ranker.rank_candidates(self.graph, [SEED], DEFECT, self.evidence)]
        self.assertEqual(order[0], GOOD)
        self.assertEqual(order[-1], BLOCKER)
        self.assertEqual(order[-2], PATHWAY_BLOCKER)
        candidates = self.by_id()
        self.assertAlmostEqual(candidates[BLOCKER].score, candidates[GOOD].score, places=10)

    def test_a_top_cut_keeps_the_rejected_class(self) -> None:
        """A reviewer must see what the ranker rejected, so a cut never hides it."""
        candidates = ranker.rank_candidates(self.graph, [SEED], DEFECT, self.evidence)
        kept, ranks = ranker.select_rows(candidates, 2)
        self.assertEqual([item.node_id for item in kept][:2], [GOOD, KNOWN])
        self.assertIn(BLOCKER, [item.node_id for item in kept])
        self.assertIn(PATHWAY_BLOCKER, [item.node_id for item in kept])
        self.assertEqual(ranks[-1], 5)
        self.assertEqual(len(kept), len(ranks))

    def test_a_top_cut_without_a_defect_keeps_the_leading_rows_only(self) -> None:
        """A run with no verdict has no rejected class, so the cut is the plain top list."""
        candidates = ranker.rank_candidates(self.graph, [SEED], evidence=self.evidence)
        kept, ranks = ranker.select_rows(candidates, 2)
        self.assertEqual(len(kept), 2)
        self.assertEqual(ranks, [1, 2])

    def test_a_run_without_a_declared_defect_does_not_reorder(self) -> None:
        """The default run gives no verdict, so the order is the evidence score alone."""
        order = [item.node_id for item in
                 ranker.rank_candidates(self.graph, [SEED], evidence=self.evidence)]
        self.assertEqual(order, [KNOWN, BLOCKER, GOOD, PATHWAY_BLOCKER, FAR])
        for candidate in ranker.rank_candidates(self.graph, [SEED], evidence=self.evidence):
            self.assertEqual(candidate.direction_verdict, "not_assessed")

    def test_a_gene_seed_is_a_defect_gene_on_its_own(self) -> None:
        """A run seeded on the gene needs no separate declaration of the defect gene."""
        options = ranker.RankOptions(defect="loss_of_function")
        candidates = {item.node_id: item for item in
                      ranker.rank_candidates(self.graph, [STRONG], options, self.evidence)}
        self.assertEqual(candidates[GOOD].direction_verdict, "restores_target")
        self.assertEqual(candidates[BLOCKER].direction_verdict, "reinforces_target")

    def test_a_gain_of_function_defect_inverts_the_verdicts(self) -> None:
        """The same table serves the opposite defect without a second rule."""
        options = ranker.RankOptions(defect="gain_of_function", defect_genes=(STRONG,))
        candidates = {item.node_id: item for item in
                      ranker.rank_candidates(self.graph, [SEED], options, self.evidence)}
        self.assertEqual(candidates[GOOD].direction_verdict, "reinforces_target")
        self.assertEqual(candidates[BLOCKER].direction_verdict, "restores_target")

    def by_id(self) -> dict[str, ranker.Candidate]:
        """Return the candidates of the defect run by identifier.

        Returns:
            Every candidate, keyed by node identifier.
        """
        return {item.node_id: item for item in
                ranker.rank_candidates(self.graph, [SEED], DEFECT, self.evidence)}

    def setUp(self) -> None:
        """Write the hand graph and load it."""
        self.folder = tempfile.TemporaryDirectory()
        nodes_path, edges_path, drugs_path = write_graph(Path(self.folder.name))
        self.graph = ranker.load_graph(nodes_path, edges_path)
        self.evidence = ranker.read_drug_evidence(drugs_path)

    def tearDown(self) -> None:
        """Remove the temporary directory."""
        self.folder.cleanup()


class PermutationTest(unittest.TestCase):
    """Check that the permutation keeps the degree sequence and answers the degree question."""

    def test_the_degree_of_every_node_is_preserved(self) -> None:
        """Ten permutations of the fixture graph leave every degree unchanged."""
        before = controls.degree_sequence(self.graph)
        generator = random.Random(7)
        for _ in range(10):
            permuted = controls.permute_graph(self.graph, generator)
            self.assertEqual(controls.degree_sequence(permuted), before)

    def test_the_edge_attributes_travel_with_the_edge(self) -> None:
        """The permutation moves endpoints only. The relations and the PMIDs do not change."""
        permuted = controls.permute_graph(self.graph, random.Random(11))
        self.assertEqual(len(permuted.edges), len(self.graph.edges))
        self.assertEqual(sorted(edge.relation for edge in permuted.edges),
                         sorted(edge.relation for edge in self.graph.edges))
        self.assertEqual(sorted(edge.pmids for edge in permuted.edges),
                         sorted(edge.pmids for edge in self.graph.edges))

    def test_the_permutation_changes_the_graph(self) -> None:
        """A permutation that changed nothing would make the control meaningless."""
        permuted = controls.permute_graph(self.graph, random.Random(3))
        before = {(edge.source_id, edge.target_id, edge.relation) for edge in self.graph.edges}
        after = {(edge.source_id, edge.target_id, edge.relation) for edge in permuted.edges}
        self.assertNotEqual(before, after)

    def test_a_band_holds_one_kind_of_node_only(self) -> None:
        """A gene never takes the place of a chemical, so the permuted graph stays readable."""
        permuted = controls.permute_graph(self.graph, random.Random(5))
        for index, edge in enumerate(permuted.edges):
            self.assertEqual(permuted.kind(edge.source_id),
                             self.graph.kind(self.graph.edges[index].source_id))
            self.assertEqual(permuted.kind(edge.target_id),
                             self.graph.kind(self.graph.edges[index].target_id))

    def test_the_default_scope_rewires_every_edge_of_the_graph(self) -> None:
        """The recorded null is the whole graph and not the indication layer."""
        permuted = controls.permute_graph(self.graph, random.Random(13))
        self.assertEqual(permuted.stats["edges_in_permutation_scope"], len(self.graph.edges))

    def test_the_indication_scope_leaves_every_other_layer_in_place(self) -> None:
        """An edge outside the indication layer keeps both of its endpoints."""
        scope = controls.scope_relations(controls.NULL_SCOPE_INDICATION_EDGES,
                                         controls.TREATMENT_RELATIONS)
        permuted = controls.permute_graph(self.graph, random.Random(13), scope)
        in_scope = [index for index, edge in enumerate(self.graph.edges)
                    if edge.relation.strip().lower() in controls.TREATMENT_RELATIONS]
        self.assertEqual(permuted.stats["edges_in_permutation_scope"], len(in_scope))
        for index, edge in enumerate(self.graph.edges):
            if index in set(in_scope):
                continue
            self.assertEqual(permuted.edges[index].source_id, edge.source_id)
            self.assertEqual(permuted.edges[index].target_id, edge.target_id)

    def test_every_scope_preserves_the_degree_sequence(self) -> None:
        """A node keeps its slot count inside every group, so its degree cannot change."""
        before = controls.degree_sequence(self.graph)
        for scope in controls.NULL_SCOPES:
            names = controls.scope_relations(scope, controls.TREATMENT_RELATIONS)
            permuted = controls.permute_graph(self.graph, random.Random(17), names)
            self.assertEqual(controls.degree_sequence(permuted), before)

    def test_an_unknown_scope_stops_the_run(self) -> None:
        """A scope that the control does not implement must not fall back on a default."""
        with self.assertRaises(ValueError):
            controls.scope_relations("every_second_edge", controls.TREATMENT_RELATIONS)

    def test_the_control_writes_a_p_value_for_every_tested_candidate(self) -> None:
        """Every candidate of the real top k carries a p-value and a verdict."""
        result = controls.permutation_control(self.graph, [fixture.SEED_DISEASE], permutations=25,
                                              top_k=5, random_seed=99, evidence=self.evidence)
        self.assertTrue(result["degree_sequence_preserved"])
        self.assertEqual(len(result["candidates"]), 5)
        for row in result["candidates"]:
            self.assertGreaterEqual(row["p_value"], 1 / 26)
            self.assertLessEqual(row["p_value"], 1.0)
            self.assertIn(row["verdict"], ("evidence_supported", "degree_artifact"))

    def test_the_control_ranks_the_permuted_graph_with_the_report_constraints(self) -> None:
        """A constraint of the report is a constraint of the control, or the two disagree."""
        options = ranker.RankOptions(defect="loss_of_function",
                                     defect_genes=(fixture.DEFECT_GENE,))
        result = controls.permutation_control(self.graph, [fixture.SEED_DISEASE], permutations=5,
                                              top_k=3, random_seed=99, options=options,
                                              evidence=self.evidence)
        self.assertTrue(all(row["terminal_relation"] for row in result["candidates"]))
        for row in result["candidates"]:
            self.assertNotEqual(row["direction_verdict"], "reinforces_target")

    def test_the_control_reports_the_self_loops_it_could_not_repair(self) -> None:
        """A self loop that survives the repair is dropped, which would change a degree."""
        result = controls.permutation_control(self.graph, [fixture.SEED_DISEASE], permutations=10,
                                              top_k=3, random_seed=5, evidence=self.evidence)
        self.assertEqual(result["largest_self_loop_count_after_repair"], 0)
        self.assertTrue(result["degree_sequence_preserved"])

    def test_too_few_permutations_cannot_pass(self) -> None:
        """A p-value below alpha is impossible with five permutations, so the run says so."""
        result = controls.permutation_control(self.graph, [fixture.SEED_DISEASE], permutations=5,
                                              top_k=3, random_seed=99, evidence=self.evidence)
        self.assertEqual(result["verdict"], "INSUFFICIENT_PERMUTATIONS")

    def setUp(self) -> None:
        """Write the synthetic fixture graph and load it."""
        self.folder = tempfile.TemporaryDirectory()
        nodes_path, edges_path, drugs_path = fixture.write_fixture(Path(self.folder.name))
        self.graph = ranker.load_graph(nodes_path, edges_path)
        self.evidence = ranker.read_drug_evidence(drugs_path)

    def tearDown(self) -> None:
        """Remove the temporary directory."""
        self.folder.cleanup()


class HeldOutTest(unittest.TestCase):
    """Check the case selection, the edge removal, the recovery rate and the regime."""

    def test_only_a_disease_with_an_approved_drug_is_a_case(self) -> None:
        """The hand graph holds two diseases with an indication edge."""
        cases = controls.treatment_cases(self.graph, controls.TREATMENT_RELATIONS)
        self.assertEqual(set(cases), {SEED, OTHER_DISEASE})
        self.assertEqual(cases[SEED]["drugs"], {KNOWN})
        self.assertEqual(len(cases[SEED]["edge_indices"]), 1)

    def test_the_indication_relation_comes_from_the_declared_table(self) -> None:
        """One place declares what a clinical edge is, and the control reads that place."""
        self.assertIn("indication", controls.TREATMENT_RELATIONS)
        self.assertNotIn("associate", controls.TREATMENT_RELATIONS)

    def test_the_removal_drops_the_right_edge_and_keeps_the_rest(self) -> None:
        """The reduced graph loses the indication edge only."""
        cases = controls.treatment_cases(self.graph, controls.TREATMENT_RELATIONS)
        reduced = controls.graph_without(self.graph, cases[SEED]["edge_indices"])
        self.assertEqual(len(reduced.edges), len(self.graph.edges) - 1)
        pairs = {(edge.source_id, edge.target_id) for edge in reduced.edges}
        self.assertNotIn((SEED, KNOWN), pairs)
        self.assertIn((WEAK, KNOWN), pairs)

    def test_the_real_graph_is_untouched_after_the_run(self) -> None:
        """Every removed edge is in place again for the next case."""
        before = [(edge.source_id, edge.target_id, edge.relation) for edge in self.graph.edges]
        degrees = controls.degree_sequence(self.graph)
        controls.heldout_control(self.graph, 3, controls.TREATMENT_RELATIONS, min_cases=1,
                                 evidence=self.evidence)
        after = [(edge.source_id, edge.target_id, edge.relation) for edge in self.graph.edges]
        self.assertEqual(before, after)
        self.assertEqual(controls.degree_sequence(self.graph), degrees)

    def test_the_true_drug_returns_by_its_gene_route(self) -> None:
        """The known drug is an agonist of the weak gene, so a two edge path finds it again.

        The route through the weak gene carries less evidence than the route to the good drug, so
        the true drug comes back at rank 4 of 5. A top list of three misses it and a top list of
        five holds it. The control therefore measures the ranking and not the reachability alone.
        """
        result = controls.heldout_control(self.graph, 5, controls.TREATMENT_RELATIONS, min_cases=1,
                                          evidence=self.evidence)
        row = next(item for item in result["results"] if item["disease"] == SEED)
        self.assertEqual(row["removed_edges"], 1)
        self.assertEqual(row["best_rank"], 4)
        self.assertTrue(row["recovered"])
        self.assertTrue(result["graph_restored"])
        strict = controls.heldout_control(self.graph, 3, controls.TREATMENT_RELATIONS, min_cases=1,
                                          evidence=self.evidence)
        strict_row = next(item for item in strict["results"] if item["disease"] == SEED)
        self.assertFalse(strict_row["recovered"])

    def test_a_drug_with_no_other_route_is_not_recovered(self) -> None:
        """When the indication edge is the only route, the true drug cannot come back."""
        folder = Path(self.folder.name)
        edges = "\n".join(line for line in EDGES.strip().splitlines()
                          if not line.startswith("GENE:SYNWEAK\tCHEMBL:SYNKNOWN")) + "\n"
        nodes_path, edges_path, drugs_path = write_graph(folder, NODES, edges)
        graph = ranker.load_graph(nodes_path, edges_path)
        result = controls.heldout_control(graph, 5, controls.TREATMENT_RELATIONS, min_cases=1,
                                          evidence=ranker.read_drug_evidence(drugs_path))
        row = next(item for item in result["results"] if item["disease"] == SEED)
        self.assertFalse(row["recovered"])
        self.assertEqual(row["best_rank"], 0)

    def test_too_few_cases_report_no_rate(self) -> None:
        """Two cases are below the floor of five, so the run reports INSUFFICIENT_CASES."""
        result = controls.heldout_control(self.graph, 3, controls.TREATMENT_RELATIONS,
                                          min_cases=5, evidence=self.evidence)
        self.assertEqual(result["verdict"], "INSUFFICIENT_CASES")
        self.assertIsNone(result["full"]["recovery_rate"])

    def test_the_fixture_recovers_its_nine_planted_drugs(self) -> None:
        """The fixture plants a pharmacological route for nine of its ten cases.

        The tenth disease keeps its drugs behind no gene route at all, so no ranking can return
        them. The rate is nine of ten and not ten of ten, which is the honest result. The run uses
        20 permutations, because a PASS needs a chance probability at or below 0.05 and the
        smallest value that B permutations can reach is 1 / (B + 1).
        """
        result = controls.heldout_control(self.fixture_graph, 3, controls.TREATMENT_RELATIONS,
                                          permutations=20, generator=random.Random(2),
                                          evidence=self.fixture_evidence)
        self.assertEqual(result["full"]["diseases"], fixture.HELD_OUT_CASES + 1)
        self.assertEqual(result["full"]["recovered"], fixture.HELD_OUT_CASES)
        self.assertEqual(result["full"]["recovery_rate"], 0.9)
        self.assertLess(result["full"]["null"]["recovery_rate"], result["full"]["recovery_rate"])
        self.assertEqual(result["verdict"], "PASS")

    def test_a_run_below_the_permutation_floor_reports_it(self) -> None:
        """A PASS needs p at or below 0.05, which 5 permutations cannot reach."""
        result = controls.heldout_control(self.fixture_graph, 3, controls.TREATMENT_RELATIONS,
                                          permutations=5, generator=random.Random(2),
                                          evidence=self.fixture_evidence)
        self.assertEqual(result["verdict"], "INSUFFICIENT_PERMUTATIONS")
        self.assertEqual(result["full"]["smallest_reachable_permutation_probability"],
                         float(f"{1 / 6:.3g}"))

    def test_the_permutation_quantile_never_falls_below_its_floor(self) -> None:
        """The quantile is (1 + reached) / (B + 1), so it is at least 1 / (B + 1)."""
        result = controls.heldout_control(self.fixture_graph, 3, controls.TREATMENT_RELATIONS,
                                          permutations=20, generator=random.Random(2),
                                          evidence=self.fixture_evidence)
        summary = result["full"]
        self.assertGreaterEqual(summary["permutation_quantile_probability"],
                                summary["smallest_reachable_permutation_probability"])
        self.assertEqual(controls.permutation_quantile(9, [0, 0, 0, 0]), 0.2)
        self.assertEqual(controls.permutation_quantile(9, [9, 0, 0, 0]), 0.4)

    def test_the_null_rate_spread_reads_the_permutations_as_the_replicates(self) -> None:
        """The cases of one permutation share one graph, so B and not B times D is the count."""
        result = controls.heldout_control(self.fixture_graph, 3, controls.TREATMENT_RELATIONS,
                                          permutations=20, generator=random.Random(2),
                                          evidence=self.fixture_evidence)
        spread = result["full"]["null_rate_spread"]
        self.assertEqual(spread["permutations_used"], 20)
        self.assertEqual(len(spread["rates"]), 20)
        self.assertLessEqual(spread["lower"], spread["mean"])
        self.assertLessEqual(spread["mean"], spread["upper"])

    def test_the_tail_at_the_null_upper_limit_is_the_larger_probability(self) -> None:
        """A larger null rate cannot lower the upper tail of the recovery count."""
        result = controls.heldout_control(self.fixture_graph, 3, controls.TREATMENT_RELATIONS,
                                          permutations=20, generator=random.Random(2),
                                          evidence=self.fixture_evidence)
        summary = result["full"]
        self.assertGreaterEqual(summary["chance_probability_at_null_upper_limit"],
                                summary["chance_probability"])

    def test_a_dispersion_of_one_reproduces_the_binomial_tail(self) -> None:
        """At no dispersion the beta-binomial is the binomial, so the correction adds nothing."""
        binomial = controls.poisson_binomial_tail([0.3] * 40, 18)
        corrected = controls.beta_binomial_tail(40, 0.3, 1.0, 18)
        self.assertAlmostEqual(binomial, corrected, places=12)

    def test_a_wider_null_raises_the_far_upper_tail(self) -> None:
        """More variance puts more mass far above the mean, so the far tail grows."""
        narrow = controls.beta_binomial_tail(115, 0.3092, 1.0, 50)
        wide = controls.beta_binomial_tail(115, 0.3092, 1.54, 50)
        self.assertGreater(wide, narrow)
        self.assertLessEqual(wide, 1.0)

    def test_a_wider_null_does_not_raise_the_tail_at_the_mean(self) -> None:
        """A count at the mean sits inside the body, where more variance moves mass both ways."""
        mean_count = round(115 * 0.3092)
        narrow = controls.beta_binomial_tail(115, 0.3092, 1.0, mean_count)
        wide = controls.beta_binomial_tail(115, 0.3092, 1.54, mean_count)
        self.assertLessEqual(wide, narrow)
        self.assertGreater(wide, 0.4)

    def test_the_corrected_distribution_carries_the_two_target_moments(self) -> None:
        """The fitted mean is n times the rate and the fitted variance is phi times the binomial."""
        for cases, rate, dispersion in ((115, 0.3092, 1.0639), (54, 0.4769, 1.2107),
                                        (61, 0.0969, 1.5355)):
            mass = [controls.beta_binomial_tail(cases, rate, dispersion, hits)
                    - controls.beta_binomial_tail(cases, rate, dispersion, hits + 1)
                    for hits in range(cases + 1)]
            self.assertAlmostEqual(sum(mass), 1.0, places=9)
            mean = sum(hits * value for hits, value in enumerate(mass))
            variance = sum((hits - mean) ** 2 * value for hits, value in enumerate(mass))
            self.assertAlmostEqual(mean, cases * rate, places=6)
            self.assertAlmostEqual(variance, dispersion * cases * rate * (1.0 - rate), places=6)

    def test_the_dispersion_factor_reads_the_permuted_counts(self) -> None:
        """The factor is the sample variance of the counts over the binomial variance."""
        counts = [10, 12, 8, 14, 6, 10, 12, 8]
        factor = controls.dispersion_factor(counts, 40, 0.25)
        self.assertAlmostEqual(factor, statistics.variance(counts) / (40 * 0.25 * 0.75), places=10)
        self.assertIsNone(controls.dispersion_factor([7], 40, 0.25))
        self.assertIsNone(controls.dispersion_factor([7, 8], 40, 0.0))

    def test_the_correction_never_lowers_the_tail_below_the_binomial(self) -> None:
        """The measured error has one direction, so a factor below 1 is not applied."""
        under = controls.corrected_tail([10, 10, 10, 11, 10], 40, 0.25, 18)
        self.assertLess(under["dispersion_measured"], 1.0)
        self.assertEqual(under["dispersion_applied"], 1.0)
        self.assertAlmostEqual(under["tail"],
                               float(f"{controls.poisson_binomial_tail([0.25] * 40, 18):.3g}"),
                               places=12)

    def test_the_decision_probability_is_the_larger_of_the_two(self) -> None:
        """A PASS needs the assumption-free value and the corrected model value under alpha."""
        result = controls.heldout_control(self.fixture_graph, 3, controls.TREATMENT_RELATIONS,
                                          permutations=20, generator=random.Random(2),
                                          evidence=self.fixture_evidence)
        summary = result["full"]
        self.assertEqual(summary["decision_probability"],
                         max(summary["permutation_quantile_probability"],
                             summary["dispersion_corrected_probability"]))
        self.assertGreaterEqual(summary["decision_probability"],
                                summary["permutation_quantile_probability"])

    def test_a_second_top_k_needs_no_second_run(self) -> None:
        """The stored best rank of every case gives the recovery at any size of the top list."""
        result = controls.heldout_control(self.fixture_graph, 3, controls.TREATMENT_RELATIONS,
                                          permutations=5, generator=random.Random(2),
                                          evidence=self.fixture_evidence, report_top_k=(1, 3))
        self.assertEqual(result["by_top_k"]["3"]["full"]["recovered"],
                         result["full"]["recovered"])
        self.assertLessEqual(result["by_top_k"]["1"]["full"]["recovered"],
                             result["by_top_k"]["3"]["full"]["recovered"])

    def test_the_worker_count_does_not_change_the_result(self) -> None:
        """The per_permutation stream seeds each permutation from the run seed and its number."""
        single = controls.heldout_control(
            self.fixture_graph, 3, controls.TREATMENT_RELATIONS, permutations=4,
            evidence=self.fixture_evidence, null_stream=controls.NULL_STREAM_PER_PERMUTATION,
            workers=1, random_seed=99)
        parallel = controls.heldout_control(
            self.fixture_graph, 3, controls.TREATMENT_RELATIONS, permutations=4,
            evidence=self.fixture_evidence, null_stream=controls.NULL_STREAM_PER_PERMUTATION,
            workers=3, random_seed=99)
        self.assertEqual(single["null_results"], parallel["null_results"])
        self.assertEqual(single["full"]["chance_probability"],
                         parallel["full"]["chance_probability"])

    def test_more_than_one_worker_refuses_the_sequential_stream(self) -> None:
        """One generator cannot be shared by several processes, so the run stops."""
        with self.assertRaises(ValueError):
            controls.heldout_control(self.fixture_graph, 3, controls.TREATMENT_RELATIONS,
                                     permutations=4, evidence=self.fixture_evidence,
                                     null_stream=controls.NULL_STREAM_SEQUENTIAL, workers=2)

    def test_the_resolution_of_a_written_receipt_matches_the_run(self) -> None:
        """A receipt carries enough to recompute its own strata and its own probabilities."""
        folder = Path(self.folder.name)
        nodes_path, edges_path, _ = write_graph(folder)
        result = controls.heldout_control(self.graph, 3, controls.TREATMENT_RELATIONS,
                                          permutations=6, generator=random.Random(2),
                                          evidence=self.evidence, min_cases=1, regime="sparse",
                                          limits={"max_genes": 3, "max_drugs": 1})
        record = {"receipt_id": "fixture", "graph": {"nodes_file": str(nodes_path),
                                                     "edges_file": str(edges_path),
                                                     "edges": len(self.graph.edges),
                                                     "edges_sha256_prefix": "0" * 16},
                  "parameters": {"top_k": 3, "min_cases": 1, "min_recovery": 0.3,
                                 "permutations": 6, "regime": "sparse",
                                 "regime_limits": {"max_genes": 3, "max_drugs": 1}},
                  "result": result}
        path = folder / "fixture_receipt.json"
        path.write_text(json.dumps(record), encoding="utf-8")
        rebuilt = controls.receipt_resolution(path, (3,))
        self.assertEqual(rebuilt["by_top_k"]["3"]["full"]["recovered"], result["full"]["recovered"])
        self.assertEqual(rebuilt["by_top_k"]["3"]["full"]["chance_probability"],
                         result["full"]["chance_probability"])
        self.assertEqual(rebuilt["by_top_k"]["3"]["restricted"]["cases"],
                         result["restricted"]["cases"])

    def test_a_truncated_receipt_cannot_be_rebuilt(self) -> None:
        """A receipt that dropped null rows must refuse a recomputation, not guess one."""
        folder = Path(self.folder.name)
        nodes_path, edges_path, _ = write_graph(folder)
        record = {"receipt_id": "fixture",
                  "graph": {"nodes_file": str(nodes_path), "edges_file": str(edges_path),
                            "edges": 1, "edges_sha256_prefix": "0" * 16},
                  "parameters": {"top_k": 3, "min_cases": 1, "min_recovery": 0.3,
                                 "permutations": 6, "regime": "all", "regime_limits": {}},
                  "result": {"results": [], "null_results": [], "null_results_truncated": True}}
        path = folder / "truncated_receipt.json"
        path.write_text(json.dumps(record), encoding="utf-8")
        with self.assertRaises(ValueError):
            controls.receipt_resolution(path, (3,))

    def test_the_receipt_keeps_the_null_rows(self) -> None:
        """Another regime and another top k must be measurable from the receipt with no re-run."""
        result = controls.heldout_control(self.fixture_graph, 3, controls.TREATMENT_RELATIONS,
                                          permutations=2, generator=random.Random(4),
                                          evidence=self.fixture_evidence, min_cases=1)
        self.assertEqual(len(result["null_results"]), 2 * result["full"]["diseases"])
        self.assertFalse(result["null_results_truncated"])
        self.assertEqual({row["permutation"] for row in result["null_results"]}, {1, 2})
        for row in result["null_results"]:
            self.assertIn("best_rank", row)
            self.assertIn("candidates_after_removal", row)

    def test_a_regime_reports_a_second_rate_beside_the_full_rate(self) -> None:
        """The full number and the restricted number both enter the receipt."""
        result = controls.heldout_control(self.fixture_graph, 3, controls.TREATMENT_RELATIONS,
                                          evidence=self.fixture_evidence, regime="sparse",
                                          limits={"max_genes": 3, "max_drugs": 1}, min_cases=1)
        self.assertIsNotNone(result["restricted"])
        self.assertEqual(result["restricted"]["regime"], "sparse")
        self.assertLessEqual(result["restricted"]["diseases"], result["full"]["diseases"])
        self.assertEqual(result["verdict"], result["restricted"]["verdict"])

    def test_the_mesh_supplementary_rule_reads_the_identifier(self) -> None:
        """MeSH files a rare disease as a supplementary concept record with a C number."""
        self.assertTrue(controls.regime_mesh_supplementary(
            {"disease_id": "DISEASE:MESH:C536987"}, {}))
        self.assertFalse(controls.regime_mesh_supplementary(
            {"disease_id": "DISEASE:MESH:D001943"}, {}))

    def test_the_sparse_rule_reads_both_thresholds(self) -> None:
        """A disease passes only when it is below the gene count and the drug count."""
        limits = {"max_genes": 5, "max_drugs": 2}
        self.assertTrue(controls.regime_sparse({"genes": 3, "drugs": 1}, limits))
        self.assertFalse(controls.regime_sparse({"genes": 30, "drugs": 1}, limits))
        self.assertFalse(controls.regime_sparse({"genes": 3, "drugs": 40}, limits))

    def test_the_drug_poor_rule_reads_the_drug_count_only(self) -> None:
        """A disease with many approved drugs gives a random ranking many chances to hit one."""
        limits = {"max_genes": 5, "max_drugs": 9}
        self.assertTrue(controls.regime_drug_poor({"genes": 300, "drugs": 2}, limits))
        self.assertFalse(controls.regime_drug_poor({"genes": 1, "drugs": 40}, limits))

    def test_a_regime_with_no_disease_reports_no_rate(self) -> None:
        """A restriction that the graph cannot satisfy must say so and report no number."""
        result = controls.heldout_control(self.fixture_graph, 3, controls.TREATMENT_RELATIONS,
                                          evidence=self.fixture_evidence,
                                          regime="mesh_supplementary", min_cases=5)
        self.assertEqual(result["restricted"]["diseases"], 0)
        self.assertIsNone(result["restricted"]["recovery_rate"])
        self.assertEqual(result["verdict"], "INSUFFICIENT_CASES")

    def test_the_disease_facts_count_genes_and_drugs(self) -> None:
        """The regime rules read counts that come from the graph and from no hand written list."""
        facts = controls.disease_facts(self.graph)[SEED]
        self.assertEqual(facts["genes"], 3)
        self.assertEqual(facts["drugs"], 1)

    def test_the_random_baseline_is_a_probability(self) -> None:
        """A pool of ten with one true drug and a top three gives three tenths."""
        self.assertAlmostEqual(controls.random_hit_probability(10, 1, 3), 0.3, places=10)
        self.assertEqual(controls.random_hit_probability(10, 0, 3), 0.0)
        self.assertEqual(controls.random_hit_probability(3, 1, 5), 1.0)

    def test_the_chance_probability_of_a_fair_coin(self) -> None:
        """Four hits from four trials at one half has probability one sixteenth."""
        self.assertAlmostEqual(controls.poisson_binomial_tail([0.5] * 4, 4), 0.0625, places=10)
        self.assertAlmostEqual(controls.poisson_binomial_tail([0.5] * 4, 0), 1.0, places=10)

    def setUp(self) -> None:
        """Write the hand graph and the fixture graph, then load both."""
        self.folder = tempfile.TemporaryDirectory()
        nodes_path, edges_path, drugs_path = write_graph(Path(self.folder.name))
        self.graph = ranker.load_graph(nodes_path, edges_path)
        self.evidence = ranker.read_drug_evidence(drugs_path)
        paths = fixture.write_fixture(Path(self.folder.name) / "fixture")
        self.fixture_graph = ranker.load_graph(paths[0], paths[1])
        self.fixture_evidence = ranker.read_drug_evidence(paths[2])

    def tearDown(self) -> None:
        """Remove the temporary directory."""
        self.folder.cleanup()


class NullCorrectionReportTest(unittest.TestCase):
    """Check the table that scripts/track2_null_correction.py builds for the report."""

    def test_a_probability_at_its_floor_prints_as_a_bound(self) -> None:
        """A comparison against B permuted graphs cannot resolve a value below 1 / (B + 1)."""
        self.assertEqual(correction.reported_probability(0.00498, 0.00498), "0.005 or less")
        self.assertEqual(correction.reported_probability(0.167, 0.167), "0.167 or less")

    def test_a_probability_above_its_floor_prints_as_a_point_value(self) -> None:
        """A value that the permutation evidence reaches is a measurement and not a bound."""
        self.assertEqual(correction.reported_probability(0.00995, 0.00498), "0.00995")
        self.assertEqual(correction.reported_probability(0.5, 0.167), "0.500")

    def test_a_stratum_with_no_probability_says_so(self) -> None:
        """A stratum below the case floor reports no rate and no probability."""
        self.assertEqual(correction.reported_probability(None, 0.00498), "not reported")


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
