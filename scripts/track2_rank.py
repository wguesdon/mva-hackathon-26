"""Rank chemical nodes of the Track 2 knowledge graph as drug repurposing candidates.

The ranker answers one question: starting from a seed disease or a seed gene set, which chemical
nodes sit at the end of a short, well supported path, and what is that path?

The score is deliberately small and readable, because a panel reads the explanation and not the
number. For one candidate chemical the ranker keeps every shortest path of at most `--max-edges`
edges and scores the candidate from the best of them:

    coverage(p)   = share of the edges of p that carry at least one PMID
    path_term(p)  = log1p(min evidence_count over the edges of p)
                  + coverage(p) * log1p(number of distinct PMIDs over the edges of p)
    score(c)      = LENGTH_DECAY ** (L - 1)
                    * (max path_term over the shortest paths
                       + BREADTH_WEIGHT * log1p(number of shortest paths - 1))

The weakest edge of a path caps the inference, so the evidence term uses the minimum and not the
sum. Distinct PMIDs count independent literature support. The coverage factor stops a well studied
first edge from carrying a last edge that no publication supports. The breadth term rewards a
candidate that several shortest paths reach, and it is the term that degree bias inflates, which is
why `scripts/track2_controls.py` must run before any candidate is reported.

The score alone ranks a laboratory reagent first. A PubTator3 `associate` edge between a gene and a
chemical states that the two appear in the same publication. The chemicals that share the most
publications with a spindle assembly checkpoint gene are the compounds that arrest mitosis in those
experiments, and they carry more PMIDs than any drug. Two constraints separate a co-mention from a
pharmacological claim.

1. The terminal constraint. The last edge into the candidate must belong to a class in
   `--terminal-classes`, which holds `pharmacological` and `clinical` by default. A co-mention edge
   then supports a step inside the path but cannot nominate a drug by itself. `--terminal-classes
   all` turns the constraint off and reproduces the unconstrained ranking.
2. The eligibility constraint. The candidate must carry a drug record in `drug_evidence.tsv` or an
   indication edge in the graph. A compound that was never a drug is not a repurposing candidate.
   `--no-require-drug-evidence` turns the constraint off.

`RELATION_CLASS` declares the class of every relation in one place. An unknown relation falls back
to `co_mention`, which is the conservative class, and the run prints its name.

The graph is undirected for path finding. A relation keeps its own name in the output, so a reader
sees the mechanism of each step and can reject a candidate on mechanism alone.

Usage:
    uv run python scripts/track2_rank.py --seed DISEASE:MESH:C536987
    uv run python scripts/track2_rank.py --seed DISEASE:MESH:C536987 --terminal-classes all \
        --no-require-drug-evidence
    uv run python scripts/track2_rank.py --seed DISEASE:MESH:C536987 --gene-layer functional

The `--gene-layer` option chooses the gene to gene edges. The default, `comention`, is the graph as
built. `functional` and `union` read the STRING layer of `scripts/build_track2_functional_layer.py`.
"""

from __future__ import annotations

import argparse
import csv
import functools
import math
import operator
import sys
from collections import deque
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GRAPH_DIR = ROOT / "reports" / "track2" / "graph"
DEFAULT_NODES = DEFAULT_GRAPH_DIR / "nodes.tsv"
DEFAULT_EDGES = DEFAULT_GRAPH_DIR / "edges.tsv"
DEFAULT_DRUG_EVIDENCE = DEFAULT_GRAPH_DIR / "drug_evidence.tsv"
DEFAULT_OUT_DIR = ROOT / "reports" / "track2" / "ranking"
DEFAULT_FUNCTIONAL_DIR = ROOT / "reports" / "track2" / "functional"

# The gene to gene layer of the graph. `comention` is the default graph as built, where every gene
# to gene edge is a PubTator3 co-mention. `functional` replaces those edges with the STRING layer of
# `scripts/build_track2_functional_layer.py`. `union` keeps both. Only the gene to gene edges
# change. A disease, compound or pathway edge of the default graph stays in every layer.
GENE_LAYER_COMENTION = "comention"
GENE_LAYER_FUNCTIONAL = "functional"
GENE_LAYER_UNION = "union"
GENE_LAYERS = (GENE_LAYER_COMENTION, GENE_LAYER_FUNCTIONAL, GENE_LAYER_UNION)
LITERATURE_SOURCES = ("pubtator3", "pubtator3_consequence")

NODE_COLUMNS = ("node_id", "kind", "label", "source")
EDGE_COLUMNS = ("source_id", "target_id", "relation", "evidence_count", "pmids", "source")
DRUG_COLUMNS = ("drug_id", "target_gene_id", "action_type", "max_clinical_stage_overall",
                "mechanism_of_action")

LENGTH_DECAY = 0.5
BREADTH_WEIGHT = 0.25
MAX_PATHS_PER_CANDIDATE = 200

CANDIDATE_KIND = "chemical"

# The kind of a node that the node table does not name. The graph builder writes CHEMICAL for a
# PubTator3 chemical and CHEMBL for a ChEMBL compound, and both are chemical nodes here.
PREFIX_KIND = {"GENE": "gene", "DISEASE": "disease", "CHEMICAL": "chemical", "CHEMBL": "chemical",
               "DRUG": "chemical", "PATHWAY": "pathway", "REACTOME": "pathway",
               "MONDO": "disease"}

PHARMACOLOGICAL = "pharmacological"
CLINICAL = "clinical"
CO_MENTION = "co_mention"
STRUCTURAL = "structural"
FUNCTIONAL = "functional"
RELATION_CLASSES = (PHARMACOLOGICAL, CLINICAL, CO_MENTION, STRUCTURAL, FUNCTIONAL)

# The class of every relation of the graph, declared in one place.
#
# pharmacological: an Open Targets mechanism of action. The record states what the drug does to the
#     target, so the edge supports a claim about a drug.
# clinical: a ChEMBL indication. The record states that the drug is used against the disease.
# co_mention: a PubTator3 relation. The record states that the two entities appear in the same
#     publication. The relation name is the direction of the sentence, not a pharmacological claim.
#     PubTator3 `inhibit` is a co-mention and Open Targets `inhibitor` is not. The two names differ
#     by one letter and the class separates them.
# structural: a membership edge such as a gene inside a pathway.
# functional: a STRING interaction from the experimental and database channels. The record states
#     that two proteins interact or share a curated complex or pathway. It carries no sign and no
#     PMID, and it is present only when the optional functional gene layer is loaded.
RELATION_CLASS = {
    "inhibitor": PHARMACOLOGICAL,
    "agonist": PHARMACOLOGICAL,
    "antagonist": PHARMACOLOGICAL,
    "modulator": PHARMACOLOGICAL,
    "positive modulator": PHARMACOLOGICAL,
    "negative modulator": PHARMACOLOGICAL,
    "negative allosteric modulator": PHARMACOLOGICAL,
    "activator": PHARMACOLOGICAL,
    "degrader": PHARMACOLOGICAL,
    "binding agent": PHARMACOLOGICAL,
    "stabiliser": PHARMACOLOGICAL,
    "hydrolytic enzyme": PHARMACOLOGICAL,
    "antisense inhibitor": PHARMACOLOGICAL,
    "rnai inhibitor": PHARMACOLOGICAL,
    "cross-linking agent": PHARMACOLOGICAL,
    "exogenous gene": PHARMACOLOGICAL,
    "vaccine antigen": PHARMACOLOGICAL,
    "other": PHARMACOLOGICAL,
    "indication": CLINICAL,
    "associate": CO_MENTION,
    "positive_correlate": CO_MENTION,
    "negative_correlate": CO_MENTION,
    "stimulate": CO_MENTION,
    "inhibit": CO_MENTION,
    "interact": CO_MENTION,
    "compare": CO_MENTION,
    "cotreat": CO_MENTION,
    "drug_interact": CO_MENTION,
    "prevent": CO_MENTION,
    "treat": CO_MENTION,
    "cause": CO_MENTION,
    "participates_in": STRUCTURAL,
    "functional_interaction": FUNCTIONAL,
}
DEFAULT_RELATION_CLASS = CO_MENTION
DEFAULT_TERMINAL_CLASSES = (PHARMACOLOGICAL, CLINICAL)

INCREASE = "increase"
DECREASE = "decrease"
UNCLEAR = "unclear"
OPPOSITE_DIRECTION = {INCREASE: DECREASE, DECREASE: INCREASE}

# What an Open Targets action type does to the function of its target, declared in one place.
# A name that is absent returns `unclear`, and the ranker then states that the direction is unknown
# instead of guessing it.
ACTION_DIRECTION = {
    "inhibitor": DECREASE,
    "antagonist": DECREASE,
    "degrader": DECREASE,
    "antisense inhibitor": DECREASE,
    "rnai inhibitor": DECREASE,
    "negative modulator": DECREASE,
    "negative allosteric modulator": DECREASE,
    "agonist": INCREASE,
    "activator": INCREASE,
    "positive modulator": INCREASE,
    "stabiliser": INCREASE,
    "exogenous gene": INCREASE,
}

# The direction that the defect itself pushes the function of the defect gene.
DEFECT_DIRECTION = {"loss_of_function": DECREASE, "gain_of_function": INCREASE, "none": ""}

# The order of the verdicts. A candidate that reinforces the defect sinks below every other
# candidate. It stays in the table with its reason, because a reviewer must see the rejected class.
DIRECTION_TIER = {"restores_target": 0, "restores_pathway": 1, "unclear_action": 2,
                  "unclear_target": 2, "unclear_composition": 2, "not_assessed": 2,
                  "reinforces_pathway": 3, "reinforces_target": 4}

OUTPUT_COLUMNS = ("rank", "chembl_id", "label", "score", "direction_verdict", "direction_reason",
                  "action_type", "action_direction", "path_length", "path", "relations",
                  "terminal_relation", "terminal_class", "terminal_target", "clinical_stage",
                  "mechanism", "eligibility", "pmids", "n_pmids", "pmid_coverage",
                  "n_shortest_paths", "min_evidence_count", "path_kinds", "intermediates", "seed")


@dataclass(frozen=True)
class Node:
    """One node of the graph.

    Attributes:
        node_id: The identifier, for example `GENE:BUB1B`.
        kind: One of gene, disease, chemical, pathway.
        label: The human readable name.
        source: The database the node comes from.
    """

    node_id: str
    kind: str
    label: str
    source: str = ""


@dataclass(frozen=True)
class Edge:
    """One undirected edge with its provenance.

    Attributes:
        source_id: The first endpoint.
        target_id: The second endpoint.
        relation: The relation name of the source database.
        evidence_count: The number of supporting records the source database reports.
        pmids: The supporting PubMed identifiers.
        source: The database the edge comes from.
    """

    source_id: str
    target_id: str
    relation: str
    evidence_count: int
    pmids: tuple[str, ...] = ()
    source: str = ""


@dataclass
class Graph:
    """A node table, an edge list and the undirected adjacency built from them.

    Attributes:
        nodes: Every node by identifier.
        edges: Every edge, in file order after the duplicate merge.
        adjacency: For each node, the pairs of neighbour identifier and edge index.
        stats: Counts collected while the graph was loaded.
    """

    nodes: dict[str, Node] = field(default_factory=dict)
    edges: list[Edge] = field(default_factory=list)
    adjacency: dict[str, list[tuple[str, int]]] = field(default_factory=dict)
    stats: dict[str, int] = field(default_factory=dict)

    def degree(self, node_id: str) -> int:
        """Return the undirected degree of one node.

        Args:
            node_id: The node identifier.

        Returns:
            The number of incident edge endpoints.
        """
        return len(self.adjacency.get(node_id, ()))

    def label(self, node_id: str) -> str:
        """Return the label of one node.

        Args:
            node_id: The node identifier.

        Returns:
            The label, or the identifier when the node table holds no label.
        """
        node = self.nodes.get(node_id)
        return node.label if node and node.label else node_id

    def kind(self, node_id: str) -> str:
        """Return the kind of one node.

        Args:
            node_id: The node identifier.

        Returns:
            The kind string, or the kind inferred from the identifier prefix.
        """
        node = self.nodes.get(node_id)
        return node.kind if node else kind_from_id(node_id)


@dataclass(frozen=True)
class GraphPath:
    """One shortest path from a seed to a candidate.

    Attributes:
        node_ids: The node identifiers, from the seed to the candidate.
        edge_indices: The edge index of every step.
    """

    node_ids: tuple[str, ...]
    edge_indices: tuple[int, ...]

    @property
    def length(self) -> int:
        """Return the number of edges of the path.

        Returns:
            The edge count.
        """
        return len(self.edge_indices)


@dataclass(frozen=True)
class RankOptions:
    """The parameters of one ranking run.

    Attributes:
        max_edges: The largest path length in edges.
        hub_degree_cap: When above 0, a node of higher degree is never expanded.
        max_paths: The largest number of shortest paths kept per candidate.
        terminal_classes: The relation classes that the last edge of a path may carry. Every class
            of `RELATION_CLASSES` turns the terminal constraint off.
        extra_pharmacological: Relation names that count as pharmacological on top of the declared
            table.
        require_drug_evidence: Keep a candidate only when it carries a drug record or an indication
            edge.
        defect: The direction of the defect of the seed, `loss_of_function`, `gain_of_function` or
            `none`. The caller declares it, because the graph does not carry it.
        defect_genes: The genes that carry the defect. A gene seed enters this set on its own.
        blocked_intermediate_kinds: The node kinds that a path may not pass through. A disease is
            blocked by default. A path that leaves the seed disease, reaches a gene and returns to a
            second disease reads its drugs from that second disease, and a co-mention does not carry
            the claim that the two diseases share a treatment.
        candidate_kind: The node kind to rank.
    """

    max_edges: int = 3
    hub_degree_cap: int = 0
    max_paths: int = MAX_PATHS_PER_CANDIDATE
    terminal_classes: tuple[str, ...] = DEFAULT_TERMINAL_CLASSES
    extra_pharmacological: tuple[str, ...] = ()
    require_drug_evidence: bool = True
    defect: str = "none"
    defect_genes: tuple[str, ...] = ()
    blocked_intermediate_kinds: tuple[str, ...] = ("disease",)
    candidate_kind: str = CANDIDATE_KIND


@dataclass
class DrugEvidence:
    """The Open Targets drug records of the graph.

    Attributes:
        stage: The clinical stage of every drug, as the drug evidence table reports it.
        actions: The action type and the mechanism of every drug and target gene pair.
    """

    stage: dict[str, str] = field(default_factory=dict)
    actions: dict[tuple[str, str], tuple[str, str]] = field(default_factory=dict)

    def has_record(self, drug_id: str) -> bool:
        """Say whether one chemical carries a drug record.

        Args:
            drug_id: The chemical node identifier.

        Returns:
            True when the drug evidence table names the chemical.
        """
        return drug_id in self.stage

    def action(self, drug_id: str, gene_id: str) -> str:
        """Return the recorded action type of one drug on one gene.

        Args:
            drug_id: The chemical node identifier.
            gene_id: The gene node identifier.

        Returns:
            The action type, or an empty string when the pair has no record.
        """
        return self.actions.get((drug_id, gene_id), ("", ""))[0]

    def mechanism(self, drug_id: str, gene_id: str) -> str:
        """Return the mechanism of action of one drug on one gene.

        Args:
            drug_id: The chemical node identifier.
            gene_id: The gene node identifier.

        Returns:
            The mechanism text, or an empty string when the pair has no record.
        """
        return self.actions.get((drug_id, gene_id), ("", ""))[1]


@dataclass
class Candidate:
    """One ranked chemical with the evidence that supports it.

    Attributes:
        node_id: The chemical node identifier.
        label: The chemical label.
        score: The final score.
        path_length: The shortest path length in edges.
        n_paths: The number of distinct shortest paths kept.
        best_path: The path with the strongest evidence term.
        pmids: The distinct PMIDs of the best path.
        min_evidence: The smallest evidence count along the best path.
        seed_id: The seed the best path starts from.
        terminal_relation: The relation of the last edge of the best path.
        terminal_class: The class of that relation.
        terminal_target: The node that the last edge joins to the candidate.
        clinical_stage: The highest clinical stage of the drug record.
        mechanism: The mechanism of action on the gene of the last step.
        eligibility: Why the candidate is a drug at all.
        action_type: The action of the last edge, for example `inhibitor`.
        action_direction: What that action does to the function of the target.
        direction_verdict: Whether the action restores or reinforces the defect of the seed.
        direction_reason: One sentence that states why the verdict holds.
    """

    node_id: str
    label: str
    score: float
    path_length: int
    n_paths: int
    best_path: GraphPath
    pmids: tuple[str, ...]
    min_evidence: int
    seed_id: str
    terminal_relation: str = ""
    terminal_class: str = ""
    terminal_target: str = ""
    clinical_stage: str = ""
    mechanism: str = ""
    eligibility: str = ""
    action_type: str = ""
    action_direction: str = ""
    direction_verdict: str = "not_assessed"
    direction_reason: str = ""


def kind_from_id(node_id: str) -> str:
    """Infer a node kind from the identifier prefix.

    Args:
        node_id: A node identifier such as `GENE:BUB1B`.

    Returns:
        The inferred kind, or "unknown".
    """
    prefix = node_id.split(":", 1)[0].upper()
    return PREFIX_KIND.get(prefix, "unknown")


def relation_class(relation: str, extra_pharmacological: tuple[str, ...] = ()) -> str:
    """Return the class of one relation.

    Args:
        relation: The relation name of the source database.
        extra_pharmacological: Relation names that count as pharmacological on top of the table.

    Returns:
        One of `RELATION_CLASSES`. An unknown relation returns the conservative class `co_mention`.
    """
    name = (relation or "").strip().lower()
    if name in {item.strip().lower() for item in extra_pharmacological}:
        return PHARMACOLOGICAL
    return RELATION_CLASS.get(name, DEFAULT_RELATION_CLASS)


def unclassified_relations(graph: Graph, extra_pharmacological: tuple[str, ...] = ()) -> list[str]:
    """List the relations of a graph that the declared table does not name.

    A new action type of a source database must not enter the co-mention class in silence. The
    command line run prints this list.

    Args:
        graph: The graph.
        extra_pharmacological: Relation names given on the command line.

    Returns:
        The unknown relation names in sorted order.
    """
    known = set(RELATION_CLASS) | {item.strip().lower() for item in extra_pharmacological}
    return sorted({edge.relation.strip().lower() for edge in graph.edges
                   if edge.relation.strip().lower() not in known})


def check_columns(header: list[str] | None, required: tuple[str, ...], path: Path) -> None:
    """Stop when a table does not carry the columns of the input contract.

    Args:
        header: The column names read from the file.
        required: The column names the contract demands.
        path: The file, used in the error message.

    Raises:
        ValueError: When the header is absent or a column is missing.
    """
    if header is None:
        raise ValueError(f"{path} holds no header row.")
    missing = [name for name in required if name not in header]
    if missing:
        raise ValueError(f"{path} misses the columns {', '.join(missing)}.")


def parse_pmids(text: str) -> tuple[str, ...]:
    """Split the comma separated PMID field of an edge.

    Args:
        text: The raw field value, possibly empty.

    Returns:
        The distinct PMIDs in sorted order.
    """
    items = {part.strip() for part in (text or "").split(",")}
    return tuple(sorted(item for item in items if item))


def parse_evidence(text: str) -> int:
    """Read the evidence count of an edge.

    Args:
        text: The raw field value.

    Returns:
        The integer value, or 0 when the field is empty or not a number.
    """
    try:
        return max(0, int(float(text)))
    except (TypeError, ValueError):
        return 0


def read_nodes(path: Path) -> dict[str, Node]:
    """Read the node table.

    Args:
        path: The nodes.tsv path.

    Returns:
        Every node by identifier.
    """
    nodes: dict[str, Node] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        check_columns(reader.fieldnames, NODE_COLUMNS, path)
        for row in reader:
            node_id = (row["node_id"] or "").strip()
            if not node_id:
                continue
            nodes[node_id] = Node(node_id=node_id, kind=(row["kind"] or "").strip().lower(),
                                  label=(row["label"] or "").strip() or node_id,
                                  source=(row["source"] or "").strip())
    return nodes


def read_edges(path: Path) -> list[Edge]:
    """Read the edge table without merging duplicates.

    Args:
        path: The edges.tsv path.

    Returns:
        One Edge for every data row that carries two endpoints.
    """
    edges: list[Edge] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        check_columns(reader.fieldnames, EDGE_COLUMNS, path)
        for row in reader:
            source_id = (row["source_id"] or "").strip()
            target_id = (row["target_id"] or "").strip()
            if not source_id or not target_id:
                continue
            edges.append(Edge(source_id=source_id, target_id=target_id,
                              relation=(row["relation"] or "").strip() or "related",
                              evidence_count=parse_evidence(row["evidence_count"]),
                              pmids=parse_pmids(row["pmids"]),
                              source=(row["source"] or "").strip()))
    return edges


def read_drug_evidence(path: Path) -> DrugEvidence:
    """Read the Open Targets drug evidence table.

    Args:
        path: The drug_evidence.tsv path.

    Returns:
        The drug records. A missing file returns an empty record set.
    """
    evidence = DrugEvidence()
    if not Path(path).is_file():
        return evidence
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        check_columns(reader.fieldnames, DRUG_COLUMNS, Path(path))
        for row in reader:
            drug_id = (row["drug_id"] or "").strip()
            if not drug_id:
                continue
            # The column is already the maximum over the whole drug, so the first value of a drug
            # is the value of every row of that drug. An empty cell must not hide a later value.
            stage = (row["max_clinical_stage_overall"] or "").strip()
            if stage or drug_id not in evidence.stage:
                evidence.stage.setdefault(drug_id, stage)
            gene_id = (row["target_gene_id"] or "").strip()
            if gene_id:
                evidence.actions[(drug_id, gene_id)] = ((row["action_type"] or "").strip(),
                                                        (row["mechanism_of_action"] or "").strip())
    return evidence


def merge_edges(edges: list[Edge]) -> tuple[list[Edge], int]:
    """Merge parallel edges that share the endpoint pair and the relation.

    Two rows that state the same relation between the same pair are one piece of evidence with two
    records. The merge sums the evidence counts and unions the PMIDs, so the same fact cannot raise
    the breadth term twice.

    Args:
        edges: The edges as read from the file.

    Returns:
        The merged edges and the number of rows removed by the merge.
    """
    merged: dict[tuple[str, str, str], Edge] = {}
    order: list[tuple[str, str, str]] = []
    for edge in edges:
        first, second = sorted((edge.source_id, edge.target_id))
        key = (first, second, edge.relation)
        current = merged.get(key)
        if current is None:
            merged[key] = edge
            order.append(key)
            continue
        sources = sorted({part for part in (current.source, edge.source) if part})
        merged[key] = Edge(source_id=current.source_id, target_id=current.target_id,
                           relation=current.relation,
                           evidence_count=current.evidence_count + edge.evidence_count,
                           pmids=tuple(sorted(set(current.pmids) | set(edge.pmids))),
                           source=",".join(sources))
    return [merged[key] for key in order], len(edges) - len(order)


def build_graph(nodes: dict[str, Node], edges: list[Edge], merge: bool = True) -> Graph:
    """Assemble a Graph and its undirected adjacency.

    An edge endpoint that the node table does not name becomes a node whose kind comes from the
    identifier prefix. The count of such endpoints enters the statistics.

    An edge that joins a node to itself is dropped. No simple path can use it, and a graph that
    holds one has a degree that the adjacency does not show, which breaks the degree check of the
    permutation control.

    Args:
        nodes: The node table.
        edges: The edge list.
        merge: Merge parallel edges that share the endpoint pair and the relation.

    Returns:
        The assembled graph.
    """
    merged, merged_away = (merge_edges(edges) if merge else (list(edges), 0))
    kept = [edge for edge in merged if edge.source_id != edge.target_id]
    self_loops = len(merged) - len(kept)
    graph = Graph(nodes=dict(nodes), edges=kept)
    implied = 0
    for index, edge in enumerate(kept):
        for node_id in (edge.source_id, edge.target_id):
            if node_id not in graph.nodes:
                graph.nodes[node_id] = Node(node_id=node_id, kind=kind_from_id(node_id),
                                            label=node_id, source="implied_by_edge")
                implied += 1
        graph.adjacency.setdefault(edge.source_id, []).append((edge.target_id, index))
        graph.adjacency.setdefault(edge.target_id, []).append((edge.source_id, index))
    for node_id in graph.nodes:
        graph.adjacency.setdefault(node_id, [])
    graph.stats = {"nodes": len(graph.nodes), "edges": len(kept), "merged_away": merged_away,
                   "implied_nodes": implied, "self_loops_dropped": self_loops}
    return graph


def load_graph(nodes_path: Path, edges_path: Path, merge: bool = True) -> Graph:
    """Read both tables and build the graph.

    Args:
        nodes_path: The nodes.tsv path.
        edges_path: The edges.tsv path.
        merge: Merge parallel edges that share the endpoint pair and the relation.

    Returns:
        The assembled graph.
    """
    return build_graph(read_nodes(Path(nodes_path)), read_edges(Path(edges_path)), merge=merge)


def is_literature_gene_edge(nodes: dict[str, Node], edge: Edge) -> bool:
    """Say whether one edge is a gene to gene co-mention of the literature layer.

    Args:
        nodes: The node table.
        edge: The edge as read from the file.

    Returns:
        True when both endpoints are genes and the edge comes from a PubTator3 layer.
    """
    kinds = []
    for node_id in (edge.source_id, edge.target_id):
        node = nodes.get(node_id)
        kinds.append(node.kind if node else kind_from_id(node_id))
    return kinds == ["gene", "gene"] and edge.source in LITERATURE_SOURCES


def apply_gene_layer(nodes: dict[str, Node], edges: list[Edge], layer: str,
                     functional_nodes: dict[str, Node],
                     functional_edges: list[Edge]) -> tuple[dict[str, Node], list[Edge]]:
    """Choose the gene to gene layer of the graph.

    Args:
        nodes: The node table of the default graph.
        edges: The edge list of the default graph.
        layer: One of `GENE_LAYERS`.
        functional_nodes: The nodes that the functional layer adds.
        functional_edges: The edges of the functional layer.

    Returns:
        The node table and the edge list of the chosen layer. The default layer returns the input
        unchanged.

    Raises:
        ValueError: When the layer name is unknown.
    """
    if layer not in GENE_LAYERS:
        raise ValueError(f"Unknown gene layer: {layer}. Use one of {', '.join(GENE_LAYERS)}.")
    if layer == GENE_LAYER_COMENTION:
        return nodes, edges
    merged_nodes = dict(nodes)
    for node_id, node in functional_nodes.items():
        merged_nodes.setdefault(node_id, node)
    kept = edges
    if layer == GENE_LAYER_FUNCTIONAL:
        kept = [edge for edge in edges if not is_literature_gene_edge(nodes, edge)]
    return merged_nodes, list(kept) + list(functional_edges)


def load_layered_graph(nodes_path: Path, edges_path: Path, layer: str = GENE_LAYER_COMENTION,
                       functional_dir: Path = DEFAULT_FUNCTIONAL_DIR, merge: bool = True) -> Graph:
    """Read the default graph, choose its gene to gene layer and build the graph.

    The default layer reads no file of the functional layer, so it builds the same graph as
    `load_graph`.

    Args:
        nodes_path: The nodes.tsv path of the default graph.
        edges_path: The edges.tsv path of the default graph.
        layer: One of `GENE_LAYERS`.
        functional_dir: The directory of the functional layer tables.
        merge: Merge parallel edges that share the endpoint pair and the relation.

    Returns:
        The assembled graph.
    """
    if layer == GENE_LAYER_COMENTION:
        return load_graph(nodes_path, edges_path, merge=merge)
    nodes, edges = apply_gene_layer(read_nodes(Path(nodes_path)), read_edges(Path(edges_path)),
                                    layer, read_nodes(Path(functional_dir) / "nodes.tsv"),
                                    read_edges(Path(functional_dir) / "edges.tsv"))
    return build_graph(nodes, edges, merge=merge)


def load_layered_drug_evidence(path: Path, layer: str = GENE_LAYER_COMENTION,
                               functional_dir: Path = DEFAULT_FUNCTIONAL_DIR) -> DrugEvidence:
    """Read the drug records of the default graph and, for another layer, of the new genes.

    Args:
        path: The drug_evidence.tsv path of the default graph.
        layer: One of `GENE_LAYERS`.
        functional_dir: The directory of the functional layer tables.

    Returns:
        The drug records. A record of the default graph is never overwritten.
    """
    evidence = read_drug_evidence(path)
    if layer == GENE_LAYER_COMENTION:
        return evidence
    extra = read_drug_evidence(Path(functional_dir) / "drug_evidence.tsv")
    for drug_id, stage in extra.stage.items():
        evidence.stage.setdefault(drug_id, stage)
    for key, value in extra.actions.items():
        evidence.actions.setdefault(key, value)
    return evidence


def traversable(graph: Graph, node_id: str, steps: int, options: RankOptions) -> bool:
    """Say whether a path may pass through one node.

    A seed always passes. A blocked kind never passes, and a node above the hub degree cap never
    passes. The search and the path reconstruction both call this function, so a node that the
    search refused to expand cannot return as an intermediate of a rebuilt path.

    Args:
        graph: The graph.
        node_id: The node identifier.
        steps: The distance of the node from the seed set.
        options: The run parameters.

    Returns:
        True when a path may pass through the node.
    """
    if steps == 0:
        return True
    if graph.kind(node_id) in options.blocked_intermediate_kinds:
        return False
    return options.hub_degree_cap <= 0 or graph.degree(node_id) <= options.hub_degree_cap


def shortest_distances(graph: Graph, seeds: list[str],
                       options: RankOptions) -> tuple[dict[str, int], dict[str, str]]:
    """Run a breadth first search from the seed set.

    Args:
        graph: The graph.
        seeds: The seed node identifiers.
        options: The run parameters. A node that `traversable` refuses is reached but not expanded,
            so it can end a path and cannot sit inside one.

    Returns:
        The distance of every reached node and the seed that reached it first.
    """
    distance: dict[str, int] = {}
    origin: dict[str, str] = {}
    queue: deque[str] = deque()
    for seed in seeds:
        if seed in graph.nodes and seed not in distance:
            distance[seed] = 0
            origin[seed] = seed
            queue.append(seed)
    while queue:
        current = queue.popleft()
        step = distance[current] + 1
        if step > options.max_edges:
            continue
        if not traversable(graph, current, distance[current], options):
            continue
        for neighbour, _ in graph.adjacency.get(current, ()):
            if neighbour in distance:
                continue
            distance[neighbour] = step
            origin[neighbour] = origin[current]
            queue.append(neighbour)
    return distance, origin


def shortest_paths(graph: Graph, target: str, distance: dict[str, int],
                   options: RankOptions) -> list[GraphPath]:
    """Collect the shortest paths from the seed set to one target.

    The walk goes backwards from the target and only accepts a neighbour whose distance is one less.
    Every path it returns therefore has the shortest length. A predecessor that `traversable`
    refuses is skipped, so a blocked node never sits inside a path.

    Args:
        graph: The graph.
        target: The target node identifier.
        distance: The distances of `shortest_distances`.
        options: The run parameters.

    Returns:
        The distinct shortest paths, sorted for a deterministic result.
    """
    if target not in distance:
        return []
    found: list[GraphPath] = []
    stack: list[tuple[str, tuple[str, ...], tuple[int, ...]]] = [(target, (target,), ())]
    while stack and len(found) < options.max_paths:
        current, node_ids, edge_indices = stack.pop()
        if distance[current] == 0:
            found.append(GraphPath(node_ids=node_ids, edge_indices=edge_indices))
            continue
        for neighbour, index in sorted(graph.adjacency.get(current, ())):
            steps = distance.get(neighbour, -1)
            if steps != distance[current] - 1 or neighbour in node_ids:
                continue
            if not traversable(graph, neighbour, steps, options):
                continue
            stack.append((neighbour, (neighbour,) + node_ids, (index,) + edge_indices))
    return sorted(distinct_paths(graph, found), key=operator.attrgetter("node_ids"))


def distinct_paths(graph: Graph, paths: list[GraphPath]) -> list[GraphPath]:
    """Remove the paths that repeat a node sequence and a relation sequence.

    Args:
        graph: The graph.
        paths: The collected paths.

    Returns:
        One path per distinct node and relation sequence.
    """
    unique: dict[tuple[tuple[str, ...], tuple[str, ...]], GraphPath] = {}
    for path in paths:
        key = (path.node_ids, tuple(graph.edges[index].relation for index in path.edge_indices))
        unique.setdefault(key, path)
    return list(unique.values())


def candidate_paths(graph: Graph, node_id: str, distance: dict[str, int], options: RankOptions,
                    prefix_cache: dict[str, list[GraphPath]] | None = None) -> list[GraphPath]:
    """Collect the shortest paths into one candidate that satisfy the terminal constraint.

    The last edge must carry a relation of an allowed class. The search therefore does not read the
    distance of the candidate itself. It reads the distance of every neighbour that an allowed edge
    joins to the candidate, and it keeps the shortest result. A candidate whose closest route is a
    co-mention can still enter the ranking through a longer pharmacological route.

    Args:
        graph: The graph.
        node_id: The candidate node identifier.
        distance: The distances of `shortest_distances`.
        options: The run parameters.
        prefix_cache: The shortest paths to each intermediate node, reused across candidates of one
            run. Several drugs bind the same gene, so the prefix is computed once for that gene.

    Returns:
        The distinct shortest allowed paths, sorted for a deterministic result.
    """
    cache = {} if prefix_cache is None else prefix_cache
    allowed = set(options.terminal_classes)
    reachable = []
    for neighbour, index in graph.adjacency.get(node_id, ()):
        steps = distance.get(neighbour)
        if steps is None or steps + 1 > options.max_edges or neighbour == node_id:
            continue
        if not traversable(graph, neighbour, steps, options):
            continue
        if relation_class(graph.edges[index].relation, options.extra_pharmacological) not in allowed:
            continue
        reachable.append((steps + 1, neighbour, index))
    if not reachable:
        return []
    best_length = min(item[0] for item in reachable)
    found: list[GraphPath] = []
    for length, neighbour, index in sorted(reachable):
        if length != best_length or len(found) >= options.max_paths:
            continue
        if neighbour not in cache:
            cache[neighbour] = shortest_paths(graph, neighbour, distance, options)
        for prefix in cache[neighbour]:
            if node_id in prefix.node_ids:
                continue
            found.append(GraphPath(node_ids=prefix.node_ids + (node_id,),
                                   edge_indices=prefix.edge_indices + (index,)))
            if len(found) >= options.max_paths:
                break
    return sorted(distinct_paths(graph, found), key=operator.attrgetter("node_ids"))


def path_pmids(graph: Graph, path: GraphPath) -> tuple[str, ...]:
    """Collect the distinct PMIDs of one path.

    Args:
        graph: The graph.
        path: The path.

    Returns:
        The distinct PMIDs in sorted order.
    """
    pmids: set[str] = set()
    for index in path.edge_indices:
        pmids.update(graph.edges[index].pmids)
    return tuple(sorted(pmids))


def path_min_evidence(graph: Graph, path: GraphPath) -> int:
    """Return the smallest evidence count along one path.

    Args:
        graph: The graph.
        path: The path.

    Returns:
        The minimum evidence count, or 0 for an empty path.
    """
    counts = [graph.edges[index].evidence_count for index in path.edge_indices]
    return min(counts) if counts else 0


def path_pmid_coverage(graph: Graph, path: GraphPath) -> float:
    """Return the share of the steps of a path that carry at least one PMID.

    Args:
        graph: The graph.
        path: The path.

    Returns:
        A value between 0 and 1. An empty path returns 0.
    """
    if not path.edge_indices:
        return 0.0
    covered = sum(1 for index in path.edge_indices if graph.edges[index].pmids)
    return covered / len(path.edge_indices)


def path_term(graph: Graph, path: GraphPath) -> float:
    """Score the evidence of one path.

    The weakest edge caps the inference, so the evidence term uses the minimum evidence count. The
    PMID term counts the distinct publications that support any step of the path, and the coverage
    factor scales it down when a step carries no publication at all.

    Args:
        graph: The graph.
        path: The path.

    Returns:
        The evidence term of the path.
    """
    literature = math.log1p(len(path_pmids(graph, path))) * path_pmid_coverage(graph, path)
    return math.log1p(path_min_evidence(graph, path)) + literature


def path_sort_key(graph: Graph, path: GraphPath) -> tuple[float, int, tuple[str, ...]]:
    """Return the sort key that puts the best supported path first.

    Args:
        graph: The graph.
        path: The path.

    Returns:
        The negative evidence term, the negative PMID count and the node identifiers. The last
        field makes the order deterministic when two paths carry the same evidence.
    """
    return -path_term(graph, path), -len(path_pmids(graph, path)), path.node_ids


def candidate_sort_key(candidate: Candidate) -> tuple[int, float, int, str]:
    """Return the sort key of the ranking.

    The direction tier comes first, so a candidate that reinforces the declared defect sinks below
    every other candidate whatever its evidence score. A run that declares no defect gives every
    candidate the same tier and the order is then the score alone.

    Args:
        candidate: The scored candidate.

    Returns:
        The direction tier, the negative score, the path length and the node identifier.
    """
    return (DIRECTION_TIER.get(candidate.direction_verdict, 2), -candidate.score,
            candidate.path_length, candidate.node_id)


def score_candidate(graph: Graph, node_id: str, paths: list[GraphPath]) -> Candidate:
    """Score one chemical from its shortest paths.

    Args:
        graph: The graph.
        node_id: The chemical node identifier.
        paths: The shortest paths from the seed set to that chemical.

    Returns:
        The scored candidate.

    Raises:
        ValueError: When the path list is empty.
    """
    if not paths:
        raise ValueError(f"{node_id} has no path to score.")
    ranked = sorted(paths, key=functools.partial(path_sort_key, graph))
    best = ranked[0]
    length = best.length
    breadth = BREADTH_WEIGHT * math.log1p(max(0, len(paths) - 1))
    score = (LENGTH_DECAY ** max(0, length - 1)) * (path_term(graph, best) + breadth)
    return Candidate(node_id=node_id, label=graph.label(node_id), score=score, path_length=length,
                     n_paths=len(paths), best_path=best, pmids=path_pmids(graph, best),
                     min_evidence=path_min_evidence(graph, best), seed_id=best.node_ids[0])


def indication_chemicals(graph: Graph) -> set[str]:
    """Collect every chemical that a clinical edge joins to a disease.

    Args:
        graph: The graph.

    Returns:
        The chemical node identifiers.
    """
    found: set[str] = set()
    for edge in graph.edges:
        if relation_class(edge.relation) != CLINICAL:
            continue
        for node_id in (edge.source_id, edge.target_id):
            if graph.kind(node_id) == CANDIDATE_KIND:
                found.add(node_id)
    return found


def eligibility_of(node_id: str, evidence: DrugEvidence, indications: set[str]) -> str:
    """Say why one chemical is a drug at all.

    A repurposing candidate must already be a drug. A compound that appears in the literature but
    never reached a drug record and never had an indication is a reagent and not a candidate.

    Args:
        node_id: The chemical node identifier.
        evidence: The drug records.
        indications: The chemicals that carry a clinical edge.

    Returns:
        `drug_record`, `indication_edge` or an empty string.
    """
    if evidence.has_record(node_id):
        return "drug_record"
    if node_id in indications:
        return "indication_edge"
    return ""


@dataclass
class DefectContext:
    """What the caller declares about the defect of the seed.

    Attributes:
        direction: The direction that the defect pushes the function of the defect gene. An empty
            string means that the caller declared no defect and that no verdict is given.
        genes: The defect genes.
        partners: For every gene that shares a pathway with a defect gene, the smallest shared
            pathway and the defect gene it shares it with.
    """

    direction: str = ""
    genes: frozenset[str] = frozenset()
    partners: dict[str, tuple[str, str]] = field(default_factory=dict)


def pathway_partners(graph: Graph, defect_genes: frozenset[str]) -> dict[str, tuple[str, str]]:
    """Find every gene that shares a pathway with a defect gene.

    The smallest shared pathway is the most specific one, so the reason text names that one.

    Args:
        graph: The graph.
        defect_genes: The defect genes.

    Returns:
        For each partner gene, the label of the smallest shared pathway and the defect gene.
    """
    best: dict[str, tuple[int, str, str]] = {}
    for defect_gene in sorted(defect_genes):
        for pathway, _ in graph.adjacency.get(defect_gene, ()):
            if graph.kind(pathway) != "pathway":
                continue
            members = [node for node, _ in graph.adjacency.get(pathway, ())
                       if graph.kind(node) == "gene"]
            for member in members:
                if member in defect_genes:
                    continue
                current = best.get(member)
                if current is None or len(members) < current[0]:
                    best[member] = (len(members), graph.label(pathway), defect_gene)
    return {gene: (label, defect_gene) for gene, (_, label, defect_gene) in best.items()}


def build_defect_context(graph: Graph, seeds: list[str], options: RankOptions) -> DefectContext:
    """Collect the defect direction, the defect genes and their pathway partners.

    A gene seed is a defect gene on its own. A disease seed carries no gene, so the caller names the
    defect genes with `--defect-gene`.

    Args:
        graph: The graph.
        seeds: The seed node identifiers.
        options: The run parameters.

    Returns:
        The defect context. A run without a declared defect returns an empty context.
    """
    direction = DEFECT_DIRECTION.get(options.defect, "")
    if not direction:
        return DefectContext()
    genes = {seed for seed in seeds if graph.kind(seed) == "gene"}
    genes.update(gene for gene in options.defect_genes if gene in graph.nodes)
    frozen = frozenset(genes)
    return DefectContext(direction=direction, genes=frozen,
                         partners=pathway_partners(graph, frozen))


def action_direction_of(action_type: str) -> str:
    """Return what one action type does to the function of its target.

    Args:
        action_type: An Open Targets action type, for example `inhibitor`.

    Returns:
        `increase`, `decrease` or `unclear`.
    """
    return ACTION_DIRECTION.get((action_type or "").strip().lower(), UNCLEAR)


def direction_verdict(graph: Graph, candidate: Candidate, defect: DefectContext) -> tuple[str, str]:
    """Judge whether the action of a candidate restores or reinforces the defect.

    The judgement covers one step only. The last edge states what the drug does to its target, and
    the verdict compares that action with the defect. A gene to gene step of this graph is a
    co-mention that carries no sign, so the direction does not compose along a longer path. A target
    that is neither a defect gene nor a pathway partner of one therefore returns
    `unclear_composition` and not a verdict.

    Args:
        graph: The graph.
        candidate: The scored candidate with its terminal fields filled.
        defect: The declared defect.

    Returns:
        The verdict and the sentence that states why it holds.
    """
    if not defect.direction:
        return "not_assessed", "The run declares no defect direction."
    target = candidate.terminal_target
    if graph.kind(target) != "gene":
        return "unclear_target", (f"The last edge is a {candidate.terminal_relation} to "
                                  f"{graph.label(target)}, which is not a gene, so it names no "
                                  f"target function.")
    action = candidate.action_direction
    if action == UNCLEAR:
        return "unclear_action", (f"The action type {candidate.action_type or 'unknown'} has no "
                                  f"declared direction.")
    restores = action == OPPOSITE_DIRECTION[defect.direction]
    if target in defect.genes:
        verdict = "restores_target" if restores else "reinforces_target"
        return verdict, (f"{candidate.action_type} gives {action} of {graph.label(target)}, the "
                         f"declared defect gene, and the defect gives {defect.direction}.")
    if target in defect.partners:
        pathway, defect_gene = defect.partners[target]
        verdict = "restores_pathway" if restores else "reinforces_pathway"
        return verdict, (f"{candidate.action_type} gives {action} of {graph.label(target)}, which "
                         f"shares the pathway {pathway} with {graph.label(defect_gene)}.")
    return "unclear_composition", (f"{graph.label(target)} is neither a defect gene nor a pathway "
                                   f"partner of one, and a co-mention step carries no direction.")


def annotate_candidate(graph: Graph, candidate: Candidate, evidence: DrugEvidence,
                       eligibility: str, options: RankOptions,
                       defect: DefectContext) -> Candidate:
    """Add the terminal relation, the clinical stage, the mechanism and the direction verdict.

    Args:
        graph: The graph.
        candidate: The scored candidate.
        evidence: The drug records.
        eligibility: The eligibility string of the candidate.
        options: The run parameters.
        defect: The declared defect.

    Returns:
        The same candidate, with the explanation fields filled.
    """
    terminal = graph.edges[candidate.best_path.edge_indices[-1]]
    previous = candidate.best_path.node_ids[-2]
    candidate.terminal_relation = terminal.relation
    candidate.terminal_class = relation_class(terminal.relation, options.extra_pharmacological)
    candidate.terminal_target = previous
    candidate.clinical_stage = evidence.stage.get(candidate.node_id, "")
    candidate.mechanism = evidence.mechanism(candidate.node_id, previous)
    candidate.eligibility = eligibility
    recorded = evidence.action(candidate.node_id, previous)
    candidate.action_type = (terminal.relation if candidate.terminal_class == PHARMACOLOGICAL
                             else recorded)
    candidate.action_direction = action_direction_of(candidate.action_type)
    candidate.direction_verdict, candidate.direction_reason = direction_verdict(graph, candidate,
                                                                               defect)
    return candidate


def rank_candidates(graph: Graph, seeds: list[str], options: RankOptions | None = None,
                    evidence: DrugEvidence | None = None) -> list[Candidate]:
    """Rank every reachable chemical of the graph.

    Args:
        graph: The graph.
        seeds: The seed node identifiers.
        options: The run parameters. The default applies the terminal constraint, the eligibility
            constraint and the disease hop constraint, and it declares no defect.
        evidence: The drug records. An absent record set leaves eligibility to the indication edges.

    Returns:
        The candidates, best first. The direction tier sorts before the score. A tie breaks on the
        shorter path and then on the identifier, so the order is deterministic.
    """
    options = options or RankOptions()
    evidence = evidence or DrugEvidence()
    defect = build_defect_context(graph, seeds, options)
    distance, _ = shortest_distances(graph, seeds, options)
    indications = indication_chemicals(graph) if options.require_drug_evidence else set()
    seed_set = set(seeds)
    prefix_cache: dict[str, list[GraphPath]] = {}
    candidates: list[Candidate] = []
    for node_id in graph.nodes:
        if node_id in seed_set or graph.kind(node_id) != options.candidate_kind:
            continue
        eligibility = eligibility_of(node_id, evidence, indications)
        if options.require_drug_evidence and not eligibility:
            continue
        paths = candidate_paths(graph, node_id, distance, options, prefix_cache)
        if paths:
            candidates.append(annotate_candidate(graph, score_candidate(graph, node_id, paths),
                                                 evidence, eligibility, options, defect))
    return sorted(candidates, key=candidate_sort_key)


def path_string(graph: Graph, path: GraphPath) -> str:
    """Render a path as an arrow joined string of node labels.

    Args:
        graph: The graph.
        path: The path.

    Returns:
        A readable path, for example `MVA -> BUB1B -> AURKA -> alisertib`.
    """
    return " -> ".join(graph.label(node_id) for node_id in path.node_ids)


def relation_string(graph: Graph, path: GraphPath) -> str:
    """Render the relation of every step of a path.

    Args:
        graph: The graph.
        path: The path.

    Returns:
        The relation names joined by an arrow.
    """
    return " -> ".join(graph.edges[index].relation for index in path.edge_indices)


def kind_string(graph: Graph, path: GraphPath) -> str:
    """Render the node kind of every step of a path.

    Args:
        graph: The graph.
        path: The path.

    Returns:
        The kinds joined by an arrow.
    """
    return " -> ".join(graph.kind(node_id) for node_id in path.node_ids)


def intermediate_string(graph: Graph, path: GraphPath) -> str:
    """Name the nodes between the seed and the candidate.

    The reader rejects a candidate on its mechanism, so the intermediate nodes carry their kind.

    Args:
        graph: The graph.
        path: The path.

    Returns:
        The intermediate nodes, or "direct" when the path has one edge.
    """
    middle = path.node_ids[1:-1]
    if not middle:
        return "direct"
    return "; ".join(f"{graph.label(node_id)} ({graph.kind(node_id)})" for node_id in middle)


def candidate_row(graph: Graph, candidate: Candidate, rank: int) -> dict[str, str]:
    """Build the output row of one candidate.

    Args:
        graph: The graph.
        candidate: The scored candidate.
        rank: The rank, starting at 1.

    Returns:
        The row, keyed by the output columns.
    """
    return {"rank": str(rank), "chembl_id": candidate.node_id, "label": candidate.label,
            "score": f"{candidate.score:.4f}", "path_length": str(candidate.path_length),
            "path": path_string(graph, candidate.best_path),
            "relations": relation_string(graph, candidate.best_path),
            "direction_verdict": candidate.direction_verdict,
            "direction_reason": candidate.direction_reason,
            "action_type": candidate.action_type, "action_direction": candidate.action_direction,
            "terminal_relation": candidate.terminal_relation,
            "terminal_class": candidate.terminal_class,
            "terminal_target": candidate.terminal_target,
            "clinical_stage": candidate.clinical_stage, "mechanism": candidate.mechanism,
            "eligibility": candidate.eligibility,
            "pmids": ",".join(candidate.pmids), "n_pmids": str(len(candidate.pmids)),
            "pmid_coverage": f"{path_pmid_coverage(graph, candidate.best_path):.2f}",
            "n_shortest_paths": str(candidate.n_paths),
            "min_evidence_count": str(candidate.min_evidence),
            "path_kinds": kind_string(graph, candidate.best_path),
            "intermediates": intermediate_string(graph, candidate.best_path),
            "seed": candidate.seed_id}


def timestamp() -> str:
    """Return the current UTC time for an output file name.

    Returns:
        A snake_case timestamp, for example `2026_09_15_142530`.
    """
    return datetime.now(timezone.utc).strftime("%Y_%m_%d_%H%M%S")


def select_rows(candidates: list[Candidate], top: int) -> tuple[list[Candidate], list[int]]:
    """Choose the rows to write and their true ranks.

    A `--top` cut must not hide the rejected class. Every candidate that reinforces the declared
    defect therefore stays in the file with its true rank, even when it sits far below the cut. A
    reviewer must see what the ranker rejected and why.

    Args:
        candidates: The ranked candidates.
        top: The number of leading rows to keep. 0 keeps every row.

    Returns:
        The candidates to write and the rank of each one.
    """
    if top <= 0:
        return candidates, list(range(1, len(candidates) + 1))
    chosen = [(rank, candidate) for rank, candidate in enumerate(candidates, start=1)
              if rank <= top or candidate.direction_verdict.startswith("reinforces")]
    return [candidate for _, candidate in chosen], [rank for rank, _ in chosen]


def write_ranking(graph: Graph, candidates: list[Candidate], out_dir: Path,
                  stamp: str | None = None, ranks: list[int] | None = None) -> Path:
    """Write the ranking table.

    Args:
        graph: The graph.
        candidates: The ranked candidates.
        out_dir: The output directory.
        stamp: The timestamp of the file name. The current time is used when it is absent.
        ranks: The rank of each candidate. The position in the list is used when it is absent.

    Returns:
        The path of the written file.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"candidates_{stamp or timestamp()}.tsv"
    numbers = ranks if ranks is not None else list(range(1, len(candidates) + 1))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(OUTPUT_COLUMNS), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for rank, candidate in zip(numbers, candidates):
            writer.writerow(candidate_row(graph, candidate, rank))
    return path


def resolve_seeds(graph: Graph, seeds: list[str]) -> list[str]:
    """Check that every seed is a node of the graph.

    Args:
        graph: The graph.
        seeds: The seed identifiers given on the command line.

    Returns:
        The seeds in the order given, without a duplicate.

    Raises:
        ValueError: When a seed is not a node of the graph.
    """
    missing = [seed for seed in seeds if seed not in graph.nodes]
    if missing:
        raise ValueError(f"The graph holds no node for: {', '.join(missing)}.")
    ordered: list[str] = []
    for seed in seeds:
        if seed not in ordered:
            ordered.append(seed)
    return ordered


def parse_terminal_classes(text: str) -> tuple[str, ...]:
    """Read the `--terminal-classes` value.

    Args:
        text: A comma separated list of class names, or the word `all`.

    Returns:
        The allowed classes.

    Raises:
        ValueError: When a name is not a relation class.
    """
    if text.strip().lower() == "all":
        return RELATION_CLASSES
    names = tuple(part.strip().lower() for part in text.split(",") if part.strip())
    unknown = [name for name in names if name not in RELATION_CLASSES]
    if unknown:
        raise ValueError(f"Unknown relation class: {', '.join(unknown)}. "
                         f"Use one of {', '.join(RELATION_CLASSES)} or the word all.")
    return names


def add_rank_arguments(parser: argparse.ArgumentParser) -> None:
    """Add the arguments that set the ranking behaviour.

    The controls share these arguments, so a control runs the ranker that the report uses.

    Args:
        parser: The parser or the subcommand parser.
    """
    parser.add_argument("--nodes", type=Path, default=DEFAULT_NODES, help="The nodes.tsv path.")
    parser.add_argument("--edges", type=Path, default=DEFAULT_EDGES, help="The edges.tsv path.")
    parser.add_argument("--drug-evidence", type=Path, default=DEFAULT_DRUG_EVIDENCE,
                        help="The drug_evidence.tsv path.")
    parser.add_argument("--max-edges", type=int, default=3, help="The largest path length.")
    parser.add_argument("--hub-degree-cap", type=int, default=0,
                        help="Never expand a node of higher degree. 0 turns the cap off.")
    parser.add_argument("--max-paths", type=int, default=MAX_PATHS_PER_CANDIDATE,
                        help="The largest number of shortest paths kept per candidate.")
    parser.add_argument("--terminal-classes", default=",".join(DEFAULT_TERMINAL_CLASSES),
                        help="The relation classes that the last edge into a candidate may carry. "
                             "Use `all` to turn the terminal constraint off.")
    parser.add_argument("--extra-pharmacological", default="",
                        help="Relation names that count as pharmacological on top of the declared "
                             "table.")
    parser.add_argument("--require-drug-evidence", action=argparse.BooleanOptionalAction,
                        default=True,
                        help="Keep a candidate only when it carries a drug record or an indication "
                             "edge.")
    parser.add_argument("--defect", choices=sorted(DEFECT_DIRECTION), default="none",
                        help="The direction of the defect of the seed. A candidate whose action "
                             "pushes the same way sinks below every other candidate.")
    parser.add_argument("--defect-gene", action="append", default=[],
                        help="A gene that carries the defect. Repeat the option. A gene seed counts "
                             "on its own.")
    parser.add_argument("--gene-layer", choices=GENE_LAYERS, default=GENE_LAYER_COMENTION,
                        help="The gene to gene layer. comention is the default graph. functional "
                             "replaces the co-mention gene edges with the STRING layer. union keeps "
                             "both.")
    parser.add_argument("--functional-dir", type=Path, default=DEFAULT_FUNCTIONAL_DIR,
                        help="The directory of the functional layer tables.")
    parser.add_argument("--allow-disease-hops", action=argparse.BooleanOptionalAction, default=False,
                        help="Let a path pass through a second disease node. The default forbids "
                             "it, because the drugs of a co-mentioned disease are not evidence "
                             "about the seed disease.")


def options_from_args(args: argparse.Namespace) -> RankOptions:
    """Build the run parameters from the command line.

    Args:
        args: The parsed arguments.

    Returns:
        The run parameters.

    Raises:
        ValueError: When a relation class name is unknown.
    """
    extra = tuple(part.strip() for part in args.extra_pharmacological.split(",") if part.strip())
    return RankOptions(max_edges=args.max_edges, hub_degree_cap=args.hub_degree_cap,
                       max_paths=args.max_paths,
                       terminal_classes=parse_terminal_classes(args.terminal_classes),
                       extra_pharmacological=extra,
                       require_drug_evidence=args.require_drug_evidence,
                       defect=args.defect, defect_genes=tuple(args.defect_gene),
                       blocked_intermediate_kinds=() if args.allow_disease_hops else ("disease",))


def build_parser() -> argparse.ArgumentParser:
    """Build the command line parser.

    Returns:
        The parser.
    """
    parser = argparse.ArgumentParser(description="Rank chemicals as repurposing candidates.")
    add_rank_arguments(parser)
    parser.add_argument("--seed", action="append", default=[],
                        help="A seed node identifier. Repeat the option for a seed set.")
    parser.add_argument("--seeds-file", type=Path,
                        help="A file with one seed node identifier per line.")
    parser.add_argument("--top", type=int, default=0,
                        help="Write this many rows, plus every candidate that reinforces the "
                             "declared defect. 0 writes every candidate.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                        help="The output directory.")
    return parser


def read_seed_file(path: Path) -> list[str]:
    """Read seed identifiers from a file.

    Args:
        path: The seed file. One identifier per line, `#` starts a comment.

    Returns:
        The seed identifiers.
    """
    seeds = []
    for line in path.read_text(encoding="utf-8").splitlines():
        text = line.split("#", 1)[0].strip()
        if text:
            seeds.append(text)
    return seeds


def main(argv: list[str] | None = None) -> int:
    """Run the ranker from the command line.

    Args:
        argv: The arguments. `sys.argv[1:]` is used when it is absent.

    Returns:
        0 on success, 2 on a usage or input error.
    """
    args = build_parser().parse_args(argv)
    seeds = list(args.seed) + (read_seed_file(args.seeds_file) if args.seeds_file else [])
    if not seeds:
        print("Give at least one --seed or a --seeds-file.", file=sys.stderr)
        return 2
    try:
        graph = load_layered_graph(args.nodes, args.edges, args.gene_layer, args.functional_dir)
        seeds = resolve_seeds(graph, seeds)
        options = options_from_args(args)
        evidence = load_layered_drug_evidence(args.drug_evidence, args.gene_layer,
                                              args.functional_dir)
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    candidates = rank_candidates(graph, seeds, options, evidence)
    kept, ranks = select_rows(candidates, args.top)
    path = write_ranking(graph, kept, args.out_dir, ranks=ranks)
    unknown = unclassified_relations(graph, options.extra_pharmacological)
    print(f"graph: {graph.stats['nodes']} nodes, {graph.stats['edges']} edges, "
          f"{graph.stats['merged_away']} parallel rows merged")
    print(f"gene layer: {args.gene_layer}")
    print(f"drug records: {len(evidence.stage)} chemicals")
    print(f"seeds: {', '.join(seeds)}")
    print(f"terminal classes: {', '.join(options.terminal_classes)}; "
          f"require drug evidence: {options.require_drug_evidence}; "
          f"blocked intermediate kinds: {', '.join(options.blocked_intermediate_kinds) or 'none'}")
    if unknown:
        print(f"unclassified relations, treated as {DEFAULT_RELATION_CLASS}: {', '.join(unknown)}")
    defect = build_defect_context(graph, seeds, options)
    if defect.direction:
        print(f"defect: {args.defect}, {len(defect.genes)} defect genes, "
              f"{len(defect.partners)} pathway partner genes")
    counts: dict[str, int] = {}
    for candidate in candidates:
        counts[candidate.direction_verdict] = counts.get(candidate.direction_verdict, 0) + 1
    print("direction verdicts: " + ", ".join(f"{name} {counts[name]}" for name in sorted(counts)))
    print(f"candidates: {len(candidates)} within {options.max_edges} edges, {len(kept)} written")
    print(f"output: {path}")
    print("Read the path and the relations before the score. A control run decides whether a "
          "score means anything: scripts/track2_controls.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
