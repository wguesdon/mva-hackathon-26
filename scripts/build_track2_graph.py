"""Assemble the Track 2 knowledge graph from the cached PubTator3 and Open Targets responses.

The script reads only the files under `reports/track2/raw/`. It performs no network call, so the
graph rebuilds without a refetch. It writes three files that form the contract with the ranking
agent, plus two side tables.

* `nodes.tsv` holds `node_id`, `kind`, `label` and `source`.
* `edges.tsv` holds `source_id`, `target_id`, `relation`, `evidence_count`, `specificity_npmi`,
  `pmids` and `source`. The `specificity_npmi` column is the normalized pointwise mutual
  information of the entity pair. It corrects the degree bias of a raw co-mention count, and it
  sits beside `evidence_count` rather than replacing it, so both weightings can be compared.
* `manifest.json` holds the creation time, the counts, the endpoints and a SHA-256 of each table.
* `drug_evidence.tsv` holds the maximum clinical phase of each compound against its target gene,
  because the edge schema has no column for it.
* `indication_evidence.tsv` holds the maximum clinical phase of each disease to compound
  indication, for the same reason.

The build also filters the literature edges of the seed disease node. PubTator3 maps the words
"mitotic abnormalities" and other descriptions of a mitotic phenotype onto the seed concept, so most
of those edges state nothing about the syndrome. `scripts/track2_seed_filter.py` applies the rule of
`literature/track2_seed_disease_surface_forms_v1.yaml` and this build removes every seed edge that
the rule rejects. `--seed-mention-filter off` restores the earlier behaviour, so a receipt written
before 2026-09-15 stays reproducible.

Every other disease node of the graph carries the same class of mapping error.
`--mention-filter-scope graph` applies the rule of `scripts/track2_graph_mention_rule.py` to every
disease node and not to the seed alone, under ADR 0010. The default is `seed`, which is the
behaviour of every build before that record, so an earlier receipt reproduces without an option.

Usage:
    uv run python scripts/build_track2_graph.py
    uv run python scripts/build_track2_graph.py --seed-mention-filter off
    uv run python scripts/build_track2_graph.py --mention-filter-scope graph \\
        --graph-dir reports/track2/graph_all_disease_nodes
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402
import track2_graph_mention_rule as graph_rule  # noqa: E402
import track2_seed_filter as seed_filter  # noqa: E402

NODE_HEADER = ["node_id", "kind", "label", "source"]
EDGE_HEADER = ["source_id", "target_id", "relation", "evidence_count", "specificity_npmi",
               "pmids", "source"]
DRUG_HEADER = ["drug_id", "drug_name", "drug_type", "target_gene_id", "action_type",
               "max_clinical_stage_on_target", "max_clinical_stage_overall", "mechanism_of_action"]
INDICATION_HEADER = ["disease_id", "disease_label", "drug_id", "drug_name", "max_phase_for_ind",
                     "efo_id", "efo_term", "reference_count"]

SOURCE_PUBTATOR = "pubtator3"
SOURCE_OPEN_TARGETS = "open_targets"
SOURCE_REACTOME = "open_targets_reactome"
SOURCE_CHEMBL = "chembl"
SOURCE_CONSEQUENCE = "pubtator3_consequence"
PATHWAY_RELATION = "participates_in"
INDICATION_RELATION = "indication"
MAX_PMIDS_PER_EDGE = 30
LITERATURE_SOURCES = (SOURCE_PUBTATOR, SOURCE_CONSEQUENCE)
SEED_NODE = "DISEASE:MESH:C536987"
SEED_FILTER_STRICT = "strict"
SEED_FILTER_OFF = "off"
SCOPE_SEED = "seed"
SCOPE_GRAPH = "graph"


def entity_identity(accession: str, entities: dict) -> tuple[str, str]:
    """Return the node identifier and the label of one PubTator3 accession.

    Args:
        accession: The accession, for example `@CHEMICAL_Paclitaxel`.
        entities: The resolved entity cache, keyed by accession.

    Returns:
        The node identifier and the label. The identifier is empty for an unsupported biotype.
    """
    entry = entities.get(accession, {})
    node_id = common.node_id_for_accession(accession, entry.get("identifier", ""),
                                           entry.get("database", ""))
    label = entry.get("name") or common.accession_label(accession)
    return node_id, label


def build_literature_nodes(accessions: set[str], entities: dict) -> dict[str, dict]:
    """Build the gene, disease and chemical nodes of the literature layer.

    Args:
        accessions: The node accessions that the selection returned.
        entities: The resolved entity cache, keyed by accession.

    Returns:
        The node rows, keyed by node identifier.
    """
    nodes: dict[str, dict] = {}
    for accession in sorted(accessions):
        kind = common.accession_kind(accession)
        if kind not in common.PREFIX_BY_KIND:
            continue
        node_id, label = entity_identity(accession, entities)
        if not node_id:
            continue
        nodes.setdefault(node_id, {"node_id": node_id, "kind": kind, "label": label,
                                   "source": SOURCE_PUBTATOR})
    return nodes


def pmids_for(record: dict, evidence: dict) -> list[str]:
    """Read the PMIDs that the search endpoint returned for one relation record.

    Args:
        record: The relation record.
        evidence: The evidence cache, keyed by `type|source|target`.

    Returns:
        The PMIDs in the order that the endpoint returned them, without a duplicate.
    """
    key = f"{record.get('type', '')}|{record.get('source', '')}|{record.get('target', '')}"
    items = (evidence.get(key) or {}).get("evidence", [])
    seen: list[str] = []
    for item in items:
        pmid = item.get("pmid", "")
        if pmid and pmid not in seen:
            seen.append(pmid)
    return seen[:MAX_PMIDS_PER_EDGE]


def specificity_for(record: dict, totals: dict[tuple[str, str], int], counts: dict[str, int],
                    corpus_size: int) -> str:
    """Compute the specificity weight of one relation record.

    The weight belongs to the entity pair, not to the relation type, so the co-mention count is
    the sum over every relation type between the two entities.

    Args:
        record: The relation record.
        totals: The summed co-mention count of each unordered pair.
        counts: The total publication count of each entity.
        corpus_size: The number of records in the corpus.

    Returns:
        The value rounded to four places as a string, or an empty string when it is undefined.
    """
    source, target = record.get("source", ""), record.get("target", "")
    pair_count = totals.get(tuple(sorted((source, target))), 0)
    value = common.npmi(pair_count, counts.get(source, 0), counts.get(target, 0), corpus_size)
    return "" if value is None else f"{value:.4f}"


def build_literature_edges(records: list[dict], nodes: dict[str, dict], entities: dict,
                           evidence: dict, totals: dict[tuple[str, str], int],
                           counts: dict[str, int], corpus_size: int,
                           source: str = SOURCE_PUBTATOR) -> list[dict]:
    """Build the edge rows of the literature layer.

    Args:
        records: The selected relation records.
        nodes: The node rows, keyed by node identifier.
        entities: The resolved entity cache, keyed by accession.
        evidence: The evidence cache, keyed by `type|source|target`.
        totals: The summed co-mention count of each unordered pair.
        counts: The total publication count of each entity.
        corpus_size: The number of records in the corpus.
        source: The value of the `source` column.

    Returns:
        The edge rows, with a PMID list where the search endpoint returned one.
    """
    edges: list[dict] = []
    seen: set[tuple[str, str, str]] = set()
    for record in records:
        source_id, _ = entity_identity(record.get("source", ""), entities)
        target_id, _ = entity_identity(record.get("target", ""), entities)
        relation = record.get("type", "")
        if not source_id or not target_id or source_id not in nodes or target_id not in nodes:
            continue
        key = common.edge_key(source_id, target_id, relation)
        if key in seen:
            continue
        seen.add(key)
        edges.append({"source_id": source_id, "target_id": target_id, "relation": relation,
                      "evidence_count": int(record.get("publications", 0) or 0),
                      "specificity_npmi": specificity_for(record, totals, counts, corpus_size),
                      "pmids": ",".join(pmids_for(record, evidence)),
                      "source": source})
    return edges


def build_drug_layer(targets: dict, nodes: dict[str, dict]) -> tuple[list[dict], list[dict]]:
    """Build the compound nodes, the compound edges and the clinical phase table.

    A compound edge points from the compound to the gene, because the action type describes what
    the compound does to the target.

    Args:
        targets: The Open Targets cache, keyed by gene symbol.
        nodes: The node rows to extend, keyed by node identifier.

    Returns:
        The compound edge rows and the clinical phase rows.
    """
    edges: list[dict] = []
    phases: list[dict] = []
    seen: set[tuple[str, str, str]] = set()
    for symbol in sorted(targets):
        gene_id = f"GENE:{symbol}"
        if gene_id not in nodes:
            continue
        for drug in targets[symbol].get("drugs", []):
            drug_id = drug.get("drug_id", "")
            if not drug_id:
                continue
            node_id = f"CHEMBL:{drug_id}"
            nodes.setdefault(node_id, {"node_id": node_id, "kind": "chemical",
                                       "label": drug.get("drug_name", drug_id),
                                       "source": SOURCE_OPEN_TARGETS})
            relation = drug.get("action_type", "") or "targets"
            key = common.edge_key(node_id, gene_id, relation)
            if key in seen:
                continue
            seen.add(key)
            edges.append({"source_id": node_id, "target_id": gene_id, "relation": relation,
                          "evidence_count": 1, "specificity_npmi": "", "pmids": "",
                          "source": SOURCE_OPEN_TARGETS})
            phases.append({"drug_id": node_id, "drug_name": drug.get("drug_name", ""),
                           "drug_type": drug.get("drug_type", ""), "target_gene_id": gene_id,
                           "action_type": relation,
                           "max_clinical_stage_on_target": drug.get("max_clinical_stage_on_target", ""),
                           "max_clinical_stage_overall": drug.get("max_clinical_stage_overall", ""),
                           "mechanism_of_action": drug.get("mechanism_of_action", "")})
    return edges, phases


def build_pathway_layer(targets: dict, nodes: dict[str, dict]) -> list[dict]:
    """Build the Reactome pathway nodes and the gene to pathway edges.

    Args:
        targets: The Open Targets cache, keyed by gene symbol.
        nodes: The node rows to extend, keyed by node identifier.

    Returns:
        The pathway edge rows.
    """
    edges: list[dict] = []
    seen: set[tuple[str, str, str]] = set()
    for symbol in sorted(targets):
        gene_id = f"GENE:{symbol}"
        if gene_id not in nodes:
            continue
        for pathway in targets[symbol].get("pathways", []):
            pathway_id = pathway.get("pathway_id", "")
            if not pathway_id:
                continue
            node_id = f"PATHWAY:{pathway_id}"
            nodes.setdefault(node_id, {"node_id": node_id, "kind": "pathway",
                                       "label": pathway.get("pathway", pathway_id),
                                       "source": SOURCE_REACTOME})
            key = common.edge_key(gene_id, node_id, PATHWAY_RELATION)
            if key in seen:
                continue
            seen.add(key)
            edges.append({"source_id": gene_id, "target_id": node_id,
                          "relation": PATHWAY_RELATION, "evidence_count": 1,
                          "specificity_npmi": "", "pmids": "", "source": SOURCE_REACTOME})
    return edges


def build_indication_layer(indications: list[dict], nodes: dict[str, dict]) -> tuple[list[dict],
                                                                                    list[dict]]:
    """Build the disease to compound indication edges and the clinical phase table.

    The edge runs from the disease to the compound, so it reads as the disease being treated by
    the compound. An indication is kept only when both nodes already exist, which means the
    compound also reaches the disease through its target gene and the held-out control has a
    case that stays connected after the direct edge is removed.

    The ChEMBL endpoint gives no PMID, so the PMID field stays empty. The evidence count holds
    the number of distinct clinical trial and regulatory references of the pair.

    Args:
        indications: The joined ChEMBL indication records.
        nodes: The node rows, keyed by node identifier.

    Returns:
        The indication edge rows and the indication phase rows.
    """
    edges: list[dict] = []
    phases: list[dict] = []
    seen: set[tuple[str, str, str]] = set()
    for entry in indications:
        drug_id = f"CHEMBL:{entry.get('molecule_chembl_id', '')}"
        disease_id = f"DISEASE:MESH:{entry.get('mesh_id', '')}"
        if drug_id not in nodes or disease_id not in nodes:
            continue
        key = common.edge_key(disease_id, drug_id, INDICATION_RELATION)
        if key in seen:
            continue
        seen.add(key)
        edges.append({"source_id": disease_id, "target_id": drug_id,
                      "relation": INDICATION_RELATION,
                      "evidence_count": int(entry.get("reference_count", 1) or 1),
                      "specificity_npmi": "", "pmids": "", "source": SOURCE_CHEMBL})
        phases.append({"disease_id": disease_id, "disease_label": nodes[disease_id]["label"],
                       "drug_id": drug_id, "drug_name": nodes[drug_id]["label"],
                       "max_phase_for_ind": entry.get("max_phase_for_ind", ""),
                       "efo_id": entry.get("efo_id", ""), "efo_term": entry.get("efo_term", ""),
                       "reference_count": entry.get("reference_count", "")})
    return edges, phases


def gene_route_cases(edges: list[dict], indication_edges: list[dict]) -> dict[str, int]:
    """Count the indication cases that survive the removal of their direct edge.

    The held-out control removes a disease to compound edge and asks whether the ranker finds the
    compound again. The case only teaches something when the compound stays reachable. This
    function counts the cases that keep a two step route, where the disease shares a gene with the
    compound through a literature edge and an Open Targets mechanism edge.

    Args:
        edges: Every edge of the graph.
        indication_edges: The indication edges.

    Returns:
        The number of cases with a gene route, the number without, and the number of diseases
        that keep at least one such case.
    """
    genes_of_disease: dict[str, set[str]] = {}
    genes_of_drug: dict[str, set[str]] = {}
    for edge in edges:
        source_id, target_id = edge["source_id"], edge["target_id"]
        if edge["source"] == SOURCE_PUBTATOR:
            if source_id.startswith("DISEASE:") and target_id.startswith("GENE:"):
                genes_of_disease.setdefault(source_id, set()).add(target_id)
            elif target_id.startswith("DISEASE:") and source_id.startswith("GENE:"):
                genes_of_disease.setdefault(target_id, set()).add(source_id)
        elif edge["source"] == SOURCE_OPEN_TARGETS and source_id.startswith("CHEMBL:"):
            genes_of_drug.setdefault(source_id, set()).add(target_id)
    with_route = 0
    diseases: set[str] = set()
    for edge in indication_edges:
        shared = genes_of_disease.get(edge["source_id"], set()) & \
            genes_of_drug.get(edge["target_id"], set())
        if shared:
            with_route += 1
            diseases.add(edge["source_id"])
    return {"cases_with_a_gene_route": with_route,
            "cases_without_a_gene_route": len(indication_edges) - with_route,
            "diseases_with_a_recoverable_case": len(diseases)}


def seed_edge_count(edges: list[dict], seed: str) -> int:
    """Count the literature edges that touch the seed node.

    Args:
        edges: The edge rows.
        seed: The seed node identifier.

    Returns:
        The number of distinct source, target and relation keys of the seed.
    """
    return len({(row["source_id"], row["target_id"], row["relation"]) for row in edges
                if seed in (row["source_id"], row["target_id"])
                and row["source"] in LITERATURE_SOURCES})


def apply_seed_filter(edges: list[dict], seed: str, entities: dict, rule_path: Path,
                      evidence_path: Path, verdicts_path: Path) -> tuple[list[dict], dict]:
    """Remove every literature edge of the seed whose sentences never name the syndrome.

    The rule lives in the list file at `rule_path` and the decision comes from
    `scripts/track2_seed_filter.py`. This function only removes the edge rows that the decision
    rejects, and it writes the verdict table so that a reader checks each edge.

    Args:
        edges: The edge rows of the whole graph.
        seed: The seed node identifier.
        entities: The resolved entity cache, keyed by accession.
        rule_path: The surface form list file.
        evidence_path: The cached PubTator3 search results of the seed edges.
        verdicts_path: The verdict table to write.

    Returns:
        The kept edge rows and the record of the filter for the manifest.

    Raises:
        FileNotFoundError: When the cached search results are missing.
    """
    cache = common.read_json(evidence_path) or {}
    if not cache.get("pairs"):
        raise FileNotFoundError(f"no cached seed evidence at {evidence_path}; run "
                                "scripts/track2_edge_provenance.py audit-seed first")
    rule = seed_filter.load_rule(rule_path)
    verdicts = seed_filter.pair_verdicts(cache["pairs"], rule)
    removed_keys = seed_filter.dropped_edge_keys(verdicts, entities)
    before = seed_edge_count(edges, seed)
    kept = [row for row in edges
            if not (row["source"] in LITERATURE_SOURCES
                    and (row["source_id"], row["target_id"], row["relation"]) in removed_keys)]
    after = seed_edge_count(kept, seed)
    rows = [{"relation": record["relation"], "source_accession": record["source_accession"],
             "target_accession": record["target_accession"], "sentences": record["sentences"],
             "accepted_sentences": record["accepted_sentences"], "verdict": record["verdict"],
             "accepted_forms": " | ".join(record["accepted_forms"]),
             "rejected_forms": " | ".join(record["rejected_forms"])}
            for record in sorted(verdicts.values(),
                                 key=lambda r: (r["verdict"], -r["accepted_sentences"],
                                                r["relation"], r["source_accession"],
                                                r["target_accession"]))]
    written = common.write_tsv(verdicts_path, seed_filter.VERDICT_HEADER, rows)
    record = seed_filter.summary(verdicts, rule, seed_filter.listed_forms(rule_path))
    record.update({
        "mode": SEED_FILTER_STRICT,
        "evidence_file": seed_filter.repo_path(evidence_path),
        "verdict_table": seed_filter.repo_path(written),
        "seed_literature_edges_before": before,
        "seed_literature_edges_after": after,
        "seed_literature_edges_removed": before - after,
        "edge_rows_removed": len(edges) - len(kept),
    })
    return kept, record


def count_by(rows: list[dict], column: str) -> dict[str, int]:
    """Count the rows that share a value in one column.

    Args:
        rows: The node rows or the edge rows.
        column: The column to group by.

    Returns:
        The count of each value, in descending order of count.
    """
    counter = Counter(str(row.get(column, "")) for row in rows)
    return dict(sorted(counter.items(), key=lambda item: (-item[1], item[0])))


def main() -> int:
    """Read the raw cache, assemble the graph and write the three contract files.

    Returns:
        0 on success, 1 when a required cache file is missing.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--raw-dir", type=Path, default=common.RAW_DIR,
                        help="The raw cache directory.")
    parser.add_argument("--graph-dir", type=Path, default=common.GRAPH_DIR,
                        help="The output directory of the graph.")
    parser.add_argument("--tier2-chemical-min-publications", type=int, default=20,
                        help="The smallest publication count of a second tier chemical.")
    parser.add_argument("--class-accessions", type=Path, default=None,
                        help="The ingest file of a named chemical class, under ADR 0016. Every "
                             "accession of that file enters the graph whatever its publication "
                             "count is. The default reads no file.")
    parser.add_argument("--seed-mention-filter", choices=(SEED_FILTER_STRICT, SEED_FILTER_OFF),
                        default=SEED_FILTER_STRICT,
                        help="Keep a literature edge of the seed only when a sentence names the "
                             "syndrome. `off` restores the behaviour of a build before 2026-09-15.")
    parser.add_argument("--mention-filter-scope", choices=(SCOPE_SEED, SCOPE_GRAPH),
                        default=SCOPE_SEED,
                        help="`seed` decides the seed node alone, which is the behaviour of every "
                             "build before ADR 0010. `graph` decides every disease node.")
    parser.add_argument("--mention-rules-dir", type=Path, default=graph_rule.DEFAULT_RULES_DIR,
                        help="The directory that holds one rule file per disease node.")
    parser.add_argument("--mention-cache-dir", type=Path, default=graph_rule.DEFAULT_CACHE_DIR,
                        help="The sentence cache directory of the graph-wide rule.")
    parser.add_argument("--mention-source-graph-dir", type=Path,
                        default=graph_rule.DEFAULT_GRAPH_DIR,
                        help="The unfiltered graph directory whose node table the graph-wide rule "
                             "reads.")
    parser.add_argument("--seed-node", default=SEED_NODE, help="The seed disease node.")
    parser.add_argument("--seed-rule", type=Path, default=seed_filter.DEFAULT_RULE,
                        help="The surface form list file of the seed.")
    parser.add_argument("--seed-evidence", type=Path, default=seed_filter.DEFAULT_EVIDENCE,
                        help="The cached PubTator3 search results of the seed edges.")
    parser.add_argument("--seed-verdicts-out", type=Path, default=seed_filter.DEFAULT_VERDICTS,
                        help="The verdict table of the seed edges.")
    args = parser.parse_args()

    index = common.read_json(args.raw_dir / "pubtator_relations_index.json")
    if not index:
        print("FAIL: run scripts/fetch_pubtator_relations.py first.")
        return 1
    seeds = index["seed_genes"]
    neighbours = index["derived_neighbours"]
    records = common.load_relation_cache(args.raw_dir, seeds + neighbours, "pubtator_relations")
    class_ingest = common.read_json(args.class_accessions) or {} if args.class_accessions else {}
    class_accessions = frozenset(class_ingest.get("accessions", []))
    accessions, selected = common.select_graph(records, seeds, neighbours,
                                               args.tier2_chemical_min_publications,
                                               class_accessions)

    entities = (common.read_json(args.raw_dir / "pubtator_entities.json") or {}).get("entities", {})
    for accession, entry in (class_ingest.get("entities") or {}).items():
        entities.setdefault(accession, entry)
    evidence_cache = common.read_json(args.raw_dir / "pubtator_evidence.json") or {}
    evidence = evidence_cache.get("pairs", {})
    drug_cache = common.read_json(args.raw_dir / "open_targets_targets.json") or {}
    targets = drug_cache.get("targets", {})
    indication_cache = common.read_json(args.raw_dir / "chembl_indications.json") or {}
    indications = indication_cache.get("indications", [])
    count_cache = common.read_json(args.raw_dir / "pubtator_entity_counts.json") or {}
    counts = count_cache.get("counts", {})
    corpus_size = int(count_cache.get("corpus_size", 0) or 0)
    consequence_cache = common.read_json(args.raw_dir / "pubtator_consequences.json") or {}
    consequences = consequence_cache.get("consequences", [])
    consequence_records = {
        entry["accession"]: common.read_json(
            args.raw_dir / "pubtator_consequences"
            / f"{common.split_accession(entry['accession'])[1]}.json") or []
        for entry in consequences}
    extra_accessions, consequence_selected = common.select_consequence_layer(
        consequences, consequence_records, accessions,
        args.tier2_chemical_min_publications)
    accessions.update(extra_accessions)

    totals = common.pair_total(selected + consequence_selected)
    nodes = build_literature_nodes(accessions, entities)
    edges = build_literature_edges(selected, nodes, entities, evidence, totals, counts,
                                   corpus_size)
    consequence_edges = build_literature_edges(consequence_selected, nodes, entities, evidence,
                                               totals, counts, corpus_size, SOURCE_CONSEQUENCE)
    drug_edges, phases = build_drug_layer(targets, nodes)
    pathway_edges = build_pathway_layer(targets, nodes)
    indication_edges, indication_phases = build_indication_layer(indications, nodes)
    edges.extend(consequence_edges)
    edges.extend(drug_edges)
    edges.extend(pathway_edges)
    edges.extend(indication_edges)
    if args.seed_mention_filter == SEED_FILTER_STRICT \
            and args.mention_filter_scope == SCOPE_GRAPH:
        before = len(edges)
        edges, seed_record = graph_rule.apply_graph_filter(edges, args.mention_rules_dir,
                                                           args.mention_cache_dir, entities,
                                                           args.mention_source_graph_dir)
        seed_record["seed_literature_edges_after"] = seed_edge_count(edges, args.seed_node)
        seed_record["graph_edge_rows_before"] = before
        seed_record["graph_edge_rows_after"] = len(edges)
    elif args.seed_mention_filter == SEED_FILTER_STRICT:
        edges, seed_record = apply_seed_filter(edges, args.seed_node, entities, args.seed_rule,
                                               args.seed_evidence, args.seed_verdicts_out)
    else:
        seed_record = {"mode": SEED_FILTER_OFF,
                       "note": "Every literature edge of the seed stays in the graph, including "
                               "the edges whose sentences never name the syndrome.",
                       "seed_literature_edges_before": seed_edge_count(edges, args.seed_node),
                       "seed_literature_edges_after": seed_edge_count(edges, args.seed_node),
                       "seed_literature_edges_removed": 0, "edge_rows_removed": 0}
    routes = gene_route_cases(edges, indication_edges)

    node_rows = sorted(nodes.values(), key=lambda row: (row["kind"], row["node_id"]))
    edge_rows = sorted(edges, key=lambda row: (row["source"], -int(row["evidence_count"]),
                                               row["source_id"], row["target_id"],
                                               row["relation"]))
    nodes_path = common.write_tsv(args.graph_dir / "nodes.tsv", NODE_HEADER, node_rows)
    edges_path = common.write_tsv(args.graph_dir / "edges.tsv", EDGE_HEADER, edge_rows)
    drugs_path = common.write_tsv(args.graph_dir / "drug_evidence.tsv", DRUG_HEADER, phases)
    indications_path = common.write_tsv(args.graph_dir / "indication_evidence.tsv",
                                        INDICATION_HEADER, indication_phases)

    with_pmids = sum(1 for row in edge_rows if row["pmids"])
    weighted = sum(1 for row in edge_rows if row["specificity_npmi"])
    manifest = {
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "seed_genes": seeds,
        "derived_neighbour_genes": neighbours,
        "tier2_chemical_min_publications": args.tier2_chemical_min_publications,
        "api_endpoints": [
            "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/relations",
            "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/entity/autocomplete/",
            "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/",
            "https://api.platform.opentargets.org/api/v4/graphql",
            "https://www.ebi.ac.uk/chembl/api/data/drug_indication.json",
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi",
        ],
        "source_licences": {"pubtator3": "public domain (NCBI)", "open_targets": "CC0",
                            "open_targets_reactome": "CC0", "chembl": "CC BY-SA 3.0"},
        "node_count": len(node_rows),
        "edge_count": len(edge_rows),
        "nodes_by_kind": count_by(node_rows, "kind"),
        "nodes_by_source": count_by(node_rows, "source"),
        "edges_by_relation": count_by(edge_rows, "relation"),
        "edges_by_source": count_by(edge_rows, "source"),
        "edges_with_pmids": with_pmids,
        "specificity": {
            "measure": "normalized pointwise mutual information",
            "column": "specificity_npmi",
            "corpus_size": corpus_size,
            "corpus_date": count_cache.get("corpus_date", ""),
            "entities_with_a_count": len(counts),
            "edges_with_a_value": weighted,
        },
        "seed_mention_filter": seed_record,
        "consequence_layer": {
            "entities": [entry["accession"] for entry in consequences],
            "min_npmi": consequence_cache.get("min_npmi"),
            "min_seeds": consequence_cache.get("min_seeds"),
            "edges_before_the_seed_filter": len(consequence_edges),
            "edges": sum(1 for row in edge_rows if row["source"] == SOURCE_CONSEQUENCE),
            "nodes_added": len(extra_accessions),
        },
        "distinct_pmids": len({p for row in edge_rows for p in row["pmids"].split(",") if p}),
        "drug_evidence_rows": len(phases),
        "indication_edges": len(indication_edges),
        "diseases_with_an_indication": len({e["source_id"] for e in indication_edges}),
        "compounds_with_an_indication": len({e["target_id"] for e in indication_edges}),
        "indication_control_eligibility": routes,
        "files": {
            "nodes.tsv": {"sha256": common.sha256_of(nodes_path),
                          "bytes": nodes_path.stat().st_size},
            "edges.tsv": {"sha256": common.sha256_of(edges_path),
                          "bytes": edges_path.stat().st_size},
            "drug_evidence.tsv": {"sha256": common.sha256_of(drugs_path),
                                  "bytes": drugs_path.stat().st_size},
            "indication_evidence.tsv": {"sha256": common.sha256_of(indications_path),
                                        "bytes": indications_path.stat().st_size},
        },
    }
    manifest_path = common.write_json(args.graph_dir / "manifest.json", manifest)
    print(f"  nodes {len(node_rows)} {manifest['nodes_by_kind']}")
    print(f"  edges {len(edge_rows)} {manifest['edges_by_source']}")
    print(f"  {with_pmids} edges carry a PMID, {manifest['distinct_pmids']} distinct PMIDs")
    print(f"  {len(indication_edges)} indication edges over "
          f"{manifest['diseases_with_an_indication']} diseases, {routes}")
    print(f"  {manifest['edges_by_source'].get(SOURCE_CONSEQUENCE, 0)} consequence edges from "
          f"{len(consequences)} consequence entities")
    print(f"  {weighted} edges carry a specificity weight")
    if seed_record["mode"] == SCOPE_GRAPH:
        print(f"  mention filter graph: {seed_record['disease_nodes']} disease nodes, "
              f"{seed_record['disease_nodes_with_a_rule']} with a rule, "
              f"{seed_record['literature_edges_decided']} literature edges decided, "
              f"{seed_record['literature_edges_dropped']} dropped, "
              f"{seed_record['edge_rows_removed']} edge rows removed")
    else:
        print(f"  seed mention filter {seed_record['mode']}: "
              f"{seed_record['seed_literature_edges_before']} seed literature edges before, "
              f"{seed_record['seed_literature_edges_after']} after, "
              f"{seed_record['seed_literature_edges_removed']} removed")
    if seed_record.get("forms_not_in_the_list_file"):
        print(f"  WARNING: surface forms not in the list file: "
              f"{seed_record['forms_not_in_the_list_file']}")
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
