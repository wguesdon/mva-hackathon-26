"""Compare the co-mention, functional and union gene layers of the Track 2 graph.

The comparison answers queue item R4: does a functional interaction layer change the gene
neighbour set and the Track 2 ranking? Each arm loads the default graph with one gene to gene layer
of `scripts/track2_rank.py` and runs the same ranker, the same held-out control and the same
permutation control. Only the gene to gene edges differ between the arms.

Every parameter of the comparison sits in `PARAMETERS`. The values copy the parameters of the
existing receipts `t2_heldout_real_v1` and `t2_permutation_real_v1`, and they were committed before
the first held-out result of a functional arm existed. No value is tuned on a held-out result.

An optional reproduction step runs the default layer on a pinned copy of the graph that the existing
receipts used, and compares the case rows and the p-values with those receipts. The step shows that
the new option leaves the default behaviour unchanged.

Usage:
    uv run python scripts/track2_functional_compare.py
    uv run python scripts/track2_functional_compare.py --reproduce-graph-dir /path/to/graph_snapshot
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_controls as controls  # noqa: E402
import track2_rank as ranker  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RECEIPT_DIR = ROOT / "tracking" / "receipts"
RECEIPT_ID = "t2_functional_layer_v1"

SEED = "DISEASE:MESH:C536987"
DEFECT_GENES = ("GENE:BUB1B", "GENE:CEP57", "GENE:TRIP13")
NEIGHBOUR_GENE = "GENE:BUB1B"

# The parameters of the comparison, fixed before any held-out result of a functional arm existed.
PARAMETERS = {
    "arms": list(ranker.GENE_LAYERS),
    "reference_arm": ranker.GENE_LAYER_COMENTION,
    "seed": SEED,
    "report_rank_options": {"defect": "loss_of_function", "defect_genes": list(DEFECT_GENES)},
    "unconstrained_rank_options": {"terminal_classes": "all", "require_drug_evidence": False,
                                   "blocked_intermediate_kinds": []},
    # The watched list names the compounds that the Track 2 report discusses. It steers no
    # parameter of the comparison and it selects no candidate. It only reports a rank.
    "watched_chemicals": {"CHEMICAL:MESH:D009537": "Nicotinamide Mononucleotide",
                          "CHEMICAL:MESH:D009536": "Niacinamide",
                          "CHEMICAL:MESH:D009243": "NAD",
                          "CHEMBL:CHEMBL1233528": "Volasertib",
                          "CHEMBL:CHEMBL483158": "Alisertib"},
    "top_n_overlap": 10,
    "heldout": {"top_k": 10, "report_top_k": [1, 3, 10], "indication_relations": ["indication"],
                "permutations": 5, "random_seed": 20260915, "min_cases": 5, "min_recovery": 0.3,
                "rank_options": "RankOptions() default, as in t2_heldout_real_v1",
                "alpha": controls.DEFAULT_ALPHA,
                "regime_limits": {"max_genes": 10, "max_drugs": 9},
                "regimes": {"all": "every disease with a curated indication",
                            "drug_rich": "10 or more drugs with a clinical edge",
                            "drug_poor": "9 or fewer drugs with a clinical edge",
                            "rare_monogenic": "MeSH supplementary concept record"}},
    "permutation": {"permutations": 200, "top_k": 10, "alpha": 0.05, "random_seed": 20260915},
    "reference_receipts": {"heldout": "tracking/receipts/t2_heldout_real_v1.json",
                           "permutation": "tracking/receipts/t2_permutation_real_v1.json"},
}


def report_options() -> ranker.RankOptions:
    """Return the ranking options of the Track 2 report.

    Returns:
        The constrained options with the declared loss of function defect.
    """
    return ranker.RankOptions(defect="loss_of_function", defect_genes=DEFECT_GENES)


def unconstrained_options() -> ranker.RankOptions:
    """Return the ranking options of the unconstrained reagent ranking.

    Returns:
        The options with every terminal class, no eligibility constraint and no blocked kind.
    """
    return ranker.RankOptions(terminal_classes=ranker.RELATION_CLASSES,
                              require_drug_evidence=False, blocked_intermediate_kinds=())


def jaccard(first: set[str], second: set[str]) -> float:
    """Return the Jaccard index of two sets.

    Args:
        first: The first set.
        second: The second set.

    Returns:
        The size of the intersection over the size of the union, or 0 for two empty sets.
    """
    union = first | second
    return round(len(first & second) / len(union), 4) if union else 0.0


def gene_neighbours(graph: ranker.Graph, gene_id: str) -> set[str]:
    """Return the genes that a gene to gene edge joins to one gene.

    Args:
        graph: The graph.
        gene_id: The gene node identifier.

    Returns:
        The neighbour gene identifiers.
    """
    return {node_id for node_id, _ in graph.adjacency.get(gene_id, ())
            if graph.kind(node_id) == "gene" and node_id != gene_id}


def average_ranks(values: list[float]) -> list[float]:
    """Rank values from 1 upward and give tied values their average rank.

    Args:
        values: The values.

    Returns:
        The rank of each value, in input order.
    """
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start
        while end + 1 < len(order) and values[order[end + 1]] == values[order[start]]:
            end += 1
        for position in range(start, end + 1):
            ranks[order[position]] = (start + end) / 2 + 1
        start = end + 1
    return ranks


def spearman(first: list[float], second: list[float]) -> float | None:
    """Return the Spearman rank correlation of two paired lists.

    Args:
        first: The first values.
        second: The second values, paired with the first.

    Returns:
        The correlation, or None when fewer than 3 pairs exist or a list is constant.
    """
    if len(first) < 3:
        return None
    rank_a, rank_b = average_ranks(first), average_ranks(second)
    mean_a, mean_b = sum(rank_a) / len(rank_a), sum(rank_b) / len(rank_b)
    covariance = sum((a - mean_a) * (b - mean_b) for a, b in zip(rank_a, rank_b))
    spread_a = sum((a - mean_a) ** 2 for a in rank_a) ** 0.5
    spread_b = sum((b - mean_b) ** 2 for b in rank_b) ** 0.5
    if spread_a == 0 or spread_b == 0:
        return None
    return round(covariance / (spread_a * spread_b), 4)


def ranking_record(graph: ranker.Graph, candidates: list[ranker.Candidate], top_n: int) -> dict:
    """Summarize one ranking.

    Args:
        graph: The graph.
        candidates: The ranked candidates.
        top_n: The number of leading rows to keep.

    Returns:
        The candidate count, the rank of every candidate and the leading rows.
    """
    return {"candidates": len(candidates),
            "ranks": controls.rank_map(candidates),
            "top": [{"rank": rank, "node_id": candidate.node_id, "label": candidate.label,
                     "score": round(candidate.score, 4),
                     "direction_verdict": candidate.direction_verdict,
                     "path": ranker.path_string(graph, candidate.best_path),
                     "relations": ranker.relation_string(graph, candidate.best_path)}
                    for rank, candidate in enumerate(candidates[:top_n], start=1)]}


def compare_rankings(reference: dict, arm: dict, top_n: int) -> dict:
    """Compare the ranking of one arm with the ranking of the reference arm.

    Args:
        reference: The ranking record of the reference arm.
        arm: The ranking record of the arm.
        top_n: The size of the top list.

    Returns:
        The shared candidate count, the Spearman correlation over the shared candidates, and the
        overlap of the two top lists.
    """
    shared = sorted(set(reference["ranks"]) & set(arm["ranks"]))
    top_ref = {row["node_id"] for row in reference["top"][:top_n]}
    top_arm = {row["node_id"] for row in arm["top"][:top_n]}
    return {"shared_candidates": len(shared),
            "only_in_reference": len(set(reference["ranks"]) - set(arm["ranks"])),
            "only_in_arm": len(set(arm["ranks"]) - set(reference["ranks"])),
            "spearman_shared": spearman([reference["ranks"][c] for c in shared],
                                        [arm["ranks"][c] for c in shared]),
            "top_overlap": len(top_ref & top_arm), "top_n": top_n,
            "top_jaccard": jaccard(top_ref, top_arm)}


def watched_ranks(ranking: dict, watched: dict[str, str]) -> dict:
    """Report the rank of each watched chemical, or why it has none.

    Args:
        ranking: The ranking record.
        watched: The watched chemical identifiers and labels.

    Returns:
        The rank of each watched chemical, or None when the chemical is not a ranked candidate.
    """
    return {label: ranking["ranks"].get(node_id) for node_id, label in watched.items()}


def regime_sets(graph: ranker.Graph, limits: dict) -> dict[str, set[str]]:
    """Select the diseases of every regime of the comparison.

    Args:
        graph: The graph.
        limits: The regime thresholds.

    Returns:
        The disease identifiers of each regime.
    """
    every = controls.regime_diseases(graph, "all", limits)
    poor = controls.regime_diseases(graph, "drug_poor", limits)
    return {"all": every, "drug_rich": every - poor, "drug_poor": poor,
            "rare_monogenic": controls.regime_diseases(graph, "mesh_supplementary", limits)}


def rows_at_k(rows: list[dict], top_k: int) -> list[dict]:
    """Recompute the recovery flag of case rows for another top k.

    Args:
        rows: The case rows. Each carries the best rank of a true drug, 0 when none is ranked.
        top_k: The size of the top list.

    Returns:
        Copies of the rows with the recovery flag of that top k.
    """
    return [dict(row, recovered=bool(row["best_rank"] and row["best_rank"] <= top_k))
            for row in rows]


def regime_summaries(result: dict, regimes: dict[str, set[str]], settings: dict) -> dict:
    """Summarize the held-out result of one arm for every regime and every reported top k.

    Args:
        result: The record of `controls.heldout_control`.
        regimes: The disease identifiers of each regime.
        settings: The held-out parameters.

    Returns:
        For each top k and each regime, the real rate, the null rate and the verdict.
    """
    out: dict[str, dict] = {}
    for top_k in settings["report_top_k"]:
        real = rows_at_k(result["results"], top_k)
        null = rows_at_k(result["null_results"], top_k)
        block = {}
        for name, diseases in regimes.items():
            summary = controls.summarize_recovery(
                [row for row in real if row["disease"] in diseases],
                [row for row in null if row["disease"] in diseases],
                top_k, settings["min_cases"], settings["min_recovery"], settings["permutations"])
            block[name] = {"diseases": summary["diseases"], "recovered": summary["recovered"],
                           "recovery_rate": summary["recovery_rate"],
                           "null_rate": summary["null"]["recovery_rate"] if summary["null"] else None,
                           "null_rate_smoothed": summary["null_rate_smoothed"],
                           "chance_probability": summary["chance_probability"],
                           "above_null": (summary["recovery_rate"] is not None
                                          and summary["null"] is not None
                                          and summary["recovery_rate"]
                                          > summary["null"]["recovery_rate"]),
                           "verdict": summary["verdict"]}
        out[f"top_{top_k}"] = block
    return out


def drug_poor_verdict(arms: dict, reference_arm: str, alpha: float) -> tuple[str, dict]:
    """Say whether a functional arm lifts the drug-poor regime above its degree-matched null.

    A point estimate above the null is not enough. The null of this control has 5 permutations per
    disease, and the recovery rate of the drug-poor regime rests on a handful of cases, so the
    chance probability of the comparison must also reach alpha.

    Args:
        arms: The arm records.
        reference_arm: The co-mention arm, which is the thing to beat.
        alpha: The threshold of the chance probability.

    Returns:
        The verdict and the drug-poor block of every arm at the primary top k.
    """
    blocks = {name: record["heldout"]["summaries"]["top_10"]["drug_poor"]
              for name, record in arms.items()}
    moved = [name for name, block in blocks.items()
             if name != reference_arm and block["above_null"]
             and block["chance_probability"] <= alpha]
    verdict = ("FUNCTIONAL_LAYER_MOVES_DRUG_POOR_ABOVE_NULL" if moved
               else "FUNCTIONAL_LAYER_DOES_NOT_MOVE_DRUG_POOR_ABOVE_NULL")
    return verdict, blocks


def run_heldout(graph: ranker.Graph, evidence: ranker.DrugEvidence) -> dict:
    """Run the held-out control with the fixed parameters.

    Args:
        graph: The graph of one arm.
        evidence: The drug records of that arm.

    Returns:
        The record of `controls.heldout_control`.
    """
    settings = PARAMETERS["heldout"]
    return controls.heldout_control(
        graph, settings["top_k"], tuple(settings["indication_relations"]),
        permutations=settings["permutations"],
        generator=random.Random(settings["random_seed"]), options=ranker.RankOptions(),
        evidence=evidence, min_cases=settings["min_cases"],
        min_recovery=settings["min_recovery"], regime="drug_poor",
        limits=dict(settings["regime_limits"]))


def run_permutation(graph: ranker.Graph, evidence: ranker.DrugEvidence) -> dict:
    """Run the permutation control with the fixed parameters.

    Args:
        graph: The graph of one arm.
        evidence: The drug records of that arm.

    Returns:
        The record of `controls.permutation_control`.
    """
    settings = PARAMETERS["permutation"]
    return controls.permutation_control(graph, [SEED], settings["permutations"], settings["top_k"],
                                        settings["random_seed"], options=report_options(),
                                        evidence=evidence, alpha=settings["alpha"])


def compact_cases(rows: list[dict]) -> list[dict]:
    """Keep the fields of a case row that a later reader needs.

    Args:
        rows: The case rows.

    Returns:
        The compact rows.
    """
    return [{"disease": row["disease"], "best_rank": row["best_rank"],
             "candidates_after_removal": row["candidates_after_removal"],
             "true_drugs_reachable": row["true_drugs_reachable"], "recovered": row["recovered"]}
            for row in rows]


def reproduce(graph_dir: Path) -> dict:
    """Run the default layer on the pinned receipt graph and compare with the existing receipts.

    Args:
        graph_dir: The directory of the pinned graph copy.

    Returns:
        The digests of the pinned copy and whether each receipt reproduces.
    """
    nodes, edges = graph_dir / "nodes.tsv", graph_dir / "edges.tsv"
    graph = ranker.load_layered_graph(nodes, edges, ranker.GENE_LAYER_COMENTION)
    evidence = ranker.load_layered_drug_evidence(graph_dir / "drug_evidence.tsv",
                                                 ranker.GENE_LAYER_COMENTION)
    heldout_receipt = json.loads((ROOT / PARAMETERS["reference_receipts"]["heldout"]).read_text())
    permutation_receipt = json.loads(
        (ROOT / PARAMETERS["reference_receipts"]["permutation"]).read_text())
    heldout = run_heldout(graph, evidence)
    permutation = run_permutation(graph, evidence)
    keys = ("disease", "best_rank", "recovered")
    same_real = ([tuple(row[k] for k in keys) for row in heldout["results"]]
                 == [tuple(row[k] for k in keys) for row in heldout_receipt["result"]["results"]])
    same_null = ([tuple(row[k] for k in keys) for row in heldout["null_results"]]
                 == [tuple(row[k] for k in keys)
                     for row in heldout_receipt["result"]["null_results"]])
    same_p = ([(row["node_id"], row["p_value"]) for row in permutation["candidates"]]
              == [(row["node_id"], row["p_value"])
                  for row in permutation_receipt["result"]["candidates"]])
    return {"graph_dir_note": "A pinned copy of the graph that the reference receipts used.",
            "nodes_sha256_prefix": controls.file_digest(nodes),
            "edges_sha256_prefix": controls.file_digest(edges),
            "matches_receipt_edges_digest": (controls.file_digest(edges)
                                             == heldout_receipt["graph"]["edges_sha256_prefix"]),
            "heldout_real_rows_identical": same_real,
            "heldout_null_rows_identical": same_null,
            "heldout_full_rate": heldout["full"]["recovery_rate"],
            "heldout_drug_poor_rate": heldout["restricted"]["recovery_rate"],
            "permutation_p_values_identical": same_p,
            "permutation_verdict": permutation["verdict"]}


def run_arm(layer: str, args: argparse.Namespace) -> tuple[ranker.Graph, dict]:
    """Run every measure of the comparison on one arm.

    Args:
        layer: The gene layer of the arm.
        args: The parsed arguments.

    Returns:
        The graph of the arm and the record of the arm.
    """
    graph = ranker.load_layered_graph(args.nodes, args.edges, layer, args.functional_dir)
    evidence = ranker.load_layered_drug_evidence(args.drug_evidence, layer, args.functional_dir)
    top_n = PARAMETERS["top_n_overlap"]
    report = ranker.rank_candidates(graph, [SEED], report_options(), evidence)
    loose = ranker.rank_candidates(graph, [SEED], unconstrained_options(), evidence)
    print(f"[{layer}] graph {graph.stats['nodes']} nodes, {graph.stats['edges']} edges; "
          f"{len(report)} report candidates, {len(loose)} unconstrained candidates", flush=True)
    heldout = run_heldout(graph, evidence)
    regimes = regime_sets(graph, PARAMETERS["heldout"]["regime_limits"])
    summaries = regime_summaries(heldout, regimes, PARAMETERS["heldout"])
    print(f"[{layer}] held-out top 10: " + ", ".join(
        f"{name} {block['recovery_rate']} vs null {block['null_rate']}"
        for name, block in summaries["top_10"].items()), flush=True)
    permutation = run_permutation(graph, evidence)
    print(f"[{layer}] permutation: {permutation['verdict']}, survivors "
          f"{permutation['survivors']} of {len(permutation['candidates'])}", flush=True)
    record = {
        "graph": dict(graph.stats),
        "bub1b_gene_neighbours": sorted(gene_neighbours(graph, NEIGHBOUR_GENE)),
        "report_ranking": ranking_record(graph, report, top_n),
        "unconstrained_ranking": ranking_record(graph, loose, top_n),
        "heldout": {"regime_sizes_in_graph": {name: len(ids) for name, ids in regimes.items()},
                    "summaries": summaries,
                    "full_verdict": heldout["full"]["verdict"],
                    "drug_poor_verdict": heldout["restricted"]["verdict"],
                    "cases": compact_cases(heldout["results"]),
                    "null_cases": compact_cases(heldout["null_results"])},
        "permutation": {key: permutation[key] for key in
                        ("permutations", "alpha", "top_k", "candidates_total",
                         "degree_sequence_preserved", "survivors", "verdict")},
    }
    record["permutation"]["candidates"] = [
        {key: row[key] for key in ("rank", "node_id", "label", "p_value", "verdict",
                                   "direction_verdict", "path")}
        for row in permutation["candidates"]]
    return graph, record


def neighbour_sets(functional_dir: Path, graph_dir: Path) -> dict:
    """Compare the derived seed neighbour sets of the arms.

    Args:
        functional_dir: The directory of the functional layer tables.
        graph_dir: The directory of the default graph.

    Returns:
        The three neighbour sets with their overlap and Jaccard index against the co-mention set.
    """
    layer_manifest = json.loads((functional_dir / "manifest.json").read_text())
    graph_manifest = json.loads((graph_dir / "manifest.json").read_text())
    comention = set(graph_manifest["derived_neighbour_genes"])
    functional = {row["symbol"] for row in layer_manifest["neighbours"]["functional"]}
    union = comention | functional
    return {"comention": {"size": len(comention), "genes": sorted(comention)},
            "functional": {"size": len(functional), "genes": sorted(functional),
                           "overlap_with_comention": len(functional & comention),
                           "jaccard_with_comention": jaccard(functional, comention),
                           "shared": sorted(functional & comention),
                           "not_in_default_graph": layer_manifest["neighbours"][
                               "functional_not_in_default_graph"]},
            "union": {"size": len(union), "overlap_with_comention": len(union & comention),
                      "jaccard_with_comention": jaccard(union, comention)}}


def main() -> int:
    """Run the three arms, compare them and write the receipt.

    Returns:
        0 on success, 2 on an input error.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--nodes", type=Path, default=ranker.DEFAULT_NODES)
    parser.add_argument("--edges", type=Path, default=ranker.DEFAULT_EDGES)
    parser.add_argument("--drug-evidence", type=Path, default=ranker.DEFAULT_DRUG_EVIDENCE)
    parser.add_argument("--functional-dir", type=Path, default=ranker.DEFAULT_FUNCTIONAL_DIR)
    parser.add_argument("--reproduce-graph-dir", type=Path, default=None,
                        help="A pinned copy of the graph of the reference receipts.")
    parser.add_argument("--receipt-dir", type=Path, default=RECEIPT_DIR)
    parser.add_argument("--receipt-id", default=RECEIPT_ID)
    args = parser.parse_args()
    if not (args.functional_dir / "edges.tsv").is_file():
        print("ERROR: run scripts/build_track2_functional_layer.py first.", file=sys.stderr)
        return 2

    started = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    reproduction = reproduce(args.reproduce_graph_dir) if args.reproduce_graph_dir else None
    if reproduction:
        print(f"reproduction: {reproduction}", flush=True)
    arms: dict[str, dict] = {}
    graphs: dict[str, ranker.Graph] = {}
    for layer in ranker.GENE_LAYERS:
        graphs[layer], arms[layer] = run_arm(layer, args)

    reference = arms[PARAMETERS["reference_arm"]]
    top_n = PARAMETERS["top_n_overlap"]
    reference_neighbours = set(reference["bub1b_gene_neighbours"])
    comparison = {}
    for layer, record in arms.items():
        neighbours = set(record["bub1b_gene_neighbours"])
        comparison[layer] = {
            "bub1b_gene_neighbours": {"size": len(neighbours),
                                      "overlap_with_comention": len(neighbours
                                                                    & reference_neighbours),
                                      "jaccard_with_comention": jaccard(neighbours,
                                                                        reference_neighbours)},
            "report_ranking": compare_rankings(reference["report_ranking"],
                                               record["report_ranking"], top_n),
            "unconstrained_ranking": compare_rankings(reference["unconstrained_ranking"],
                                                      record["unconstrained_ranking"], top_n),
            "watched_report_ranks": watched_ranks(record["report_ranking"],
                                                  PARAMETERS["watched_chemicals"]),
            "watched_unconstrained_ranks": watched_ranks(record["unconstrained_ranking"],
                                                         PARAMETERS["watched_chemicals"]),
        }
    for record in arms.values():
        for key in ("report_ranking", "unconstrained_ranking"):
            record[key].pop("ranks")

    layer_manifest = json.loads((args.functional_dir / "manifest.json").read_text())
    verdict, poor = drug_poor_verdict(arms, PARAMETERS["reference_arm"],
                                      PARAMETERS["heldout"]["alpha"])
    record = {
        "receipt_id": args.receipt_id,
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "started_utc": started,
        "track": "t2_rank", "control": "functional_gene_layer_comparison",
        "queue_item": "R4",
        "data_class_of_run": "A", "data_class_of_receipt": "A",
        "receipt_scope": "Public graph identifiers, STRING identifiers and PMIDs only. No patient "
                         "data enters this run.",
        "scripts": ["scripts/build_track2_functional_layer.py",
                    "scripts/track2_functional_compare.py", "scripts/track2_rank.py",
                    "scripts/track2_controls.py"],
        "source": layer_manifest["source"],
        "layer_parameters": layer_manifest["parameters"],
        "mapping": layer_manifest["mapping"],
        "layer_counts": layer_manifest["counts"],
        "default_graph": {"nodes_sha256_prefix": controls.file_digest(args.nodes),
                          "edges_sha256_prefix": controls.file_digest(args.edges)},
        "functional_layer_files": layer_manifest["files"],
        "parameters": PARAMETERS,
        "reproduction_of_reference_receipts": reproduction,
        "neighbour_sets": neighbour_sets(args.functional_dir, args.nodes.parent),
        "comparison": comparison,
        "arms": arms,
        "drug_poor_at_top_10": poor,
        "verdict": verdict,
        "verdict_rule": "A functional arm moves the drug-poor regime only when its recovery rate "
                        "is above the degree-matched null and the chance probability of that "
                        "comparison is at or below alpha.",
    }
    path = controls.write_receipt(record, args.receipt_id, args.receipt_dir)
    for layer in arms:
        block = comparison[layer]
        print(f"{layer}: BUB1B neighbours {block['bub1b_gene_neighbours']}, report ranking "
              f"{block['report_ranking']}, watched {block['watched_report_ranks']} / "
              f"{block['watched_unconstrained_ranks']}")
    print(f"verdict: {record['verdict']}")
    print(f"receipt: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
