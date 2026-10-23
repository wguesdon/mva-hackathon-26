"""Run the two controls that decide whether a Track 2 ranking means anything.

A knowledge graph ranking fails in two known ways, and each control tests one of them.

Control a, the degree-preserving permutation. A drug with many edges wins a naive ranking whatever
the evidence says. The control shuffles the edge endpoints inside a band of node kind and node
degree, so every node keeps its exact degree and every relation keeps its endpoint kinds, then it
ranks again. A candidate that holds its real rank in the permuted graphs owes that rank to its
degree. The control writes an empirical p-value for every candidate of the real top k.

Control b, the held-out indication recovery. The graph can state the answer in a single edge. The
control takes each disease that already has an approved drug, removes every direct edge between that
disease and any chemical, ranks again from the disease, and asks whether the true drug returns
inside the top k. The recovery rate is the only number that says whether the method works at all.
A run with too few eligible diseases reports INSUFFICIENT_CASES and no rate.

The same evaluation on a degree-preserving permutation gives the null rate of control b. A uniform
draw from the candidate pool is the wrong null, because the drug of a permuted indication edge keeps
the degree of a real drug and a node of high degree is easier to reach. The receipt reports both.

The permutation has two scopes and the receipt names the scope of the run. `all_edges` rewires every
edge of the ranking graph, so no layer of the graph carries real structure and the null answers
whether the recovery beats a degree-matched random graph. `indication_edges` rewires the indication
layer only and leaves the gene, pathway and drug target layers in place, so the null answers whether
the disease to drug pairing carries information beyond those layers. Both scopes keep the exact
number of slots of every node inside every group, so every disease keeps its indication degree and
the eligible disease set of a permuted graph is the eligible disease set of the real graph.

The held-out control reports several probabilities for the same comparison. The permutation quantile
is `(1 + the number of permutations whose recovery count reaches the real count) / (permutations +
1)`. It assumes nothing and its smallest reachable value is `1 / (permutations + 1)`. The Poisson
binomial tail treats the smoothed null rate as a known per-case success probability. It reaches a
smaller value than the permutation quantile, and a value below the floor is a property of that
model and not a measurement of the permutation distribution. The permutations also measure the
variance of the recovery count, which is larger than that tail allows, so the tail is
anti-conservative. The corrected tail is a beta-binomial whose first two moments are the first two
moments of the permuted count, and the decision probability that the verdict reads is the larger of
the permutation quantile and that corrected tail. The receipt carries every value, the standard
error of the null rate across permutations, and both tails recomputed at the upper confidence limit
of the null rate.

A rate over every disease of the graph answers a question about common diseases with approved drugs,
because those diseases dominate the set. `--regime` restricts the same rows to one regime and reports
a second rate beside the first. Each rule of `REGIMES` reads the graph, and no rule holds a hand
written disease list. The receipt keeps the case rows of the real run and of the null, so a second
regime and a second top k need no second run.

Both controls rank with the parameters of the report, so a constraint of the report is a constraint
of the control. Both write a receipt JSON into tracking/receipts/.

Usage:
    uv run python scripts/track2_controls.py permutation --seed DISEASE:MESH:C536987 -n 200
    uv run python scripts/track2_controls.py heldout --top-k 10 -n 200 --regime drug_poor \
        --workers 12 --null-stream per_permutation
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import multiprocessing
import random
import re
import statistics
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_rank as ranker  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RECEIPT_DIR = ROOT / "tracking" / "receipts"

# The relation names that state an approved indication come from the declared class table of the
# ranker, so one place defines what a clinical edge is.
TREATMENT_RELATIONS = tuple(sorted(name for name, group in ranker.RELATION_CLASS.items()
                                   if group == ranker.CLINICAL))
DEFAULT_PERMUTATIONS = 200
DEFAULT_ALPHA = 0.05
DEFAULT_TOP_K = 10
DEFAULT_MIN_CASES = 5
DEFAULT_MIN_RECOVERY = 0.3
# The null rows enter the receipt so that another regime and another top k can be measured from the
# receipt alone, with no second run. A very large null is summarized instead.
MAX_NULL_ROWS_IN_RECEIPT = 5000
# The Poisson binomial tail of a uniform draw costs work in proportion to the square of the number
# of rows, and it answers a question about one graph. A set of permuted rows pools several graphs,
# so a set above this size reports no uniform tail. The run of 2026-09-18 at 200 permutations
# reached 23,000 null rows, and the tail of that set needed minutes of compute for a field that no
# verdict reads.
MAX_ROWS_FOR_UNIFORM_TAIL = 2000
# The scope of the permutation. `all_edges` rewires every edge of the ranking graph. The null then
# holds no real structure in any layer. `indication_edges` rewires the indication layer only and
# leaves the gene, pathway and drug target layers in place.
NULL_SCOPE_ALL_EDGES = "all_edges"
NULL_SCOPE_INDICATION_EDGES = "indication_edges"
NULL_SCOPES = (NULL_SCOPE_ALL_EDGES, NULL_SCOPE_INDICATION_EDGES)
# The random stream of the null. `sequential` draws every permutation from one generator, which is
# the stream of every run recorded before 2026-09-18 and which no worker can reproduce in parallel.
# `per_permutation` seeds permutation i from the run seed and i, so the result of a run does not
# depend on the number of workers.
NULL_STREAM_SEQUENTIAL = "sequential"
NULL_STREAM_PER_PERMUTATION = "per_permutation"
NULL_STREAMS = (NULL_STREAM_SEQUENTIAL, NULL_STREAM_PER_PERMUTATION)
# The values of k that the receipt summarizes beside the k of the run, so a table of several top
# lists needs no second run.
DEFAULT_REPORT_TOP_K = (1, 3, 10, 25)
DEFAULT_CONFIDENCE = 0.95
# The smallest dispersion factor that the correction of section 15 of docs/benchmarks/
# t2_null_model_v1.md applies. The measured error of the binomial tail has one direction, and the
# sampling standard deviation of the logarithm of the factor is sqrt(2 / (B - 1)), which is 0.100 at
# 200 permutations and 0.707 at 5. A factor below 1 is not separated from 1 at either count and a
# tail that reads it is smaller than the binomial tail, so the correction never moves that way.
MIN_DISPERSION_FACTOR = 1.0
# The state of one worker process. The pool initializer fills it once per process, so the graph
# crosses the process boundary once and not once per permutation.
WORKER_STATE: dict = {}


def degree_band(graph: ranker.Graph, node_id: str) -> tuple[str, int]:
    """Return the permutation band of one node.

    A band holds nodes of the same kind and of a similar degree. The band keeps a gene in the place
    of a gene and a hub in the place of a hub, so the permuted graph stays comparable with the real
    one.

    Args:
        graph: The graph.
        node_id: The node identifier.

    Returns:
        The node kind and the base two degree bucket.
    """
    degree = graph.degree(node_id)
    return graph.kind(node_id), (0 if degree <= 1 else int(math.log2(degree)))


def slot_groups(graph: ranker.Graph, relations: tuple[str, ...] | None = None
                ) -> dict[tuple[str, int, tuple[str, int]], list[tuple[int, int]]]:
    """Group every edge endpoint slot by relation, by position and by band.

    Args:
        graph: The graph.
        relations: The relation names in the scope of the permutation, in lower case. `None` puts
            every edge of the graph in scope.

    Returns:
        For each group key, the list of slots. A slot is an edge index and a position, 0 for the
        source endpoint and 1 for the target endpoint.
    """
    groups: dict[tuple[str, int, tuple[str, int]], list[tuple[int, int]]] = defaultdict(list)
    for index, edge in enumerate(graph.edges):
        if relations is not None and edge.relation.strip().lower() not in relations:
            continue
        for position, node_id in ((0, edge.source_id), (1, edge.target_id)):
            groups[(edge.relation, position, degree_band(graph, node_id))].append((index, position))
    return groups


def scope_relations(scope: str, relations: tuple[str, ...]) -> tuple[str, ...] | None:
    """Return the relation names that one permutation scope rewires.

    Args:
        scope: A value of `NULL_SCOPES`.
        relations: The relation names that state an approved indication.

    Returns:
        `None` for the whole graph, or the indication relation names.

    Raises:
        ValueError: When the scope is not a value of `NULL_SCOPES`.
    """
    if scope == NULL_SCOPE_ALL_EDGES:
        return None
    if scope == NULL_SCOPE_INDICATION_EDGES:
        return relations
    raise ValueError(f"unknown null scope {scope}.")


def count_self_loops(endpoints: dict[tuple[int, int], str]) -> int:
    """Count the edges whose two endpoints are the same node.

    Args:
        endpoints: The node of every slot.

    Returns:
        The number of self loops.
    """
    return sum(1 for (index, position) in endpoints
               if position == 0 and endpoints[(index, 1)] == endpoints[(index, 0)])


def swap_out_self_loops(endpoints: dict[tuple[int, int], str],
                        groups: dict[tuple[str, int, tuple[str, int]], list[tuple[int, int]]],
                        generator: random.Random, passes: int = 6) -> int:
    """Remove self loops by swapping two slots of one group.

    Args:
        endpoints: The node of every slot. The mapping is edited in place.
        groups: The slots of every group.
        generator: The random generator.
        passes: The number of repair passes.

    Returns:
        The number of self loops that remain after the last pass.
    """
    slot_group = {slot: key for key, slots in groups.items() for slot in slots}
    for _ in range(passes):
        loops = sorted({index for (index, position) in endpoints
                        if position == 0 and endpoints[(index, 0)] == endpoints[(index, 1)]})
        if not loops:
            return 0
        for index in loops:
            if endpoints[(index, 0)] != endpoints[(index, 1)]:
                continue
            partners = list(groups[slot_group[(index, 1)]])
            generator.shuffle(partners)
            for partner in partners:
                if partner == (index, 1):
                    continue
                other_index = partner[0]
                mine, theirs = endpoints[(index, 1)], endpoints[partner]
                if theirs == endpoints[(index, 0)]:
                    continue
                if mine == endpoints[(other_index, 1 - partner[1])]:
                    continue
                endpoints[(index, 1)], endpoints[partner] = theirs, mine
                break
    return count_self_loops(endpoints)


def permute_graph(graph: ranker.Graph, generator: random.Random,
                  relations: tuple[str, ...] | None = None) -> ranker.Graph:
    """Build a degree-preserving permutation of the graph.

    Every edge keeps its relation, its evidence count and its PMIDs. Only the endpoints move, and
    they move only inside a group of the same relation, the same endpoint position and the same
    band of node kind and node degree. A permutation of a group is a permutation of the nodes of
    that group, so every node keeps the exact number of endpoint slots it had in that group. The
    degree of a node is the sum of those counts, so the degree sequence is preserved exactly and so
    is the number of edges of each relation at each node.

    Args:
        graph: The real graph.
        generator: The random generator.
        relations: The relation names in the scope of the permutation, in lower case. `None`
            rewires every edge of the graph. A relation outside the scope keeps its endpoints.

    Returns:
        The permuted graph. Its statistics carry the number of self loops that remain and the
        number of edges that the permutation moved.
    """
    groups = slot_groups(graph, relations)
    endpoints: dict[tuple[int, int], str] = {}
    for key in sorted(groups):
        slots = groups[key]
        nodes = [graph.edges[index].source_id if position == 0 else graph.edges[index].target_id
                 for index, position in slots]
        generator.shuffle(nodes)
        for slot, node_id in zip(slots, nodes):
            endpoints[slot] = node_id
    for index, edge in enumerate(graph.edges):
        endpoints.setdefault((index, 0), edge.source_id)
        endpoints.setdefault((index, 1), edge.target_id)
    remaining = swap_out_self_loops(endpoints, groups, generator)
    permuted = [ranker.Edge(source_id=endpoints[(index, 0)], target_id=endpoints[(index, 1)],
                            relation=edge.relation, evidence_count=edge.evidence_count,
                            pmids=edge.pmids, source=edge.source)
                for index, edge in enumerate(graph.edges)]
    new_graph = ranker.build_graph(graph.nodes, permuted, merge=False)
    new_graph.stats["self_loops_after_repair"] = remaining
    new_graph.stats["edges_in_permutation_scope"] = len(
        {index for slots in groups.values() for index, _ in slots})
    return new_graph


def degree_sequence(graph: ranker.Graph) -> dict[str, int]:
    """Return the degree of every node.

    Args:
        graph: The graph.

    Returns:
        The degree of every node identifier.
    """
    return {node_id: graph.degree(node_id) for node_id in graph.nodes}


def rank_map(candidates: list[ranker.Candidate]) -> dict[str, int]:
    """Map every candidate identifier to its rank.

    Args:
        candidates: The ranked candidates.

    Returns:
        The rank of every candidate, starting at 1.
    """
    return {candidate.node_id: rank for rank, candidate in enumerate(candidates, start=1)}


def permutation_control(graph: ranker.Graph, seeds: list[str], permutations: int, top_k: int,
                        random_seed: int, options: ranker.RankOptions | None = None,
                        evidence: ranker.DrugEvidence | None = None,
                        alpha: float = DEFAULT_ALPHA) -> dict:
    """Run the degree-preserving permutation control.

    Args:
        graph: The real graph.
        seeds: The seed node identifiers.
        permutations: The number of permuted graphs.
        top_k: The number of real candidates to test.
        random_seed: The seed of the random generator, so the run repeats.
        options: The ranking parameters. The control ranks the permuted graph with the parameters
            of the report, so a constraint of the report is also a constraint of the control.
        evidence: The drug records.
        alpha: The p-value threshold of a candidate that is not a degree artifact.

    Returns:
        The result record that the receipt carries.
    """
    options = options or ranker.RankOptions()
    evidence = evidence or ranker.DrugEvidence()
    generator = random.Random(random_seed)
    real = ranker.rank_candidates(graph, seeds, options, evidence)
    tested = real[:top_k]
    real_ranks = rank_map(real)
    holds = {candidate.node_id: 0 for candidate in tested}
    score_sum = {candidate.node_id: 0.0 for candidate in tested}
    reached = {candidate.node_id: 0 for candidate in tested}
    degree_kept = True
    loops_left = 0
    real_degrees = degree_sequence(graph)
    for _ in range(permutations):
        permuted = permute_graph(graph, generator)
        loops_left = max(loops_left, permuted.stats.get("self_loops_after_repair", 0))
        if degree_sequence(permuted) != real_degrees:
            degree_kept = False
        ranked = ranker.rank_candidates(permuted, seeds, options, evidence)
        ranks = rank_map(ranked)
        scores = {candidate.node_id: candidate.score for candidate in ranked}
        for candidate in tested:
            position = ranks.get(candidate.node_id)
            if position is None:
                continue
            reached[candidate.node_id] += 1
            score_sum[candidate.node_id] += scores[candidate.node_id]
            if position <= real_ranks[candidate.node_id]:
                holds[candidate.node_id] += 1
    rows = []
    for candidate in tested:
        count = holds[candidate.node_id]
        p_value = (1 + count) / (permutations + 1)
        rows.append({"rank": real_ranks[candidate.node_id], "node_id": candidate.node_id,
                     "label": candidate.label, "score": round(candidate.score, 4),
                     "path": ranker.path_string(graph, candidate.best_path),
                     "relations": ranker.relation_string(graph, candidate.best_path),
                     "n_pmids": len(candidate.pmids),
                     "direction_verdict": candidate.direction_verdict,
                     "terminal_relation": candidate.terminal_relation,
                     "degree": graph.degree(candidate.node_id),
                     "permutations_at_or_above_real_rank": count,
                     "permutations_reached": reached[candidate.node_id],
                     "mean_permuted_score": round(
                         score_sum[candidate.node_id] / reached[candidate.node_id], 4)
                     if reached[candidate.node_id] else 0.0,
                     "p_value": round(p_value, 4),
                     "verdict": "evidence_supported" if p_value <= alpha else "degree_artifact"})
    survivors = [row for row in rows if row["verdict"] == "evidence_supported"]
    if permutations + 1 < 1 / alpha:
        verdict = "INSUFFICIENT_PERMUTATIONS"
    else:
        verdict = "PASS" if survivors else "FAIL"
    return {"permutations": permutations, "alpha": alpha, "top_k": top_k,
            "candidates_total": len(real), "degree_sequence_preserved": degree_kept,
            "largest_self_loop_count_after_repair": loops_left,
            "survivors": len(survivors), "verdict": verdict, "candidates": rows}


def is_treatment_edge(graph: ranker.Graph, edge: ranker.Edge, relations: tuple[str, ...]) -> bool:
    """Say whether one edge states an approved indication between a disease and a chemical.

    Args:
        graph: The graph.
        edge: The edge.
        relations: The relation names that state an indication.

    Returns:
        True when the edge joins a disease and a chemical with an indication relation.
    """
    kinds = {graph.kind(edge.source_id), graph.kind(edge.target_id)}
    if kinds != {"disease", "chemical"}:
        return False
    return edge.relation.strip().lower() in relations


def treatment_cases(graph: ranker.Graph, relations: tuple[str, ...]) -> dict[str, dict]:
    """Collect every disease of the graph that has an approved drug.

    Args:
        graph: The graph.
        relations: The relation names that state an indication.

    Returns:
        For each disease identifier, the true drugs and every direct disease to chemical edge index.
        Every direct edge is removed, not only the indication edge, because another relation name
        between the same pair leaks the same answer.
    """
    cases: dict[str, dict] = {}
    for index, edge in enumerate(graph.edges):
        kinds = {graph.kind(edge.source_id), graph.kind(edge.target_id)}
        if kinds != {"disease", "chemical"}:
            continue
        disease = (edge.source_id if graph.kind(edge.source_id) == "disease" else edge.target_id)
        chemical = (edge.target_id if disease == edge.source_id else edge.source_id)
        case = cases.setdefault(disease, {"drugs": set(), "edge_indices": set(),
                                          "treatment_edges": 0})
        case["edge_indices"].add(index)
        if is_treatment_edge(graph, edge, relations):
            case["drugs"].add(chemical)
            case["treatment_edges"] += 1
    return {disease: case for disease, case in cases.items() if case["drugs"]}


def graph_without(graph: ranker.Graph, drop: set[int]) -> ranker.Graph:
    """Build a copy of the graph without the given edges.

    The real graph is never edited, so every case starts from the same graph and the removed edges
    are back in place for the next case.

    Args:
        graph: The graph.
        drop: The edge indices to remove.

    Returns:
        The reduced graph.
    """
    kept = [edge for index, edge in enumerate(graph.edges) if index not in drop]
    return ranker.build_graph(graph.nodes, kept, merge=False)


def random_hit_probability(pool: int, true_drugs: int, top_k: int) -> float:
    """Return the chance that a random order puts a true drug inside the top k.

    Args:
        pool: The number of ranked candidates.
        true_drugs: The number of true drugs inside that pool.
        top_k: The size of the top list.

    Returns:
        The probability, between 0 and 1.
    """
    if pool <= 0 or true_drugs <= 0:
        return 0.0
    miss = 1.0
    for step in range(min(top_k, pool)):
        remaining = pool - step
        if remaining <= 0:
            break
        miss *= max(0.0, (remaining - true_drugs)) / remaining
    return 1.0 - miss


def poisson_binomial_tail(probabilities: list[float], at_least: int) -> float:
    """Return the chance that independent trials give at least this many hits.

    Each held-out case has its own chance of a random hit, because each case ranks a pool of its
    own size. The exact distribution of the total is a Poisson binomial, and this function sums its
    upper tail.

    Args:
        probabilities: The random hit probability of every case.
        at_least: The observed number of recoveries.

    Returns:
        The upper tail probability, between 0 and 1.
    """
    if at_least <= 0:
        return 1.0
    distribution = [1.0]
    for probability in probabilities:
        shifted = [0.0] * (len(distribution) + 1)
        for hits, mass in enumerate(distribution):
            shifted[hits] += mass * (1.0 - probability)
            shifted[hits + 1] += mass * probability
        distribution = shifted
    return sum(distribution[min(at_least, len(distribution) - 1):])


def dispersion_factor(counts: list[int], cases: int, rate: float) -> float | None:
    """Return the variance of the permuted count divided by the binomial variance.

    The Poisson binomial tail of this control treats the cases of one graph as independent Bernoulli
    trials with one common rate, so its variance is `cases * rate * (1 - rate)`. Each permutation
    gives one draw of the recovery count over the same diseases, so the sample variance of those
    counts estimates the true variance of the null. The ratio of the two is the dispersion factor.

    Args:
        counts: The recovery count of every permutation.
        cases: The number of held-out cases of the stratum.
        rate: The null rate that the tail reads.

    Returns:
        The ratio, or `None` when fewer than two permutations exist or the binomial variance is 0.
    """
    if len(counts) < 2 or cases <= 0:
        return None
    expected = cases * rate * (1.0 - rate)
    if expected <= 0.0:
        return None
    return statistics.variance(counts) / expected


def beta_binomial_tail(cases: int, rate: float, dispersion: float, at_least: int) -> float:
    """Return the upper tail of a beta-binomial matched to the permuted count.

    The distribution has mean `cases * rate` and variance `dispersion * cases * rate * (1 - rate)`,
    so it carries the first two moments of the permutation distribution. The intra-class correlation
    is `rho = (dispersion - 1) / (cases - 1)` and the two shape parameters are
    `a = rate * (1 - rho) / rho` and `b = (1 - rate) * (1 - rho) / rho`. At `dispersion` equal to 1
    the correlation is 0 and the distribution is the binomial of every recorded run, so the
    correction reduces to the uncorrected value when the permutations show no dispersion.

    Args:
        cases: The number of held-out cases of the stratum.
        rate: The null rate, between 0 and 1.
        dispersion: The dispersion factor, at or above 1.
        at_least: The recovery count of the real graph.

    Returns:
        The upper tail probability, between 0 and 1.
    """
    if at_least <= 0:
        return 1.0
    if cases <= 0 or at_least > cases:
        return 0.0
    if rate <= 0.0:
        return 0.0
    if rate >= 1.0:
        return 1.0
    correlation = (dispersion - 1.0) / (cases - 1.0) if cases > 1 else 0.0
    if correlation <= 0.0:
        return poisson_binomial_tail([rate] * cases, at_least)
    if correlation >= 1.0:
        return rate
    scale = (1.0 - correlation) / correlation
    alpha = rate * scale
    beta = (1.0 - rate) * scale
    constant = math.lgamma(alpha + beta) - math.lgamma(alpha) - math.lgamma(beta)
    total = 0.0
    for hits in range(at_least, cases + 1):
        log_mass = (math.lgamma(cases + 1) - math.lgamma(hits + 1) - math.lgamma(cases - hits + 1)
                    + math.lgamma(hits + alpha) + math.lgamma(cases - hits + beta)
                    - math.lgamma(cases + alpha + beta) + constant)
        total += math.exp(log_mass)
    return min(1.0, total)


def corrected_tail(counts: list[int], cases: int, rate: float, at_least: int) -> dict:
    """Correct the binomial tail for the dispersion that the permutations measure.

    Args:
        counts: The recovery count of every permutation.
        cases: The number of held-out cases of the stratum.
        rate: The null rate that the tail reads.
        at_least: The recovery count of the real graph.

    Returns:
        The measured dispersion factor, the factor that the correction applies, the intra-class
        correlation and the corrected upper tail.
    """
    measured = dispersion_factor(counts, cases, rate)
    applied = MIN_DISPERSION_FACTOR if measured is None else max(MIN_DISPERSION_FACTOR, measured)
    correlation = (applied - 1.0) / (cases - 1.0) if cases > 1 else 0.0
    return {"dispersion_measured": None if measured is None else round(measured, 4),
            "dispersion_applied": round(applied, 4),
            "intraclass_correlation": round(correlation, 6),
            "tail": float(f"{beta_binomial_tail(cases, rate, applied, at_least):.3g}")}


def heldout_rows(graph: ranker.Graph, top_k: int, relations: tuple[str, ...],
                 options: ranker.RankOptions, evidence: ranker.DrugEvidence,
                 tag: int = 0) -> list[dict]:
    """Evaluate every eligible disease of one graph.

    The graph is never edited. Each case ranks a copy that lacks every direct edge between the
    disease and any chemical, so the removed edges are in place again for the next case.

    Args:
        graph: The graph to evaluate.
        top_k: The size of the top list that counts as a recovery.
        relations: The relation names that state an indication.
        options: The ranking parameters.
        evidence: The drug records.
        tag: The permutation number, 0 for the real graph.

    Returns:
        One row per eligible disease.
    """
    cases = treatment_cases(graph, relations)
    rows = []
    for disease in sorted(cases):
        case = cases[disease]
        reduced = graph_without(graph, case["edge_indices"])
        ranked = ranker.rank_candidates(reduced, [disease], options, evidence)
        ranks = rank_map(ranked)
        indications = (ranker.indication_chemicals(reduced) if options.require_drug_evidence
                       else set())
        eligible = [drug for drug in sorted(case["drugs"])
                    if not options.require_drug_evidence
                    or ranker.eligibility_of(drug, evidence, indications)]
        found = {drug: ranks[drug] for drug in sorted(case["drugs"]) if drug in ranks}
        best = min(found.values()) if found else 0
        baseline = random_hit_probability(len(ranked), len(found), top_k)
        rows.append({"permutation": tag, "disease": disease, "label": graph.label(disease),
                     "true_drugs": sorted(case["drugs"]),
                     "removed_edges": len(case["edge_indices"]),
                     "removed_indication_edges": case["treatment_edges"],
                     "candidates_after_removal": len(ranked),
                     "true_drugs_eligible": len(eligible),
                     "true_drugs_reachable": len(found),
                     "best_rank": best,
                     "recovered": bool(best and best <= top_k),
                     "random_hit_probability": round(baseline, 4)})
    return rows


def rows_at_top_k(rows: list[dict], top_k: int) -> list[dict]:
    """Rewrite case rows for another size of the top list.

    Every row carries the best rank of a true drug, the size of the candidate pool and the number
    of reachable true drugs, so a second value of k needs no second ranking.

    Args:
        rows: The case rows of one run.
        top_k: The new size of the top list.

    Returns:
        The same rows with `recovered` and `random_hit_probability` recomputed.
    """
    rewritten = []
    for row in rows:
        best = row["best_rank"]
        baseline = random_hit_probability(row["candidates_after_removal"],
                                          row["true_drugs_reachable"], top_k)
        rewritten.append({**row, "recovered": bool(best and best <= top_k),
                          "random_hit_probability": round(baseline, 4)})
    return rewritten


def worker_initializer(graph: ranker.Graph, top_k: int, relations: tuple[str, ...],
                       options: ranker.RankOptions, evidence: ranker.DrugEvidence,
                       random_seed: int, scope_names: tuple[str, ...] | None) -> None:
    """Fill the state of one worker process.

    Args:
        graph: The real graph.
        top_k: The size of the top list that counts as a recovery.
        relations: The relation names that state an approved indication.
        options: The ranking parameters.
        evidence: The drug records.
        random_seed: The seed of the run.
        scope_names: The relation names in the scope of the permutation, or `None` for every edge.
    """
    WORKER_STATE.update({"graph": graph, "top_k": top_k, "relations": relations,
                         "options": options, "evidence": evidence, "random_seed": random_seed,
                         "scope_names": scope_names})


def permutation_generator(random_seed: int, step: int) -> random.Random:
    """Return the generator of one permutation of the `per_permutation` stream.

    The seed of a permutation depends on the run seed and on the permutation number only, so the
    result of a run does not depend on the number of workers.

    Args:
        random_seed: The seed of the run.
        step: The permutation number, starting at 1.

    Returns:
        The generator of that permutation.
    """
    return random.Random(f"track2-null:{random_seed}:{step}")


def worker_null_rows(step: int) -> list[dict]:
    """Evaluate one permutation inside a worker process.

    Args:
        step: The permutation number, starting at 1.

    Returns:
        The case rows of that permutation.
    """
    state = WORKER_STATE
    permuted = permute_graph(state["graph"], permutation_generator(state["random_seed"], step),
                             state["scope_names"])
    return heldout_rows(permuted, state["top_k"], state["relations"], state["options"],
                        state["evidence"], tag=step)


def build_null_rows(graph: ranker.Graph, top_k: int, relations: tuple[str, ...],
                    options: ranker.RankOptions, evidence: ranker.DrugEvidence, permutations: int,
                    generator: random.Random, random_seed: int, scope: str, stream: str,
                    workers: int) -> list[dict]:
    """Evaluate every permutation of the null.

    Args:
        graph: The real graph.
        top_k: The size of the top list that counts as a recovery.
        relations: The relation names that state an approved indication.
        options: The ranking parameters.
        evidence: The drug records.
        permutations: The number of permuted graphs.
        generator: The generator of the `sequential` stream.
        random_seed: The seed of the run, which the `per_permutation` stream reads.
        scope: A value of `NULL_SCOPES`.
        stream: A value of `NULL_STREAMS`.
        workers: The number of worker processes. 1 runs in this process.

    Returns:
        The case rows of every permutation, tagged with the permutation number.

    Raises:
        ValueError: When more than one worker is asked for with the `sequential` stream.
    """
    scope_names = scope_relations(scope, relations)
    if permutations <= 0:
        return []
    if workers > 1 and stream != NULL_STREAM_PER_PERMUTATION:
        raise ValueError("more than one worker needs the per_permutation null stream.")
    if workers > 1:
        context = multiprocessing.get_context("fork")
        with context.Pool(processes=workers, initializer=worker_initializer,
                          initargs=(graph, top_k, relations, options, evidence, random_seed,
                                    scope_names)) as pool:
            batches = pool.map(worker_null_rows, range(1, permutations + 1), chunksize=1)
        return [row for batch in batches for row in batch]
    rows: list[dict] = []
    for step in range(1, permutations + 1):
        source = (generator if stream == NULL_STREAM_SEQUENTIAL
                  else permutation_generator(random_seed, step))
        rows += heldout_rows(permute_graph(graph, source, scope_names), top_k, relations, options,
                             evidence, tag=step)
    return rows


def recovery_summary(rows: list[dict], top_k: int) -> dict:
    """Summarize one set of case rows.

    The Poisson binomial tail against a uniform draw costs work in proportion to the square of the
    number of rows, and it answers a question about one graph. A pooled set of permuted rows is not
    one graph, so the field is reported for a set at or below `MAX_ROWS_FOR_UNIFORM_TAIL` rows and
    it is `None` above that size.

    Args:
        rows: The case rows.
        top_k: The size of the top list that counts as a recovery.

    Returns:
        The counts, the rate and the uniform random baseline of that set.
    """
    hits = sum(1 for row in rows if row["recovered"])
    rate = hits / len(rows) if rows else 0.0
    uniform = (sum(row["random_hit_probability"] for row in rows) / len(rows)) if rows else 0.0
    tail = None
    if len(rows) <= MAX_ROWS_FOR_UNIFORM_TAIL:
        tail = float(
            f"{poisson_binomial_tail([row['random_hit_probability'] for row in rows], hits):.3g}")
    return {"top_k": top_k, "cases": len(rows),
            "diseases": len({row["disease"] for row in rows}),
            "recovered": hits, "recovery_rate": round(rate, 4),
            "uniform_baseline_rate": round(uniform, 4),
            "chance_probability_uniform": tail}


def regularized_incomplete_beta(x: float, a: float, b: float) -> float:
    """Return the regularized incomplete beta function.

    The function is the continued fraction of Lentz, which is the standard evaluation. It supports
    the Student t distribution below, so the control needs no third party numerical package.

    Args:
        x: The upper limit, between 0 and 1.
        a: The first shape parameter, above 0.
        b: The second shape parameter, above 0.

    Returns:
        `I_x(a, b)`, between 0 and 1.
    """
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    if x > (a + 1.0) / (a + b + 2.0):
        return 1.0 - regularized_incomplete_beta(1.0 - x, b, a)
    front = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                     + a * math.log(x) + b * math.log(1.0 - x)) / a
    return front * beta_continued_fraction(x, a, b)


def beta_continued_fraction(x: float, a: float, b: float) -> float:
    """Return the continued fraction of the incomplete beta function.

    The evaluation is the modified Lentz algorithm.

    Args:
        x: The upper limit, between 0 and 1.
        a: The first shape parameter, above 0.
        b: The second shape parameter, above 0.

    Returns:
        The value of the continued fraction.
    """
    tiny = 1e-30
    c = 1.0
    d = 1.0 - (a + b) * x / (a + 1.0)
    d = tiny if abs(d) < tiny else d
    d = 1.0 / d
    result = d
    for step in range(1, 301):
        even = step * (b - step) * x / ((a + 2.0 * step - 1.0) * (a + 2.0 * step))
        d = 1.0 + even * d
        d = tiny if abs(d) < tiny else d
        c = 1.0 + even / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        result *= d * c
        odd = -(a + step) * (a + b + step) * x / ((a + 2.0 * step) * (a + 2.0 * step + 1.0))
        d = 1.0 + odd * d
        d = tiny if abs(d) < tiny else d
        c = 1.0 + odd / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        change = d * c
        result *= change
        if abs(change - 1.0) < 1e-14:
            break
    return result


def student_t_cdf(value: float, degrees: int) -> float:
    """Return the cumulative distribution function of the Student t distribution.

    Args:
        value: The quantile.
        degrees: The degrees of freedom, at least 1.

    Returns:
        The probability that a t variable is at or below the value.
    """
    tail = regularized_incomplete_beta(degrees / (degrees + value * value), degrees / 2.0, 0.5)
    return 1.0 - tail / 2.0 if value > 0 else tail / 2.0


def student_t_quantile(probability: float, degrees: int) -> float:
    """Return the quantile of the Student t distribution by bisection.

    Args:
        probability: The cumulative probability, strictly between 0 and 1.
        degrees: The degrees of freedom, at least 1.

    Returns:
        The quantile.
    """
    low, high = -1000.0, 1000.0
    for _ in range(200):
        middle = (low + high) / 2.0
        if student_t_cdf(middle, degrees) < probability:
            low = middle
        else:
            high = middle
    return (low + high) / 2.0


def null_counts_by_permutation(rows: list[dict], permutations: int) -> list[int]:
    """Return the recovery count of every permutation.

    Every permutation keeps the exact indication degree of every node, so a permuted graph holds
    the same eligible diseases as the real graph and the counts are paired with the real count.

    Args:
        rows: The null case rows, each with a `permutation` tag and a `recovered` flag.
        permutations: The number of permutations of the run.

    Returns:
        One count per permutation, in permutation order. A permutation with no row counts 0.
    """
    counts = {step: 0 for step in range(1, permutations + 1)}
    for row in rows:
        if row["recovered"]:
            counts[row["permutation"]] = counts.get(row["permutation"], 0) + 1
    return [counts[step] for step in sorted(counts)]


def null_case_counts_by_permutation(rows: list[dict], permutations: int) -> list[int]:
    """Return the number of evaluated cases of every permutation.

    Args:
        rows: The null case rows, each with a `permutation` tag.
        permutations: The number of permutations of the run.

    Returns:
        One case count per permutation, in permutation order.
    """
    counts = {step: 0 for step in range(1, permutations + 1)}
    for row in rows:
        counts[row["permutation"]] = counts.get(row["permutation"], 0) + 1
    return [counts[step] for step in sorted(counts)]


def permutation_quantile(real_recovered: int, counts: list[int]) -> float:
    """Return the permutation quantile probability of the recovery count.

    The value is `(1 + the number of permutations whose count reaches the real count) /
    (permutations + 1)`. It assumes nothing about the distribution of the count and its smallest
    reachable value is `1 / (permutations + 1)`.

    Args:
        real_recovered: The recovery count of the real graph.
        counts: The recovery count of every permutation.

    Returns:
        The probability, between 0 and 1.
    """
    if not counts:
        return 1.0
    reached = sum(1 for count in counts if count >= real_recovered)
    return (1 + reached) / (len(counts) + 1)


def null_rate_statistics(counts: list[int], case_counts: list[int],
                         confidence: float = DEFAULT_CONFIDENCE) -> dict:
    """Describe the uncertainty of the null rate estimate.

    The cases of one permutation share one permuted graph, so the permutations and not the case
    rows are the independent replicates. The standard error is therefore the standard error of the
    mean of the per-permutation rates, and its degrees of freedom are the number of permutations
    less one.

    Args:
        counts: The recovery count of every permutation.
        case_counts: The number of evaluated cases of every permutation.
        confidence: The confidence level of the interval.

    Returns:
        The per-permutation rates, their mean, their sample standard deviation, the standard error
        of the mean and the two limits of the interval. A run with fewer than two permutations
        reports no standard deviation and no interval.
    """
    rates = [count / cases for count, cases in zip(counts, case_counts) if cases]
    if not rates:
        return {"permutations_used": 0, "rates": [], "mean": None, "sample_sd": None,
                "standard_error": None, "confidence": confidence, "lower": None, "upper": None}
    mean = statistics.fmean(rates)
    if len(rates) < 2:
        return {"permutations_used": len(rates), "rates": [round(rate, 4) for rate in rates],
                "mean": round(mean, 4), "sample_sd": None, "standard_error": None,
                "confidence": confidence, "lower": None, "upper": None}
    sample_sd = statistics.stdev(rates)
    standard_error = sample_sd / math.sqrt(len(rates))
    width = student_t_quantile(1.0 - (1.0 - confidence) / 2.0, len(rates) - 1) * standard_error
    return {"permutations_used": len(rates), "rates": [round(rate, 4) for rate in rates],
            "mean": round(mean, 6), "sample_sd": round(sample_sd, 6),
            "standard_error": round(standard_error, 6), "confidence": confidence,
            "lower": round(max(0.0, mean - width), 6), "upper": round(min(1.0, mean + width), 6)}


def disease_facts(graph: ranker.Graph) -> dict[str, dict]:
    """Collect the properties of every disease that a regime rule reads.

    Args:
        graph: The real graph.

    Returns:
        For each disease, the number of associated genes, the number of drugs with a clinical edge,
        the number of distinct PMIDs on its edges and its identifier.
    """
    facts: dict[str, dict] = {}
    for node_id in graph.nodes:
        if graph.kind(node_id) != "disease":
            continue
        genes = set()
        drugs = set()
        pmids: set[str] = set()
        for neighbour, index in graph.adjacency.get(node_id, ()):
            edge = graph.edges[index]
            pmids.update(edge.pmids)
            if graph.kind(neighbour) == "gene":
                genes.add(neighbour)
            elif (graph.kind(neighbour) == ranker.CANDIDATE_KIND
                  and ranker.relation_class(edge.relation) == ranker.CLINICAL):
                drugs.add(neighbour)
        facts[node_id] = {"disease_id": node_id, "genes": len(genes), "drugs": len(drugs),
                          "pmids": len(pmids)}
    return facts


def regime_all(facts: dict, limits: dict) -> bool:
    """Accept every disease.

    Args:
        facts: The properties of one disease.
        limits: The thresholds of the run.

    Returns:
        True.
    """
    return True


def regime_mesh_supplementary(facts: dict, limits: dict) -> bool:
    """Accept a disease that MeSH files as a supplementary concept record.

    MeSH gives a main heading a `D` number and a supplementary concept record a `C` number, and it
    files most rare diseases as supplementary concept records. The rule reads the identifier that
    the graph already carries, so no hand written disease list enters the control.

    Args:
        facts: The properties of one disease.
        limits: The thresholds of the run.

    Returns:
        True when the MeSH identifier is a supplementary concept record.
    """
    return bool(re.search(r":C\d+$", facts["disease_id"]))


def regime_sparse(facts: dict, limits: dict) -> bool:
    """Accept a disease with few associated genes and few approved drugs.

    The rule measures the coverage of the graph and not the biology. A disease of low coverage is
    not the same object as a rare monogenic disease, and the notes say so.

    Args:
        facts: The properties of one disease.
        limits: The thresholds of the run.

    Returns:
        True when both counts are at or below the thresholds.
    """
    return facts["genes"] <= limits["max_genes"] and facts["drugs"] <= limits["max_drugs"]


def regime_drug_poor(facts: dict, limits: dict) -> bool:
    """Accept a disease with few approved drugs.

    The number of approved drugs of a disease is the property that decides whether the held-out
    control can measure anything. A disease with 168 approved drugs gives a random ranking many
    chances to hit one. A disease with one approved drug gives it one chance. The seed of this work
    has none at all, so the regime of the seed is the low end of this scale.

    Args:
        facts: The properties of one disease.
        limits: The thresholds of the run.

    Returns:
        True when the drug count is at or below the threshold.
    """
    return facts["drugs"] <= limits["max_drugs"]


def regime_drug_rich(facts: dict, limits: dict) -> bool:
    """Accept a disease with many approved drugs.

    The rule is the complement of `regime_drug_poor` inside the eligible diseases. A report that
    prints the two strata beside each other needs both from one run, because the null rows of a run
    do not reach the receipt in full.

    Args:
        facts: The properties of one disease.
        limits: The thresholds of the run.

    Returns:
        True when the drug count is above the threshold.
    """
    return facts["drugs"] > limits["max_drugs"]


REGIMES = {"all": regime_all, "mesh_supplementary": regime_mesh_supplementary,
           "sparse": regime_sparse, "drug_poor": regime_drug_poor,
           "drug_rich": regime_drug_rich}


def regime_diseases(graph: ranker.Graph, regime: str, limits: dict) -> set[str]:
    """Select the diseases of one regime.

    Args:
        graph: The real graph.
        regime: A key of `REGIMES`.
        limits: The thresholds of the run.

    Returns:
        The disease identifiers that the regime accepts.
    """
    rule = REGIMES[regime]
    return {node_id for node_id, facts in disease_facts(graph).items() if rule(facts, limits)}


def summarize_recovery(real_rows: list[dict], null_rows: list[dict], top_k: int, min_cases: int,
                       min_recovery: float, permutations: int,
                       confidence: float = DEFAULT_CONFIDENCE) -> dict:
    """Compare a real recovery with the permuted null of the same diseases.

    The summary carries five probabilities for one comparison. `permutation_quantile_probability`
    compares the real recovery count with the recovery count of each permutation and assumes
    nothing. `chance_probability` is the Poisson binomial tail that treats the smoothed null rate as
    a known per-case success probability, and the permutations show that the recovery count of a
    permuted graph is more variable than that tail allows, so the value is anti-conservative.
    `dispersion_corrected_probability` is the beta-binomial tail whose first two moments are the
    first two moments of the permuted count, which corrects that error.
    `chance_probability_at_null_upper_limit` and `dispersion_corrected_probability_at_null_upper_limit`
    repeat the two tails at the upper confidence limit of the null rate, which propagates the
    uncertainty of an estimate that rests on `permutations` graphs. `decision_probability` is the
    larger of the permutation quantile and the corrected tail, and the verdict reads it, so a PASS
    needs the assumption-free value and the corrected model value both at or below alpha.
    `smallest_reachable_permutation_probability` is `1 / (permutations + 1)`, and
    `probability_resolution` says whether the decision probability is finer than the permutation
    evidence.

    Args:
        real_rows: The case rows of the real graph.
        null_rows: The case rows of the permuted graphs.
        top_k: The size of the top list that counts as a recovery.
        min_cases: The smallest number of eligible diseases for a reported rate.
        min_recovery: The recovery rate floor of a PASS.
        permutations: The number of permuted graphs.
        confidence: The confidence level of the null rate interval.

    Returns:
        The summary and its verdict.
    """
    real = recovery_summary(real_rows, top_k)
    null = recovery_summary(null_rows, top_k) if null_rows else None
    quantile = None
    spread = None
    chance_upper = None
    corrected = None
    corrected_upper = None
    if null and null["cases"]:
        counts = null_counts_by_permutation(null_rows, permutations)
        cases = null_case_counts_by_permutation(null_rows, permutations)
        spread = null_rate_statistics(counts, cases, confidence)
        used = [count for count, size in zip(counts, cases) if size]
        quantile = permutation_quantile(real["recovered"], used)
        null_rate = (null["recovered"] + 1) / (null["cases"] + 1)
        chance = poisson_binomial_tail([null_rate] * real["cases"], real["recovered"])
        corrected = corrected_tail(used, real["cases"], null_rate, real["recovered"])
        if spread["upper"] is not None:
            chance_upper = poisson_binomial_tail([spread["upper"]] * real["cases"],
                                                 real["recovered"])
            corrected_upper = beta_binomial_tail(real["cases"], spread["upper"],
                                                 corrected["dispersion_applied"],
                                                 real["recovered"])
        comparison = "degree_preserving_permutation"
    else:
        null_rate = real["uniform_baseline_rate"]
        tail = real["chance_probability_uniform"]
        chance = tail if tail is not None else poisson_binomial_tail(
            [null_rate] * real["cases"], real["recovered"])
        comparison = "uniform_random_order"
    decision = chance if corrected is None else max(quantile, corrected["tail"])
    floor = 1.0 / (permutations + 1) if permutations > 0 else None
    resolution = "no_permutation_null"
    if floor is not None:
        resolution = "resolved" if decision >= floor else "below_permutation_resolution"
    if real["diseases"] < min_cases:
        verdict = "INSUFFICIENT_CASES"
    elif real["recovery_rate"] < min_recovery:
        verdict = "FAIL"
    elif comparison == "degree_preserving_permutation" and permutations + 1 < 1 / DEFAULT_ALPHA:
        verdict = "INSUFFICIENT_PERMUTATIONS"
    elif decision <= DEFAULT_ALPHA:
        verdict = "PASS"
    else:
        verdict = "FAIL"
    return {"top_k": top_k, "cases": real["cases"], "diseases": real["diseases"],
            "recovered": real["recovered"],
            "recovery_rate": real["recovery_rate"] if real["diseases"] >= min_cases else None,
            "uniform_baseline_rate": real["uniform_baseline_rate"],
            "null": null, "null_rate_smoothed": round(null_rate, 4), "comparison": comparison,
            "chance_probability": float(f"{chance:.3g}"), "alpha": DEFAULT_ALPHA,
            "permutations": permutations, "verdict": verdict,
            "null_rate_spread": spread,
            "permutation_quantile_probability": (None if quantile is None
                                                 else float(f"{quantile:.3g}")),
            "smallest_reachable_permutation_probability": (None if floor is None
                                                           else float(f"{floor:.3g}")),
            "chance_probability_at_null_upper_limit": (None if chance_upper is None
                                                       else float(f"{chance_upper:.3g}")),
            "null_count_dispersion_factor": (None if corrected is None
                                             else corrected["dispersion_measured"]),
            "dispersion_factor_applied": (None if corrected is None
                                          else corrected["dispersion_applied"]),
            "dispersion_corrected_probability": (None if corrected is None
                                                 else corrected["tail"]),
            "dispersion_corrected_probability_at_null_upper_limit":
                (None if corrected_upper is None else float(f"{corrected_upper:.3g}")),
            "decision_probability": float(f"{decision:.3g}"),
            "probability_resolution": resolution}


def stratum_summaries(real_rows: list[dict], null_rows: list[dict], selected: set[str] | None,
                      top_k: int, min_cases: int, min_recovery: float, permutations: int,
                      confidence: float) -> dict:
    """Summarize one stratum at one size of the top list.

    Args:
        real_rows: The case rows of the real graph.
        null_rows: The case rows of the permuted graphs.
        selected: The diseases of the stratum, or `None` for every disease.
        top_k: The size of the top list that counts as a recovery.
        min_cases: The smallest number of eligible diseases for a reported rate.
        min_recovery: The recovery rate floor of a PASS.
        permutations: The number of permuted graphs.
        confidence: The confidence level of the null rate interval.

    Returns:
        The summary of that stratum at that k.
    """
    real = rows_at_top_k(real_rows, top_k)
    null = rows_at_top_k(null_rows, top_k)
    if selected is not None:
        real = [row for row in real if row["disease"] in selected]
        null = [row for row in null if row["disease"] in selected]
    return summarize_recovery(real, null, top_k, min_cases, min_recovery, permutations, confidence)


def heldout_control(graph: ranker.Graph, top_k: int, relations: tuple[str, ...],
                    permutations: int = 0, generator: random.Random | None = None,
                    options: ranker.RankOptions | None = None,
                    evidence: ranker.DrugEvidence | None = None,
                    min_cases: int = DEFAULT_MIN_CASES,
                    min_recovery: float = DEFAULT_MIN_RECOVERY,
                    regime: str = "all", limits: dict | None = None,
                    null_scope: str = NULL_SCOPE_ALL_EDGES,
                    null_stream: str = NULL_STREAM_SEQUENTIAL, workers: int = 1,
                    random_seed: int = 0,
                    report_top_k: tuple[int, ...] = DEFAULT_REPORT_TOP_K,
                    confidence: float = DEFAULT_CONFIDENCE,
                    extra_regimes: tuple[str, ...] = ()) -> dict:
    """Run the held-out indication recovery control.

    The real graph gives the recovery rate. A degree-preserving permutation of the same graph gives
    the null rate, and the null rate is the right comparison. A uniform draw from the candidate pool
    is not, because the drug of a permuted indication edge keeps the degree of a real drug and a
    high degree node is easier to reach.

    A recovery rate over every disease of the graph answers a question about common diseases with
    approved drugs, because those diseases dominate the set. The regime restricts the same rows to
    the diseases of one regime and reports a second rate. Both rates enter the receipt.

    Args:
        graph: The real graph.
        top_k: The size of the top list that counts as a recovery.
        relations: The relation names that state an indication.
        permutations: The number of permuted graphs of the null. 0 turns the null off.
        generator: The random generator of the permutations.
        options: The ranking parameters.
        evidence: The drug records.
        min_cases: The smallest number of eligible diseases for a reported rate.
        min_recovery: The recovery rate floor of a PASS.
        regime: A key of `REGIMES`. The restricted block covers the diseases of that regime.
        limits: The thresholds that a regime rule reads.
        null_scope: A value of `NULL_SCOPES`. `all_edges` rewires every edge of the graph.
            `indication_edges` rewires the indication layer only.
        null_stream: A value of `NULL_STREAMS`.
        workers: The number of worker processes of the null.
        random_seed: The seed of the run, which the `per_permutation` stream reads.
        report_top_k: The further values of k that the result summarizes from the same rows.
        confidence: The confidence level of the null rate interval.
        extra_regimes: Further keys of `REGIMES`. Each one adds a stratum to `by_top_k`, so a
            report that prints several strata of one graph needs one run and not several.

    Returns:
        The result record that the receipt carries.
    """
    options = options or ranker.RankOptions()
    evidence = evidence or ranker.DrugEvidence()
    limits = limits or {"max_genes": 10, "max_drugs": 3}
    before = len(graph.edges)
    rows = heldout_rows(graph, top_k, relations, options, evidence)
    null_rows = build_null_rows(graph, top_k, relations, options, evidence, permutations,
                                generator or random.Random(0), random_seed, null_scope,
                                null_stream, workers)
    full = summarize_recovery(rows, null_rows, top_k, min_cases, min_recovery, permutations,
                              confidence)
    compact = [{"permutation": row["permutation"], "disease": row["disease"],
                "best_rank": row["best_rank"],
                "candidates_after_removal": row["candidates_after_removal"],
                "true_drugs_reachable": row["true_drugs_reachable"],
                "recovered": row["recovered"],
                "random_hit_probability": row["random_hit_probability"]}
               for row in null_rows[:MAX_NULL_ROWS_IN_RECEIPT]]
    restricted = None
    selected = None
    if regime != "all":
        selected = regime_diseases(graph, regime, limits)
        restricted = summarize_recovery([row for row in rows if row["disease"] in selected],
                                        [row for row in null_rows if row["disease"] in selected],
                                        top_k, min_cases, min_recovery, permutations, confidence)
        restricted["regime"] = regime
        restricted["regime_diseases_in_graph"] = len(selected)
    strata = {name: regime_diseases(graph, name, limits) for name in extra_regimes}
    by_top_k = {}
    for value in sorted(set(report_top_k) | {top_k}):
        block = {"full": stratum_summaries(rows, null_rows, None, value, min_cases, min_recovery,
                                           permutations, confidence)}
        if selected is not None:
            block["restricted"] = stratum_summaries(rows, null_rows, selected, value, min_cases,
                                                    min_recovery, permutations, confidence)
        for name, diseases in strata.items():
            block[name] = stratum_summaries(rows, null_rows, diseases, value, min_cases,
                                            min_recovery, permutations, confidence)
        by_top_k[str(value)] = block
    return {"top_k": top_k, "indication_relations": list(relations), "regime": regime,
            "regime_limits": limits, "min_cases": min_cases, "min_recovery": min_recovery,
            "null_scope": null_scope, "null_stream": null_stream, "workers": workers,
            "extra_regimes": list(extra_regimes),
            "extra_regime_diseases_in_graph": {name: len(diseases)
                                               for name, diseases in strata.items()},
            "full": full, "restricted": restricted, "by_top_k": by_top_k,
            "edges_before": before, "edges_after": len(graph.edges),
            "graph_restored": before == len(graph.edges),
            "verdict": (restricted or full)["verdict"], "results": rows,
            "null_results": compact,
            "null_results_truncated": len(null_rows) > MAX_NULL_ROWS_IN_RECEIPT}


def count_dispersion(summary: dict) -> dict | None:
    """Compare the spread of the permuted recovery count with the binomial spread.

    The Poisson binomial tail assumes that the cases of one graph are independent Bernoulli trials.
    Each permutation gives one draw of the recovery count over the same diseases, so the
    permutations test that assumption. A ratio above 1 means that the count is more variable than
    the assumption allows, and the tail then understates the probability of a large count. The
    sampling standard deviation of the natural logarithm of the ratio is `sqrt(2 / (B - 1))`, which
    the field `z` divides into that logarithm.

    Args:
        summary: The record of `summarize_recovery`.

    Returns:
        The two variances, their ratio and the standardized logarithm of the ratio, or `None` when
        the summary holds fewer than two permutations or no case.
    """
    spread = summary.get("null_rate_spread") or {}
    rates = spread.get("rates") or []
    cases = summary.get("cases") or 0
    if len(rates) < 2 or cases <= 0:
        return None
    mean_rate = statistics.fmean(rates)
    expected = cases * mean_rate * (1.0 - mean_rate)
    if expected <= 0:
        return None
    observed = statistics.variance([rate * cases for rate in rates])
    ratio = observed / expected
    standard = math.sqrt(2.0 / (len(rates) - 1))
    return {"observed_variance": round(observed, 4), "binomial_variance": round(expected, 4),
            "ratio": round(ratio, 4),
            "z": round(math.log(ratio) / standard, 3) if ratio > 0 else None,
            "permutations_used": len(rates)}


def counts_from_stored_rates(summary: dict) -> list[int]:
    """Rebuild the recovery count of every permutation from a stored summary.

    A receipt written at 200 permutations truncates its null case rows, so the per-permutation
    rates of `null_rate_spread` are the only record of the permutation distribution. A rate is
    stored to four decimal places and a stratum holds at most 397 cases, so the product of the two
    recovers the integer count with no ambiguity.

    Args:
        summary: The record of `summarize_recovery` as a receipt stores it.

    Returns:
        One count per permutation, or an empty list when the summary holds no rate.
    """
    spread = summary.get("null_rate_spread") or {}
    rates = spread.get("rates") or []
    cases = summary.get("cases") or 0
    if not rates or cases <= 0:
        return []
    return [round(rate * cases) for rate in rates]


def correct_summary(summary: dict, min_cases: int, min_recovery: float) -> dict:
    """Add the dispersion-corrected probability and the corrected verdict to a stored summary.

    The correction is section 15 of `docs/benchmarks/t2_null_model_v1.md`. The function reranks
    nothing and it reads one stored summary. It leaves `chance_probability` unchanged, because every
    recorded run reports that field.

    Args:
        summary: The record of `summarize_recovery` as a receipt stores it.
        min_cases: The smallest number of eligible diseases for a reported rate.
        min_recovery: The recovery rate floor of a PASS.

    Returns:
        A copy of the summary with the dispersion fields, the decision probability, the resolution
        and the corrected verdict, and with the verdict of the receipt under `verdict_as_recorded`.
    """
    corrected = dict(summary)
    corrected["verdict_as_recorded"] = summary.get("verdict")
    counts = counts_from_stored_rates(summary)
    cases = summary.get("cases") or 0
    permutations = summary.get("permutations") or 0
    if len(counts) < 2 or cases <= 0:
        return corrected
    rate = summary["null_rate_smoothed"]
    block = corrected_tail(counts, cases, rate, summary["recovered"])
    quantile = permutation_quantile(summary["recovered"], counts)
    spread = summary.get("null_rate_spread") or {}
    upper = spread.get("upper")
    decision = max(quantile, block["tail"])
    floor = 1.0 / (permutations + 1) if permutations > 0 else None
    corrected.update({
        "null_count_dispersion_factor": block["dispersion_measured"],
        "dispersion_factor_applied": block["dispersion_applied"],
        "intraclass_correlation": block["intraclass_correlation"],
        "dispersion_corrected_probability": block["tail"],
        "dispersion_corrected_probability_at_null_upper_limit":
            (None if upper is None
             else float(f"{beta_binomial_tail(cases, upper, block['dispersion_applied'], summary['recovered']):.3g}")),
        "permutation_quantile_recomputed": float(f"{quantile:.3g}"),
        "permutation_quantile_matches_receipt":
            summary.get("permutation_quantile_probability") == float(f"{quantile:.3g}"),
        "decision_probability": float(f"{decision:.3g}"),
        "probability_resolution": ("no_permutation_null" if floor is None
                                   else ("resolved" if decision >= floor
                                         else "below_permutation_resolution"))})
    if (summary.get("diseases") or 0) < min_cases:
        corrected["verdict"] = "INSUFFICIENT_CASES"
    elif (summary.get("recovery_rate") or 0.0) < min_recovery:
        corrected["verdict"] = "FAIL"
    elif permutations + 1 < 1 / DEFAULT_ALPHA:
        corrected["verdict"] = "INSUFFICIENT_PERMUTATIONS"
    elif decision <= DEFAULT_ALPHA:
        corrected["verdict"] = "PASS"
    else:
        corrected["verdict"] = "FAIL"
    return corrected


def resolution_block(summary: dict) -> dict:
    """Keep the fields of one summary that describe the resolution of its probability.

    Args:
        summary: The record of `summarize_recovery`.

    Returns:
        The counts, the two rates, every probability, the dispersion of the permuted count and the
        resolution verdict.
    """
    return {"null_count_dispersion": count_dispersion(summary),
            "top_k": summary["top_k"], "cases": summary["cases"],
            "recovered": summary["recovered"], "recovery_rate": summary["recovery_rate"],
            "null_rate_smoothed": summary["null_rate_smoothed"],
            "null_rate_spread": summary["null_rate_spread"],
            "chance_probability": summary["chance_probability"],
            "chance_probability_at_null_upper_limit":
                summary["chance_probability_at_null_upper_limit"],
            "null_count_dispersion_factor": summary.get("null_count_dispersion_factor"),
            "dispersion_factor_applied": summary.get("dispersion_factor_applied"),
            "dispersion_corrected_probability": summary.get("dispersion_corrected_probability"),
            "dispersion_corrected_probability_at_null_upper_limit":
                summary.get("dispersion_corrected_probability_at_null_upper_limit"),
            "decision_probability": summary.get("decision_probability"),
            "verdict_as_recorded": summary.get("verdict_as_recorded", summary["verdict"]),
            "permutation_quantile_recomputed": summary.get("permutation_quantile_recomputed"),
            "permutation_quantile_matches_receipt":
                summary.get("permutation_quantile_matches_receipt"),
            "permutation_quantile_probability": summary["permutation_quantile_probability"],
            "smallest_reachable_permutation_probability":
                summary["smallest_reachable_permutation_probability"],
            "probability_resolution": summary["probability_resolution"],
            "verdict": summary["verdict"]}


def receipt_resolution(path: Path, report_top_k: tuple[int, ...],
                       confidence: float = DEFAULT_CONFIDENCE,
                       extra_regimes: tuple[str, ...] = ()) -> dict:
    """Recompute the resolution of every probability of one held-out receipt.

    The receipt names its graph, its regime and its permutation count, so the record needs no other
    input. The function reranks nothing. It reads the stored case rows, which carry the best rank of
    every case, and it recomputes the strata and the probabilities from them.

    Args:
        path: The held-out receipt.
        report_top_k: The values of k to recompute.
        confidence: The confidence level of the null rate interval.
        extra_regimes: Further keys of `REGIMES` to report beside the regime of the receipt.

    Returns:
        The record of that receipt, with one block per value of k and per stratum.

    Raises:
        ValueError: When the receipt truncated its null rows, because the strata cannot be rebuilt.
    """
    record = json.loads(path.read_text(encoding="utf-8"))
    result = record["result"]
    parameters = record["parameters"]
    if result.get("by_top_k"):
        min_cases = parameters["min_cases"]
        min_recovery = parameters["min_recovery"]
        return {"receipt": str(path), "receipt_id": record["receipt_id"],
                "permutations": parameters["permutations"],
                "regime": parameters.get("regime", "all"),
                "null_scope": result.get("null_scope", NULL_SCOPE_ALL_EDGES),
                "graph_edges": record["graph"]["edges"],
                "edges_sha256_prefix": record["graph"]["edges_sha256_prefix"],
                "source": "the by_top_k block of the receipt",
                "by_top_k": {value: {name: resolution_block(
                                         correct_summary(summary, min_cases, min_recovery))
                                     for name, summary in block.items()}
                             for value, block in result["by_top_k"].items()}}
    if result.get("null_results_truncated"):
        raise ValueError(f"{path} truncated its null rows, so no stratum can be rebuilt.")
    permutations = parameters["permutations"]
    graph = ranker.load_graph(Path(record["graph"]["nodes_file"]),
                              Path(record["graph"]["edges_file"]))
    regime = parameters.get("regime", "all")
    selected = (regime_diseases(graph, regime, parameters["regime_limits"])
                if regime != "all" else None)
    strata = {name: regime_diseases(graph, name, parameters["regime_limits"])
              for name in extra_regimes}
    blocks: dict[str, dict] = {}
    for value in sorted(set(report_top_k) | {parameters["top_k"]}):
        block = {"full": resolution_block(stratum_summaries(
            result["results"], result["null_results"], None, value, parameters["min_cases"],
            parameters["min_recovery"], permutations, confidence))}
        if selected is not None:
            block["restricted"] = resolution_block(stratum_summaries(
                result["results"], result["null_results"], selected, value,
                parameters["min_cases"], parameters["min_recovery"], permutations, confidence))
        for name, diseases in strata.items():
            block[name] = resolution_block(stratum_summaries(
                result["results"], result["null_results"], diseases, value,
                parameters["min_cases"], parameters["min_recovery"], permutations, confidence))
        blocks[str(value)] = block
    return {"receipt": str(path), "receipt_id": record["receipt_id"],
            "permutations": permutations, "regime": regime,
            "null_scope": result.get("null_scope", NULL_SCOPE_ALL_EDGES),
            "graph_edges": record["graph"]["edges"],
            "edges_sha256_prefix": record["graph"]["edges_sha256_prefix"],
            "source": "the stored case rows of the receipt",
            "by_top_k": blocks}


def file_digest(path: Path) -> str:
    """Return the SHA-256 prefix of a file.

    Args:
        path: The file.

    Returns:
        The first 16 characters of the hexadecimal digest.
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def write_receipt(record: dict, receipt_id: str, out_dir: Path) -> Path:
    """Write a receipt JSON.

    Args:
        record: The receipt content.
        receipt_id: The receipt identifier and the file stem.
        out_dir: The receipt directory.

    Returns:
        The path of the written receipt.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{receipt_id}.json"
    path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return path


def receipt_header(receipt_id: str, control: str, nodes_path: Path, edges_path: Path,
                   graph: ranker.Graph, seeds: list[str], note: str = "",
                   gene_layer: str = ranker.GENE_LAYER_COMENTION,
                   functional_dir: Path = ranker.DEFAULT_FUNCTIONAL_DIR) -> dict:
    """Build the fields that both receipts share.

    Args:
        receipt_id: The receipt identifier.
        control: The control name.
        nodes_path: The nodes.tsv path.
        edges_path: The edges.tsv path.
        graph: The loaded graph.
        seeds: The seed identifiers, empty for the held-out control.
        note: One sentence about the run, for example the origin of the graph.
        gene_layer: The gene to gene layer of the graph.
        functional_dir: The directory of the functional layer tables.

    Returns:
        The shared receipt fields.
    """
    layer = {"gene_layer": gene_layer}
    if gene_layer != ranker.GENE_LAYER_COMENTION:
        layer["functional_edges_sha256_prefix"] = file_digest(Path(functional_dir) / "edges.tsv")
        layer["functional_nodes_sha256_prefix"] = file_digest(Path(functional_dir) / "nodes.tsv")
    return {"receipt_id": receipt_id, "note": note,
            "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "track": "t2_rank", "control": control,
            "data_class_of_run": "A", "data_class_of_receipt": "A",
            "receipt_scope": "Public graph identifiers and PMIDs only. No patient data enters this "
                             "run.",
            "ranker": {"script": "scripts/track2_rank.py", "length_decay": ranker.LENGTH_DECAY,
                       "breadth_weight": ranker.BREADTH_WEIGHT,
                       "max_paths_per_candidate": ranker.MAX_PATHS_PER_CANDIDATE},
            "graph": {"nodes_file": str(nodes_path), "edges_file": str(edges_path),
                      "nodes_sha256_prefix": file_digest(nodes_path),
                      "edges_sha256_prefix": file_digest(edges_path),
                      "nodes": graph.stats["nodes"], "edges": graph.stats["edges"],
                      "merged_away": graph.stats["merged_away"],
                      "self_loops_dropped": graph.stats["self_loops_dropped"],
                      "implied_nodes": graph.stats["implied_nodes"], **layer},
            "seeds": seeds}


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    """Add the arguments that both controls share.

    Args:
        parser: The subcommand parser.
    """
    ranker.add_rank_arguments(parser)
    parser.add_argument("--receipt-dir", type=Path, default=DEFAULT_RECEIPT_DIR,
                        help="The receipt directory.")
    parser.add_argument("--receipt-id", default="", help="The receipt identifier and file stem.")
    parser.add_argument("--note", default="",
                        help="One sentence that the receipt carries, for example the origin of "
                             "the graph.")


def build_parser() -> argparse.ArgumentParser:
    """Build the command line parser.

    Returns:
        The parser with one subcommand per control.
    """
    parser = argparse.ArgumentParser(description="Run the Track 2 ranking controls.")
    sub = parser.add_subparsers(dest="control", required=True)

    permutation = sub.add_parser("permutation", help="Degree-preserving permutation control.")
    add_common_arguments(permutation)
    permutation.add_argument("--seed", action="append", default=[],
                             help="A seed node identifier. Repeat the option for a seed set.")
    permutation.add_argument("-n", "--permutations", type=int, default=DEFAULT_PERMUTATIONS,
                             help="The number of permuted graphs.")
    permutation.add_argument("--top-k", type=int, default=DEFAULT_TOP_K,
                             help="The number of real candidates to test.")
    permutation.add_argument("--alpha", type=float, default=DEFAULT_ALPHA,
                             help="The p-value threshold of a candidate that is not an artifact.")
    permutation.add_argument("--random-seed", type=int, default=20260915,
                             help="The seed of the random generator.")

    heldout = sub.add_parser("heldout", help="Held-out indication recovery control.")
    add_common_arguments(heldout)
    heldout.add_argument("--top-k", type=int, default=DEFAULT_TOP_K,
                         help="The size of the top list that counts as a recovery.")
    heldout.add_argument("--relations", default=",".join(TREATMENT_RELATIONS),
                         help="The relation names that state an approved indication.")
    heldout.add_argument("--min-cases", type=int, default=DEFAULT_MIN_CASES,
                         help="The smallest number of eligible diseases for a reported rate.")
    heldout.add_argument("--min-recovery", type=float, default=DEFAULT_MIN_RECOVERY,
                         help="The recovery rate floor of a PASS.")
    heldout.add_argument("-n", "--permutations", type=int, default=20,
                         help="The number of permuted graphs of the null. 0 turns the null off "
                              "and the verdict falls back on a uniform random order.")
    heldout.add_argument("--random-seed", type=int, default=20260915,
                         help="The seed of the random generator of the null.")
    heldout.add_argument("--regime", choices=sorted(REGIMES), default="all",
                         help="Report a second recovery rate over the diseases of this regime. "
                              "The rule reads the graph and no hand written disease list.")
    heldout.add_argument("--max-disease-genes", type=int, default=10,
                         help="The gene count threshold of the sparse regime.")
    heldout.add_argument("--max-disease-drugs", type=int, default=3,
                         help="The drug count threshold of the sparse regime and of the drug_poor "
                              "regime.")
    heldout.add_argument("--null-scope", choices=NULL_SCOPES, default=NULL_SCOPE_ALL_EDGES,
                         help="all_edges rewires every edge of the ranking graph. "
                              "indication_edges rewires the indication layer only and leaves the "
                              "gene, pathway and drug target layers in place.")
    heldout.add_argument("--null-stream", choices=NULL_STREAMS, default=NULL_STREAM_SEQUENTIAL,
                         help="sequential draws every permutation from one generator and "
                              "reproduces every run recorded before 2026-09-18. per_permutation "
                              "seeds each permutation from the run seed and the permutation "
                              "number, so the result does not depend on the worker count.")
    heldout.add_argument("--workers", type=int, default=1,
                         help="The number of worker processes of the null. More than one needs "
                              "the per_permutation null stream.")
    heldout.add_argument("--report-top-k", default=",".join(str(k) for k in DEFAULT_REPORT_TOP_K),
                         help="The further values of k that the receipt summarizes from the same "
                              "rows, with no second run.")
    heldout.add_argument("--confidence", type=float, default=DEFAULT_CONFIDENCE,
                         help="The confidence level of the null rate interval.")
    heldout.add_argument("--extra-regime", action="append", default=[], choices=sorted(REGIMES),
                         help="Add a further stratum to the by_top_k block of the receipt. Repeat "
                              "the option for a set of strata.")

    resolution = sub.add_parser("resolution",
                                help="Measure what the probabilities of a held-out receipt can "
                                     "resolve.")
    resolution.add_argument("--receipt", action="append", default=[], type=Path,
                            help="A held-out receipt. Repeat the option for a set of receipts.")
    resolution.add_argument("--report-top-k",
                            default=",".join(str(k) for k in DEFAULT_REPORT_TOP_K),
                            help="The values of k to recompute from the stored rows.")
    resolution.add_argument("--confidence", type=float, default=DEFAULT_CONFIDENCE,
                            help="The confidence level of the null rate interval.")
    resolution.add_argument("--receipt-dir", type=Path, default=DEFAULT_RECEIPT_DIR,
                            help="The receipt directory.")
    resolution.add_argument("--receipt-id", default="",
                            help="The receipt identifier and file stem.")
    resolution.add_argument("--note", default="", help="One sentence that the receipt carries.")
    resolution.add_argument("--extra-regime", action="append", default=[],
                            choices=sorted(REGIMES),
                            help="Report a further stratum beside the regime of the receipt. "
                                 "Repeat the option for a set of strata.")
    return parser


def rank_option_record(options: ranker.RankOptions) -> dict:
    """Write the ranking parameters into a receipt.

    Args:
        options: The ranking parameters.

    Returns:
        The parameters as plain values.
    """
    return {"max_edges": options.max_edges, "hub_degree_cap": options.hub_degree_cap,
            "terminal_classes": list(options.terminal_classes),
            "require_drug_evidence": options.require_drug_evidence,
            "blocked_intermediate_kinds": list(options.blocked_intermediate_kinds),
            "defect": options.defect, "defect_genes": list(options.defect_genes),
            "extra_pharmacological": list(options.extra_pharmacological)}


def print_recovery(title: str, summary: dict) -> None:
    """Print one recovery summary.

    Args:
        title: The name of the disease set.
        summary: The summary record.
    """
    if summary["recovery_rate"] is None:
        print(f"{title}: {summary['diseases']} diseases, below the floor. No rate is reported. "
              f"{summary['verdict']}")
        return
    line = (f"{title}: {summary['recovered']}/{summary['cases']} = {summary['recovery_rate']} "
            f"at top {summary['top_k']}")
    if summary["null"]:
        line += (f", null {summary['null']['recovered']}/{summary['null']['cases']} = "
                 f"{summary['null']['recovery_rate']}")
    print(line + f", decision p {summary.get('decision_probability')}, {summary['verdict']}")
    if summary.get("permutation_quantile_probability") is not None:
        spread = summary.get("null_rate_spread") or {}
        print(f"    permutation quantile {summary['permutation_quantile_probability']}, "
              f"floor {summary['smallest_reachable_permutation_probability']}, "
              f"corrected tail {summary.get('dispersion_corrected_probability')} at dispersion "
              f"{summary.get('dispersion_factor_applied')}, "
              f"uncorrected tail {summary['chance_probability']}, "
              f"corrected tail at the null upper limit "
              f"{summary.get('dispersion_corrected_probability_at_null_upper_limit')}, "
              f"null rate {spread.get('mean')} "
              f"[{spread.get('lower')}, {spread.get('upper')}], "
              f"{summary['probability_resolution']}")


def run_permutation(args: argparse.Namespace) -> int:
    """Run the permutation control from the command line.

    Args:
        args: The parsed arguments.

    Returns:
        0 when the verdict is PASS, 1 otherwise, 2 on an input error.
    """
    try:
        graph = ranker.load_layered_graph(args.nodes, args.edges, args.gene_layer,
                                          args.functional_dir)
        seeds = ranker.resolve_seeds(graph, list(args.seed))
        options = ranker.options_from_args(args)
        evidence = ranker.load_layered_drug_evidence(args.drug_evidence, args.gene_layer,
                                                     args.functional_dir)
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    if not seeds:
        print("Give at least one --seed.", file=sys.stderr)
        return 2
    result = permutation_control(graph, seeds, args.permutations, args.top_k, args.random_seed,
                                 options=options, evidence=evidence, alpha=args.alpha)
    receipt_id = args.receipt_id or f"t2_permutation_{ranker.timestamp()}"
    record = receipt_header(receipt_id, "degree_preserving_permutation", args.nodes, args.edges,
                            graph, seeds, args.note, args.gene_layer, args.functional_dir)
    record["parameters"] = {"permutations": args.permutations, "top_k": args.top_k,
                            "alpha": args.alpha, "random_seed": args.random_seed}
    record["rank_options"] = rank_option_record(options)
    record["result"] = result
    record["verdict"] = result["verdict"]
    record["reading"] = ("The p-value is the share of permuted graphs that reach the real rank of "
                         "the candidate or a better one. A p-value above alpha means that a graph "
                         "with the same degree sequence and no real evidence reaches that rank "
                         "often, so the rank is not separated from node degree. Report such a "
                         "candidate only with this warning.")
    path = write_receipt(record, receipt_id, args.receipt_dir)
    print(f"permutations: {args.permutations}, degree sequence preserved: "
          f"{result['degree_sequence_preserved']}")
    for row in result["candidates"]:
        print(f"  {row['rank']:>3}  {row['label'][:32]:<32} p={row['p_value']:<7} "
              f"degree={row['degree']:<5} {row['verdict']}")
    print(f"verdict: {result['verdict']}, survivors: {result['survivors']} of "
          f"{len(result['candidates'])}")
    print(f"receipt: {path}")
    return 0 if result["verdict"] == "PASS" else 1


def run_heldout(args: argparse.Namespace) -> int:
    """Run the held-out recovery control from the command line.

    Args:
        args: The parsed arguments.

    Returns:
        0 when the verdict is PASS, 1 otherwise, 2 on an input error.
    """
    try:
        graph = ranker.load_layered_graph(args.nodes, args.edges, args.gene_layer,
                                          args.functional_dir)
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    try:
        options = ranker.options_from_args(args)
        evidence = ranker.load_layered_drug_evidence(args.drug_evidence, args.gene_layer,
                                                     args.functional_dir)
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    relations = tuple(part.strip().lower() for part in args.relations.split(",") if part.strip())
    limits = {"max_genes": args.max_disease_genes, "max_drugs": args.max_disease_drugs}
    report_top_k = tuple(int(part) for part in args.report_top_k.split(",") if part.strip())
    try:
        result = heldout_control(graph, args.top_k, relations, permutations=args.permutations,
                                 generator=random.Random(args.random_seed), options=options,
                                 evidence=evidence, min_cases=args.min_cases,
                                 min_recovery=args.min_recovery, regime=args.regime, limits=limits,
                                 null_scope=args.null_scope, null_stream=args.null_stream,
                                 workers=args.workers, random_seed=args.random_seed,
                                 report_top_k=report_top_k, confidence=args.confidence,
                                 extra_regimes=tuple(args.extra_regime))
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    receipt_id = args.receipt_id or f"t2_heldout_{ranker.timestamp()}"
    record = receipt_header(receipt_id, "held_out_indication_recovery", args.nodes, args.edges,
                            graph, [], args.note, args.gene_layer, args.functional_dir)
    record["parameters"] = {"top_k": args.top_k, "indication_relations": list(relations),
                            "min_cases": args.min_cases, "min_recovery": args.min_recovery,
                            "permutations": args.permutations, "random_seed": args.random_seed,
                            "regime": args.regime, "regime_limits": limits,
                            "null_scope": args.null_scope, "null_stream": args.null_stream,
                            "workers": args.workers, "report_top_k": list(report_top_k),
                            "confidence": args.confidence,
                            "extra_regimes": list(args.extra_regime)}
    record["rank_options"] = rank_option_record(options)
    record["result"] = result
    record["verdict"] = result["verdict"]
    record["reading"] = ("The recovery rate is the share of held-out diseases whose true drug "
                         "returns inside the top k after every direct disease to chemical edge of "
                         "that disease is removed. The null rate comes from the same evaluation on "
                         "degree-preserving permutations of the graph, and it is the rate that a "
                         "ranking with no real evidence reaches. A run below min_cases reports no "
                         "rate at all. Read `permutation_quantile_probability` first, because it "
                         "assumes nothing and its smallest reachable value is "
                         "`1 / (permutations + 1)`. `chance_probability` is a Poisson binomial "
                         "tail that treats the smoothed null rate as known, and a value below "
                         "that floor is a property of the model and not a measurement. "
                         "`chance_probability_at_null_upper_limit` repeats the tail at the upper "
                         "confidence limit of the null rate.")
    path = write_receipt(record, receipt_id, args.receipt_dir)
    if len(result["results"]) <= 20:
        for row in result["results"]:
            print(f"  {row['label'][:34]:<34} removed={row['removed_edges']:<3} "
                  f"best_rank={row['best_rank']:<5} recovered={row['recovered']}")
    print_recovery("full graph", result["full"])
    if result["restricted"]:
        print_recovery(f"regime {result['restricted']['regime']}", result["restricted"])
    print(f"verdict: {result['verdict']}")
    print(f"receipt: {path}")
    return 0 if result["verdict"] == "PASS" else 1


def run_resolution(args: argparse.Namespace) -> int:
    """Measure what the probabilities of one or more held-out receipts can resolve.

    Args:
        args: The parsed arguments.

    Returns:
        0 on success, 2 on an input error.
    """
    if not args.receipt:
        print("Give at least one --receipt.", file=sys.stderr)
        return 2
    report_top_k = tuple(int(part) for part in args.report_top_k.split(",") if part.strip())
    records = []
    try:
        for path in args.receipt:
            records.append(receipt_resolution(path, report_top_k, args.confidence,
                                              tuple(args.extra_regime)))
    except (OSError, ValueError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    receipt_id = args.receipt_id or f"t2_null_resolution_{ranker.timestamp()}"
    record = {"receipt_id": receipt_id, "note": args.note,
              "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
              "track": "t2_rank", "control": "held_out_probability_resolution",
              "data_class_of_run": "A", "data_class_of_receipt": "A",
              "receipt_scope": "Public graph identifiers only. No patient data enters this run.",
              "scripts": {"analysis": "scripts/track2_controls.py resolution"},
              "parameters": {"report_top_k": list(report_top_k), "confidence": args.confidence,
                             "extra_regimes": list(args.extra_regime)},
              "reading": ("A comparison against B permuted graphs cannot resolve a probability "
                          "below 1 / (B + 1). `chance_probability` is a Poisson binomial tail that "
                          "treats the smoothed null rate as a known constant, and a value below "
                          "that floor is a property of the model. "
                          "`permutation_quantile_probability` is the value that the permutation "
                          "evidence supports. "
                          "`chance_probability_at_null_upper_limit` repeats the tail at the upper "
                          "confidence limit of the null rate, and the limit reads the permutations "
                          "and not the case rows as the independent replicates. "
                          "`dispersion_corrected_probability` is the beta-binomial tail whose "
                          "first two moments are the first two moments of the permuted recovery "
                          "count, which corrects the over-dispersion that "
                          "`null_count_dispersion_factor` measures. `decision_probability` is the "
                          "larger of the permutation quantile and that corrected tail, and the "
                          "field `verdict` reads it, while `verdict_as_recorded` keeps the verdict "
                          "of the source receipt. The correction is section 15 of "
                          "docs/benchmarks/t2_null_model_v1.md."),
              "receipts": records}
    path = write_receipt(record, receipt_id, args.receipt_dir)
    for item in records:
        for value, block in item["by_top_k"].items():
            for name, summary in block.items():
                spread = summary["null_rate_spread"] or {}
                print(f"{item['receipt_id']} top {value:>2} {name:<10} "
                      f"rate {summary['recovery_rate']} null {summary['null_rate_smoothed']} "
                      f"[{spread.get('lower')}, {spread.get('upper')}] "
                      f"uncorrected {summary['chance_probability']} "
                      f"corrected {summary.get('dispersion_corrected_probability')} "
                      f"corrected upper "
                      f"{summary.get('dispersion_corrected_probability_at_null_upper_limit')} "
                      f"quantile {summary['permutation_quantile_probability']} "
                      f"decision {summary.get('decision_probability')} "
                      f"dispersion {(summary.get('null_count_dispersion') or {}).get('ratio')} "
                      f"floor {summary['smallest_reachable_permutation_probability']} "
                      f"{summary['probability_resolution']} "
                      f"{summary.get('verdict_as_recorded')} -> {summary['verdict']}")
    print(f"receipt: {path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Run one control from the command line.

    Args:
        argv: The arguments. `sys.argv[1:]` is used when it is absent.

    Returns:
        The exit code of the control.
    """
    args = build_parser().parse_args(argv)
    if args.control == "permutation":
        return run_permutation(args)
    if args.control == "resolution":
        return run_resolution(args)
    return run_heldout(args)


if __name__ == "__main__":
    sys.exit(main())
