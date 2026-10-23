#!/usr/bin/env python3
"""Run the mismatched phenotype negative control of the B0 phenotype layer. Data class A.

Defect R2-3 of the second jury round states that the phenotype layer shipped under a gate that its
own negative control would have failed, and that the gate had no mismatched phenotype arm. A layer
that raises the score of the right gene with the wrong phenotype is not reading the phenotype.

Four arms run on the same candidate records and the same scorer.

    variant_only                the shipped ranker with the layer off
    phenotype_matched           the case profile of the disease of the scenario gene
    phenotype_mismatched        the case profile of an unrelated disease, under a constrained
                                derangement of the profile assignment
    phenotype_source_held_out   the matched profile, with every annotation row that produced a term
                                of the profile removed from the gene profiles

The derangement is constrained. A scenario receives the profile of a donor scenario only when the
gene of the donor and the gene of the receiver share no Mendelian disease in `genes_to_disease.txt`.
The receipt records every scenario that the constraint could not place.

The design and the amended deployment gate are in `docs/benchmarks/t1_phenotype_mismatch_v1.md`,
recorded before the mismatched arm ran.

Every input is public. The case profiles come from the HPO release alone. No term of the child
enters any profile and the script opens no patient file.

Usage:
    uv run python scripts/track1_phenotype_mismatch_control.py --background bg.vep.vcf.gz \\
        --spikes spikes.vep.vcf.gz --manifest reports/benchmarks/spikein_v1/manifest.json \\
        --literature reports/known_alleles/known_alleles_grch38_2026_09_15.tsv \\
        --hpo-dir tracking/cache/hpo_v2026_09_01 --out-dir <dir> \\
        --receipt tracking/receipts/t1_phenotype_mismatch_v1.json
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import random
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "pipelines" / "track1_b0" / "bin"
sys.path.insert(0, str(BIN))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import b0_rank  # noqa: E402
import phenotype_layer  # noqa: E402
import track1_interval_stats as stats  # noqa: E402
import track1_phenotype_suite as suite  # noqa: E402

DERANGEMENT_SEED = 20260917
DERANGEMENT_ATTEMPTS = 200
ARMS = ("variant_only", "phenotype_matched", "phenotype_mismatched", "phenotype_source_held_out")


def disease_sets(gene_diseases: dict[str, list[str]]) -> dict[str, set[str]]:
    """Return the Mendelian disease set of every gene.

    Args:
        gene_diseases: The Mendelian diseases of every gene.

    Returns:
        The gene symbol mapped to its disease set.
    """
    return {gene: set(diseases) for gene, diseases in gene_diseases.items()}


def unrelated(first: str, second: str, sets: dict[str, set[str]]) -> bool:
    """State whether two genes share no Mendelian disease.

    Args:
        first: The first gene symbol.
        second: The second gene symbol.
        sets: The disease set of every gene.

    Returns:
        True when the two genes differ and share no disease.
    """
    if first == second:
        return False
    return not (sets.get(first, set()) & sets.get(second, set()))


def constrained_derangement(order: list[str], genes: dict[str, str], sets: dict[str, set[str]],
                            rng: random.Random) -> dict[str, str]:
    """Assign to every scenario the case profile of an unrelated scenario.

    The function shuffles the donors and keeps the first assignment in which every scenario receives
    the profile of a scenario whose gene shares no Mendelian disease with its own gene. When no
    attempt succeeds, the function returns the best attempt and the caller records the scenarios
    that the constraint could not place.

    Args:
        order: The scenario identifiers in suite order.
        genes: The gene of every scenario.
        sets: The disease set of every gene.
        rng: The seeded random generator.

    Returns:
        The scenario identifier mapped to the donor scenario identifier.
    """
    best: dict[str, str] = {}
    best_violations = len(order) + 1
    for _ in range(DERANGEMENT_ATTEMPTS):
        donors = order[:]
        rng.shuffle(donors)
        assignment = dict(zip(order, donors))
        violations = sum(1 for receiver, donor in assignment.items()
                         if not unrelated(genes[receiver], genes[donor], sets))
        if violations < best_violations:
            best, best_violations = assignment, violations
        if violations == 0:
            return assignment
    return best


def unplaced(assignment: dict[str, str], genes: dict[str, str],
             sets: dict[str, set[str]]) -> list[str]:
    """List the scenarios whose donor is not an unrelated gene.

    Args:
        assignment: The donor of every scenario.
        genes: The gene of every scenario.
        sets: The disease set of every gene.

    Returns:
        The scenario identifiers, sorted.
    """
    return sorted(receiver for receiver, donor in assignment.items()
                  if not unrelated(genes[receiver], genes[donor], sets))


def gate_verdict(metrics: dict) -> dict:
    """Evaluate the four conditions of the amended deployment gate.

    Args:
        metrics: The arm metrics of the four arms.

    Returns:
        The four conditions, the outcome and the weight that follows.
    """
    conditions = {
        "A1_matched_reaches_variant_only_rank_points": {
            "matched": metrics["phenotype_matched"]["mean_rank_points"],
            "variant_only": metrics["variant_only"]["mean_rank_points"],
            "holds": metrics["phenotype_matched"]["mean_rank_points"]
                     >= metrics["variant_only"]["mean_rank_points"],
            "prospective": True},
        "A2_matched_reaches_variant_only_f_max": {
            "matched": metrics["phenotype_matched"]["mean_f_max"],
            "variant_only": metrics["variant_only"]["mean_f_max"],
            "holds": metrics["phenotype_matched"]["mean_f_max"]
                     >= metrics["variant_only"]["mean_f_max"],
            "prospective": True},
        "A3_source_held_out_reaches_variant_only_rank_points": {
            "source_held_out": metrics["phenotype_source_held_out"]["mean_rank_points"],
            "variant_only": metrics["variant_only"]["mean_rank_points"],
            "holds": metrics["phenotype_source_held_out"]["mean_rank_points"]
                     >= metrics["variant_only"]["mean_rank_points"],
            "prospective": False,
            "note": "The value of this arm was known on 2026-09-17 before the gate was amended, so "
                    "this condition is not an independent test."},
        "A4_matched_beats_mismatched_rank_points": {
            "matched": metrics["phenotype_matched"]["mean_rank_points"],
            "mismatched": metrics["phenotype_mismatched"]["mean_rank_points"],
            "holds": metrics["phenotype_matched"]["mean_rank_points"]
                     > metrics["phenotype_mismatched"]["mean_rank_points"],
            "prospective": True},
    }
    holds = all(item["holds"] for item in conditions.values())
    return {"record": "docs/benchmarks/t1_phenotype_mismatch_v1.md",
            "conditions": conditions,
            "outcome": "GATE MET" if holds else "GATE NOT MET",
            "shipped_phenotype_weight": 2.0 if holds else 0.0}


def main() -> int:
    """Run the four arms and write the receipt.

    Returns:
        Zero after the receipt is written.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--background", required=True, type=Path)
    parser.add_argument("--spikes", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--literature", required=True, type=Path)
    parser.add_argument("--hpo-dir", required=True, type=Path)
    parser.add_argument("--hpo-release", default="v2026-09-01")
    parser.add_argument("--weight", type=float, default=suite.PHENOTYPE_WEIGHT)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()

    obo = args.hpo_dir / "hp.obo"
    genes_to_phenotype = args.hpo_dir / "genes_to_phenotype.txt"
    genes_to_disease = args.hpo_dir / "genes_to_disease.txt"
    hpoa = args.hpo_dir / "phenotype.hpoa"

    parents, alternates, _ = phenotype_layer.parse_obo(obo)
    closure = phenotype_layer.ancestor_closure(parents)
    descendants = phenotype_layer.descendant_map(closure)
    profiles = phenotype_layer.load_gene_phenotypes(genes_to_phenotype, closure, alternates)
    parsed = {"closure": closure, "descendants": descendants, "alternates": alternates,
              "profiles": profiles}
    abnormality_terms = sorted(term for term in descendants[suite.PHENOTYPIC_ABNORMALITY]
                               if term != suite.PHENOTYPIC_ABNORMALITY)

    gene_diseases = suite.read_gene_diseases(genes_to_disease)
    annotations = suite.read_disease_annotations(hpoa)
    disease_genes: dict[str, set[str]] = defaultdict(set)
    for gene, diseases in gene_diseases.items():
        for disease in diseases:
            disease_genes[disease].add(gene)
    sets = disease_sets(gene_diseases)

    manifest = json.loads(args.manifest.read_text())
    literature = b0_rank.load_literature(args.literature)
    spike_ids = {f"{v['chrom']}_{v['pos']}_{v['ref']}_{v['alt']}"
                 for s in manifest["scenarios"] for v in s["variants"]}
    background = b0_rank.parse_vcf(args.background, literature)
    spikes = b0_rank.parse_vcf(args.spikes, literature, spike_ids)

    # The same seed and the same order as track1_phenotype_suite, so the profiles are identical.
    rng = random.Random(suite.SEED)
    case_profiles = {}
    for scenario in manifest["scenarios"]:
        case_profiles[scenario["scenario_id"]] = suite.build_case_profile(
            scenario["gene"], rng, gene_diseases, annotations, parents, abnormality_terms)

    all_references = {reference for profile in case_profiles.values()
                      for reference in profile["references"]}
    by_reference, support = suite.annotation_support(all_references, annotations, disease_genes)

    order = [scenario["scenario_id"] for scenario in manifest["scenarios"]]
    genes = {scenario["scenario_id"]: scenario["gene"] for scenario in manifest["scenarios"]}
    assignment = constrained_derangement(order, genes, sets, random.Random(DERANGEMENT_SEED))
    not_placed = unplaced(assignment, genes, sets)

    matched_model = suite.build_model(parsed, None)
    results: dict[str, list[dict]] = {arm: [] for arm in ARMS}
    per_scenario_layer = []
    for arm in ARMS:
        (args.out_dir / arm).mkdir(parents=True, exist_ok=True)

    for scenario in manifest["scenarios"]:
        scenario_id = scenario["scenario_id"]
        keys = {(v["chrom"], int(v["pos"]), v["ref"].upper(), v["alt"].upper())
                for v in scenario["variants"]}
        records = background + [record for record in spikes if record["key"] in keys]
        candidate_genes = sorted({record["gene"] for record in records})
        profile = case_profiles[scenario_id]
        donor_id = assignment[scenario_id]
        donor_profile = case_profiles[donor_id]

        matched_terms = matched_model.normalize(profile["terms"])
        matched_scores = {gene: matched_model.score_gene(gene, matched_terms)
                          for gene in candidate_genes}
        donor_terms = matched_model.normalize(donor_profile["terms"])
        donor_scores = {gene: matched_model.score_gene(gene, donor_terms)
                        for gene in candidate_genes}
        drop = suite.held_out_pairs(set(profile["references"]), by_reference, support)
        held_model = suite.build_model(parsed, drop)
        held_terms = held_model.normalize(profile["terms"])
        held_scores = {gene: held_model.score_gene(gene, held_terms) for gene in candidate_genes}

        common = {"scenario_id": scenario_id, "group": scenario["group"], "gene": scenario["gene"]}
        arm_bonus = {
            "variant_only": None,
            "phenotype_matched": {g: args.weight * v for g, v in matched_scores.items()},
            "phenotype_mismatched": {g: args.weight * v for g, v in donor_scores.items()},
            "phenotype_source_held_out": {g: args.weight * v for g, v in held_scores.items()},
        }
        for arm in ARMS:
            outcome = suite.score_scenario(records, keys, arm_bonus[arm],
                                           args.out_dir / arm / f"{scenario_id}.csv")
            results[arm].append(common | outcome)

        per_scenario_layer.append({
            "scenario_id": scenario_id, "gene": scenario["gene"], "group": scenario["group"],
            "own_disease": profile["disease"], "donor_scenario": donor_id,
            "donor_gene": genes[donor_id], "donor_disease": donor_profile["disease"],
            "donor_is_unrelated": unrelated(scenario["gene"], genes[donor_id], sets),
            "matched_score_of_the_scenario_gene": round(matched_scores.get(scenario["gene"], 0.0), 4),
            "mismatched_score_of_the_scenario_gene": round(donor_scores.get(scenario["gene"], 0.0), 4),
            "held_out_score_of_the_scenario_gene": round(held_scores.get(scenario["gene"], 0.0), 4),
        })
        print(f"{scenario_id} done", file=sys.stderr)

    metrics = {arm: suite.arm_metrics(results[arm]) for arm in ARMS}
    intervals = {arm: stats.arm_intervals(results[arm]) for arm in ARMS}
    gate = gate_verdict(metrics)
    receipt = {
        "receipt_id": "t1_phenotype_mismatch_v1",
        "created_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "track": "t1_rank", "queue_item": "R21", "data_class": "A",
        "defect": "R2-3 of reports/jury/jury_review_track1_2026-09-17T154227Z.md",
        "gate_record": "docs/benchmarks/t1_phenotype_mismatch_v1.md",
        "receipt_scope": "A GIAB background, synthetic spike-in variants and simulated case profiles "
                         "from the HPO release. No patient term and no patient file enters this run.",
        "hpo": {"release": args.hpo_release,
                "hp_obo_sha256": phenotype_layer.sha256_of(obo),
                "genes_to_phenotype_sha256": phenotype_layer.sha256_of(genes_to_phenotype),
                "genes_to_disease_sha256": phenotype_layer.sha256_of(genes_to_disease),
                "phenotype_hpoa_sha256": phenotype_layer.sha256_of(hpoa)},
        "suite": {"name": manifest["suite"], "scenarios": len(manifest["scenarios"]),
                  "manifest_sha256": phenotype_layer.sha256_of(args.manifest)},
        "scorer": "hackathon_evaluation.py, Space commit 1c761cc23d90aebe6a011fd5b0b99517df42408c",
        "weight_under_test": args.weight,
        "derangement": {
            "seed": DERANGEMENT_SEED, "attempts": DERANGEMENT_ATTEMPTS,
            "constraint": "The donor gene and the receiver gene share no Mendelian disease in "
                          "genes_to_disease.txt of the same HPO release.",
            "scenarios": len(order),
            "scenarios_placed_under_the_constraint": len(order) - len(not_placed),
            "scenarios_not_placed": not_placed,
            "scenarios_that_kept_their_own_profile": sorted(
                key for key, value in assignment.items() if key == value)},
        "phenotype_scores": {
            "mean_matched_score_of_the_scenario_gene": round(statistics.mean(
                item["matched_score_of_the_scenario_gene"] for item in per_scenario_layer), 4),
            "mean_mismatched_score_of_the_scenario_gene": round(statistics.mean(
                item["mismatched_score_of_the_scenario_gene"] for item in per_scenario_layer), 4),
            "mean_held_out_score_of_the_scenario_gene": round(statistics.mean(
                item["held_out_score_of_the_scenario_gene"] for item in per_scenario_layer), 4),
            "scenarios_where_the_mismatched_score_is_at_least_the_matched_score": sum(
                1 for item in per_scenario_layer
                if item["mismatched_score_of_the_scenario_gene"]
                >= item["matched_score_of_the_scenario_gene"])},
        "arms": metrics,
        "arm_intervals": intervals,
        "paired_tests": {
            "phenotype_matched_against_variant_only":
                suite.paired_test(results["phenotype_matched"], results["variant_only"]),
            "phenotype_mismatched_against_variant_only":
                suite.paired_test(results["phenotype_mismatched"], results["variant_only"]),
            "phenotype_matched_against_phenotype_mismatched":
                suite.paired_test(results["phenotype_matched"], results["phenotype_mismatched"]),
            "phenotype_source_held_out_against_variant_only":
                suite.paired_test(results["phenotype_source_held_out"], results["variant_only"]),
        },
        "decision_gate": gate,
        "per_scenario_phenotype": per_scenario_layer,
        "per_scenario": results,
        "verdict": "MEASURED",
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"arms": metrics, "phenotype_scores": receipt["phenotype_scores"],
                      "paired_tests": receipt["paired_tests"], "gate": gate,
                      "derangement": {key: value for key, value in receipt["derangement"].items()
                                      if key != "constraint"}}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
