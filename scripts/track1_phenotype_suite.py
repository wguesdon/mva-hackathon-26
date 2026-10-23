"""Validate the phenotype layer of the B0 ranker on the public spike-in suite. Data class A.

Defect T1-3 of the jury review of 2026-09-17 states that the Track 1 method does not use the
phenotype, which the task names as an input. This script measures the phenotype layer of
`pipelines/track1_b0/bin/phenotype_layer.py` on the 50 scenarios of the suite `spikein_v1`, before
the layer touches the patient path. The decision gate is in `docs/benchmarks/t1_analysis_gaps_v1.md`.

Every input is public. The background is GIAB HG002, the spiked alleles come from ClinVar and from
the resolved literature table, and every case phenotype profile is built from the Human Phenotype
Ontology release alone. No term of the child enters any profile and the script opens no patient file.

Three arms run on the same candidate records.

    variant_only                    the shipped ranker, phenotype weight 0.0
    phenotype_matched_corpus        the gene profile comes from the whole HPO release
    phenotype_source_held_out       every annotation row that produced a term of the case profile is
                                    removed from the gene profiles before the model is built

Usage:
    uv run python scripts/track1_phenotype_suite.py --background bg.vep.vcf.gz \
        --spikes spikes.vep.vcf.gz --manifest reports/benchmarks/spikein_v1/manifest.json \
        --literature reports/known_alleles/known_alleles_grch38_2026_09_15.tsv \
        --hpo-dir tracking/cache/hpo_v2026_09_01 --out-dir reports/benchmarks/t1_phenotype_v1 \
        --receipt tracking/receipts/t1_phenotype_layer_v1.json
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import random
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "pipelines" / "track1_b0" / "bin"
sys.path.insert(0, str(BIN))

import b0_rank  # noqa: E402
import hackathon_evaluation as evaluation  # noqa: E402
import phenotype_layer  # noqa: E402

SEED = 20260917
TERMS_PER_PROFILE = 8
IMPRECISE_TERMS = 2
NOISE_TERMS = 2
PHENOTYPE_WEIGHT = 2.0
PHENOTYPIC_ABNORMALITY = "HP:0000118"
BASELINE_MEAN_RANK_POINTS = 94.0
BASELINE_MEAN_F_MAX = 0.96


def read_gene_diseases(path: Path) -> dict[str, list[str]]:
    """Read the Mendelian disease identifiers of every gene.

    Args:
        path: The `genes_to_disease.txt` file of the HPO release.

    Returns:
        The gene symbol mapped to its Mendelian disease identifiers.
    """
    diseases: dict[str, list[str]] = defaultdict(list)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["association_type"] == "MENDELIAN":
                diseases[row["gene_symbol"]].append(row["disease_id"])
    return dict(diseases)


def read_disease_annotations(path: Path) -> dict[str, list[tuple[str, str]]]:
    """Read the phenotypic abnormality rows of every disease.

    Args:
        path: The `phenotype.hpoa` file of the HPO release.

    Returns:
        The disease identifier mapped to a list of (HPO term, reference) pairs.
    """
    rows: dict[str, list[tuple[str, str]]] = defaultdict(list)
    with path.open(newline="") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "database_id":
                continue
            if len(fields) < 11 or fields[2] == "NOT" or fields[10] != "P":
                continue
            rows[fields[0]].append((fields[3], fields[4]))
    return dict(rows)


def disease_of_gene(gene: str, gene_diseases: dict[str, list[str]],
                    annotations: dict[str, list[tuple[str, str]]]) -> str | None:
    """Pick the Mendelian disease of a gene that carries the most annotation rows.

    Args:
        gene: The gene symbol.
        gene_diseases: The Mendelian diseases of every gene.
        annotations: The phenotypic abnormality rows of every disease.

    Returns:
        The disease identifier, or None when the gene has no annotated Mendelian disease.
    """
    candidates = [disease for disease in gene_diseases.get(gene, []) if annotations.get(disease)]
    if not candidates:
        return None
    return max(sorted(set(candidates)), key=lambda disease: len(annotations[disease]))


def parent_of(term: str, parents: dict[str, set[str]]) -> str:
    """Return one parent of a term inside the phenotypic abnormality branch.

    Args:
        term: The term identifier.
        parents: The direct parents of every term.

    Returns:
        The first parent in sorted order, or the term itself when it has no usable parent.
    """
    live = sorted(parent for parent in parents.get(term, set()) if parent != "HP:0000001")
    return live[0] if live else term


def build_case_profile(gene: str, rng: random.Random, gene_diseases: dict[str, list[str]],
                       annotations: dict[str, list[tuple[str, str]]],
                       parents: dict[str, set[str]],
                       abnormality_terms: list[str]) -> dict:
    """Build one simulated case phenotype profile from public annotation rows only.

    The protocol follows the simulated patient design of the Exomiser papers. The script samples
    terms of one disease of the gene, generalizes a fixed number of them to a parent term and adds a
    fixed number of random terms as noise.

    Args:
        gene: The gene of the scenario.
        rng: The random generator, seeded by the caller.
        gene_diseases: The Mendelian diseases of every gene.
        annotations: The phenotypic abnormality rows of every disease.
        parents: The direct parents of every term.
        abnormality_terms: The terms of the phenotypic abnormality branch, for the noise draw.

    Returns:
        The profile, the disease that produced it and the annotation references that it used.
    """
    disease = disease_of_gene(gene, gene_diseases, annotations)
    if disease is None:
        return {"gene": gene, "disease": None, "terms": [], "references": [], "sampled": 0,
                "imprecise": 0, "noise": 0}
    rows = sorted(set(annotations[disease]))
    sample = rng.sample(rows, min(TERMS_PER_PROFILE, len(rows)))
    terms = [term for term, _ in sample]
    references = sorted({reference for _, reference in sample if reference})
    for index in range(min(IMPRECISE_TERMS, len(terms))):
        terms[index] = parent_of(terms[index], parents)
    noise = rng.sample(abnormality_terms, NOISE_TERMS)
    profile = []
    for term in terms + noise:
        if term not in profile:
            profile.append(term)
    return {"gene": gene, "disease": disease, "terms": profile, "references": references,
            "sampled": len(sample), "imprecise": min(IMPRECISE_TERMS, len(terms)),
            "noise": NOISE_TERMS}


def annotation_support(references: set[str], annotations: dict[str, list[tuple[str, str]]],
                       disease_genes: dict[str, set[str]]) -> tuple[dict[str, set[tuple[str, str]]],
                                                                    dict[tuple[str, str], set[str]]]:
    """Index the gene and term pairs that a set of annotation references supports.

    The index is restricted to the pairs that at least one of the given references touches, so the
    two passes over `phenotype.hpoa` run once for the whole suite.

    Args:
        references: Every annotation reference that any case profile of the suite used.
        annotations: The phenotypic abnormality rows of every disease.
        disease_genes: The genes of every disease.

    Returns:
        The pairs of each reference, and the whole reference support of each of those pairs.
    """
    by_reference: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for disease, rows in annotations.items():
        genes = disease_genes.get(disease)
        if not genes:
            continue
        for term, reference in rows:
            if reference in references:
                for gene in genes:
                    by_reference[reference].add((gene, term))
    touched = {pair for pairs in by_reference.values() for pair in pairs}
    support: dict[tuple[str, str], set[str]] = defaultdict(set)
    for disease, rows in annotations.items():
        genes = disease_genes.get(disease)
        if not genes:
            continue
        for term, reference in rows:
            for gene in genes:
                if (gene, term) in touched:
                    support[(gene, term)].add(reference)
    return dict(by_reference), dict(support)


def held_out_pairs(references: set[str], by_reference: dict[str, set[tuple[str, str]]],
                   support: dict[tuple[str, str], set[str]]) -> set[tuple[str, str]]:
    """Find the gene and term pairs whose whole annotation support lies in the held out references.

    Args:
        references: The annotation references that produced one case profile.
        by_reference: The pairs of each reference.
        support: The whole reference support of each touched pair.

    Returns:
        The (gene, term) pairs to remove from the gene profiles.
    """
    touched: set[tuple[str, str]] = set()
    for reference in references:
        touched.update(by_reference.get(reference, ()))
    return {pair for pair in touched if support.get(pair, set()) <= references}


def build_model(parsed: dict, drop: set[tuple[str, str]] | None) -> phenotype_layer.PhenotypeModel:
    """Build a phenotype model from the parsed release, with an optional removal of rows.

    Args:
        parsed: The parsed ontology and the gene profiles of the release.
        drop: The (gene, term) pairs to remove, or None.

    Returns:
        The model.
    """
    profiles = {gene: set(terms) for gene, terms in parsed["profiles"].items()}
    if drop:
        for gene, term in drop:
            if gene in profiles:
                profiles[gene].discard(term)
        profiles = {gene: terms for gene, terms in profiles.items() if terms}
    content = phenotype_layer.information_content(profiles, parsed["closure"])
    return phenotype_layer.PhenotypeModel(parsed["closure"], parsed["descendants"], content,
                                          profiles, parsed["alternates"])


def score_scenario(records: list[dict], keys: set[tuple], bonus: dict[str, float] | None,
                   out_path: Path) -> dict:
    """Rank one scenario, write its submission file and score it with the official scorer.

    Args:
        records: The candidate variant records of the scenario.
        keys: The truth variants of the scenario.
        bonus: The phenotype bonus of every gene, or None when the layer is off.
        out_path: Where to write the submission file.

    Returns:
        The rank points, the F-max and the two match ranks.
    """
    rows = b0_rank.rows_for(b0_rank.hypotheses(records, bonus), "PROBAND01")
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=b0_rank.HEADER)
        writer.writeheader()
        writer.writerows(rows)
    submission = evaluation.load_submission(str(out_path))
    score = evaluation.score_proband("PROBAND01", submission["PROBAND01"], frozenset(keys))
    return {"rows": len(rows), "full_match_rank": score.full_match_rank,
            "partial_match_rank": score.partial_match_rank,
            "rank_points": score.rank_points, "f_max": round(score.f_max, 4)}


def arm_metrics(results: list[dict]) -> dict:
    """Summarize the per-scenario results of one arm.

    Args:
        results: The per-scenario results.

    Returns:
        The mean rank points, the two hit rates and the mean F-max, overall and by group.
    """
    by_group: dict[str, list[dict]] = defaultdict(list)
    for result in results:
        by_group[result["group"]].append(result)
    return {
        "scenarios": len(results),
        "mean_rank_points": round(statistics.mean(r["rank_points"] for r in results), 2),
        "pair_top1": round(sum(r["full_match_rank"] == 1 for r in results) / len(results), 3),
        "pair_top10": round(sum(r["full_match_rank"] is not None for r in results) / len(results), 3),
        "mean_f_max": round(statistics.mean(r["f_max"] for r in results), 3),
        "by_group": {group: {"n": len(items),
                             "mean_rank_points": round(statistics.mean(i["rank_points"] for i in items), 2),
                             "mean_f_max": round(statistics.mean(i["f_max"] for i in items), 3)}
                     for group, items in sorted(by_group.items())},
    }


def mcnemar_exact(first_only: int, second_only: int) -> dict:
    """Run the two-sided exact McNemar test on two discordant counts.

    Args:
        first_only: The scenarios that only the first arm placed at rank 1.
        second_only: The scenarios that only the second arm placed at rank 1.

    Returns:
        The discordant counts, the two-sided p value and the null hypothesis.
    """
    discordant = first_only + second_only
    if discordant == 0:
        probability = 1.0
    else:
        smaller = min(first_only, second_only)
        tail = sum(math.comb(discordant, index) for index in range(smaller + 1)) * 0.5 ** discordant
        probability = min(1.0, 2.0 * tail)
    return {"first_only": first_only, "second_only": second_only, "discordant": discordant,
            "p_two_sided": round(probability, 4),
            "null": "A discordant scenario is placed at rank 1 by either arm with probability 0.5.",
            "smallest_reachable_p": round(min(1.0, 2.0 * 0.5 ** discordant), 4) if discordant else None}


def paired_test(first: list[dict], second: list[dict]) -> dict:
    """Compare two arms on the rank 1 outcome of every scenario.

    Args:
        first: The per-scenario results of the first arm.
        second: The per-scenario results of the second arm.

    Returns:
        The exact McNemar result.
    """
    first_map = {result["scenario_id"]: result["full_match_rank"] == 1 for result in first}
    second_map = {result["scenario_id"]: result["full_match_rank"] == 1 for result in second}
    keys = sorted(set(first_map) & set(second_map))
    only_first = sum(1 for key in keys if first_map[key] and not second_map[key])
    only_second = sum(1 for key in keys if second_map[key] and not first_map[key])
    return mcnemar_exact(only_first, only_second)


def main() -> int:
    """Run the three arms and write the receipt.

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
    parser.add_argument("--weight", type=float, default=PHENOTYPE_WEIGHT)
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
    abnormality_terms = sorted(term for term in descendants[PHENOTYPIC_ABNORMALITY]
                               if term != PHENOTYPIC_ABNORMALITY)

    gene_diseases = read_gene_diseases(genes_to_disease)
    annotations = read_disease_annotations(hpoa)
    disease_genes: dict[str, set[str]] = defaultdict(set)
    for gene, diseases in gene_diseases.items():
        for disease in diseases:
            disease_genes[disease].add(gene)

    manifest = json.loads(args.manifest.read_text())
    literature = b0_rank.load_literature(args.literature)
    spike_ids = {f"{v['chrom']}_{v['pos']}_{v['ref']}_{v['alt']}"
                 for s in manifest["scenarios"] for v in s["variants"]}
    background = b0_rank.parse_vcf(args.background, literature)
    spikes = b0_rank.parse_vcf(args.spikes, literature, spike_ids)

    rng = random.Random(SEED)
    case_profiles = {}
    for scenario in manifest["scenarios"]:
        case_profiles[scenario["scenario_id"]] = build_case_profile(
            scenario["gene"], rng, gene_diseases, annotations, parents, abnormality_terms)

    all_references = {reference for profile in case_profiles.values()
                      for reference in profile["references"]}
    by_reference, support = annotation_support(all_references, annotations, disease_genes)

    matched_model = build_model(parsed, None)
    arms: dict[str, list[dict]] = {"variant_only": [], "phenotype_matched_corpus": [],
                                   "phenotype_source_held_out": []}
    per_scenario_layer = []
    for arm in arms:
        (args.out_dir / arm).mkdir(parents=True, exist_ok=True)

    for scenario in manifest["scenarios"]:
        scenario_id = scenario["scenario_id"]
        keys = {(v["chrom"], int(v["pos"]), v["ref"].upper(), v["alt"].upper())
                for v in scenario["variants"]}
        records = background + [record for record in spikes if record["key"] in keys]
        genes = sorted({record["gene"] for record in records})
        profile = case_profiles[scenario_id]
        terms = matched_model.normalize(profile["terms"])

        base = score_scenario(records, keys, None, args.out_dir / "variant_only" / f"{scenario_id}.csv")
        matched_scores = {gene: matched_model.score_gene(gene, terms) for gene in genes}
        matched = score_scenario(records, keys,
                                 {gene: args.weight * value for gene, value in matched_scores.items()},
                                 args.out_dir / "phenotype_matched_corpus" / f"{scenario_id}.csv")
        drop = held_out_pairs(set(profile["references"]), by_reference, support)
        held_model = build_model(parsed, drop)
        held_terms = held_model.normalize(profile["terms"])
        held_scores = {gene: held_model.score_gene(gene, held_terms) for gene in genes}
        held = score_scenario(records, keys,
                              {gene: args.weight * value for gene, value in held_scores.items()},
                              args.out_dir / "phenotype_source_held_out" / f"{scenario_id}.csv")

        common = {"scenario_id": scenario_id, "group": scenario["group"], "gene": scenario["gene"]}
        arms["variant_only"].append(common | base)
        arms["phenotype_matched_corpus"].append(common | matched)
        arms["phenotype_source_held_out"].append(common | held)
        per_scenario_layer.append({
            "scenario_id": scenario_id, "gene": scenario["gene"], "group": scenario["group"],
            "disease": profile["disease"], "case_terms": len(terms),
            "annotation_references": len(profile["references"]),
            "gene_term_pairs_removed": len(drop),
            "phenotype_score_of_the_scenario_gene": round(matched_scores.get(scenario["gene"], 0.0), 4),
            "phenotype_score_source_held_out": round(held_scores.get(scenario["gene"], 0.0), 4),
            "max_phenotype_score_of_any_candidate_gene": round(max(matched_scores.values(), default=0.0), 4),
        })
        print(f"{scenario_id} done", file=sys.stderr)

    metrics = {arm: arm_metrics(results) for arm, results in arms.items()}
    baseline_reproduced = (metrics["variant_only"]["mean_rank_points"] == BASELINE_MEAN_RANK_POINTS
                           and metrics["variant_only"]["mean_f_max"] == BASELINE_MEAN_F_MAX)
    gate = {
        "record": "docs/benchmarks/t1_analysis_gaps_v1.md, section T1-3",
        "condition_1_mean_rank_points": {
            "with_the_layer": metrics["phenotype_matched_corpus"]["mean_rank_points"],
            "without_the_layer": metrics["variant_only"]["mean_rank_points"],
            "holds": metrics["phenotype_matched_corpus"]["mean_rank_points"]
                     >= metrics["variant_only"]["mean_rank_points"]},
        "condition_2_mean_f_max": {
            "with_the_layer": metrics["phenotype_matched_corpus"]["mean_f_max"],
            "without_the_layer": metrics["variant_only"]["mean_f_max"],
            "holds": metrics["phenotype_matched_corpus"]["mean_f_max"]
                     >= metrics["variant_only"]["mean_f_max"]},
    }
    gate["outcome"] = "GATE MET" if gate["condition_1_mean_rank_points"]["holds"] and \
        gate["condition_2_mean_f_max"]["holds"] else "GATE NOT MET"
    receipt = {
        "receipt_id": "t1_phenotype_layer_v1",
        "created_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "track": "t1_rank", "queue_item": "T1-3", "data_class": "A",
        "receipt_scope": "A GIAB background, synthetic spike-in variants and simulated case profiles "
                         "from the HPO release. No patient term and no patient file enters this run.",
        "hpo": {"release": args.hpo_release,
                "url": "https://github.com/obophenotype/human-phenotype-ontology/releases/tag/"
                       + args.hpo_release,
                "hp_obo_sha256": phenotype_layer.sha256_of(obo),
                "genes_to_phenotype_sha256": phenotype_layer.sha256_of(genes_to_phenotype),
                "genes_to_disease_sha256": phenotype_layer.sha256_of(genes_to_disease),
                "phenotype_hpoa_sha256": phenotype_layer.sha256_of(hpoa),
                "terms_in_the_release": len(closure),
                "genes_with_a_phenotype_profile": len(profiles),
                "max_information_content_nats": round(matched_model.max_content, 4)},
        "similarity": {
            "term_measure": "Resnik information content of the most informative common ancestor",
            "term_measure_source": "Resnik 1995, IJCAI proceedings, pages 448 to 453",
            "set_measure": "best match average",
            "set_measure_source": "Pesquita and colleagues 2009, PMID 19649320",
            "used_by": ["Phenomizer, PMID 19800049", "Exomiser, PMID 26562621"],
            "normalization": "The best match average is divided by the largest information content "
                             "of the corpus, so the score lies between 0 and 1.",
            "combination_with_the_variant_score":
                "hypothesis score = variant score + weight * phenotype score of the gene",
            "weight": args.weight,
            "default_weight_of_the_ranker": 0.0},
        "case_profiles": {
            "source": "phenotype.hpoa and genes_to_disease.txt of the same release",
            "protocol": "Sample terms of the Mendelian disease of the gene that carries the most "
                        "annotation rows, generalize a fixed number to a parent term, add a fixed "
                        "number of random terms of the phenotypic abnormality branch.",
            "protocol_source": "the simulated patient design of the Exomiser papers, PMID 26562621",
            "seed": SEED, "terms_sampled": TERMS_PER_PROFILE, "imprecise_terms": IMPRECISE_TERMS,
            "noise_terms": NOISE_TERMS,
            "no_patient_term": True,
            "scenarios_with_a_profile": sum(1 for item in per_scenario_layer if item["case_terms"]),
        },
        "suite": {"name": manifest["suite"], "scenarios": len(manifest["scenarios"]),
                  "manifest_sha256": phenotype_layer.sha256_of(args.manifest),
                  "background_candidates": len(background)},
        "scorer": "hackathon_evaluation.py, Space commit 1c761cc23d90aebe6a011fd5b0b99517df42408c",
        "baseline_reproduced": {"expected_mean_rank_points": BASELINE_MEAN_RANK_POINTS,
                                "expected_mean_f_max": BASELINE_MEAN_F_MAX,
                                "observed_mean_rank_points": metrics["variant_only"]["mean_rank_points"],
                                "observed_mean_f_max": metrics["variant_only"]["mean_f_max"],
                                "agrees": baseline_reproduced,
                                "source": "tracking/receipts/s4_b0_suite_v1.json, field arms.b0"},
        "arms": metrics,
        "paired_tests": {
            "phenotype_matched_corpus_against_variant_only":
                paired_test(arms["phenotype_matched_corpus"], arms["variant_only"]),
            "phenotype_source_held_out_against_variant_only":
                paired_test(arms["phenotype_source_held_out"], arms["variant_only"]),
        },
        "decision_gate": gate,
        "per_scenario_phenotype": per_scenario_layer,
        "per_scenario": arms,
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"receipt": str(args.receipt), "arms": metrics,
                      "gate": gate["outcome"], "baseline_reproduced": baseline_reproduced}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
