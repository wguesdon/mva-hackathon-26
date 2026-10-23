#!/usr/bin/env python3
"""Score the held-out time split suite with the frozen Track 1 ranker. Data class A.

The suite `timesplit_v1` holds compound heterozygous pairs of ClinVar alleles that entered the
database after the cutoff of the frozen resources. The frozen annotation resource is the Ensembl VEP
116.2 offline cache, whose header states `ClinVar="202509"`, so the cache cannot carry a
significance term for a record that ClinVar received after September 2025.

Three arms run on the same background and the same scorer.

    frozen_full_features    the shipped ranker, with the ClinVar term at 2.0 points
    frozen_clinvar_ablated  the same ranker with the ClinVar term at 0.0 points
    reference_spikein_v1    the shipped ranker on the suite spikein_v1, as a reproduction check

The design and the decision gates are in `docs/benchmarks/t1_time_split_v1.md`, recorded before the
result.

Every input is public. The background is GIAB HG002 and the script opens no patient file.

Usage:
    uv run python scripts/track1_time_split_benchmark.py --background bg.vep.vcf.gz \\
        --heldout-spikes heldout.vep.vcf.gz --heldout-manifest reports/benchmarks/timesplit_v1/manifest.json \\
        --reference-spikes spikes.vep.vcf.gz --reference-manifest reports/benchmarks/spikein_v1/manifest.json \\
        --literature reports/known_alleles/known_alleles_grch38_2026_09_15.tsv \\
        --out-dir <dir> --receipt tracking/receipts/t1_time_split_v1.json
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "pipelines" / "track1_b0" / "bin"
sys.path.insert(0, str(BIN))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import b0_rank  # noqa: E402
import track1_interval_stats as stats  # noqa: E402
import track1_leave_one_publication_out as lopo  # noqa: E402

ARM_WEIGHTS = {"frozen_full_features": 2.0, "frozen_clinvar_ablated": 0.0}
FLOOR_RANK_POINTS = 25.0
FLOOR_PAIR_TOP10 = 0.5
SHIPPED_SUITE = {"mean_rank_points": 94.0, "mean_f_max": 0.96, "pair_top1": 0.88}
ABLATED_SUITE = {"mean_rank_points": 48.8, "mean_f_max": 0.648, "pair_top1": 0.08}


def spike_identifiers(manifest: dict) -> set[str]:
    """Return the VCF identifiers of every spiked allele of a suite.

    Args:
        manifest: The suite manifest.

    Returns:
        The identifiers.
    """
    return {f"{v['chrom']}_{v['pos']}_{v['ref']}_{v['alt']}"
            for scenario in manifest["scenarios"] for v in scenario["variants"]}


def clinvar_terms(path: Path, keep_ids: set[str]) -> dict[tuple, str]:
    """Read the CLIN_SIG value that the frozen cache gave each spiked allele.

    Args:
        path: The annotated spike VCF.
        keep_ids: The identifiers of the spiked alleles.

    Returns:
        The key (chrom, pos, ref, alt) mapped to its CLIN_SIG text.
    """
    fields: list[str] = []
    terms: dict[tuple, str] = {}
    with b0_rank.open_text(path) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                fields = b0_rank.csq_format(line)
                continue
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt = columns[:5]
            if vid not in keep_ids:
                continue
            csq = next((item[4:] for item in columns[7].split(";") if item.startswith("CSQ=")), "")
            values = [dict(zip(fields, entry.split("|"))).get("CLIN_SIG", "")
                      for entry in csq.split(",")]
            key = (chrom if chrom.startswith("chr") else "chr" + chrom, int(pos), ref.upper(),
                   alt.upper())
            terms[key] = ";".join(sorted({value for value in values if value}))
    return terms


def label_audit(terms: dict[tuple, str]) -> dict:
    """Count the held-out alleles that the frozen cache still calls pathogenic.

    Args:
        terms: The CLIN_SIG text of every spiked allele.

    Returns:
        The audit of the time split.
    """
    pathogenic = {key: value for key, value in terms.items()
                  if b0_rank.clinvar_points(value) > 0.0}
    return {"spiked_alleles": len(terms),
            "alleles_with_any_clinvar_term_in_the_frozen_cache":
                sum(1 for value in terms.values() if value),
            "alleles_that_the_frozen_cache_calls_pathogenic": len(pathogenic),
            "examples": sorted(f"{key[0]}:{key[1]}" for key in pathogenic)[:10],
            "reading": "The time split holds when the frozen cache calls no held-out allele "
                       "pathogenic. Every allele that it does call pathogenic is a record that "
                       "ClinVar held under another identifier before the cutoff."}


def run_arm(weight: float, background_path: Path, spikes_path: Path, literature: dict,
            spike_ids: set[str], manifest: dict, out_dir: Path) -> list[dict]:
    """Score every scenario of one suite under one ClinVar weight.

    Args:
        weight: The points that a pathogenic ClinVar term adds.
        background_path: The annotated HG002 background VCF.
        spikes_path: The annotated spike VCF.
        literature: The literature allele key map.
        spike_ids: The identifiers of the spike records.
        manifest: The suite manifest.
        out_dir: The directory for the submission file of each scenario.

    Returns:
        The per-scenario result of the arm.
    """
    b0_rank.CLINVAR_WEIGHT = weight
    background = b0_rank.parse_vcf(background_path, literature)
    spikes = b0_rank.parse_vcf(spikes_path, literature, spike_ids)
    out_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for scenario in manifest["scenarios"]:
        keys = {(v["chrom"], int(v["pos"]), v["ref"].upper(), v["alt"].upper())
                for v in scenario["variants"]}
        records = background + [record for record in spikes if record["key"] in keys]
        outcome = lopo.score_scenario(records, keys, out_dir / f"{scenario['scenario_id']}.csv")
        results.append({"scenario_id": scenario["scenario_id"], "group": scenario["group"],
                        "gene": scenario["gene"]} | outcome)
    b0_rank.CLINVAR_WEIGHT = ARM_WEIGHTS["frozen_full_features"]
    return results


def subset(results: list[dict], group: str) -> list[dict]:
    """Return the results of one scenario group.

    Args:
        results: The per-scenario results of one arm.
        group: The group name.

    Returns:
        The results of that group.
    """
    return [item for item in results if item["group"] == group]


def gate_h1(metrics: dict) -> dict:
    """Evaluate the two floors of the shipped B0 gate on the held-out arm.

    Args:
        metrics: The arm metrics of the held-out full feature arm.

    Returns:
        The gate verdict.
    """
    rank_holds = metrics["mean_rank_points"] >= FLOOR_RANK_POINTS
    top10_holds = metrics["pair_top10"] >= FLOOR_PAIR_TOP10
    return {"floor_mean_rank_points": FLOOR_RANK_POINTS,
            "observed_mean_rank_points": metrics["mean_rank_points"],
            "floor_pair_top10": FLOOR_PAIR_TOP10,
            "observed_pair_top10": metrics["pair_top10"],
            "outcome": "GATE MET" if rank_holds and top10_holds else "GATE NOT MET",
            "source_of_the_floors": "tracking/receipts/s4_b0_suite_v1.json, field floors"}


def gate_h2(test: dict, audit: dict) -> dict:
    """Evaluate the label contribution on the held-out suite.

    Args:
        test: The exact McNemar result between the two held-out arms.
        audit: The label audit of the held-out spikes.

    Returns:
        The gate verdict.
    """
    no_contribution = test["discordant"] == 0
    return {"discordant_scenarios": test["discordant"], "p_two_sided": test["p_two_sided"],
            "held_out_alleles_that_the_frozen_cache_calls_pathogenic":
                audit["alleles_that_the_frozen_cache_calls_pathogenic"],
            "outcome": "NO LABEL CONTRIBUTION" if no_contribution else "LABEL STILL CONTRIBUTES",
            "reading": "NO LABEL CONTRIBUTION means that the ClinVar feature changes no rank on an "
                       "allele that the frozen resources do not know, so the held-out number is the "
                       "number that applies to an allele with no database label."}


def main() -> int:
    """Run the arms and write the receipt.

    Returns:
        Zero after the receipt is written.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--background", required=True, type=Path)
    parser.add_argument("--heldout-spikes", required=True, type=Path)
    parser.add_argument("--heldout-manifest", required=True, type=Path)
    parser.add_argument("--reference-spikes", required=True, type=Path)
    parser.add_argument("--reference-manifest", required=True, type=Path)
    parser.add_argument("--literature", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()

    heldout = json.loads(args.heldout_manifest.read_text())
    reference = json.loads(args.reference_manifest.read_text())
    rows = lopo.read_table(args.literature)
    literature = lopo.table_keys(rows)
    heldout_ids = spike_identifiers(heldout)
    reference_ids = spike_identifiers(reference)

    results = {arm: run_arm(weight, args.background, args.heldout_spikes, literature, heldout_ids,
                            heldout, args.out_dir / arm)
               for arm, weight in ARM_WEIGHTS.items()}
    results["reference_spikein_v1"] = run_arm(
        ARM_WEIGHTS["frozen_full_features"], args.background, args.reference_spikes, literature,
        reference_ids, reference, args.out_dir / "reference_spikein_v1")

    audit = label_audit(clinvar_terms(args.heldout_spikes, heldout_ids))
    clinvar_group = {arm: subset(items, "heldout_clinvar") + subset(items, "heldout_mva")
                     for arm, items in results.items() if arm != "reference_spikein_v1"}
    literature_group = {arm: subset(items, "heldout_literature")
                        for arm, items in results.items() if arm != "reference_spikein_v1"}

    metrics = {arm: lopo.arm_metrics(items) for arm, items in results.items()}
    clinvar_metrics = {arm: lopo.arm_metrics(items) for arm, items in clinvar_group.items()}
    intervals = {arm: stats.arm_intervals(items) for arm, items in results.items()}
    clinvar_intervals = {arm: stats.arm_intervals(items) for arm, items in clinvar_group.items()}
    literature_metrics = {arm: (lopo.arm_metrics(items) if items else {"scenarios": 0})
                          for arm, items in literature_group.items()}
    paired = lopo.paired_test(results["frozen_clinvar_ablated"], results["frozen_full_features"])
    reproduction = {
        "expected": SHIPPED_SUITE,
        "observed": {"mean_rank_points": metrics["reference_spikein_v1"]["mean_rank_points"],
                     "mean_f_max": metrics["reference_spikein_v1"]["mean_f_max"],
                     "pair_top1": metrics["reference_spikein_v1"]["pair_top1"]},
        "source": "tracking/receipts/s4_b0_suite_v1.json, field arms.b0"}
    reproduction["agrees"] = reproduction["expected"] == reproduction["observed"]

    drop = round(SHIPPED_SUITE["mean_rank_points"]
                 - clinvar_metrics["frozen_full_features"]["mean_rank_points"], 2)
    receipt = {
        "receipt_id": "t1_time_split_v1",
        "created_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "track": "t1_rank", "queue_item": "R22", "data_class": "A",
        "defect": "R2-1 of reports/jury/jury_review_track1_2026-09-17T154227Z.md",
        "gate_record": "docs/benchmarks/t1_time_split_v1.md",
        "receipt_scope": "A GIAB HG002 background and synthetic spike-in variants from public "
                         "ClinVar releases. No patient file enters this run.",
        "frozen_resources": {
            "ranker": "pipelines/track1_b0/bin/b0_rank.py, CLINVAR_WEIGHT 2.0, MAX_AF 0.01, "
                      "MIN_SCORE 2.0, MAX_HET_PER_GENE 4, MAX_ROWS 10, phenotype weight 0.0",
            "annotation": "Ensembl VEP 116.2 offline cache homo_sapiens_vep_116_GRCh38, SHA-256 "
                          "014c7dd9bb5ad06665866d62eb80f31ca761197bcb9b59280300676f996e600d",
            "cache_source_releases": "ClinVar 202509, dbSNP 156, GENCODE 50, gnomADe v4.1, "
                                     "gnomADg v4.1, assembly GRCh38.p14",
            "literature_table": str(args.literature.name)},
        "time_split": heldout["time_split"],
        "pool": heldout["pool"],
        "composition": heldout["composition"],
        "literature_arm_manifest": heldout["literature_arm"],
        "suites": {
            "heldout": {"name": heldout["suite"], "scenarios": len(heldout["scenarios"]),
                        "manifest_sha256": lopo.sha256_of(args.heldout_manifest)},
            "reference": {"name": reference["suite"], "scenarios": len(reference["scenarios"]),
                          "manifest_sha256": lopo.sha256_of(args.reference_manifest)}},
        "scorer": "hackathon_evaluation.py, Space commit 1c761cc23d90aebe6a011fd5b0b99517df42408c",
        "label_audit_of_the_held_out_alleles": audit,
        "reproduction_of_the_shipped_suite": reproduction,
        "arms_all_scenarios": metrics,
        "arms_heldout_clinvar_only": clinvar_metrics,
        "arms_heldout_literature_only": literature_metrics,
        "intervals_all_scenarios": intervals,
        "intervals_heldout_clinvar_only": clinvar_intervals,
        "paired_test_clinvar_ablated_against_full_features": paired,
        "comparison_with_the_shipped_suite": {
            "shipped_suite_spikein_v1": SHIPPED_SUITE,
            "clinvar_ablated_suite_spikein_v1": ABLATED_SUITE,
            "held_out_full_features": {
                "mean_rank_points": clinvar_metrics["frozen_full_features"]["mean_rank_points"],
                "mean_f_max": clinvar_metrics["frozen_full_features"]["mean_f_max"],
                "pair_top1": clinvar_metrics["frozen_full_features"]["pair_top1"]},
            "fall_in_mean_rank_points_against_the_shipped_suite": drop},
        "decision_gates": {
            "h1": gate_h1(clinvar_metrics["frozen_full_features"]),
            "h2": gate_h2(paired, audit),
            "h3": {
                "expected_on_an_allele_with_no_clinvar_pathogenic_term":
                    clinvar_metrics["frozen_full_features"]["mean_rank_points"],
                "expected_on_an_allele_with_a_clinvar_pathogenic_term":
                    SHIPPED_SUITE["mean_rank_points"],
                "record": "docs/benchmarks/t1_time_split_v1.md, gate H3"}},
        "per_scenario": results,
        "verdict": "MEASURED",
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({key: receipt[key] for key in
                      ("time_split", "pool", "composition", "label_audit_of_the_held_out_alleles",
                       "reproduction_of_the_shipped_suite", "arms_heldout_clinvar_only",
                       "arms_heldout_literature_only", "intervals_heldout_clinvar_only",
                       "paired_test_clinvar_ablated_against_full_features",
                       "comparison_with_the_shipped_suite", "decision_gates")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
