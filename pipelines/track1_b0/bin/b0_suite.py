#!/usr/bin/env python3
"""Run B0 on every spike-in scenario, score each file with the official scorer and apply the gate.

Two arms run: `b0` with the literature allele list (the arm that ships) and `b0_no_literature`.
The gate applies to `b0` only. The floors come from the command line and enter the receipt.

Usage:
    b0_suite.py --background bg.vep.vcf.gz --spikes spikes.vep.vcf.gz --manifest manifest.json
                --literature known_alleles.tsv --min-mean-rank-points 25 --min-pair-top10 0.5 --out-dir out
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
import sys
from pathlib import Path

import b0_rank
import hackathon_evaluation as evaluation


def spike_id(variant: dict) -> str:
    """Return the VCF ID that build_spikein_suite.py gives a spiked variant.

    Args:
        variant: A manifest variant.

    Returns:
        The identifier.
    """
    return f"{variant['chrom']}_{variant['pos']}_{variant['ref']}_{variant['alt']}"


def run_arm(name: str, background: list[dict], spikes: list[dict], manifest: dict, out_dir: Path) -> dict:
    """Rank and score every scenario for one arm.

    Args:
        name: The arm name.
        background: Candidate records of the background.
        spikes: Candidate records of all spiked variants.
        manifest: The suite manifest.
        out_dir: The output directory.

    Returns:
        The arm metrics and per-scenario results.
    """
    arm_dir = out_dir / name
    arm_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for scenario in manifest["scenarios"]:
        keys = {(v["chrom"], int(v["pos"]), v["ref"].upper(), v["alt"].upper()) for v in scenario["variants"]}
        records = background + [record for record in spikes if record["key"] in keys]
        rows = b0_rank.rows_for(b0_rank.hypotheses(records), "PROBAND01")
        path = arm_dir / f"{scenario['scenario_id']}.csv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=b0_rank.HEADER)
            writer.writeheader()
            writer.writerows(rows)
        submission = evaluation.load_submission(str(path))
        score = evaluation.score_proband("PROBAND01", submission["PROBAND01"], frozenset(keys))
        results.append({"scenario_id": scenario["scenario_id"], "group": scenario["group"], "gene": scenario["gene"],
                        "rows": len(rows), "full_match_rank": score.full_match_rank,
                        "partial_match_rank": score.partial_match_rank, "rank_points": score.rank_points,
                        "f_max": round(score.f_max, 4)})
    by_group: dict[str, list[dict]] = {}
    for result in results:
        by_group.setdefault(result["group"], []).append(result)
    metrics = {
        "scenarios": len(results),
        "mean_rank_points": round(statistics.mean(r["rank_points"] for r in results), 2),
        "pair_top1": round(sum(r["full_match_rank"] == 1 for r in results) / len(results), 3),
        "pair_top10": round(sum(r["full_match_rank"] is not None for r in results) / len(results), 3),
        "any_variant_top10": round(sum(r["full_match_rank"] is not None or r["partial_match_rank"] is not None
                                       for r in results) / len(results), 3),
        "mean_f_max": round(statistics.mean(r["f_max"] for r in results), 3),
        "by_group": {group: {"n": len(items), "mean_rank_points": round(statistics.mean(i["rank_points"] for i in items), 2),
                             "pair_top10": round(sum(i["full_match_rank"] is not None for i in items) / len(items), 3)}
                     for group, items in sorted(by_group.items())},
    }
    return {"metrics": metrics, "results": results}


def main() -> int:
    """Run both arms and write the gate receipt.

    Returns:
        Zero after the receipt is written. The verdict is in the receipt, so a failed gate still publishes it.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--background", type=Path, required=True)
    parser.add_argument("--spikes", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--literature", type=Path, required=True)
    parser.add_argument("--min-mean-rank-points", type=float, required=True)
    parser.add_argument("--min-pair-top10", type=float, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    ids = {spike_id(v) for s in manifest["scenarios"] for v in s["variants"]}
    literature = b0_rank.load_literature(args.literature)
    arms = {}
    for name, lit in (("b0", literature), ("b0_no_literature", {})):
        background = b0_rank.parse_vcf(args.background, lit)
        spikes = b0_rank.parse_vcf(args.spikes, lit, ids)
        arms[name] = run_arm(name, background, spikes, manifest, args.out_dir)
        arms[name]["background_candidates"] = len(background)
    shipped = arms["b0"]["metrics"]
    checks = {
        "serializer": all(r["rows"] <= 10 for r in arms["b0"]["results"]),
        "mean_rank_points": shipped["mean_rank_points"] >= args.min_mean_rank_points,
        "pair_top10": shipped["pair_top10"] >= args.min_pair_top10,
    }
    receipt = {
        "receipt_id": "s4_b0_spikein_v1", "data_class": "A", "suite": manifest["suite"],
        "manifest_sha256": hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
        "scorer": "hackathon_evaluation.py, Space commit 1c761cc23d90aebe6a011fd5b0b99517df42408c",
        "floors": {"min_mean_rank_points": args.min_mean_rank_points, "min_pair_top10": args.min_pair_top10},
        "checks": checks, "verdict": "PASS" if all(checks.values()) else "FAIL",
        "arms": {name: {"metrics": arm["metrics"], "background_candidates": arm["background_candidates"]}
                 for name, arm in arms.items()},
        "per_scenario": {name: arm["results"] for name, arm in arms.items()},
    }
    (args.out_dir / "s4_b0_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"verdict": receipt["verdict"], "checks": checks, "b0": shipped,
                      "b0_no_literature": arms["b0_no_literature"]["metrics"]}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
