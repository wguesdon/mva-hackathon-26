"""Measure the B0 ranker on the spike-in suite with each spike-in publication removed. Data class A.

Defect T1-6 of the jury review of 2026-09-17 states that the spike-in suite inserts alleles that the
literature table also holds, so the ablation of the table measures the retrieval of a known entry and
not generalization to an allele that the table does not hold. This script repeats the suite with the
publication that reports each spiked allele removed from the table for that scenario. The decision
gate is in `docs/benchmarks/t1_analysis_gaps_v1.md`.

Three arms run on the same candidate records.

    full_table                   the shipped ranker, with the whole literature table
    leave_one_publication_out    for each scenario, every allele of the table whose PMID matches the
                                 PMID of an allele that the scenario spikes is removed
    no_literature                no literature table at all

Every input is public. The background is GIAB HG002 and the script opens no patient file.

Usage:
    uv run python scripts/track1_leave_one_publication_out.py --background bg.vep.vcf.gz \
        --spikes spikes.vep.vcf.gz --manifest reports/benchmarks/spikein_v1/manifest.json \
        --literature reports/known_alleles/known_alleles_grch38_2026_09_15.tsv \
        --out-dir <dir> --receipt tracking/receipts/t1_leave_one_publication_out_v1.json
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "pipelines" / "track1_b0" / "bin"
sys.path.insert(0, str(BIN))

import b0_rank  # noqa: E402
import hackathon_evaluation as evaluation  # noqa: E402

ARMS = ("full_table", "leave_one_publication_out", "no_literature")


def sha256_of(path: Path) -> str:
    """Hash a file.

    Args:
        path: The file.

    Returns:
        The hexadecimal SHA-256.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def read_table(path: Path) -> list[dict]:
    """Read the resolved rows of the literature allele table.

    Args:
        path: The table.

    Returns:
        The rows that carry a position.
    """
    with path.open(newline="") as handle:
        return [row for row in csv.DictReader(handle, delimiter="\t")
                if row["status"] == "resolved" and row["pos"]]


def table_keys(rows: list[dict]) -> dict[tuple[str, int, str, str], str]:
    """Convert table rows to the scorer tuples that the ranker uses.

    Args:
        rows: The resolved table rows.

    Returns:
        (chrom, pos, ref, alt) with a chr prefix, mapped to the gene symbol.
    """
    keys = {}
    for row in rows:
        key = ("chr" + row["chrom"].removeprefix("chr"), int(row["pos"]),
               row["ref"].upper(), row["alt"].upper())
        keys[key] = row["gene"]
    return keys


def scenario_publications(scenario: dict, rows: list[dict]) -> list[str]:
    """Find the publications of the table that report an allele of one scenario.

    Args:
        scenario: A suite scenario.
        rows: The resolved table rows.

    Returns:
        The PMIDs to remove for that scenario, sorted.
    """
    spiked = {("chr" + v["chrom"].removeprefix("chr"), int(v["pos"]),
               v["ref"].upper(), v["alt"].upper()) for v in scenario["variants"]}
    pmids = set()
    for row in rows:
        key = ("chr" + row["chrom"].removeprefix("chr"), int(row["pos"]),
               row["ref"].upper(), row["alt"].upper())
        if key in spiked and row["pmid"]:
            pmids.add(row["pmid"])
    return sorted(pmids)


def reduced_keys(rows: list[dict], pmids: set[str]) -> dict[tuple[str, int, str, str], str]:
    """Build the literature key map with every allele of the given publications removed.

    Args:
        rows: The resolved table rows.
        pmids: The PMIDs to remove.

    Returns:
        The reduced key map.
    """
    return table_keys([row for row in rows if row["pmid"] not in pmids])


def score_scenario(records: list[dict], keys: set[tuple], out_path: Path) -> dict:
    """Rank one scenario, write its submission file and score it with the official scorer.

    Args:
        records: The candidate variant records of the scenario.
        keys: The truth variants of the scenario.
        out_path: Where to write the submission file.

    Returns:
        The rank points, the F-max and the two match ranks.
    """
    rows = b0_rank.rows_for(b0_rank.hypotheses(records), "PROBAND01")
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


def paired_test(first: list[dict], second: list[dict], scope: set[str] | None = None) -> dict:
    """Compare two arms on the rank 1 outcome of every scenario of a scope.

    Args:
        first: The per-scenario results of the first arm.
        second: The per-scenario results of the second arm.
        scope: The scenario identifiers to include, or None for every scenario.

    Returns:
        The exact McNemar result with the number of scenarios compared.
    """
    first_map = {r["scenario_id"]: r["full_match_rank"] == 1 for r in first}
    second_map = {r["scenario_id"]: r["full_match_rank"] == 1 for r in second}
    keys = sorted(set(first_map) & set(second_map) & (scope if scope is not None else set(first_map)))
    only_first = sum(1 for key in keys if first_map[key] and not second_map[key])
    only_second = sum(1 for key in keys if second_map[key] and not first_map[key])
    return mcnemar_exact(only_first, only_second) | {"scenarios_compared": len(keys)}


def subset_metrics(results: list[dict], scope: set[str]) -> dict:
    """Summarize the results of a subset of scenarios.

    Args:
        results: The per-scenario results of one arm.
        scope: The scenario identifiers to include.

    Returns:
        The mean rank points and the mean F-max of the subset.
    """
    items = [r for r in results if r["scenario_id"] in scope]
    if not items:
        return {"scenarios": 0, "mean_rank_points": None, "mean_f_max": None}
    return {"scenarios": len(items),
            "mean_rank_points": round(statistics.mean(r["rank_points"] for r in items), 2),
            "mean_f_max": round(statistics.mean(r["f_max"] for r in items), 3)}


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
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    rows = read_table(args.literature)
    full = table_keys(rows)
    spike_ids = {f"{v['chrom']}_{v['pos']}_{v['ref']}_{v['alt']}"
                 for s in manifest["scenarios"] for v in s["variants"]}

    background_full = b0_rank.parse_vcf(args.background, full)
    background_none = b0_rank.parse_vcf(args.background, {})
    background_insensitive = background_full == background_none
    spikes_by_table: dict[tuple[str, ...], list[dict]] = {}

    per_scenario_tables = {}
    for scenario in manifest["scenarios"]:
        per_scenario_tables[scenario["scenario_id"]] = scenario_publications(scenario, rows)
    affected = {key for key, pmids in per_scenario_tables.items() if pmids}

    for arm in ARMS:
        (args.out_dir / arm).mkdir(parents=True, exist_ok=True)
    results: dict[str, list[dict]] = {arm: [] for arm in ARMS}

    for scenario in manifest["scenarios"]:
        scenario_id = scenario["scenario_id"]
        keys = {(v["chrom"], int(v["pos"]), v["ref"].upper(), v["alt"].upper())
                for v in scenario["variants"]}
        common = {"scenario_id": scenario_id, "group": scenario["group"], "gene": scenario["gene"],
                  "publications_removed": per_scenario_tables[scenario_id]}
        for arm in ARMS:
            if arm == "full_table":
                literature, background = full, background_full
            elif arm == "no_literature":
                literature, background = {}, background_none
            else:
                removed = set(per_scenario_tables[scenario_id])
                literature = reduced_keys(rows, removed)
                if not removed or background_insensitive:
                    background = background_full
                else:
                    background = b0_rank.parse_vcf(args.background, literature)
            signature = tuple(sorted(f"{k[0]}:{k[1]}:{k[2]}:{k[3]}" for k in literature))
            if signature not in spikes_by_table:
                spikes_by_table[signature] = b0_rank.parse_vcf(args.spikes, literature, spike_ids)
            spikes = spikes_by_table[signature]
            records = background + [record for record in spikes if record["key"] in keys]
            result = score_scenario(records, keys, args.out_dir / arm / f"{scenario_id}.csv")
            results[arm].append(common | result)
        print(f"{scenario_id} done", file=sys.stderr)

    metrics = {arm: arm_metrics(results[arm]) for arm in ARMS}
    affected_metrics = {arm: subset_metrics(results[arm], affected) for arm in ARMS}
    gate_holds = (affected_metrics["leave_one_publication_out"]["mean_rank_points"]
                  > affected_metrics["no_literature"]["mean_rank_points"])
    receipt = {
        "receipt_id": "t1_leave_one_publication_out_v1",
        "created_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "track": "t1_rank", "queue_item": "T1-6", "data_class": "A",
        "receipt_scope": "A GIAB background and synthetic spike-in variants. No patient file enters "
                         "this run.",
        "gate_record": "docs/benchmarks/t1_analysis_gaps_v1.md, section T1-6",
        "suite": {"name": manifest["suite"], "scenarios": len(manifest["scenarios"]),
                  "manifest_sha256": sha256_of(args.manifest)},
        "literature_table": {"path": str(args.literature),
                             "sha256": sha256_of(args.literature),
                             "resolved_alleles": len(rows),
                             "publications": sorted({row["pmid"] for row in rows if row["pmid"]})},
        "scorer": "hackathon_evaluation.py, Space commit 1c761cc23d90aebe6a011fd5b0b99517df42408c",
        "background_is_insensitive_to_the_table": background_insensitive,
        "scenarios_that_spike_a_table_allele": sorted(affected),
        "publications_removed_per_scenario": {key: value for key, value in per_scenario_tables.items() if value},
        "arms_all_scenarios": metrics,
        "arms_affected_scenarios": affected_metrics,
        "paired_tests": {
            "leave_one_publication_out_against_full_table_all":
                paired_test(results["leave_one_publication_out"], results["full_table"]),
            "leave_one_publication_out_against_full_table_affected":
                paired_test(results["leave_one_publication_out"], results["full_table"], affected),
            "leave_one_publication_out_against_no_literature_affected":
                paired_test(results["leave_one_publication_out"], results["no_literature"], affected),
            "full_table_against_no_literature_all":
                paired_test(results["full_table"], results["no_literature"]),
        },
        "decision_gate": {
            "condition": "The leave-one-publication-out arm scores strictly more mean rank points "
                         "than the arm with no literature table, over the affected scenarios.",
            "leave_one_publication_out": affected_metrics["leave_one_publication_out"]["mean_rank_points"],
            "no_literature": affected_metrics["no_literature"]["mean_rank_points"],
            "outcome": "GATE MET" if gate_holds else "GATE NOT MET",
            "reading": "GATE MET means that the table generalizes to an allele that it does not "
                       "hold. GATE NOT MET means that the table acts as a lookup of a known entry."},
        "per_scenario": results,
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"receipt": str(args.receipt), "arms_all": metrics,
                      "arms_affected": affected_metrics,
                      "gate": receipt["decision_gate"]["outcome"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
