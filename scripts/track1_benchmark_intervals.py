#!/usr/bin/env python3
"""Interval estimates and paired tests for the Track 1 benchmarks. Class A.

The jury review of 2026-09-17 recorded that no Track 1 comparison carried a statistical test, a null
or an interval. This script supplies both where the committed record holds the counts that a test
needs, and it names the comparisons where the record does not.

Every count is read from a receipt under `tracking/receipts/`. No new run happens here, no patient
file is opened and no AWS job is started. The inputs are public GIAB HG002 benchmark counts and
public spike-in scenario counts, so the output is data class A.

Two estimators run.

    Wilson score interval
        A two-sided 95 percent interval for a recall. The unit is one truth record, the trial count
        is the number of truth records of the class and the success count is the number of truth
        records that the arm recovered. The interval assumes that the truth records are independent
        Bernoulli trials with one shared probability. Truth records inside one repeat array are not
        independent, so the interval is a lower bound on the width.

    McNemar exact test
        A two-sided exact binomial test on the discordant truth records of two arms that ran on the
        same truth set. The null is that a discordant truth record is equally likely to be recovered
        by either arm, which is a probability of 0.5 on the discordant records. The test conditions
        on the discordant count and it uses no normal approximation.

Usage:
    uv run python scripts/track1_benchmark_intervals.py --out tracking/receipts/t1_benchmark_intervals_v1.json
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECEIPTS = ROOT / "tracking" / "receipts"
Z_95 = 1.959963984540054


def wilson_interval(successes: int, trials: int, z: float = Z_95) -> dict:
    """Compute the Wilson score interval of a binomial proportion.

    Args:
        successes: The number of recovered truth records.
        trials: The number of truth records.
        z: The standard normal quantile of the interval. The default is the 95 percent quantile.

    Returns:
        The point estimate, the lower limit, the upper limit and the trial count. The limits are
        None when the trial count is zero.
    """
    if trials == 0:
        return {"estimate": None, "low": None, "high": None, "trials": 0, "successes": 0}
    proportion = successes / trials
    denominator = 1.0 + z * z / trials
    center = (proportion + z * z / (2 * trials)) / denominator
    spread = z / denominator * math.sqrt(
        proportion * (1.0 - proportion) / trials + z * z / (4 * trials * trials))
    return {"estimate": round(proportion, 4), "low": round(max(0.0, center - spread), 4),
            "high": round(min(1.0, center + spread), 4), "trials": trials, "successes": successes}


def mcnemar_exact(first_only: int, second_only: int) -> dict:
    """Run the two-sided exact McNemar test on two discordant counts.

    Args:
        first_only: The truth records that the first arm recovered and the second arm missed.
        second_only: The truth records that the second arm recovered and the first arm missed.

    Returns:
        The two discordant counts, the discordant total, the two-sided p value and the null.
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
            "null": "A discordant truth record is recovered by either arm with probability 0.5.",
            "smallest_reachable_p": round(min(1.0, 2.0 * 0.5 ** discordant), 4) if discordant else None}


def read_receipt(name: str) -> dict:
    """Read one receipt by its identifier.

    Args:
        name: The receipt identifier without the suffix.

    Returns:
        The receipt content.
    """
    return json.loads((RECEIPTS / f"{name}.json").read_text())


def panel_counts(receipt: dict, arm: str) -> dict:
    """Read the recovered truth record counts of one control panel arm.

    Args:
        receipt: The R13 receipt content.
        arm: The arm key inside `scopes.control_panel.arms`.

    Returns:
        The recovered count of all classes, of deletions and of insertions.
    """
    panel = receipt["scopes"]["control_panel"]["arms"][arm]
    return {"all": panel["truvari_summary"]["TP-base"],
            "DEL": panel["by_class"]["DEL"]["tp_base"],
            "INS": panel["by_class"]["INS"]["tp_base"]}


def panel_truth(receipt: dict) -> dict:
    """Read the truth record counts of the control panel by class.

    Args:
        receipt: The R13 receipt content.

    Returns:
        The truth record count of all classes, of deletions and of insertions.
    """
    panel = receipt["scopes"]["control_panel"]
    arm = panel["arms"]["pipeline_tiddit_manta_cnvnator"]
    return {"all": panel["truth_records"],
            "DEL": arm["by_class"]["DEL"]["truth_records"],
            "INS": arm["by_class"]["INS"]["truth_records"]}


def discordant_from_union(first: int, second: int, union: int) -> tuple[int, int]:
    """Derive the two discordant counts of a paired comparison from the union count.

    Two arms ran on the same truth set. The union of the recovered truth records is the count of
    truth records that at least one arm recovered, so the intersection is the sum of the two arms
    less the union, and each discordant count is the arm count less the intersection.

    Args:
        first: The truth records that the first arm recovered.
        second: The truth records that the second arm recovered.
        union: The truth records that at least one arm recovered.

    Returns:
        The records recovered by the first arm alone and by the second arm alone.

    Raises:
        ValueError: The union count is not consistent with the two arm counts.
    """
    if union > first + second or union < max(first, second):
        raise ValueError(f"A union of {union} cannot come from arms of {first} and {second}.")
    intersection = first + second - union
    return first - intersection, second - intersection


def panel_report(r13: dict) -> dict:
    """Build the interval and test report of the 72 window control panel.

    Args:
        r13: The R13 receipt content.

    Returns:
        The Wilson intervals of the three arms and the paired tests of the two single arms.
    """
    truth = panel_truth(r13)
    union = r13["scopes"]["control_panel"]["union_of_both_call_sets"]
    union_counts = {"all": union["tp_base"], "DEL": union["DEL"]["tp_base"],
                    "INS": union["INS"]["tp_base"]}
    arms = {"pipeline": panel_counts(r13, "pipeline_tiddit_manta_cnvnator"),
            "reassembly": panel_counts(r13, "r12_local_reassembly"),
            "merged": panel_counts(r13, "r13_merged")}
    intervals = {name: {group: wilson_interval(counts[group], truth[group]) for group in truth}
                 for name, counts in arms.items()}
    tests = {}
    for group in truth:
        first_only, second_only = discordant_from_union(
            arms["pipeline"][group], arms["reassembly"][group], union_counts[group])
        tests[group] = mcnemar_exact(first_only, second_only)
    return {"truth_records": truth, "recovered": arms, "union_recovered": union_counts,
            "wilson_95": intervals, "mcnemar_pipeline_against_reassembly": tests,
            "no_paired_test_for_the_merged_arm": (
                "The committed record holds the recovered count of the merged arm and the union "
                "count of the two single arms. It does not hold the match status of each truth "
                "record under the merged arm, so the discordant counts of a merged against single "
                "arm comparison cannot be derived and no paired test is reported for that pair.")}


def genome_report(t1: dict) -> dict:
    """Build the interval report of the genome-wide structural variant benchmark.

    The report of this work quotes the arm that Truvari ran with `--dup-to-ins`, which is the field
    `results.v06_dup_as_insertion` of the receipt. The arm without that option is reported beside it,
    because the two arms differ in the insertion counts.

    Args:
        t1: The R2 benchmark receipt content.

    Returns:
        The Wilson intervals of the deletion, insertion and overall recall of both arms, and the
        interval of the CMRG control.
    """
    report = {}
    for arm in ("v06_dup_as_insertion", "v06_as_called"):
        by_class = t1["results"][arm]["by_class"]
        summary = t1["results"][arm]["truvari_summary"]
        report[arm] = {
            "deletion": wilson_interval(by_class["DEL"]["tp_base"],
                                        by_class["DEL"]["truth_records"]),
            "insertion": wilson_interval(by_class["INS"]["tp_base"],
                                         by_class["INS"]["truth_records"]),
            "all": wilson_interval(summary["TP-base"], summary["base cnt"])}
    cmrg = t1["results"]["cmrg_v1_00"]
    report["cmrg_v1_00"] = wilson_interval(cmrg["TP-base"], cmrg["base cnt"])
    report["size_bands_dup_as_insertion"] = {
        f"{svtype}_{band}": wilson_interval(counts["tp_base"], counts["truth_records"])
        for svtype, bands in t1["results"]["v06_dup_as_insertion"]["by_class_and_band"].items()
        for band, counts in bands.items() if counts["truth_records"]}
    return report


def ranker_report(s4: dict) -> dict:
    """Build the paired test of the literature allele ablation on the spike-in suite.

    The two arms score the same 50 scenarios, so each scenario is a paired observation. The receipt
    records that the literature table moves the two literature scenarios and changes nothing in the
    two ClinVar groups, so the discordant scenarios are 2 in favour of the literature arm and 0
    against it.

    Args:
        s4: The S4 suite receipt content.

    Returns:
        The scenario counts, the paired test and the Wilson interval of the top-one pair rate.
    """
    scenarios = s4["suite"]["scenarios"]
    literature_scenarios = s4["suite"]["groups"]["mva_literature"]
    with_table = s4["arms"]["b0"]
    without_table = s4["arms"]["b0_no_literature"]
    test = mcnemar_exact(0, literature_scenarios)
    return {"scenarios": scenarios,
            "mean_rank_points_with_the_table": with_table["mean_rank_points"],
            "mean_rank_points_without_the_table": without_table["mean_rank_points"],
            "scenarios_that_improve": literature_scenarios,
            "scenarios_that_worsen": 0,
            "mcnemar_literature_ablation": test,
            "pair_top1_with_the_table": wilson_interval(
                round(with_table["pair_top1"] * scenarios), scenarios),
            "pair_top1_without_the_table": wilson_interval(
                round(without_table["pair_top1"] * scenarios), scenarios)}


def main() -> int:
    """Compute every interval and every test and write the receipt.

    Returns:
        Zero when the report was written.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    r13 = read_receipt("r13_merged_call_set_v1")
    t1 = read_receipt("t1_sv_mei_benchmark_v1")
    s4 = read_receipt("s4_b0_suite_v1")

    report = {
        "receipt_id": "t1_benchmark_intervals_v1",
        "created_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "track": "t1_call",
        "data_class_of_run": "A",
        "data_class_of_receipt": "A",
        "receipt_scope": ("Interval estimates and paired tests computed from counts that are "
                          "already in the receipts named below. No run happened, no patient file "
                          "was opened and no AWS job started."),
        "code": "scripts/track1_benchmark_intervals.py",
        "tools": {"python": "uv run python, standard library only",
                  "runtime": "local, no AWS job, spend 0.00 USD"},
        "inputs": {"control_panel": "tracking/receipts/r13_merged_call_set_v1.json",
                   "genome_wide": "tracking/receipts/t1_sv_mei_benchmark_v1.json",
                   "spike_in_suite": "tracking/receipts/s4_b0_suite_v1.json"},
        "estimators": {
            "wilson_95": ("Wilson score interval at 95 percent. The unit is one truth record. The "
                          "interval assumes independent Bernoulli truth records with one shared "
                          "probability, which two truth records inside one repeat array do not "
                          "satisfy, so the width is a lower bound."),
            "mcnemar_exact": ("Two-sided exact binomial test on the discordant truth records of two "
                              "arms that ran on the same truth set. The null is a recovery "
                              "probability of 0.5 on a discordant record. No normal approximation "
                              "and no continuity correction.")},
        "control_panel_72_windows": panel_report(r13),
        "genome_wide_v06": genome_report(t1),
        "spike_in_suite": ranker_report(s4),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"control_panel": report["control_panel_72_windows"]["wilson_95"],
                      "mcnemar": report["control_panel_72_windows"]["mcnemar_pipeline_against_reassembly"],
                      "genome_wide": report["genome_wide_v06"],
                      "ranker": report["spike_in_suite"]["mcnemar_literature_ablation"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
