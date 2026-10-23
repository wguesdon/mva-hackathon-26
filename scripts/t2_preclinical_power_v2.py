#!/usr/bin/env python3
"""Simulate the power of the two stages of the preclinical protocol under their own designs.

Amendment 2 of `reports/track2/preclinical_protocol_gate_v1.md` adds conditions C8, C9 and C10. The
first power calculation of `scripts/t2_preclinical_protocol.py` simulated a two sample t test of the
top dose against the vehicle arm, and the protocol states a different analysis for each stage. This
script simulates the analysis that each stage states.

Stage 1. The primary contrast is the slope of log2 BUBR1 abundance on log10 concentration in the
*BUB1B*-reduced lines. The generating model carries a donor random intercept, a passage random
intercept inside a donor and a donor random slope:

    y[d, p, c] = mu + u[d] + v[d, p] + (beta + s[d]) * x[c] + e[d, p, c]

    u[d]    ~ N(0, tau ** 2)      one value per donor line
    v[d, p] ~ N(0, gamma ** 2)    one value per passage inside a donor
    s[d]    ~ N(0, omega ** 2)    one value per donor line, the random slope
    e       ~ N(0, sigma ** 2)    one value per well

The two intercept terms cancel from a slope that is estimated inside one donor and one passage, so
the estimator is the ordinary least squares slope of that passage. The passage slopes of one donor
average to a donor slope. The donor slopes are independent and they carry the donor random slope, so
the donor line is the unit of replication and the test is a one sample two sided t test on the donor
slopes with `donors - 1` degrees of freedom. Condition C9 asks for that statement.

Stage 2. The design is two genotypes by two arms with one animal per litter per arm. The litter is a
random block that the four animals of a litter share:

    y[k, g, a] = mu + l[k] + genotype[g] + arm[a] + interaction[g, a] + e[k, g, a]

    l[k] ~ N(0, lam ** 2)         one value per litter
    e    ~ N(0, sigma ** 2)       one value per animal

The litter term cancels from the contrast of the treated animal against the vehicle animal of the
same litter and the same genotype. The test is a one sample two sided t test on the litter
differences of the *BubR1*-insufficient genotype with `litters - 1` degrees of freedom.

Both stages test at the same threshold, which is 0.05 shared over the three compounds of stage 1.
The effect and the standard deviation carry no published basis, which condition C5 of the gate
records, so the script sweeps the standard deviation and the variance components instead of fixing
one value.

Usage:
    uv run python scripts/t2_preclinical_power_v2.py --out tracking/receipts/t2_preclinical_power_v2.json
"""

from __future__ import annotations

import argparse
import json
import math
import random
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from t2_preclinical_protocol import student_t_quantile  # noqa: E402  the helper lives there

# The five concentrations of stage 1 in millimolar, from Table 26 of the report.
CONCENTRATIONS_MM = [0.01, 0.0316, 0.1, 0.316, 1.0]
EFFECT_LOG2 = 0.585
ALPHA = 0.05 / 3
POWER_TARGET = 0.80
SEED = 20260918


def centered_log_doses(concentrations: list[float]) -> list[float]:
    """Give the centered log10 concentrations of the dose ladder.

    Args:
        concentrations: The concentrations of the ladder in millimolar.

    Returns:
        The log10 of each concentration, minus their mean.
    """
    logs = [math.log10(value) for value in concentrations]
    mean = sum(logs) / len(logs)
    return [value - mean for value in logs]


def slope_per_log10(effect_log2: float, doses: list[float]) -> float:
    """Convert a fold change at the top dose into a slope per log10 of concentration.

    Args:
        effect_log2: The rise at the top dose on the log2 scale, relative to the lowest dose.
        doses: The centered log10 concentrations.

    Returns:
        The slope of log2 abundance per unit of log10 concentration.
    """
    return effect_log2 / (max(doses) - min(doses))


def passage_slope(rng: random.Random, doses: list[float], intercept: float, slope: float,
                  sigma: float) -> float:
    """Draw one passage and return its ordinary least squares slope.

    Args:
        rng: The random source.
        doses: The centered log10 concentrations.
        intercept: The intercept of this passage, which the estimator removes.
        slope: The true slope of this donor.
        sigma: The residual standard deviation of one well.

    Returns:
        The least squares slope of the passage.
    """
    sxx = sum(x * x for x in doses)
    sxy = 0.0
    for x in doses:
        y = intercept + slope * x + rng.gauss(0.0, sigma)
        sxy += x * y
    return sxy / sxx


def one_sample_t(values: list[float]) -> float:
    """Give the two sided one sample t statistic of a list against zero.

    Args:
        values: The observed values.

    Returns:
        The t statistic, or zero when the sample variance is zero.
    """
    n = len(values)
    mean = sum(values) / n
    variance = sum((v - mean) ** 2 for v in values) / (n - 1)
    if variance <= 0.0:
        return 0.0
    return mean / math.sqrt(variance / n)


def stage1_power(donors: int, passages: int, sigma: float, omega: float, tau: float, gamma: float,
                 iterations: int, rng: random.Random) -> float:
    """Estimate the power of the stage 1 slope test under the two level model.

    Args:
        donors: The number of donor lines.
        passages: The number of independent passages per donor line.
        sigma: The residual standard deviation of one well, on the log2 scale.
        omega: The standard deviation of the donor random slope, per log10 of concentration.
        tau: The standard deviation of the donor random intercept.
        gamma: The standard deviation of the passage random intercept inside a donor.
        iterations: The number of simulated experiments.
        rng: The random source.

    Returns:
        The fraction of experiments that reject at the threshold.
    """
    doses = centered_log_doses(CONCENTRATIONS_MM)
    true_slope = slope_per_log10(EFFECT_LOG2, doses)
    critical = student_t_quantile(1.0 - ALPHA / 2.0, donors - 1)
    rejections = 0
    for _ in range(iterations):
        donor_slopes = []
        for _ in range(donors):
            donor_intercept = rng.gauss(0.0, tau)
            donor_slope = true_slope + rng.gauss(0.0, omega)
            estimates = []
            for _ in range(passages):
                intercept = donor_intercept + rng.gauss(0.0, gamma)
                estimates.append(passage_slope(rng, doses, intercept, donor_slope, sigma))
            donor_slopes.append(sum(estimates) / passages)
        if abs(one_sample_t(donor_slopes)) >= critical:
            rejections += 1
    return rejections / iterations


def stage2_power(litters: int, sigma: float, iterations: int, rng: random.Random) -> float:
    """Estimate the power of the stage 2 within-litter contrast.

    Args:
        litters: The number of litters, which is the group size of each arm.
        sigma: The residual standard deviation of one animal, on the log2 scale.
        iterations: The number of simulated experiments.
        rng: The random source.

    Returns:
        The fraction of experiments that reject at the threshold.
    """
    critical = student_t_quantile(1.0 - ALPHA / 2.0, litters - 1)
    rejections = 0
    for _ in range(iterations):
        differences = []
        for _ in range(litters):
            treated = EFFECT_LOG2 + rng.gauss(0.0, sigma)
            vehicle = rng.gauss(0.0, sigma)
            differences.append(treated - vehicle)
        if abs(one_sample_t(differences)) >= critical:
            rejections += 1
    return rejections / iterations


def stage1_grid(iterations: int) -> list[dict]:
    """Sweep the donor count, the passage count and the slope variance of stage 1.

    Args:
        iterations: The number of simulated experiments at each point.

    Returns:
        One row per grid point, with the power and the smallest passage count that reaches the
        target at that donor count.
    """
    rows = []
    for donors in (3, 4, 5, 6, 8):
        for omega_ratio in (0.0, 0.25, 0.50):
            doses = centered_log_doses(CONCENTRATIONS_MM)
            omega = omega_ratio * slope_per_log10(EFFECT_LOG2, doses)
            for passages in (4, 6, 8, 12):
                rng = random.Random(SEED + donors * 1000 + passages * 10
                                    + int(omega_ratio * 100))
                power = stage1_power(donors, passages, 0.30, omega, 0.30, 0.20, iterations, rng)
                rows.append({
                    "donor_lines": donors,
                    "passages_per_line": passages,
                    "residual_sd_log2": 0.30,
                    "donor_slope_sd_over_slope": omega_ratio,
                    "power": round(power, 4),
                })
    return rows


def stage2_grid(iterations: int) -> list[dict]:
    """Sweep the litter count and the residual standard deviation of stage 2.

    Args:
        iterations: The number of simulated experiments at each point.

    Returns:
        One row per grid point with the simulated power.
    """
    rows = []
    for sigma in (0.20, 0.30, 0.40):
        for litters in (4, 6, 8, 10, 12, 16):
            rng = random.Random(SEED + litters * 100 + int(sigma * 100))
            power = stage2_power(litters, sigma, iterations, rng)
            rows.append({
                "litters": litters,
                "residual_sd_log2": sigma,
                "power": round(power, 4),
            })
    return rows


def smallest_reaching(rows: list[dict], key: str, filters: dict) -> int | None:
    """Give the smallest value of one column whose row reaches the target power.

    Args:
        rows: The grid rows.
        key: The column to minimize.
        filters: The column values that select the rows of one design.

    Returns:
        The smallest value that reaches the target, or None when no row of the grid reaches it.
    """
    matches = [row for row in rows
               if all(row[name] == value for name, value in filters.items())
               and row["power"] >= POWER_TARGET]
    if not matches:
        return None
    return min(row[key] for row in matches)


def repository_commit() -> str:
    """Give the current commit of the repository.

    Returns:
        The short commit hash.
    """
    out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                         capture_output=True, text=True, check=True)
    return out.stdout.strip()


def main() -> int:
    """Run both simulations and write the receipt.

    Returns:
        Zero when the receipt is written.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path,
                        default=ROOT / "tracking" / "receipts" / "t2_preclinical_power_v2.json")
    parser.add_argument("--iterations", type=int, default=20000)
    args = parser.parse_args()

    doses = centered_log_doses(CONCENTRATIONS_MM)
    stage1 = stage1_grid(args.iterations)
    stage2 = stage2_grid(args.iterations)

    receipt = {
        "receipt_id": "t2_preclinical_power_v2",
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "track": "t2_translation",
        "control": "power_of_the_stated_analysis_of_each_preclinical_stage",
        "queue_item": "R42",
        "data_class_of_run": "A",
        "data_class_of_receipt": "A",
        "receipt_scope": "A simulation of two experimental designs. No patient data and no measured "
                         "biological value enters this run.",
        "decision_record": "reports/track2/preclinical_protocol_gate_v1.md, amendment 2",
        "conditions_answered": ["C8", "C9", "C10"],
        "script": "scripts/t2_preclinical_power_v2.py",
        "repository_commit": repository_commit(),
        "why_this_receipt_exists": "Defect 7 of reports/jury/jury_review_track2_2026-09-18T134500Z.md "
                                   "states that the power calculation of Table 27 simulates a two "
                                   "sample t test of the top dose against the vehicle arm while the "
                                   "protocol states a mixed model slope for stage 1 and a litter "
                                   "structured two by two design for stage 2, and that stage 2 "
                                   "copies the group size of stage 1.",
        "shared_parameters": {
            "effect_at_the_top_dose_log2": EFFECT_LOG2,
            "concentrations_mm": CONCENTRATIONS_MM,
            "centered_log10_doses": [round(x, 4) for x in doses],
            "slope_per_log10_of_concentration": round(slope_per_log10(EFFECT_LOG2, doses), 4),
            "alpha": round(ALPHA, 6),
            "alpha_rule": "0.05 shared over the three compounds of stage 1 at the one primary "
                          "contrast.",
            "power_target": POWER_TARGET,
            "iterations": args.iterations,
            "seed": SEED,
            "no_published_basis": "The effect of 1.5 fold and every standard deviation of this "
                                  "receipt carry no published basis. No publication reports the "
                                  "between-replicate variance of a BUBR1 immunoblot in a human "
                                  "fibroblast or in mouse skeletal muscle.",
        },
        "stage1": {
            "analysis": "The slope of log2 BUBR1 abundance on log10 concentration, estimated inside "
                        "one donor line and one passage, averaged over the passages of a donor "
                        "line, and tested by a one sample two sided t test on the donor slopes.",
            "unit_of_replication": "the donor line",
            "degrees_of_freedom": "donor lines minus one",
            "random_effects": {
                "donor_intercept_sd_log2": 0.30,
                "passage_intercept_sd_log2": 0.20,
                "donor_slope_sd_over_the_slope": [0.0, 0.25, 0.50],
            },
            "residual_sd_log2": 0.30,
            "grid": stage1,
            "protocol_design_of_table_26": {
                "donor_lines": 3,
                "passages_per_line": 8,
                "power_at_no_donor_slope_variance": next(
                    row["power"] for row in stage1
                    if row["donor_lines"] == 3 and row["passages_per_line"] == 8
                    and row["donor_slope_sd_over_slope"] == 0.0),
                "power_at_a_donor_slope_sd_of_one_quarter_of_the_slope": next(
                    row["power"] for row in stage1
                    if row["donor_lines"] == 3 and row["passages_per_line"] == 8
                    and row["donor_slope_sd_over_slope"] == 0.25),
                "power_at_a_donor_slope_sd_of_one_half_of_the_slope": next(
                    row["power"] for row in stage1
                    if row["donor_lines"] == 3 and row["passages_per_line"] == 8
                    and row["donor_slope_sd_over_slope"] == 0.50),
            },
            "smallest_donor_count_that_reaches_the_target": {
                "at_no_donor_slope_variance": min(
                    (row["donor_lines"] for row in stage1
                     if row["donor_slope_sd_over_slope"] == 0.0 and row["power"] >= POWER_TARGET),
                    default=None),
                "at_a_donor_slope_sd_of_one_quarter": min(
                    (row["donor_lines"] for row in stage1
                     if row["donor_slope_sd_over_slope"] == 0.25 and row["power"] >= POWER_TARGET),
                    default=None),
                "at_a_donor_slope_sd_of_one_half": min(
                    (row["donor_lines"] for row in stage1
                     if row["donor_slope_sd_over_slope"] == 0.50 and row["power"] >= POWER_TARGET),
                    default=None),
            },
        },
        "stage2": {
            "analysis": "The difference of the treated animal and the vehicle animal of the same "
                        "litter, inside the BubR1-insufficient genotype, tested by a one sample two "
                        "sided t test on the litter differences.",
            "unit_of_replication": "the litter",
            "degrees_of_freedom": "litters minus one",
            "random_effects": {"litter_intercept": "cancels from the within-litter difference"},
            "grid": stage2,
            "smallest_litter_count_that_reaches_the_target": {
                str(sigma): smallest_reaching(stage2, "litters", {"residual_sd_log2": sigma})
                for sigma in (0.20, 0.30, 0.40)
            },
        },
    }
    args.out.write_text(json.dumps(receipt, indent=1) + "\n")
    print(f"wrote {args.out}")
    print(json.dumps(receipt["stage1"]["protocol_design_of_table_26"], indent=1))
    print(json.dumps(receipt["stage1"]["smallest_donor_count_that_reaches_the_target"], indent=1))
    print(json.dumps(receipt["stage2"]["smallest_litter_count_that_reaches_the_target"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
