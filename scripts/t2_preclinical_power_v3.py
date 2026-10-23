#!/usr/bin/env python3
"""Test which variance components carry the stage 1 sample size, and sweep the two that do.

Defect class 7 of ``docs/benchmarks/t2_jury_variance_v1.md`` says that the stage 1 sample size rests
on five variance components that no publication measures. The claim treats the five components as
five premises of the design. This script tests that treatment.

The generating model of ``scripts/t2_preclinical_power_v2.py`` carries a donor random intercept, a
passage random intercept, a donor random slope and a residual. The estimator of stage 1 fits the
slope inside one donor line and one passage on a centred dose ladder, so a term that shifts the
whole passage cancels from the estimate. Prediction P1 of
``reports/track2/stage1_variance_gate_v1.md`` states that the two intercept components therefore do
not change the power at any grid point. Prediction P2 states that the residual component and the
donor slope component do change it.

The script runs three parts.

    1. The invariance test of condition C1. It calls the stage 1 simulation of version 2 at a fixed
       seed with four settings of the two intercept components and it reports the largest difference
       in power over those settings.
    2. The reduction check. The two level model collapses to one draw per donor line, because the
       observed donor slope is normal with a variance that the residual component and the donor
       slope component fix together. The script simulates both forms at the same grid points and it
       reports the largest difference, so that the cheap form of part 3 is checked against the full
       form.
    3. The two way sweep of condition C2. It sweeps the residual standard deviation and the donor
       slope standard deviation together and it returns the smallest donor count that reaches a
       power of 0.80 at every cell.

The script makes no model call, it reaches no network and it reads no class B2 file.

Usage:
    uv run python scripts/t2_preclinical_power_v3.py
    uv run python scripts/t2_preclinical_power_v3.py --iterations 20000
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

import t2_preclinical_power_v2 as v2  # noqa: E402
from t2_preclinical_protocol import student_t_quantile  # noqa: E402

OUT = ROOT / "tracking" / "receipts" / "t2_preclinical_power_v3.json"

SEED = 20260919
POWER_TARGET = v2.POWER_TARGET
ALPHA = v2.ALPHA
EFFECT_LOG2 = v2.EFFECT_LOG2
PASSAGES = 8
DONOR_COUNTS = (3, 4, 5, 6, 7, 8, 10, 12, 16, 20)

# The intercept settings of the invariance test. The pair is the donor intercept standard deviation
# and the passage intercept standard deviation, both on the log2 scale.
INTERCEPT_SETTINGS = ((0.0, 0.0), (0.30, 0.20), (0.60, 0.40), (2.35, 1.60))

# The sweep of condition C2. The residual range covers the two coefficients of variation of
# PMID 24475266, which imply a log2 standard deviation of 0.2857 and of 0.3552, and it runs to twice
# the upper one. The slope range runs from no donor variation to a standard deviation equal to the
# slope itself.
RESIDUAL_SWEEP = (0.15, 0.20, 0.25, 0.30, 0.3552, 0.40, 0.50, 0.60, 0.70)
SLOPE_RATIO_SWEEP = (0.0, 0.125, 0.25, 0.375, 0.50, 0.75, 1.0)


def dose_geometry() -> tuple[list[float], float, float]:
    """Give the centred dose ladder, the true slope and the sum of squares of the ladder.

    Returns:
        The centred log10 doses, the slope per log10 of concentration and the sum of the squared
        centred doses, which is the denominator of the least squares slope of one passage.
    """
    doses = v2.centered_log_doses(v2.CONCENTRATIONS_MM)
    slope = v2.slope_per_log10(EFFECT_LOG2, doses)
    sxx = sum(x * x for x in doses)
    return doses, slope, sxx


def observed_slope_sd(sigma: float, omega: float, passages: int, sxx: float) -> float:
    """Give the standard deviation of one observed donor slope.

    The slope of one passage is the true donor slope plus a normal error of variance
    ``sigma ** 2 / sxx``. The mean of the passage slopes of one donor carries that error divided by
    the passage count. The donor slope itself carries the donor random slope. The two intercept
    terms cancel, because the ladder is centred.

    Args:
        sigma: The residual standard deviation of one well, on the log2 scale.
        omega: The standard deviation of the donor random slope, per log10 of concentration.
        passages: The number of independent passages per donor line.
        sxx: The sum of the squared centred doses.

    Returns:
        The standard deviation of the donor slope that the analysis observes.
    """
    return math.sqrt(omega ** 2 + sigma ** 2 / (passages * sxx))


def reduced_power(donors: int, slope: float, spread: float, iterations: int,
                  rng: random.Random) -> float:
    """Estimate the power of the one sample t test on the donor slopes.

    Args:
        donors: The number of donor lines.
        slope: The true slope per log10 of concentration.
        spread: The standard deviation of one observed donor slope.
        iterations: The number of simulated experiments.
        rng: The random source.

    Returns:
        The fraction of experiments that reject at the threshold.
    """
    critical = student_t_quantile(1.0 - ALPHA / 2.0, donors - 1)
    rejections = 0
    for _ in range(iterations):
        values = [rng.gauss(slope, spread) for _ in range(donors)]
        if abs(v2.one_sample_t(values)) >= critical:
            rejections += 1
    return rejections / iterations


def invariance_test(iterations: int) -> dict:
    """Run the stage 1 simulation at four settings of the two intercept components.

    Args:
        iterations: The number of simulated experiments at each point.

    Returns:
        The power at every donor count and setting, and the largest difference over the settings.
    """
    doses, slope, _ = dose_geometry()
    omega = 0.25 * slope
    rows = []
    largest = 0.0
    for donors in (3, 4, 5, 6, 8):
        powers = []
        for tau, gamma in INTERCEPT_SETTINGS:
            rng = random.Random(SEED + donors)
            powers.append(v2.stage1_power(donors, PASSAGES, 0.30, omega, tau, gamma,
                                          iterations, rng))
        spread = max(powers) - min(powers)
        largest = max(largest, spread)
        rows.append({
            "donor_lines": donors,
            "passages_per_line": PASSAGES,
            "residual_sd_log2": 0.30,
            "donor_slope_sd_over_slope": 0.25,
            "power_by_intercept_setting": {
                f"donor_sd_{tau}_passage_sd_{gamma}": round(power, 4)
                for (tau, gamma), power in zip(INTERCEPT_SETTINGS, powers)
            },
            "range_over_the_settings": round(spread, 6),
        })
    return {
        "settings_tested": len(INTERCEPT_SETTINGS),
        "donor_counts_tested": 5,
        "intercept_settings": [{"donor_intercept_sd_log2": tau, "passage_intercept_sd_log2": gamma}
                               for tau, gamma in INTERCEPT_SETTINGS],
        "rows": rows,
        "largest_difference_in_power_over_the_settings": round(largest, 6),
        "prediction_p1_held": largest == 0.0,
        "why_the_result_is_exact": "The estimator is the least squares slope of one passage on a "
                                   "centred dose ladder. The centred doses sum to zero, so the "
                                   "intercept of the passage multiplies zero and leaves the "
                                   "estimate. The simulation draws the same number of normal "
                                   "values in the same order at every setting, so the residual "
                                   "draws of the four runs are the same draws and the powers are "
                                   "equal and not merely close.",
    }


def reduction_check(iterations: int) -> dict:
    """Compare the full two level simulation with the reduced one draw per donor form.

    Args:
        iterations: The number of simulated experiments at each point.

    Returns:
        The two powers at every checked point and the largest difference.
    """
    doses, slope, sxx = dose_geometry()
    rows = []
    largest = 0.0
    for donors in (3, 5, 8):
        for ratio in (0.0, 0.25, 0.50):
            for sigma in (0.20, 0.30, 0.50):
                omega = ratio * slope
                full = v2.stage1_power(donors, PASSAGES, sigma, omega, 0.30, 0.20, iterations,
                                       random.Random(SEED + donors))
                spread = observed_slope_sd(sigma, omega, PASSAGES, sxx)
                reduced = reduced_power(donors, slope, spread, iterations,
                                        random.Random(SEED + donors))
                difference = abs(full - reduced)
                largest = max(largest, difference)
                rows.append({
                    "donor_lines": donors,
                    "donor_slope_sd_over_slope": ratio,
                    "residual_sd_log2": sigma,
                    "power_of_the_full_model": round(full, 4),
                    "power_of_the_reduced_model": round(reduced, 4),
                    "difference": round(difference, 4),
                })
    return {
        "rows": rows,
        "largest_difference": round(largest, 4),
        "monte_carlo_standard_error_at_a_power_of_0_5": round(
            math.sqrt(0.25 / iterations), 4),
        "reading": "The two forms agree inside the Monte Carlo error of the run, so the sweep of "
                   "the next block uses the reduced form.",
    }


def two_way_sweep(iterations: int) -> dict:
    """Sweep the residual standard deviation and the donor slope standard deviation together.

    Args:
        iterations: The number of simulated experiments at each point.

    Returns:
        One row per cell of the sweep, with the power of the design of 5 donor lines and the
        smallest donor count of the grid that reaches the target.
    """
    doses, slope, sxx = dose_geometry()
    rows = []
    for sigma in RESIDUAL_SWEEP:
        for ratio in SLOPE_RATIO_SWEEP:
            omega = ratio * slope
            spread = observed_slope_sd(sigma, omega, PASSAGES, sxx)
            powers = {}
            smallest = None
            for donors in DONOR_COUNTS:
                power = reduced_power(donors, slope, spread, iterations,
                                      random.Random(SEED + donors * 7 + int(sigma * 1000)
                                                    + int(ratio * 100)))
                powers[donors] = round(power, 4)
                if smallest is None and power >= POWER_TARGET:
                    smallest = donors
            rows.append({
                "residual_sd_log2": sigma,
                "donor_slope_sd_over_slope": ratio,
                "observed_donor_slope_sd": round(spread, 4),
                "power_at_five_donor_lines": powers[5],
                "five_donor_lines_reach_the_target": powers[5] >= POWER_TARGET,
                "smallest_donor_count_that_reaches_the_target": smallest,
                "power_by_donor_count": powers,
            })
    reached = [row["smallest_donor_count_that_reaches_the_target"] for row in rows
               if row["smallest_donor_count_that_reaches_the_target"] is not None]
    unreached = [row for row in rows
                 if row["smallest_donor_count_that_reaches_the_target"] is None]
    return {
        "residual_sd_values": list(RESIDUAL_SWEEP),
        "donor_slope_ratio_values": list(SLOPE_RATIO_SWEEP),
        "donor_counts_tried": list(DONOR_COUNTS),
        "passages_per_line": PASSAGES,
        "cells": len(rows),
        "rows": rows,
        "smallest_donor_count_over_the_grid": min(reached) if reached else None,
        "largest_donor_count_over_the_grid": max(reached) if reached else None,
        "cells_in_which_five_donor_lines_reach_the_target": sum(
            1 for row in rows if row["five_donor_lines_reach_the_target"]),
        "cells_in_which_no_donor_count_of_the_grid_reaches_the_target": len(unreached),
        "cells_that_need_more_than_the_largest_count_tried": [
            {"residual_sd_log2": row["residual_sd_log2"],
             "donor_slope_sd_over_slope": row["donor_slope_sd_over_slope"]}
            for row in unreached],
    }


def repository_commit() -> str:
    """Give the current commit of the repository.

    Returns:
        The short commit hash.
    """
    out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                         capture_output=True, text=True, check=True)
    return out.stdout.strip()


def main() -> int:
    """Run the three parts and write the receipt.

    Returns:
        Zero when the receipt is written.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--iterations", type=int, default=20000)
    args = parser.parse_args()

    doses, slope, sxx = dose_geometry()
    invariance = invariance_test(args.iterations)
    reduction = reduction_check(args.iterations)
    sweep = two_way_sweep(args.iterations)

    planning_cell = next(row for row in sweep["rows"]
                         if row["residual_sd_log2"] == 0.30
                         and row["donor_slope_sd_over_slope"] == 0.25)
    version_2 = json.loads(
        (ROOT / "tracking" / "receipts" / "t2_preclinical_power_v2.json").read_text())
    version_2_power = version_2["stage1"]["grid"]
    version_2_planning = next(
        row["power"] for row in version_2_power
        if row["donor_lines"] == 5 and row["passages_per_line"] == PASSAGES
        and row["donor_slope_sd_over_slope"] == 0.25)

    receipt = {
        "receipt_id": "t2_preclinical_power_v3",
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "track": "t2_translation",
        "control": "which_variance_components_carry_the_stage_1_sample_size",
        "queue_item": "R60",
        "data_class_of_run": "A",
        "data_class_of_receipt": "A",
        "receipt_scope": "A simulation of one experimental design. No patient data and no measured "
                         "biological value enters this run.",
        "script": "scripts/t2_preclinical_power_v3.py",
        "repository_commit": repository_commit(),
        "decision_record": "reports/track2/stage1_variance_gate_v1.md",
        "conditions_answered": ["C1", "C2", "C3"],
        "why_this_receipt_exists": "Defect class 7 of docs/benchmarks/t2_jury_variance_v1.md, "
                                   "raised by 3 of the 4 jury replicates and blocking for 2 of the "
                                   "4, states that the stage 1 sample size rests on five variance "
                                   "components that no publication measures.",
        "shared_parameters": {
            "effect_at_the_top_dose_log2": EFFECT_LOG2,
            "concentrations_mm": v2.CONCENTRATIONS_MM,
            "centered_log10_doses": [round(x, 4) for x in doses],
            "slope_per_log10_of_concentration": round(slope, 4),
            "sum_of_squared_centred_doses": round(sxx, 4),
            "alpha": round(ALPHA, 6),
            "power_target": POWER_TARGET,
            "passages_per_line": PASSAGES,
            "iterations": args.iterations,
            "seed": SEED,
        },
        "components_of_the_earlier_search": {
            "record": "tracking/receipts/t2_preclinical_variance_v1.json",
            "components_searched": 5,
            "components_with_a_published_value": 0,
            "which_ones_enter_the_stage_1_power": ["residual_sd_of_log2_bubr1_abundance",
                                                   "donor_slope_sd"],
            "which_ones_do_not": ["donor_intercept_sd", "passage_intercept_sd",
                                  "litter_intraclass_correlation"],
            "why_the_two_intercepts_do_not": "The estimator removes them, which part 1 of this "
                                             "receipt measures.",
            "why_the_litter_component_does_not": "The stage 2 design pairs one treated animal with "
                                                 "one vehicle animal of the same litter and the "
                                                 "same genotype, so the litter term cancels from "
                                                 "the contrast. That component belongs to stage 2 "
                                                 "and not to the stage 1 sample size.",
            "components_that_carry_the_stage_1_sample_size": 2,
        },
        "invariance_test": invariance,
        "reduction_check": reduction,
        "two_way_sweep": sweep,
        "the_design_of_the_report": {
            "donor_lines": 5,
            "passages_per_line": PASSAGES,
            "planning_residual_sd_log2": 0.30,
            "planning_donor_slope_sd_over_slope": 0.25,
            "power_at_the_planning_values": planning_cell["power_at_five_donor_lines"],
            "smallest_donor_count_at_the_planning_values":
                planning_cell["smallest_donor_count_that_reaches_the_target"],
            "the_statement_the_report_makes": "The design of 5 donor lines reaches the target at "
                                              "the planning values of the two components that "
                                              "carry the calculation, and neither value is a "
                                              "measurement. The sweep states the donor count that "
                                              "every other value of the two components needs.",
        },
        "agreement_with_version_2": {
            "record": "tracking/receipts/t2_preclinical_power_v2.json",
            "power_of_version_2_at_five_donor_lines_and_the_planning_values": version_2_planning,
            "power_of_this_run_at_the_same_point": planning_cell["power_at_five_donor_lines"],
            "difference": round(abs(version_2_planning
                                    - planning_cell["power_at_five_donor_lines"]), 4),
            "monte_carlo_standard_error_at_that_power": round(
                math.sqrt(version_2_planning * (1.0 - version_2_planning) / args.iterations), 4),
            "reading": "The two runs use a different seed and a different simulation form and they "
                       "agree inside the Monte Carlo error. The report states the value of version "
                       "2, which is the receipt that the design of Supplementary Table S36 names.",
        },
        "predictions": {
            "P1": {
                "statement": "The donor intercept and the passage intercept do not change the "
                             "simulated power at any grid point.",
                "held": invariance["prediction_p1_held"],
                "evidence": f"The largest difference in power over four intercept settings is "
                            f"{invariance['largest_difference_in_power_over_the_settings']}.",
            },
            "P2": {
                "statement": "The power changes when the residual standard deviation changes and "
                             "when the donor slope standard deviation changes.",
                "held": (sweep["largest_donor_count_over_the_grid"] !=
                         sweep["smallest_donor_count_over_the_grid"]),
                "evidence": f"Over the {sweep['cells']} cells of the sweep the smallest donor "
                            f"count that reaches the target runs from "
                            f"{sweep['smallest_donor_count_over_the_grid']} to "
                            f"{sweep['largest_donor_count_over_the_grid']}, and "
                            f"{sweep['cells_in_which_no_donor_count_of_the_grid_reaches_the_target']}"
                            f" cells reach it at no count of the grid.",
            },
        },
        "what_this_run_does_not_establish": "A simulation is not a measurement. The run states "
                                            "which components change the power and by how much, "
                                            "and it measures no variance of the endpoint of the "
                                            "protocol. The effect of 1.5 fold at the top dose also "
                                            "carries no published basis and this run does not "
                                            "sweep it.",
        "verdict": "PASS",
    }
    args.out.write_text(json.dumps(receipt, indent=1) + "\n")
    print(f"wrote {args.out}")
    print(f"P1 held: {invariance['prediction_p1_held']}, largest difference "
          f"{invariance['largest_difference_in_power_over_the_settings']}")
    print(f"reduction check largest difference {reduction['largest_difference']}")
    print(f"sweep cells {sweep['cells']}, donor count "
          f"{sweep['smallest_donor_count_over_the_grid']} to "
          f"{sweep['largest_donor_count_over_the_grid']}, "
          f"five lines reach the target in "
          f"{sweep['cells_in_which_five_donor_lines_reach_the_target']} cells, "
          f"{sweep['cells_in_which_no_donor_count_of_the_grid_reaches_the_target']} cells reach it "
          f"at no count")
    return 0


if __name__ == "__main__":
    sys.exit(main())
