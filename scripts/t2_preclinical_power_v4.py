#!/usr/bin/env python3
"""Compute the stage 1 power of the primary nicotinic acid test at alpha 0.05. Data class A.

Amendment 5 of ``reports/track2/preclinical_protocol_gate_v1.md`` makes nicotinic acid the primary
compound of stage 1. It is tested first at alpha 0.05, and the three earlier compounds are tested
after it at the 0.05 shared over three that versions 2 and 3 of the power calculation already used.
Condition C17 of that amendment asks for the power of the primary test at its own threshold.

This script changes one quantity and nothing else. It sets the alpha of the two earlier power
modules to 0.05, then repeats the stage 1 grid of version 2 and the two-way sweep of version 3 under
the same seeds, the same dose geometry, the same variance components and the same iteration counts.
The secondary compounds keep the results of versions 2 and 3, because their threshold does not
change.

Usage:
    uv run python scripts/t2_preclinical_power_v4.py
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import t2_preclinical_power_v2 as v2  # noqa: E402
import t2_preclinical_power_v3 as v3  # noqa: E402

GATE = ROOT / "reports" / "track2" / "preclinical_protocol_gate_v1.md"
OUT = ROOT / "tracking" / "receipts" / "t2_preclinical_power_v4.json"
PRIMARY_ALPHA = 0.05
# The iteration counts of the two earlier runs, read from their receipts.
GRID_ITERATIONS = 50000
SWEEP_ITERATIONS = 20000
DESIGN_DONORS = 5
DESIGN_PASSAGES = 8
PLANNING_SLOPE_RATIO = 0.25


def set_alpha(alpha: float) -> None:
    """Set the alpha that the two earlier power modules read at call time.

    Args:
        alpha: The two sided significance threshold of the test.
    """
    v2.ALPHA = alpha
    v3.ALPHA = alpha


def design_row(grid: list[dict]) -> dict:
    """Pick the row of the stage 1 grid that holds the design of the report.

    Args:
        grid: The rows that the version 2 stage 1 grid returns.

    Returns:
        The row at 5 donor lines, 8 passages and the planning donor slope ratio.
    """
    return next(row for row in grid if row["donor_lines"] == DESIGN_DONORS
                and row["passages_per_line"] == DESIGN_PASSAGES
                and row["donor_slope_sd_over_slope"] == PLANNING_SLOPE_RATIO)


def repository_commit() -> str:
    """Give the current commit of the repository.

    Returns:
        The short commit hash.
    """
    out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                         capture_output=True, text=True, check=True)
    return out.stdout.strip()


def main() -> int:
    """Run the grid and the sweep at the primary threshold and write the receipt.

    Returns:
        Zero when the receipt is written.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()

    set_alpha(PRIMARY_ALPHA)
    grid = v2.stage1_grid(GRID_ITERATIONS)
    sweep = v3.two_way_sweep(SWEEP_ITERATIONS)
    design = design_row(grid)

    receipt = {
        "receipt_id": "t2_preclinical_power_v4",
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "data_class_of_run": "A",
        "data_class_of_receipt": "A",
        "gate": str(GATE.relative_to(ROOT)),
        "gate_conditions": ["C17"],
        "repository_commit": repository_commit(),
        "script": "scripts/t2_preclinical_power_v4.py",
        "why_this_receipt_exists": (
            "Amendment 5 of the gate makes nicotinic acid the primary compound of stage 1 at alpha "
            "0.05. Versions 2 and 3 computed power at 0.05 shared over three compounds."),
        "shared_parameters": {
            "alpha": PRIMARY_ALPHA,
            "alpha_rule": ("0.05 for the primary nicotinic acid test. The three secondary "
                           "compounds keep 0.05 shared over three and the results of versions 2 "
                           "and 3."),
            "concentrations_mM": v2.CONCENTRATIONS_MM,
            "effect_log2": v2.EFFECT_LOG2,
            "power_target": v2.POWER_TARGET,
            "seed": v2.SEED,
            "grid_iterations": GRID_ITERATIONS,
            "sweep_iterations": SWEEP_ITERATIONS,
            "what_changes_against_versions_2_and_3": "the alpha and nothing else",
        },
        "the_design_of_the_report": {
            "donor_lines": DESIGN_DONORS,
            "passages_per_line": DESIGN_PASSAGES,
            "donor_slope_sd_over_slope": PLANNING_SLOPE_RATIO,
            "power": design["power"],
            "reaches_the_target": design["power"] >= v2.POWER_TARGET,
        },
        "stage1_grid": grid,
        "two_way_sweep": sweep,
    }
    args.out.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"primary test at alpha {PRIMARY_ALPHA}: power {design['power']} at "
          f"{DESIGN_DONORS} donor lines and {DESIGN_PASSAGES} passages")
    print(f"sweep: donor counts that reach the target run from "
          f"{sweep['smallest_donor_count_over_the_grid']} to "
          f"{sweep['largest_donor_count_over_the_grid']}, and 5 lines reach it in "
          f"{sweep['cells_in_which_five_donor_lines_reach_the_target']} of {sweep['cells']} cells")
    print(f"receipt: {args.out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
