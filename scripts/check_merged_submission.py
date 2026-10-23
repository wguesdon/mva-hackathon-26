#!/usr/bin/env python3
"""Check the shape of the merged Track 1 submission file and write a class A gate file.

The organizer file holds at most ten rows and two components propose them. The B0 small variant
ranker writes a decreasing EPCR ladder and the structural and mobile element writer maps its call
score onto an EPCR. The organizer scorer builds one threshold per distinct EPCR, so two rows that
share one value share one threshold. `tracking/receipts/s4_row_count_calibration_v1.json`, field
`tie_counterfactual`, measures that cost: one shared value across ten rows drops mean F-max from
0.960 to 0.238 with no change in rank points.

The format check inside `pipelines/track1_b0/main.nf` reads the B0 file alone, so it cannot see a
collision between the two writers or a repeat inside the structural rows. This script reads the
merged file that the owner assembles and it tests four conditions over that file.

    rows_present                        The file holds at least one row.
    rows_within_budget                  The row count is at most the organizer budget.
    epcr_values_are_distinct            No two rows of the file share one EPCR.
    primary_above_every_secondary       Every small variant row sits above every structural row.

The merged file is class B2, because its row count and its EPCR values state how many candidates the
file of the child holds. This script therefore prints no row count, no EPCR value and no coordinate.
The gate file carries the verdicts and the checksum of the input. Run
`make privacy-scan-publish FILES=<gate file>` after the run, which is control 7 of
`docs/llm_data_privacy.md`.

Usage:
    uv run python scripts/check_merged_submission.py --submission <merged csv> \\
        --out <gate json> [--row-budget 10]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROW_BUDGET = 10
PRIMARY = "primary"
SECONDARY = "secondary"


def sha256_of(path: Path) -> str:
    """Compute the SHA-256 checksum of a file.

    Args:
        path: The file to read.

    Returns:
        The checksum in lower case hexadecimal.
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    """Read the submission rows of one file.

    Args:
        path: The merged submission file.

    Returns:
        The rows as dictionaries.
    """
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def verdicts(rows: list[dict[str, str]], row_budget: int) -> dict[str, object]:
    """Judge the four conditions over the merged rows.

    Args:
        rows: The submission rows.
        row_budget: The largest row count that the organizer file accepts.

    Returns:
        The verdict fields of the gate file. No field carries a row count or an EPCR value.
    """
    values = [float(row["epcr"]) for row in rows]
    primary = [float(row["epcr"]) for row in rows if row.get("finding_type") == PRIMARY]
    secondary = [float(row["epcr"]) for row in rows if row.get("finding_type") == SECONDARY]
    separated = None
    if primary and secondary:
        separated = min(primary) > max(secondary)
    return {
        "rows_present": "PASS" if rows else "FAIL",
        "rows_within_budget": "PASS" if len(rows) <= row_budget else "FAIL",
        "epcr_values_are_distinct": "PASS" if len(set(values)) == len(values) else "FAIL",
        "epcr_inside_the_unit_interval": "PASS" if all(0 < value <= 1 for value in values)
        else "FAIL",
        "both_writers_present": "PASS" if primary and secondary else "SINGLE WRITER",
        "primary_above_every_secondary": {True: "PASS", False: "FAIL", None: "NOT APPLICABLE"}[
            separated],
        "proband_ids": sorted({row["proband_id"] for row in rows}),
    }


def main() -> int:
    """Write the class A gate file for one merged submission file.

    Returns:
        Zero when every applicable condition reads PASS, and one otherwise.
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--submission", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--row-budget", type=int, default=ROW_BUDGET)
    args = parser.parse_args()
    gate: dict[str, object] = {
        "checker": "scripts/check_merged_submission.py",
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "file": args.submission.name,
        "sha256": sha256_of(args.submission),
        "data_class": "A",
        "scope": "The shape of the merged file, not its content.",
        "row_budget": args.row_budget,
    }
    try:
        rows = read_rows(args.submission)
        gate.update(verdicts(rows, args.row_budget))
        gate["format_check"] = "PASS"
    except (OSError, KeyError, ValueError) as error:
        gate["format_check"] = "FAIL"
        gate["error"] = f"{type(error).__name__}: {str(error)[:200]}"
    checks = [value for key, value in gate.items()
              if key.startswith(("rows_", "epcr_", "primary_", "format_"))]
    gate["verdict"] = "PASS" if all(value in ("PASS", "NOT APPLICABLE") for value in checks) \
        else "FAIL"
    args.out.write_text(json.dumps(gate, indent=2) + "\n")
    print(json.dumps({"out": str(args.out), "verdict": gate["verdict"]}, indent=2))
    return 0 if gate["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
