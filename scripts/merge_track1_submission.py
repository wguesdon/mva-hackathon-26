#!/usr/bin/env python3
"""Merge the two Track 1 writer outputs into one organizer file. Queue row 13.

The B0 small variant ranker writes the primary rows and the structural and mobile element writer
writes the secondary rows. The merged file keeps every small variant row above every structural row,
in the order that each writer gives, and it holds at most ten rows.

The structural writer maps its integer call score onto an EPCR, so two calls with one score share one
EPCR. The organizer scorer in `evaluation.py` sorts rows by EPCR in descending order and keeps ties in
file order, so a tie does not change a rank. It does merge the tied rows into one F-max threshold,
and `scripts/check_merged_submission.py` refuses a shared value for that reason. This script lowers
each tied value by a fixed step in file order. The rank of every row under the organizer sort stays
the same, and each row receives its own threshold. The output reports every changed value, so the
receipt of the run can record the deviation from the writer output.

Usage:
    uv run python scripts/merge_track1_submission.py --small <b0 csv> --structural <sv csv> \\
        --out <merged csv> [--tie-step 0.001] [--row-budget 10]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from decimal import Decimal
from pathlib import Path

HEADER = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2",
          "alt_2", "epcr", "finding_type", "notes"]
ROW_BUDGET = 10
TIE_STEP = Decimal("0.001")


def read_writer_rows(path: Path) -> list[dict[str, str]]:
    """Read the rows of one writer output and check its header.

    Args:
        path: The CSV file that one writer produced.

    Returns:
        The rows as dictionaries, in file order.

    Raises:
        ValueError: If the header differs from the organizer template.
    """
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != HEADER:
            raise ValueError(f"{path.name}: the header differs from the organizer template")
        return list(reader)


def check_order(rows: list[dict[str, str]]) -> None:
    """Check that the EPCR values never rise from one row to the next.

    A rise would make the organizer sort move a row, and then a tie step in file order would
    change a rank.

    Args:
        rows: The merged rows, small variant rows first.

    Raises:
        ValueError: If a row carries a larger EPCR than the row above it.
    """
    for index in range(1, len(rows)):
        if Decimal(rows[index]["epcr"]) > Decimal(rows[index - 1]["epcr"]):
            raise ValueError(f"Row {index + 1} has a larger EPCR than the row above it")


def separate_ties(rows: list[dict[str, str]], step: Decimal) -> list[dict[str, str]]:
    """Lower each tied EPCR by a fixed step in file order.

    Args:
        rows: The merged rows in non-increasing EPCR order.
        step: The amount by which a tied value falls below the value of the row above it.

    Returns:
        One record per changed row with the row number, the original value and the new value.

    Raises:
        ValueError: If a new value falls to zero or below.
    """
    changes = []
    for index in range(1, len(rows)):
        above = Decimal(rows[index - 1]["epcr"])
        current = Decimal(rows[index]["epcr"])
        if current < above:
            continue
        new = above - step
        if new <= 0:
            raise ValueError(f"Row {index + 1} falls to an EPCR of zero or below")
        changes.append({"row": index + 1, "from": rows[index]["epcr"], "to": str(new)})
        rows[index]["epcr"] = str(new)
    return changes


def merge(small: list[dict[str, str]], structural: list[dict[str, str]], row_budget: int,
          step: Decimal) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Build the merged rows and separate their ties.

    Args:
        small: The rows of the B0 small variant ranker.
        structural: The rows of the structural and mobile element writer.
        row_budget: The largest row count that the organizer file accepts.
        step: The tie step.

    Returns:
        The merged rows and the list of changed EPCR values.

    Raises:
        ValueError: If the merged file exceeds the budget or a row order check fails.
    """
    rows = [dict(row) for row in small] + [dict(row) for row in structural]
    if len(rows) > row_budget:
        raise ValueError(f"The merged file holds {len(rows)} rows and the budget is {row_budget}")
    check_order(rows)
    changes = separate_ties(rows, step)
    return rows, changes


def write_rows(path: Path, rows: list[dict[str, str]]) -> str:
    """Write the merged rows and return the checksum of the file.

    Args:
        path: The output CSV file.
        rows: The merged rows.

    Returns:
        The SHA-256 checksum of the written file.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    """Merge the two writer outputs and print the checksums and the changed values.

    Returns:
        Zero on success and one on a refused input.
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--small", required=True, type=Path)
    parser.add_argument("--structural", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--tie-step", type=Decimal, default=TIE_STEP)
    parser.add_argument("--row-budget", type=int, default=ROW_BUDGET)
    args = parser.parse_args()
    try:
        rows, changes = merge(read_writer_rows(args.small), read_writer_rows(args.structural),
                              args.row_budget, args.tie_step)
    except (OSError, ValueError) as error:
        print(f"REFUSED: {error}", file=sys.stderr)
        return 1
    digest = write_rows(args.out, rows)
    print(json.dumps({
        "out": str(args.out),
        "sha256": digest,
        "inputs": {
            "small": hashlib.sha256(args.small.read_bytes()).hexdigest(),
            "structural": hashlib.sha256(args.structural.read_bytes()).hexdigest(),
        },
        "tie_step": str(args.tie_step),
        "changed_epcr": changes,
    }, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
