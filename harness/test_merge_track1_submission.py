"""The Track 1 merge keeps the organizer rank of every row and gives each row its own threshold.

The fixtures of this test are synthetic. The chromosome chrT does not exist, and no patient row and
no real submission file enters the test.
"""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MERGER = ROOT / "scripts" / "merge_track1_submission.py"
HEADER = "proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type,notes"
SMALL = [
    "P,chrT,100,A,G,chrT,200,C,T,0.99,primary,pair",
    "P,chrT,300,G,A,,,,,0.98,primary,single",
]
STRUCTURAL = [
    "P,chrT,1000,N,<INV>,,,,,0.675,secondary,sv one",
    "P,chrT,2000,N,<INV>,,,,,0.675,secondary,sv two",
    "P,chrT,3000,N,<DEL>,,,,,0.66,secondary,sv three",
    "P,chrT,4000,N,<DEL>,,,,,0.66,secondary,sv four",
]


def write_csv(directory: Path, name: str, lines: list[str]) -> Path:
    """Write a synthetic writer output.

    Args:
        directory: The temporary directory.
        name: The file name.
        lines: The data lines without the header.

    Returns:
        The path of the written file.
    """
    path = directory / name
    path.write_text("\n".join([HEADER, *lines]) + "\n")
    return path


def run_merge(directory: Path, small: list[str], structural: list[str]) -> subprocess.CompletedProcess:
    """Run the merge script on two synthetic writer outputs.

    Args:
        directory: The temporary directory.
        small: The small variant rows.
        structural: The structural rows.

    Returns:
        The completed process.
    """
    return subprocess.run(
        [sys.executable, str(MERGER), "--small", str(write_csv(directory, "s.csv", small)),
         "--structural", str(write_csv(directory, "v.csv", structural)),
         "--out", str(directory / "merged.csv")],
        capture_output=True, text=True, check=False)


def organizer_key(item: tuple[int, dict[str, str]]) -> tuple[float, int]:
    """Return the sort key of the organizer scorer, EPCR descending and then file order.

    Args:
        item: The file index of a row and the row.

    Returns:
        The sort key.
    """
    index, row = item
    return (-float(row["epcr"]), index)


def organizer_order(path: Path) -> list[str]:
    """Sort the rows as the organizer scorer does and return their notes.

    Args:
        path: The merged file.

    Returns:
        The notes field of each row in organizer rank order.
    """
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    ranked = sorted(enumerate(rows), key=organizer_key)
    return [row["notes"] for _, row in ranked]


class MergeTest(unittest.TestCase):
    """Fixtures for the merge of queue row 13."""

    def test_ties_are_separated_and_ranks_kept(self) -> None:
        """Tied structural rows get distinct values and keep their organizer rank."""
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            result = run_merge(directory, SMALL, STRUCTURAL)
            self.assertEqual(result.returncode, 0, result.stderr)
            with (directory / "merged.csv").open(newline="") as handle:
                values = [Decimal(row["epcr"]) for row in csv.DictReader(handle)]
            self.assertEqual(len(set(values)), len(values))
            self.assertEqual(values, sorted(values, reverse=True))
            self.assertEqual(values[3], Decimal("0.674"))
            self.assertEqual(values[5], Decimal("0.659"))
            expected = ["pair", "single", "sv one", "sv two", "sv three", "sv four"]
            self.assertEqual(organizer_order(directory / "merged.csv"), expected)

    def test_rising_value_is_refused(self) -> None:
        """A structural row above a small variant row is refused, because a step would move it."""
        with tempfile.TemporaryDirectory() as tmp:
            rising = ["P,chrT,5000,N,<INV>,,,,,0.995,secondary,too high"]
            result = run_merge(Path(tmp), SMALL, rising)
            self.assertEqual(result.returncode, 1)
            self.assertIn("larger EPCR", result.stderr)

    def test_budget_is_enforced(self) -> None:
        """A merged file above ten rows is refused."""
        with tempfile.TemporaryDirectory() as tmp:
            many = [f"P,chrT,{9000 + i},N,<DEL>,,,,,{0.5 - i / 100:.2f},secondary,x" for i in range(9)]
            result = run_merge(Path(tmp), SMALL, many)
            self.assertEqual(result.returncode, 1)
            self.assertIn("budget", result.stderr)


if __name__ == "__main__":
    unittest.main()
