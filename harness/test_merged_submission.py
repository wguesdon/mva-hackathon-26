"""The merged submission check refuses a shared EPCR and it publishes a class A gate file.

Item B13 of `reports/audit/submission_audit_plan.md`. Two writers claim the ten rows of the
organizer file. The format check inside `pipelines/track1_b0/main.nf` reads the B0 file alone, so a
collision between the two writers and a repeat inside the structural rows both pass it. The measured
cost of one shared value is in `tracking/receipts/s4_row_count_calibration_v1.json`, field
`tie_counterfactual`: mean F-max falls from 0.960 to 0.238 over the 50 spike-in scenarios.

The fixtures of this test are synthetic. No patient row and no real submission file enters it.
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import privacy_guard

ROOT = Path(__file__).resolve().parent.parent
CHECKER = ROOT / "scripts" / "check_merged_submission.py"
RULES = privacy_guard.load_publish_rules(ROOT / "literature" / "class_a_publish_rules.json")
HEADER = ("proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type,notes")

# The two collisions that tracking/receipts/t1_epcr_ladder_v1.json enumerates. The B0 suite ladder
# starts at 0.95 and falls by 0.09, and the structural map reaches 0.77 at score 94 and 0.68 at
# score 76.
B0_SUITE_LADDER = (0.95, 0.86, 0.77)
STRUCTURAL_COLLISION = 0.77


def small_variant_row(value: float) -> str:
    """Build one synthetic B0 row at one EPCR.

    Args:
        value: The EPCR of the row.

    Returns:
        The CSV line.
    """
    return (f"PROBAND01,chr1,1000,A,G,chr1,2000,C,T,{value:.3f},primary,"
            "synthetic pair of the test fixture")


def structural_row(value: float) -> str:
    """Build one synthetic structural row at one EPCR.

    Args:
        value: The EPCR of the row.

    Returns:
        The CSV line.
    """
    return (f"PROBAND01,chr2,3000,N,<DEL>,,,,,{value:.3f},secondary,"
            "synthetic structural call of the test fixture")


def run_checker(lines: list[str], work: Path) -> tuple[int, dict]:
    """Run the merged submission checker over one synthetic file.

    Args:
        lines: The CSV rows without the header.
        work: The directory that holds the fixture and the gate file.

    Returns:
        The exit status and the gate file content.
    """
    submission = work / "merged.csv"
    submission.write_text("\n".join([HEADER, *lines]) + "\n")
    gate = work / "merged_format_check.json"
    result = subprocess.run([sys.executable, str(CHECKER), "--submission", str(submission),
                             "--out", str(gate)], capture_output=True, text=True)
    return result.returncode, json.loads(gate.read_text())


class MergedSubmissionTest(unittest.TestCase):
    """The checker separates a clean ladder from each collision class."""

    def test_a_distinct_decreasing_ladder_passes(self) -> None:
        """Ten distinct values with every primary above every secondary read PASS."""
        with tempfile.TemporaryDirectory() as name:
            work = Path(name)
            rows = [small_variant_row(value) for value in (0.99, 0.98, 0.97)]
            rows += [structural_row(value) for value in (0.835, 0.80, 0.77, 0.70, 0.68, 0.65, 0.60)]
            status, gate = run_checker(rows, work)
        self.assertEqual(status, 0)
        self.assertEqual(gate["verdict"], "PASS")
        self.assertEqual(gate["epcr_values_are_distinct"], "PASS")
        self.assertEqual(gate["primary_above_every_secondary"], "PASS")
        self.assertEqual(gate["both_writers_present"], "PASS")

    def test_a_shared_value_between_the_two_writers_fails(self) -> None:
        """The reachable collision at 0.77 reads FAIL on the distinctness condition."""
        with tempfile.TemporaryDirectory() as name:
            work = Path(name)
            rows = [small_variant_row(value) for value in B0_SUITE_LADDER]
            rows += [structural_row(STRUCTURAL_COLLISION)]
            status, gate = run_checker(rows, work)
        self.assertEqual(status, 1)
        self.assertEqual(gate["epcr_values_are_distinct"], "FAIL")
        self.assertEqual(gate["verdict"], "FAIL")

    def test_a_repeat_inside_the_structural_rows_fails(self) -> None:
        """Two structural calls of one score write one EPCR twice and the gate refuses it."""
        with tempfile.TemporaryDirectory() as name:
            work = Path(name)
            rows = [small_variant_row(value) for value in (0.99, 0.98, 0.97)]
            rows += [structural_row(0.79), structural_row(0.79)]
            status, gate = run_checker(rows, work)
        self.assertEqual(status, 1)
        self.assertEqual(gate["epcr_values_are_distinct"], "FAIL")

    def test_a_structural_row_above_a_b0_row_fails(self) -> None:
        """A secondary row above a primary row reads FAIL on the separation condition."""
        with tempfile.TemporaryDirectory() as name:
            work = Path(name)
            rows = [small_variant_row(0.70), structural_row(0.835)]
            status, gate = run_checker(rows, work)
        self.assertEqual(status, 1)
        self.assertEqual(gate["primary_above_every_secondary"], "FAIL")

    def test_the_gate_file_passes_the_class_a_publish_gate(self) -> None:
        """Gate G1 of ADR 0008. The gate file carries no row count and no EPCR value."""
        with tempfile.TemporaryDirectory() as name:
            work = Path(name)
            rows = [small_variant_row(value) for value in (0.99, 0.98, 0.97)]
            rows += [structural_row(0.835)]
            _, gate = run_checker(rows, work)
        problems = privacy_guard.class_a_publish_problems(
            Path("merged_format_check.json"), json.dumps(gate, indent=2), RULES)
        self.assertEqual(problems, [])
        text = json.dumps(gate)
        self.assertNotIn("0.835", text)
        self.assertNotIn("0.99", text)


if __name__ == "__main__":
    unittest.main(verbosity=2)
