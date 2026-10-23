"""Verify the cost guard decisions and the cost log without a Cost Explorer call."""

from __future__ import annotations

import datetime as dt
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import cost_guard

LIMITS = {"monthly_cost_cap_usd": 300, "monthly_cost_warn_usd": 200}


class CostGuardTest(unittest.TestCase):
    """Check verdicts, the recent-estimate window and the Cost Explorer parser."""

    def setUp(self) -> None:
        """Create a temporary root with the project config."""
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        shutil.copy(cost_guard.ROOT / "config.yaml", self.root / "config.yaml")
        self.now = dt.datetime(2026, 9, 13, 12, tzinfo=dt.UTC)

    def test_decide(self) -> None:
        """Spend under the warning allows, above it warns, above the cap refuses."""
        self.assertEqual(cost_guard.decide(10, 0, 20, LIMITS)[0], "ALLOWED")
        self.assertEqual(cost_guard.decide(190, 0, 20, LIMITS)[0], "WARNED")
        self.assertEqual(cost_guard.decide(250, 40, 20, LIMITS)[0], "REFUSED")

    def test_recent_estimates_count_until_cost_explorer_catches_up(self) -> None:
        """Two launches of 140 USD refuse a third launch of 30 USD on the same day."""
        self.assertEqual(cost_guard.check(self.root, "a", 140, spend=0, now=self.now)[0], "ALLOWED")
        self.assertEqual(cost_guard.check(self.root, "b", 140, spend=0, now=self.now)[0], "WARNED")
        self.assertEqual(cost_guard.check(self.root, "c", 30, spend=0, now=self.now)[0], "REFUSED")
        next_day = self.now + dt.timedelta(days=1)
        self.assertEqual(cost_guard.recent_estimates(self.root, next_day), 280.0)
        two_days_later = self.now + dt.timedelta(days=2)
        self.assertEqual(cost_guard.recent_estimates(self.root, two_days_later), 0.0)

    def test_parse_spend(self) -> None:
        """The parser sums the monthly totals."""
        output = json.dumps({"ResultsByTime": [{"Total": {"UnblendedCost": {"Amount": "12.5", "Unit": "USD"}}}]})
        self.assertEqual(cost_guard.parse_spend(output), 12.5)

    def test_month_bounds(self) -> None:
        """The period starts on the first day and ends tomorrow."""
        self.assertEqual(cost_guard.month_bounds(dt.date(2026, 9, 30)), ("2026-09-01", "2026-10-01"))


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
