"""Refuse an AWS launch when project spend plus the job estimate exceeds the monthly cap.

The account holds other projects. The guard therefore reads only the Cost Explorer spend of
resources tagged Project=<project tag> (the `Project` cost allocation tag is active). Cost
Explorer data lags by up to 24 hours, so the guard also adds the estimates of launches that the
cost log records since the start of yesterday (UTC). Launches from earlier days are already in the
billed spend.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
PROJECT_TAG = "mva-hackathon-2026"
COST_LOG = "tracking/cost_log.tsv"
LOG_FIELDS = ("utc", "job", "estimate_usd", "project_mtd_usd", "recent_estimates_usd", "cap_usd", "verdict")

def load_limits(root: Path = ROOT) -> dict:
    """Read the cost cap and the warning level from config.yaml.

    Args:
        root: The project root.

    Returns:
        The `aws` section of the configuration.
    """
    return yaml.safe_load((root / "config.yaml").read_text())["aws"]


def month_bounds(today: dt.date) -> tuple[str, str]:
    """Return the Cost Explorer period from the first day of the month to tomorrow.

    Args:
        today: The current UTC date.

    Returns:
        The start and the exclusive end in ISO format.
    """
    start = today.replace(day=1)
    return start.isoformat(), (today + dt.timedelta(days=1)).isoformat()


def project_spend(today: dt.date, project: str = PROJECT_TAG) -> float:
    """Read the month-to-date unblended cost of one project tag.

    Args:
        today: The current UTC date.
        project: The value of the Project tag.

    Returns:
        The spend in USD.
    """
    start, end = month_bounds(today)
    command = ["aws", "ce", "get-cost-and-usage", "--time-period", f"Start={start},End={end}",
               "--granularity", "MONTHLY", "--metrics", "UnblendedCost",
               "--filter", json.dumps({"Tags": {"Key": "Project", "Values": [project]}}), "--output", "json"]
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    return parse_spend(result.stdout)


def parse_spend(output: str) -> float:
    """Sum the unblended cost of a Cost Explorer response.

    Args:
        output: The JSON output of get-cost-and-usage.

    Returns:
        The total in USD.
    """
    payload = json.loads(output)
    return sum(float(item["Total"]["UnblendedCost"]["Amount"]) for item in payload["ResultsByTime"])


def recent_estimates(root: Path, now: dt.datetime) -> float:
    """Sum the estimates of allowed launches since the start of yesterday (UTC).

    Args:
        root: The project root.
        now: The current UTC time.

    Returns:
        The total estimate in USD.
    """
    path = root / COST_LOG
    if not path.exists():
        return 0.0
    total = 0.0
    window_start = dt.datetime.combine(now.date() - dt.timedelta(days=1), dt.time(), tzinfo=dt.UTC)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["verdict"] != "ALLOWED":
                continue
            if dt.datetime.fromisoformat(row["utc"]) >= window_start:
                total += float(row["estimate_usd"])
    return total


def decide(spend: float, recent: float, estimate: float, limits: dict) -> tuple[str, str]:
    """Compare the projected spend with the cap.

    Args:
        spend: The month-to-date project spend.
        recent: The estimates of recent launches that Cost Explorer may not show yet.
        estimate: The estimate of the new job.
        limits: The cap and warning level.

    Returns:
        The verdict ALLOWED, WARNED or REFUSED, and a message.
    """
    projected = spend + recent + estimate
    cap = float(limits["monthly_cost_cap_usd"])
    if projected > cap:
        return "REFUSED", f"Projected spend {projected:.2f} USD exceeds the monthly cap of {cap:.2f} USD."
    if projected > float(limits["monthly_cost_warn_usd"]):
        return "WARNED", f"Projected spend {projected:.2f} USD is above the warning level. The cap is {cap:.2f} USD."
    return "ALLOWED", f"Projected spend {projected:.2f} USD is inside the monthly cap of {cap:.2f} USD."


def check(root: Path, job: str, estimate: float, spend: float | None = None,
          now: dt.datetime | None = None) -> tuple[str, str]:
    """Run the guard for one launch and append the result to the cost log.

    Args:
        root: The project root.
        job: The job name.
        estimate: The job estimate in USD.
        spend: A month-to-date spend for tests. Cost Explorer is read when None.
        now: The current time for tests.

    Returns:
        The verdict and the message. A WARNED launch is recorded as ALLOWED.
    """
    now = now or dt.datetime.now(dt.UTC)
    limits = load_limits(root)
    spend = project_spend(now.date()) if spend is None else spend
    recent = recent_estimates(root, now)
    verdict, message = decide(spend, recent, estimate, limits)
    stored = "ALLOWED" if verdict == "WARNED" else verdict
    path = root / COST_LOG
    path.parent.mkdir(parents=True, exist_ok=True)
    new_file = not path.exists()
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LOG_FIELDS, delimiter="\t")
        if new_file:
            writer.writeheader()
        writer.writerow({"utc": now.isoformat(timespec="seconds"), "job": job, "estimate_usd": f"{estimate:.2f}",
                         "project_mtd_usd": f"{spend:.2f}", "recent_estimates_usd": f"{recent:.2f}",
                         "cap_usd": f"{float(limits['monthly_cost_cap_usd']):.2f}", "verdict": stored})
    return verdict, message


def main() -> int:
    """Check one launch or print the project spend.

    Returns:
        Zero when allowed, two when refused.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    launch = sub.add_parser("check")
    launch.add_argument("--job", required=True)
    launch.add_argument("--estimate-usd", type=float, required=True)
    sub.add_parser("status")
    args = parser.parse_args()
    if args.command == "status":
        limits = load_limits()
        spend = project_spend(dt.datetime.now(dt.UTC).date())
        print(f"Project {PROJECT_TAG} month-to-date spend is {spend:.2f} USD. The cap is "
              f"{float(limits['monthly_cost_cap_usd']):.2f} USD.")
        return 0
    verdict, message = check(ROOT, args.job, args.estimate_usd)
    print(f"{verdict}: {message}")
    return 2 if verdict == "REFUSED" else 0


if __name__ == "__main__":
    sys.exit(main())
