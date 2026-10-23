"""Run every harness test script and report the result.

Each test file runs as a plain script. Each one prints its own
result and exits non-zero on failure. This runner executes each script in a subprocess, so
one crash cannot hide the rest.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HARNESS = Path(__file__).resolve().parent


SKIP_EXIT = 77


def run_one(path: Path) -> tuple[str, str, str]:
    """Run one test script.

    Args:
        path: The test script path.

    Returns:
        The script name, the status PASS, FAIL or SKIP, and the last output line.
    """
    result = subprocess.run(
        [sys.executable, str(path)],
        capture_output=True,
        text=True,
        cwd=str(HARNESS),
        timeout=900,
    )
    output = (result.stdout + result.stderr).strip().splitlines()
    tail = output[-1] if output else "(no output)"
    if result.returncode == SKIP_EXIT:
        skip_line = next((line for line in output if line.startswith("SKIP:")), tail)
        return path.name, "SKIP", skip_line
    return path.name, "PASS" if result.returncode == 0 else "FAIL", tail


def main() -> int:
    """Run every test script and print a summary.

    Returns:
        0 when every script passed, 1 otherwise.
    """
    scripts = sorted(HARNESS.glob("test_*.py"))
    if not scripts:
        print("No test script found.")
        return 1
    failures = []
    skipped = []
    for script in scripts:
        name, status, tail = run_one(script)
        print(f"  {status}  {name:<48} {tail[:70]}")
        if status == "FAIL":
            failures.append(name)
        elif status == "SKIP":
            skipped.append(name)
    passed = len(scripts) - len(failures) - len(skipped)
    print()
    print(f"{passed} passed, {len(skipped)} skipped, {len(failures)} failed "
          f"of {len(scripts)} test scripts.")
    if skipped:
        print("skipped: " + ", ".join(skipped))
        print("A skipped guard activates once the harness holds the data it needs.")
    if failures:
        print("failed: " + ", ".join(failures))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
