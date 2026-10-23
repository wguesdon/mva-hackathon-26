#!/usr/bin/env python3
"""Read a class A gate file of a class B2 run, after the class A publish gate passes.

Control 7 of `docs/llm_data_privacy.md` says that a file which a class B2 pipeline publishes to a
class A path is built from a fixed key list. The pipeline runs the gate inside the task that writes
the file. This script is the second place where the gate runs, at the read boundary: an agent that
reads such a file reads it here, and the content reaches the terminal only after a PASS.

On a refusal the script prints the path and the reason and nothing from the file, and it exits with
code 2. The reason names a key path and never a value.

Usage:
    uv run python scripts/read_class_a_gate_file.py <path> [<path> ...]
    uv run python scripts/read_class_a_gate_file.py --s3 s3://bucket/prefix/class_a_gates.json
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "harness"))

import privacy_guard  # noqa: E402  the harness path is set above

RULES = ROOT / "literature" / "class_a_publish_rules.json"


def download(uri: str, target: Path) -> Path:
    """Copy one S3 object to a local file.

    Args:
        uri: The `s3://` URI of the gate file.
        target: The directory that receives the file.

    Returns:
        The local path of the downloaded file, under the object base name.

    Raises:
        subprocess.CalledProcessError: If the copy fails.
    """
    local = target / uri.rsplit("/", 1)[-1]
    subprocess.run(["aws", "s3", "cp", uri, str(local)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    return local


def read_one(path: Path, rules: dict, name: str) -> list[str]:
    """Run the gate over one file and print its content on a pass.

    Args:
        path: The local file to read.
        rules: The parsed class A publish rules.
        name: The name to print, which is the S3 URI for a downloaded file.

    Returns:
        The problems that the gate reported. An empty list means the content was printed.
    """
    text = path.read_text(errors="replace")
    problems = privacy_guard.class_a_publish_problems(path, text, rules)
    if problems:
        print(f"REFUSED\t{name}", file=sys.stderr)
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        return problems
    print(f"PASS\t{name}")
    print(text, end="" if text.endswith("\n") else "\n")
    return []


def main() -> int:
    """Read every named class A gate file through the gate.

    Returns:
        Zero when every file passed, two when the gate refused a file.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+")
    parser.add_argument("--s3", action="store_true",
                        help="Treat each path as an s3:// URI and download it first.")
    parser.add_argument("--rules", type=Path, default=RULES)
    args = parser.parse_args()
    rules = privacy_guard.load_publish_rules(args.rules)
    refused = False
    with tempfile.TemporaryDirectory() as name:
        scratch = Path(name)
        for item in args.paths:
            local = download(item, scratch) if args.s3 else Path(item)
            refused = bool(read_one(local, rules, item)) or refused
    return 2 if refused else 0


if __name__ == "__main__":
    sys.exit(main())
