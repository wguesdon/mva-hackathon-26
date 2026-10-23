"""Refuse model calls that could send patient data, and record every call.

The rules come from docs/llm_data_privacy.md.

* Class A is public content: code, docs, literature, benchmark scenarios.
* Class B is an owner-approved export of findings: a short candidate list, gene rankings, HPO terms.
* Class C is the genome or the full clinical record. No model arm may receive class C.

The exported dispatcher uses conservative ceilings separate from competition permissions.
Subscription CLI arms (Claude Code Max, Codex Pro) have configured ceiling A.
Bedrock arms through OpenCode have ceiling B. A task declares its data class, and the guard
refuses a task above the ceiling of its arm. The guard also scans every outgoing prompt for raw
genomic records and private identifiers, whatever the declared class.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import fcntl
import hashlib
import json
import re
import sys
from pathlib import Path

# The class A publish gate of control 7 runs inside a pipeline container, where the project package
# and PyYAML are not installed. Both imports are therefore optional, and the gate uses neither.
try:
    if __package__:
        from . import paths
    else:
        import paths
except ImportError:  # pragma: no cover - the container case
    paths = None

CLASSES = ("A", "B", "C")
BEDROCK_PREFIX = "amazon-bedrock/"
CLASS_A_PUBLISH_RULES = "literature/class_a_publish_rules.json"
CALL_LOG = "tracking/model_call_log.tsv"
CALL_FIELDS = ("utc", "arm", "provider", "plan", "data_class", "prompt_sha256", "prompt_bytes",
                   "task_id", "run_id", "verdict")

# A VCF data line: chromosome, position, identifier, reference allele, alternate allele. The line can start
# inside a JSON string or after a space.
VCF_LINE = re.compile(r"(?:^|(?<=[\s\"'`]))(?:chr)?(?:[0-9]{1,2}|X|Y|M|MT)\t\d+\t[^\t\n]+\t[ACGTNacgtn]+\t[ACGTNacgtn.,<>*\[\]:0-9]+\t",
                      re.MULTILINE)
# A sample genotype column in VCF or a genotype table: GT followed by at least one colon field.
GENOTYPE_FIELD = re.compile(r"(?<![\w/|])(?:[0-9.][/|][0-9.])(?::[0-9.,]+){1,}")
# A FASTQ record: header, bases, plus line.
FASTQ_RECORD = re.compile(r"^@[^\n]+\n[ACGTN]{30,}\n\+[^\n]*\n", re.MULTILINE)
# A SAM alignment line: name, flag, reference, position, mapping quality, CIGAR.
SAM_LINE = re.compile(r"^[^\t@\n][^\t\n]*\t\d+\t[^\t\n]+\t\d+\t\d+\t(?:\*|[0-9MIDNSHP=X]+)\t", re.MULTILINE)
# VCF header lines that name a sample column.
VCF_HEADER = re.compile(r"^#CHROM\tPOS\tID\tREF\tALT", re.MULTILINE)
# One line of `aws s3 ls --recursive`: the date, the time, the byte count and the object key.
S3_LS_LINE = re.compile(r"^\s*\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2}\s+(\d+)\s+(\S+)\s*$", re.MULTILINE)
# One line of `aws s3 ls --human-readable`: the date, the time, the size, the unit and the key.
S3_LS_HUMAN_LINE = re.compile(
    r"^\s*\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2}\s+\d+(?:\.\d+)?\s+(?:Bytes?|[KMGT]iB|[KMGT]B)\s+(\S+)\s*$",
    re.MULTILINE | re.IGNORECASE)
# One line of `ls -l`: the mode, the link count, the owner, the group, the byte count and the name.
LS_LONG_LINE = re.compile(r"^[-dlbcps][rwxsStT-]{9}[.+]?\s+\d+\s+\S+\s+\S+\s+(\d+)\s+\S+\s+\S+\s+\S+\s+(\S+)\s*$",
                          re.MULTILINE)
# A size with a unit anywhere in a line, for example `12.4 MiB` or `3 KB`.
HUMAN_SIZE = re.compile(r"\b\d+(?:\.\d+)?\s?(?:Bytes?|[KMGT]iB|[KMGT]B)\b", re.IGNORECASE)
# The `--summarize` footer of `aws s3 ls`.
S3_TOTAL_SIZE = re.compile(r"^\s*Total Size:\s*\d+", re.MULTILINE | re.IGNORECASE)


def class_rank(data_class: str) -> int:
    """Return the order of a data class.

    Args:
        data_class: One of A, B or C.

    Returns:
        0 for A, 1 for B, 2 for C.

    Raises:
        ValueError: If the class is unknown.
    """
    if data_class not in CLASSES:
        raise ValueError("The data class must be A, B or C.")
    return CLASSES.index(data_class)


def arm_profile(backend: str, model: str) -> dict:
    """Describe the provider, plan and class ceiling of one arm.

    Args:
        backend: The CLI name.
        model: The model identifier given to the CLI.

    Returns:
        The provider, the plan and the ceiling.
    """
    if backend == "claude":
        return {"provider": "anthropic", "plan": "claude_max_training_off", "ceiling": "A"}
    if backend == "codex":
        return {"provider": "openai", "plan": "chatgpt_pro_training_off", "ceiling": "A"}
    if backend == "opencode" and model.startswith(BEDROCK_PREFIX):
        return {"provider": "aws_bedrock", "plan": "bedrock_commercial_retention_none", "ceiling": "B"}
    if backend == "opencode":
        return {"provider": model.split("/")[0], "plan": "unreviewed", "ceiling": "A"}
    raise ValueError(f"The arm {backend} has no privacy profile.")


def deny_tokens() -> list[str]:
    """Read the private identifiers that must never appear in a prompt.

    The project package is absent inside a pipeline container, and the class A publish gate runs
    there. The list is empty in that case, and the pattern scan still runs.

    Returns:
        The nonempty tokens of PRIVACY_DENY_TOKENS in `.env`, or an empty list with no project package.
    """
    if paths is None:
        return []
    raw = paths.local_env().get("PRIVACY_DENY_TOKENS", "")
    return [token.strip() for token in raw.split(",") if len(token.strip()) >= 4]


def scan(text: str, tokens: list[str] | None = None) -> list[str]:
    """Find raw genomic records and private identifiers in a text.

    Args:
        text: The outgoing prompt or file content.
        tokens: The private identifiers. The `.env` tokens are used when None.

    Returns:
        One reason for each problem found. An empty list means the text passed.
    """
    problems = []
    # Forum prompts and reviewer answers are JSON. JSON escapes tabs and newlines, so the scan
    # also reads the text with those escapes decoded.
    if "\\t" in text or "\\n" in text:
        decoded = text.replace("\\t", "\t").replace("\\n", "\n")
        problems.extend(scan_records(decoded))
    problems.extend(scan_records(text))
    lowered = text.lower()
    for token in deny_tokens() if tokens is None else tokens:
        if token.lower() in lowered:
            problems.append("The text contains a private identifier from PRIVACY_DENY_TOKENS.")
            break
    return sorted(set(problems))


def scan_records(text: str) -> list[str]:
    """Find raw genomic record formats in a text.

    Args:
        text: The text to scan.

    Returns:
        One reason for each record format found.
    """
    problems = []
    if VCF_HEADER.search(text) or VCF_LINE.search(text):
        problems.append("The text contains a VCF record.")
    if GENOTYPE_FIELD.search(text):
        problems.append("The text contains a genotype field.")
    if FASTQ_RECORD.search(text):
        problems.append("The text contains a FASTQ record.")
    if SAM_LINE.search(text):
        problems.append("The text contains a SAM alignment record.")
    return problems


def load_publish_rules(path: Path) -> dict:
    """Read the class A publish rules of control 7.

    Args:
        path: The rules file, normally `literature/class_a_publish_rules.json`.

    Returns:
        The parsed rules.

    Raises:
        FileNotFoundError: If the rules file is absent. The gate fails closed.
    """
    return json.loads(path.read_text())


def flatten_json(value, prefix: str = "") -> list[tuple[str, object]]:
    """List every leaf of a parsed JSON document with its dotted key path.

    A list contributes the key path of the list itself for each of its elements, so a list of
    strings and a list of objects both keep the key path of their parent.

    Args:
        value: The parsed JSON value.
        prefix: The key path of this value.

    Returns:
        One pair of key path and leaf value for each leaf.
    """
    if isinstance(value, dict):
        leaves: list[tuple[str, object]] = []
        for key, item in value.items():
            leaves.extend(flatten_json(item, f"{prefix}.{key}" if prefix else str(key)))
        return leaves
    if isinstance(value, list):
        leaves = []
        for item in value:
            leaves.extend(flatten_json(item, prefix))
        return leaves
    return [(prefix, value)]


def path_matches(key_path: str, pattern: str) -> bool:
    """Test one key path against one allowed pattern.

    A pattern segment `*` matches exactly one segment. A trailing `**` matches one or more segments.

    Args:
        key_path: The dotted key path of a leaf.
        pattern: The allowed pattern from the rules file.

    Returns:
        True when the pattern covers the key path.
    """
    keys = key_path.split(".")
    parts = pattern.split(".")
    if parts and parts[-1] == "**":
        if len(keys) <= len(parts) - 1:
            return False
        keys = keys[:len(parts) - 1] + ["**"]
    if len(keys) != len(parts):
        return False
    return all(part in ("*", "**") or part == key for part, key in zip(parts, keys))


def quantity_problems(leaves: list[tuple[str, object]], entry: dict, pattern: re.Pattern) -> list[str]:
    """Find an undeclared count, size or length in a class A file.

    Args:
        leaves: The key path and value of every leaf of the file.
        entry: The rules entry of this file name.
        pattern: The compiled quantity key pattern.

    Returns:
        One reason for each undeclared quantity.
    """
    declared = entry.get("declared_quantities", {})
    problems = []
    for key_path, value in leaves:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            continue
        leaf = key_path.split(".")[-1]
        if pattern.search(leaf) and key_path not in declared:
            problems.append(f"The key {key_path} is an undeclared quantity in a class A file.")
    return problems


def key_path_problems(leaves: list[tuple[str, object]], entry: dict) -> list[str]:
    """Find a key path that the fixed key list of a class A file does not allow.

    Args:
        leaves: The key path and value of every leaf of the file.
        entry: The rules entry of this file name.

    Returns:
        One reason for each key path outside the list.
    """
    allowed = entry.get("allowed_key_paths", [])
    seen = []
    for key_path, _ in leaves:
        if key_path in seen:
            continue
        if not any(path_matches(key_path, pattern) for pattern in allowed):
            seen.append(key_path)
    return [f"The key {key_path} is not on the fixed key list of this class A file." for key_path in seen]


def listing_problems(text: str, prefixes: list[str], size_pattern: re.Pattern) -> list[str]:
    """Find an object size or a byte count printed for a class B2 prefix.

    Control 8 forbids such a listing. The function reads the output of a listing command, or any
    text that carries one.

    Args:
        text: The listing output or the file content.
        prefixes: The class B2 path prefixes of the rules file.
        size_pattern: The compiled size key pattern.

    Returns:
        One reason for each line that carries a size beside a class B2 path.
    """
    problems = []
    for _, key in S3_LS_LINE.findall(text):
        if any(prefix in key for prefix in prefixes):
            problems.append(f"A listing line prints the byte count of the class B2 object {key}.")
    for key in S3_LS_HUMAN_LINE.findall(text):
        if any(prefix in key for prefix in prefixes):
            problems.append(f"A listing line prints a human readable size for the class B2 object {key}.")
    for _, name in LS_LONG_LINE.findall(text):
        if any(prefix in name for prefix in prefixes):
            problems.append(f"A long listing prints the byte count of the class B2 object {name}.")
    for line in text.splitlines():
        if not any(prefix in line for prefix in prefixes):
            continue
        if size_pattern.search(line) or HUMAN_SIZE.search(line):
            problems.append("A line carries a size beside a class B2 path.")
    if S3_TOTAL_SIZE.search(text) and any(prefix in text for prefix in prefixes):
        problems.append("A listing prints the total size of a class B2 prefix.")
    return sorted(set(problems))


def class_a_publish_problems(path: Path, text: str, rules: dict) -> list[str]:
    """Run the class A publish gate of control 7 over one published file.

    Args:
        path: The file that a class B2 pipeline publishes to a class A path.
        text: The content of that file.
        rules: The parsed rules of `literature/class_a_publish_rules.json`.

    Returns:
        One reason for each problem. An empty list means the file passed.
    """
    prefixes = rules.get("class_b2_prefixes", [])
    size_pattern = re.compile(rules.get("size_key_pattern", "(?!)"), re.IGNORECASE)
    problems = scan(text)
    problems.extend(listing_problems(text, prefixes, size_pattern))
    entry = rules.get("files", {}).get(path.name)
    if entry is None:
        problems.append(f"The file name {path.name} has no entry in the class A publish rules, so the "
                        "gate fails closed.")
        return sorted(set(problems))
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        problems.append(f"The file {path.name} is declared in the class A publish rules and it is not "
                        "valid JSON, so its key list cannot be checked.")
        return sorted(set(problems))
    leaves = flatten_json(parsed)
    problems.extend(key_path_problems(leaves, entry))
    problems.extend(quantity_problems(leaves, entry,
                                      re.compile(rules.get("quantity_key_pattern", "(?!)"), re.IGNORECASE)))
    return sorted(set(problems))


def check_call(backend: str, model: str, data_class: str, prompt: str) -> dict:
    """Refuse a model call above the arm ceiling or with a private record.

    Args:
        backend: The CLI name.
        model: The model identifier.
        data_class: The declared data class of the task.
        prompt: The exact outgoing prompt.

    Returns:
        The arm profile when the call passes.

    Raises:
        ValueError: If the class is C, the class is above the ceiling, or the scan finds a problem.
    """
    profile = arm_profile(backend, model)
    if class_rank(data_class) >= class_rank("C"):
        raise ValueError("No model arm may receive class C content.")
    if class_rank(data_class) > class_rank(profile["ceiling"]):
        raise ValueError(f"The arm {backend} has ceiling {profile['ceiling']} and cannot take class {data_class}.")
    problems = scan(prompt)
    if problems:
        raise ValueError("The privacy scan refused the prompt. " + " ".join(sorted(set(problems))))
    return profile


def record_call(root: Path, backend: str, model: str, data_class: str, prompt: str,
                task_id: str, run_id: str, verdict: str) -> dict:
    """Append one row to the model call log.

    Args:
        root: The project root.
        backend: The CLI name.
        model: The model identifier.
        data_class: The declared data class.
        prompt: The exact outgoing prompt, or the refused prompt.
        task_id: The task identifier.
        run_id: The dispatch identifier.
        verdict: SENT or REFUSED.

    Returns:
        The written row.
    """
    try:
        profile = arm_profile(backend, model)
    except ValueError:
        profile = {"provider": "unknown", "plan": "unknown"}
    row = {
        "utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
        "arm": f"{backend}:{model}", "provider": profile["provider"], "plan": profile["plan"],
        "data_class": data_class, "prompt_sha256": hashlib.sha256(prompt.encode()).hexdigest(),
        "prompt_bytes": str(len(prompt.encode())), "task_id": task_id, "run_id": run_id, "verdict": verdict,
    }
    path = root / CALL_LOG
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", newline="") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        writer = csv.DictWriter(handle, fieldnames=CALL_FIELDS, delimiter="\t")
        if handle.tell() == 0:
            writer.writeheader()
        writer.writerow(row)
    return row


def guard_call(root: Path, backend: str, model: str, data_class: str, prompt: str,
               task_id: str, run_id: str) -> dict:
    """Check a call and record the decision in the model call log.

    Args:
        root: The project root.
        backend: The CLI name.
        model: The model identifier.
        data_class: The declared data class.
        prompt: The exact outgoing prompt.
        task_id: The task identifier.
        run_id: The dispatch identifier.

    Returns:
        The arm profile when the call passes.

    Raises:
        ValueError: If the guard refuses the call. The refusal is recorded first.
    """
    try:
        profile = check_call(backend, model, data_class, prompt)
    except ValueError:
        record_call(root, backend, model, data_class, prompt, task_id, run_id, "REFUSED")
        raise
    record_call(root, backend, model, data_class, prompt, task_id, run_id, "SENT")
    return profile


def rules_path(given: Path | None) -> Path:
    """Locate the class A publish rules.

    Args:
        given: The path from the command line, or None.

    Returns:
        The given path, the project path when the project tree is present, or the rules file beside
        this script.
    """
    if given is not None:
        return given
    here = Path(__file__).resolve().parent
    candidate = here.parent / CLASS_A_PUBLISH_RULES
    if candidate.exists():
        return candidate
    return here / Path(CLASS_A_PUBLISH_RULES).name


def main_gate(args: argparse.Namespace) -> int:
    """Run the class A publish gate or the listing gate over each file.

    The gate prints the verdict, the path and the reason of each problem. It never prints a value
    from a file that it refuses.

    Args:
        args: The parsed command line.

    Returns:
        Zero when every file passes, two when a file fails.
    """
    rules = load_publish_rules(rules_path(args.rules))
    prefixes = rules.get("class_b2_prefixes", [])
    size_pattern = re.compile(rules.get("size_key_pattern", "(?!)"), re.IGNORECASE)
    failed = False
    for path in args.files:
        text = sys.stdin.read() if str(path) == "-" else path.read_text(errors="replace")
        if args.listing:
            problems = listing_problems(text, prefixes, size_pattern)
        else:
            problems = class_a_publish_problems(path, text, rules)
        print(f"{'FAIL' if problems else 'PASS'}\t{path}\t{' '.join(problems)}".rstrip())
        failed = failed or bool(problems)
    return 2 if failed else 0


def main() -> int:
    """Scan files from the command line.

    Returns:
        Zero when every file passes, one when a file fails.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument("--allowlist", type=Path, default=None,
                        help="A YAML file mapping a path to the reason its content is not patient data.")
    parser.add_argument("--quiet-pass", action="store_true", help="Print a failed or allowed file only.")
    parser.add_argument("--class-a-publish", action="store_true",
                        help="Run the class A publish gate of control 7 over each file. A file name "
                             "with no entry in the rules file fails.")
    parser.add_argument("--listing", action="store_true",
                        help="Run the listing gate of control 8 over each file, or over standard "
                             "input when the file is a single dash.")
    parser.add_argument("--rules", type=Path, default=None,
                        help="The class A publish rules. The default is "
                             f"{CLASS_A_PUBLISH_RULES} beside this harness.")
    args = parser.parse_args()
    if args.class_a_publish or args.listing:
        return main_gate(args)
    allowed = {}
    if args.allowlist and args.allowlist.exists():
        import yaml
        allowed = yaml.safe_load(args.allowlist.read_text()) or {}
    failed = False
    for path in args.files:
        problems = scan(path.read_text(errors="replace"))
        key = str(path)
        if problems and key in allowed:
            print(f"ALLOWED\t{path}\t{allowed[key]}")
            continue
        if problems or not args.quiet_pass:
            print(f"{'FAIL' if problems else 'PASS'}\t{path}\t{' '.join(problems)}")
        failed = failed or bool(problems)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
