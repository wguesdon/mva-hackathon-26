#!/usr/bin/env python3
"""Verify every quantitative claim of a report against the receipt it names. Data class A.

The script answers defect R2-12 and the class of defect that the third jury round named: a
statement added to satisfy a reviewer creates a new defect when it contradicts another part of the
report. A number that a reader meets in the report has to come from a record, and the same quantity
has to carry the same value everywhere.

The checker runs five checks over the build source of one report.

    1. Receipt coverage. A section of the report that carries a claim value has to name a record that
       a reader can open. A record is a receipt under ``tracking/receipts/``, a benchmark record under
       ``docs/benchmarks/``, a named table of the report whose caption names one, or an entry of the
       exemption list with a written reason. The unit is the section, because a reader who meets a
       number in a section resolves the citation of that section.
    2. Receipt verification. Every claim value resolves against a leaf of a receipt that its section
       names. The checker reports the receipt and the JSON path of the leaf it matched. A match says
       that a cited record carries the value; it does not say that the field is the semantically
       intended one, which is why the declared quantities of check 4 carry the numbers that repeat.
    3. Derivation. A claim value that no leaf carries is verified against a declared formula over
       receipt leaves, for example a total that is the sum of three component costs.
    4. Consistency. Two claim values that resolve to the same receipt field, and two occurrences of a
       declared quantity, have to print the same value. The declared quantities are also read in the
       records outside the report that carry the same number.
    5. Prohibited phrase. A phrase that ``docs/report_style.md`` forbids does not appear in the
       report source.

A configuration may declare an open defect. An open defect is a claim value that the record refutes
and that another repair item owns. The entry carries the report value, the record value, the item
that owns the repair and the reason. The checker prints each one beside its record value and counts
it, and it does not exit non-zero on it, because the repair belongs to the item that the entry
names. An open defect with no owning item and no reason is not accepted.

A finding exits non-zero. This export supplies the checker but no completed report, claim
configuration or analysis receipts. Supply those inputs separately. The source-only test suite
does not run this checker against a manuscript.

Usage:
    uv run python scripts/check_report_claims.py <qmd> [<qmd> ...] --json <path> --quiet
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]

NUMBER = re.compile(r"(?<![\w.])(\d{1,3}(?:,\d{3})+|\d+)(?:\.(\d+))?(?:e-?(\d+))?(?![\w])")
RECEIPT_REFERENCE = re.compile(r"tracking/receipts/([A-Za-z0-9_.]+)\.json")
RECORD_REFERENCE = re.compile(r"docs/benchmarks/([A-Za-z0-9_.]+)\.md")
TABLE_DATA_REFERENCE = re.compile(r"tracking/([A-Za-z0-9_.]+)\.tsv")
BUILD_RECORD_REFERENCE = re.compile(
    r"\b((?:reports|literature)/[A-Za-z0-9_][A-Za-z0-9_./-]*\.(?:json|tsv|yaml|yml))")
TABLE_REFERENCE = re.compile(r"@(tbl|fig)-[A-Za-z0-9-]+")
TABLE_LABEL = re.compile(r"\{#(tbl|fig)-([A-Za-z0-9-]+)\}")
INCLUDE = re.compile(r"\{\{<\s*include\s+(\S+)\s*>\}\}")
BACKTICK = re.compile(r"`[^`]*`")
INLINE_MATH = re.compile(r"\$[^$]*\$")


def load_config(path: Path) -> dict:
    """Read the checker configuration.

    Args:
        path: The YAML configuration path.

    Returns:
        The configuration mapping.
    """
    with path.open() as handle:
        return yaml.safe_load(handle)


def resolve_include(qmd: Path, name: str) -> Path | None:
    """Give the file that one include directive contributes to the checked text.

    Control 6 of ``docs/llm_data_privacy.md`` and ADR 0008 forbid a class A tool from reading a
    class B2 file. A file whose name starts with ``_class_b2_`` is class B2, so this function
    returns the public stub of the report in its place, which is what
    ``writeups/build_public.sh`` renders and what ``scripts/stage_jury_review.sh`` stages. The
    function fails closed: when the stub is absent the include contributes nothing.

    Args:
        qmd: The build source path.
        name: The include target as the directive spells it.

    Returns:
        The path to read, or None when the include must not be read.
    """
    if Path(name).name.startswith("_class_b2_"):
        stub = qmd.parent / "_results_public_stub.qmd"
        return stub if stub.exists() else None
    target = qmd.parent / name
    return target if target.exists() else None


def read_source(qmd: Path) -> tuple[str, list[tuple[int, str]]]:
    """Read one report source and expand every include that a class A tool may read.

    Args:
        qmd: The build source path.

    Returns:
        The abstract text of the YAML header and the body as a list of line number and line text.
        An included file contributes its own lines under the line number of the include directive.
        A class B2 include contributes the public stub instead of its own text.
    """
    raw = qmd.read_text().splitlines()
    body: list[tuple[int, str]] = []
    header: list[str] = []
    in_header = False
    header_done = False
    for index, line in enumerate(raw, start=1):
        if index == 1 and line.strip() == "---":
            in_header = True
            continue
        if in_header and line.strip() == "---":
            in_header = False
            header_done = True
            continue
        if in_header:
            header.append(line)
            continue
        if not header_done and not in_header and index == 1:
            pass
        match = INCLUDE.search(line)
        if match:
            target = resolve_include(qmd, match.group(1))
            if target is not None:
                for included in target.read_text().splitlines():
                    body.append((index, included))
            continue
        body.append((index, line))
    parsed = yaml.safe_load("\n".join(header)) if header else {}
    abstract = parsed.get("abstract", "") if isinstance(parsed, dict) else ""
    return abstract, body


def abstract_line(qmd: Path) -> int:
    """Find the line of the abstract key in the YAML header.

    Args:
        qmd: The build source path.

    Returns:
        The one-based line number of the ``abstract:`` key, or 1 when it is absent.
    """
    for index, line in enumerate(qmd.read_text().splitlines(), start=1):
        if line.startswith("abstract:"):
            return index
    return 1


def split_blocks(abstract: str, body: list[tuple[int, str]],
                 abstract_at: int) -> list[dict]:
    """Group the report into blocks and name the section of each one.

    A block is one paragraph of prose, one table with its caption, or one sentence of the abstract.
    The abstract is one paragraph in the source and it carries many records, so each of its
    sentences is a block of its own and a section of its own, and each one has to name its own
    record.

    Args:
        abstract: The abstract text of the YAML header.
        body: The body lines with their line numbers.
        abstract_at: The line number of the abstract key.

    Returns:
        A list of blocks. Each block carries a kind, a line number, its text and its section key.
    """
    blocks: list[dict] = []
    for index, sentence in enumerate(split_sentences(abstract)):
        if sentence.strip():
            blocks.append({"kind": "abstract", "line": abstract_at, "text": sentence,
                           "section": f"Abstract, sentence {index + 1}"})
    current: list[str] = []
    start = 0
    kind = "paragraph"
    section = "Front matter"
    for line, text in body:
        stripped = text.strip()
        if not stripped or stripped.startswith("#"):
            if current:
                blocks.append({"kind": kind, "line": start, "text": "\n".join(current),
                               "section": section})
            current = []
            kind = "paragraph"
            if stripped.startswith("#"):
                section = f"{stripped.lstrip('# ').strip()} (line {line})"
                blocks.append({"kind": "heading", "line": line, "text": stripped,
                               "section": section})
            continue
        if not current:
            start = line
        if stripped.startswith("|") or stripped.startswith(": ") or TABLE_LABEL.search(stripped):
            kind = "table"
        current.append(text)
    if current:
        blocks.append({"kind": kind, "line": start, "text": "\n".join(current),
                       "section": section})
    return blocks


def split_sentences(text: str) -> list[str]:
    """Split a paragraph into sentences without cutting a decimal or an abbreviation.

    Args:
        text: The paragraph text.

    Returns:
        The sentences of the paragraph.
    """
    flat = " ".join(text.split())
    parts = re.split(r"(?<=[.;])\s+(?=[A-Z(`*“])", flat)
    return [part for part in parts if part.strip()]


def mask_exempt_spans(text: str) -> str:
    """Replace the spans that never hold a claim value with spaces.

    Code in backticks holds a path, a parameter, a checksum or a command, and inline math holds a
    formula. Neither is a claim about a result.

    Args:
        text: The block text.

    Returns:
        The text with every backtick span and every math span replaced by spaces of equal length.
    """
    def blank(match: re.Match[str]) -> str:
        return " " * len(match.group(0))

    return INLINE_MATH.sub(blank, BACKTICK.sub(blank, text))


def flatten(value: object, prefix: str = "") -> dict[str, object]:
    """Flatten a receipt into a map of JSON path to leaf value.

    Args:
        value: The receipt or a part of it.
        prefix: The path of the current node.

    Returns:
        A map of JSON path to leaf value. A list contributes one entry per element and one entry for
        its own length under the suffix ``#len``.
    """
    flat: dict[str, object] = {}
    if isinstance(value, dict):
        for key, item in value.items():
            flat.update(flatten(item, f"{prefix}.{key}" if prefix else key))
    elif isinstance(value, list):
        flat[f"{prefix}#len"] = len(value)
        for index, item in enumerate(value):
            flat.update(flatten(item, f"{prefix}[{index}]"))
    else:
        flat[prefix] = value
    return flat


def load_receipts(receipt_dir: Path,
                  raw: dict[str, object] | None = None) -> dict[str, dict[str, object]]:
    """Read every receipt and every result table and flatten each one.

    A result table under ``tracking/`` is a record that the report cites beside a receipt, so its
    cells are checkable values. Each cell becomes one leaf under the path ``column[row]``.

    Args:
        receipt_dir: The receipt directory.
        raw: A map that collects the parsed receipt of each identifier before it is flattened. A
            row derivation reads the case rows from it.

    Returns:
        A map of record identifier to its flattened leaves.
    """
    receipts: dict[str, dict[str, object]] = {}
    for path in sorted(receipt_dir.glob("*.json")):
        try:
            with path.open() as handle:
                parsed = json.load(handle)
        except (json.JSONDecodeError, OSError) as error:
            print(f"  WARN  {path.name} did not parse: {error}", file=sys.stderr)
            continue
        receipts[path.stem] = flatten(parsed)
        if raw is not None:
            raw[path.stem] = parsed
    for path in sorted(receipt_dir.parent.glob("*.tsv")):
        receipts[path.stem] = read_result_table(path)
    return receipts


def load_declared_records(config: dict, root: Path | None = None
                          ) -> dict[str, dict[str, object]]:
    """Read every build record that the configuration declares.

    A receipt under ``tracking/receipts/`` is read by identifier. A report of the build, for example
    ``reports/track2/graph/manifest.json``, and a rule file, for example a surface form list under
    ``literature/``, live outside that directory and the report cites each one by path. The
    configuration lists each such path, so the set of records that the value check may use stays
    explicit and a large table does not enter the check by accident. A JSON record, a YAML record
    and a tab separated record are all read.

    Args:
        config: The checker configuration.
        root: The directory that a declared path resolves against. The default is the repository
            root, and a test passes the directory of its fixtures.

    Returns:
        A map of record identifier to its flattened leaves. The identifier of a build record is its
        repository path with the extension removed.
    """
    base = root or ROOT
    records: dict[str, dict[str, object]] = {}
    for relative in config.get("records", []):
        path = base / relative
        if not path.exists():
            print(f"  WARN  the declared record {relative} is absent", file=sys.stderr)
            continue
        key = relative.rsplit(".", 1)[0]
        if path.suffix in (".json", ".yaml", ".yml"):
            try:
                with path.open() as handle:
                    parsed = (json.load(handle) if path.suffix == ".json"
                              else yaml.safe_load(handle))
            except (json.JSONDecodeError, yaml.YAMLError, OSError) as error:
                print(f"  WARN  {relative} did not parse: {error}", file=sys.stderr)
                continue
            records[key] = flatten(parsed)
        elif path.suffix == ".tsv":
            records[key] = read_result_table(path)
    return records


def read_result_table(path: Path) -> dict[str, object]:
    """Read one tab separated result table as a flat map of cell values.

    Args:
        path: The table path.

    Returns:
        A map of ``column[row]`` to the numeric value of that cell. A cell that is not a number is
        not returned, because the checker compares numbers.
    """
    cells: dict[str, object] = {}
    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return cells
        for row_index, row in enumerate(reader):
            for column, value in zip(header, row):
                text = value.strip().replace(",", "")
                if not NUMBER.fullmatch(text):
                    continue
                try:
                    cells[f"{column}[{row_index}]"] = float(text)
                except ValueError:
                    continue
    return cells


def claim_forms(leaf: object, decimals: int, scale_words: bool = True) -> set[str]:
    """Give every printed form of one receipt leaf at one precision.

    Args:
        leaf: The leaf value.
        decimals: The number of decimal places that the report prints.
        scale_words: Offer the fractional form of a value written in thousands, millions or
            billions. A report whose configuration sets ``scale_word_required`` passes True only
            where the sentence writes a scale word, because a count of 7,837 printed as 7.8, and a
            corpus of 41,148,566 records printed as 41, are statements only there. A report that
            does not set it passes True everywhere, which is the earlier behaviour.

    Returns:
        The set of strings that the report may carry for this leaf.
    """
    forms: set[str] = set()
    if isinstance(leaf, bool) or leaf is None:
        return forms
    if isinstance(leaf, str):
        for match in NUMBER.finditer(leaf):
            forms.add(match.group(0).replace(",", ""))
        return forms
    if not isinstance(leaf, (int, float)):
        return forms
    number = float(leaf)
    if math.isnan(number) or math.isinf(number):
        return forms
    forms.add(f"{number:.{decimals}f}")
    forms.add(f"{number:.{max(decimals, 0)}f}".rstrip("0").rstrip(".") or "0")
    if abs(number - round(number)) < 1e-12:
        forms.add(str(int(round(number))))
    exponent = f"{number:.{decimals}e}"
    forms.add(exponent)
    forms.add(exponent.replace("e-0", "e-").replace("e+0", "e"))
    if scale_words:
        for scale in (1e3, 1e6, 1e9):
            if abs(number) >= scale:
                forms.add(f"{number / scale:.{decimals}f}")
    return {form for form in forms if form}


def writes_a_scale_word(text: str) -> bool:
    """Say whether a claim sentence writes a number in thousands, millions or billions.

    Args:
        text: The block text that carries the claim.

    Returns:
        True when the text holds one of the three scale words.
    """
    lowered = text.lower()
    if any(word in lowered for word in ("thousand", "million", "billion")):
        return True
    return re.search(r"\d\s*(?:kb|mb|gb|tb|kilobase|megabase|gigabase)\b", lowered) is not None


def filter_leaf_paths(receipts: dict[str, dict[str, object]],
                      patterns: list[dict]) -> dict[str, dict[str, object]]:
    """Remove the receipt leaves that the configuration declares unreadable by the value check.

    A control receipt keeps one row per evaluated case so that another top k recomputes without a
    second run. Those rows are diagnostics and not reported quantities, and a receipt with several
    hundred of them matches almost any two digit or three digit number by accident. A report that
    declares such a path removes it from the value check and states the reason, and a claim that
    reads those rows resolves through a declared row derivation instead.

    Args:
        receipts: The flattened records.
        patterns: The ``ignored_leaf_paths`` entries of the configuration.

    Returns:
        The records with every matching leaf removed.
    """
    if not patterns:
        return receipts
    compiled = [re.compile(entry["pattern"]) for entry in patterns]
    kept: dict[str, dict[str, object]] = {}
    for record, leaves in receipts.items():
        kept[record] = {path: leaf for path, leaf in leaves.items()
                        if not any(rule.search(path) for rule in compiled)}
    return kept


def numbers_of(text: str) -> list[tuple[str, int, int]]:
    """Extract every numeric token of a block.

    Args:
        text: The block text with the exempt spans already masked.

    Returns:
        A list of the token string, the number of decimal places and the character offset.
    """
    found: list[tuple[str, int, int]] = []
    for match in NUMBER.finditer(text):
        decimals = len(match.group(2) or "")
        found.append((match.group(0), decimals, match.start()))
    return found


def is_exempt(token: str, text: str, offset: int, rules: list[dict]) -> str | None:
    """Decide whether one numeric token is a claim value or an exempt token.

    Args:
        token: The numeric token as the report prints it.
        text: The block text.
        offset: The character offset of the token inside the block text.
        rules: The exemption rules of the configuration.

    Returns:
        The reason of the first rule that matches, or None when the token is a claim value.
    """
    start = max(0, offset - 60)
    window = text[start:offset + len(token) + 30]
    inside = offset - start
    for rule in rules:
        pattern = rule.get("pattern")
        if not pattern:
            continue
        for candidate in re.finditer(pattern, window):
            if candidate.start() <= inside < candidate.end():
                return rule["reason"]
    return None


def wilson_bounds(successes: float, total: float, z: float = 1.959963984540054
                  ) -> tuple[float, float] | None:
    """Give the Wilson score interval of one proportion at 95 percent.

    Args:
        successes: The number of successes.
        total: The number of trials.
        z: The normal quantile of the confidence level. The default is 95 percent.

    Returns:
        The lower and the upper bound, or None when the inputs are not a proportion.
    """
    if total <= 0 or successes < 0 or successes > total:
        return None
    if successes != int(successes) or total != int(total):
        return None
    proportion = successes / total
    centre = proportion + z * z / (2.0 * total)
    spread = z * math.sqrt(proportion * (1.0 - proportion) / total
                           + z * z / (4.0 * total * total))
    scale = 1.0 + z * z / total
    return max(0.0, (centre - spread) / scale), min(1.0, (centre + spread) / scale)


def object_derivations(leaves: dict[str, object], parent: str) -> dict[str, str]:
    """Give every value that the numeric leaves of one receipt object derive.

    The rule is narrow on purpose. A candidate comes from two numeric leaves of the same object and
    nothing else, because a table row of this report reports a proportion, a Wilson interval and an
    F-measure of the counts printed in the same row.

    Args:
        leaves: The flattened receipt.
        parent: The path of the object.

    Returns:
        A map of the derived value, printed at four decimals, to the derivation that produced it.
    """
    numeric: dict[str, float] = {}
    for path, leaf in leaves.items():
        if not path.startswith(parent):
            continue
        tail = path[len(parent):].lstrip(".")
        if "." in tail or not tail or "#" in tail or "[" in tail:
            continue
        if isinstance(leaf, (int, float)) and not isinstance(leaf, bool):
            numeric[tail] = float(leaf)
    derived: dict[str, str] = {}
    for first, successes in numeric.items():
        for second, total in numeric.items():
            if first == second or total <= 0:
                continue
            bounds = wilson_bounds(successes, total)
            if bounds is not None:
                derived.setdefault(f"{bounds[0]:.4f}",
                                   f"wilson low of {parent}.{first} over {parent}.{second}")
                derived.setdefault(f"{bounds[1]:.4f}",
                                   f"wilson high of {parent}.{first} over {parent}.{second}")
            derived.setdefault(f"{successes / total:.4f}",
                               f"{parent}.{first} over {parent}.{second}")
            derived.setdefault(f"{100.0 * successes / total:.4f}",
                               f"{parent}.{first} over {parent}.{second} as a percentage")
            if successes + total > 0:
                derived.setdefault(f"{2.0 * successes / (successes + total):.4f}",
                                   f"harmonic form of {parent}.{first} and {parent}.{second}")
    for first, left in numeric.items():
        for second, right in numeric.items():
            if first == second:
                continue
            if left + right > 0:
                derived.setdefault(f"{2.0 * left * right / (left + right):.4f}",
                                   f"f measure of {parent}.{first} and {parent}.{second}")
            if left - right > 0:
                derived.setdefault(f"{left - right:.4f}",
                                   f"{parent}.{first} less {parent}.{second}")
                derived.setdefault(f"{left - right:.0f}",
                                   f"{parent}.{first} less {parent}.{second}")
    return derived


def derive_within_objects(value: str, decimals: int, named: set[str],
                          receipts: dict[str, dict[str, object]]) -> str | None:
    """Verify one claim value against a derivation inside a single receipt object.

    Args:
        value: The claim value as the report prints it.
        decimals: The number of decimal places that the report prints.
        named: The receipt identifiers that the section names.
        receipts: The flattened receipts.

    A claim printed with no decimal place is checked against a difference of two leaves of one
    object and against nothing else, because a ratio rounded to a whole number matches too many
    integers to carry evidence. A report whose receipts hold an object with dozens of numeric
    leaves, for example a count per observed surface form, sets ``automatic_object_derivations`` to
    false in its configuration, because a difference of two such leaves reaches almost any small
    integer.

    Returns:
        The derivation that produced the value, or None.
    """
    try:
        target = float(value)
    except ValueError:
        return None
    if decimals > 4:
        return None
    for receipt in sorted(named):
        leaves = receipts.get(receipt)
        if leaves is None:
            continue
        for parent in object_parents(leaves):
            for candidate, how in object_derivations(leaves, parent).items():
                if decimals == 0 and " less " not in how:
                    continue
                if abs(float(candidate) - target) < tolerance(decimals):
                    return f"{receipt}: {how}"
    return None


def object_parents(leaves: dict[str, object]) -> list[str]:
    """Give the path of every object of a flattened receipt.

    Args:
        leaves: The flattened receipt.

    Returns:
        The distinct parent paths, each one with its trailing separator removed.
    """
    parents = {path.rsplit(".", 1)[0] for path in leaves if "." in path}
    return sorted(parents)


def rows_of(field: str, raw: dict[str, object]) -> list[dict] | None:
    """Read one list of records out of a raw receipt.

    Args:
        field: The reference in the form ``receipt#dotted.path``.
        raw: The raw receipts, keyed by identifier.

    Returns:
        The list of mappings at that path, or None when the path is absent or holds another type.
    """
    receipt, path = field.split("#", 1)
    node: object = raw.get(receipt)
    for segment in path.split("."):
        if not isinstance(node, dict) or segment not in node:
            return None
        node = node[segment]
    if not isinstance(node, list) or any(not isinstance(item, dict) for item in node):
        return None
    return node


def count_at_or_below(rows: list[dict], field: str, limit: float) -> int:
    """Count the records whose rank field is a reached rank at or below one limit.

    Args:
        rows: The records.
        field: The field name to read.
        limit: The inclusive upper limit.

    Returns:
        The number of records whose field holds a rank of at least 1 and at most the limit. A
        record whose field is absent, null, not a number, or zero counts as a miss, because a rank
        of zero is how a control records a case that reached no rank at all.
    """
    hits = 0
    for row in rows:
        value = row.get(field)
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            continue
        if 1 <= value <= limit:
            hits += 1
    return hits


def binomial_upper_tail(rate: float, cases: int, observed: int) -> float:
    """Give the upper tail of a Poisson binomial whose cases share one success rate.

    Args:
        rate: The per-case success probability.
        cases: The number of cases.
        observed: The observed number of successes.

    Returns:
        The probability of at least ``observed`` successes.
    """
    distribution = [1.0]
    for _ in range(cases):
        nxt = [0.0] * (len(distribution) + 1)
        for index, mass in enumerate(distribution):
            nxt[index] += mass * (1.0 - rate)
            nxt[index + 1] += mass * rate
        distribution = nxt
    return sum(distribution[min(observed, cases):])


def row_denominator(spec: dict, receipts: dict[str, dict[str, object]]) -> float | None:
    """Give the denominator that one row derivation declares, or None when it declares none.

    A rate over a stratum divides the recoveries of the whole set by the size of that stratum, and
    the receipt carries the size as a field rather than as a row count.

    Args:
        spec: The derivation entry.
        receipts: The flattened records.

    Returns:
        The denominator, or None when the entry declares none or a named field is absent.
    """
    fields = spec.get("denominator_fields")
    if not fields:
        return None
    return resolve_derivation(
        {"operation": spec.get("denominator_operation", "sum"), "fields": fields}, receipts, {})


def resolve_row_derivation(spec: dict, raw: dict[str, object],
                           receipts: dict[str, dict[str, object]]) -> float | None:
    """Compute one declared derivation that reads the case rows of a control receipt.

    The held-out control receipt keeps one row per case and one row per permuted case, so a rate at
    another top k recomputes from the rows without a second run. Three operations read those rows.
    ``rank_share`` gives the share of cases recovered inside k. ``smoothed_rank_share`` gives the
    Laplace smoothed share that the control uses as the null rate. ``poisson_binomial_tail`` gives
    the reported chance probability of the real count under that smoothed rate.

    Args:
        spec: The derivation entry.
        raw: The raw receipts, keyed by identifier.
        receipts: The flattened records, which a declared denominator reads.

    Returns:
        The derived value, or None when a named list is absent.
    """
    operation = spec["operation"]
    field = spec.get("field", "best_rank")
    limit = float(spec["max_rank"])
    rows = rows_of(spec["rows"], raw)
    if rows is None or not rows:
        return None
    hits = count_at_or_below(rows, field, limit)
    declared = row_denominator(spec, receipts)
    if operation == "rank_share":
        return hits / (declared if declared else len(rows))
    if operation == "smoothed_rank_share":
        return (hits + 1) / ((declared if declared else len(rows)) + 1)
    if operation != "poisson_binomial_tail":
        return None
    null_rows = rows_of(spec["null_rows"], raw)
    if null_rows is None or not null_rows:
        return None
    null_hits = count_at_or_below(null_rows, field, limit)
    return binomial_upper_tail((null_hits + 1) / (len(null_rows) + 1), len(rows), hits)


def resolve_derivation(spec: dict, receipts: dict[str, dict[str, object]],
                       raw: dict[str, object]) -> float | None:
    """Compute one declared derivation from receipt leaves or from receipt rows.

    Args:
        spec: The derivation entry with an operation and either a list of receipt fields or a
            reference to a list of case rows.
        receipts: The flattened receipts.
        raw: The raw receipts, keyed by identifier.

    Returns:
        The derived value, or None when a named field is absent.
    """
    operation = spec["operation"]
    factor = float(spec.get("factor", 1.0))
    if operation in ("rank_share", "smoothed_rank_share", "poisson_binomial_tail"):
        rows = resolve_row_derivation(spec, raw, receipts)
        return None if rows is None else rows * factor
    values: list[float] = []
    for field in spec["fields"]:
        receipt, path = field.split("#", 1)
        leaves = receipts.get(receipt)
        if leaves is None or path not in leaves:
            return None
        leaf = leaves[path]
        if not isinstance(leaf, (int, float)) or isinstance(leaf, bool):
            return None
        values.append(float(leaf))
    if operation == "sum":
        return factor * sum(values)
    if operation == "ratio" and len(values) == 2 and values[1]:
        return factor * values[0] / values[1]
    if operation == "share_as_percentage" and len(values) == 2 and values[1]:
        return factor * 100.0 * values[0] / values[1]
    if operation == "f1_from_counts" and len(values) == 3 and (values[1] + values[2]):
        return factor * 2.0 * values[0] / (values[1] + values[2])
    if operation == "difference" and len(values) == 2:
        return factor * (values[0] - values[1])
    if operation == "product":
        product = factor
        for value in values:
            product *= value
        return product
    return None


def check_report(qmd: Path, config: dict, receipts: dict[str, dict[str, object]],
                 raw: dict[str, object] | None = None) -> dict:
    """Run every check over one report.

    Args:
        qmd: The build source path.
        config: The checker configuration.
        receipts: The flattened receipts.
        raw: The raw receipts, keyed by identifier. A row derivation reads them.

    Returns:
        The result mapping with the counts and the findings.
    """
    raw = raw or {}
    abstract, body = read_source(qmd)
    templates = config.get("record_templates", [])
    blocks = split_blocks(abstract, body, abstract_line(qmd))
    known = set(receipts)
    table_records = table_sources(blocks, templates, known)
    exemptions = config.get("exempt_tokens", [])
    literal = [(str(entry["value"]), re.compile(entry["context"]), entry["reason"])
               for entry in config.get("exempt_claims", [])]
    derivations: dict[str, list[dict]] = {}
    for entry in config.get("derivations", []):
        derivations.setdefault(str(entry["value"]), []).append(entry)

    findings: dict[str, list[dict]] = {"no_receipt": [], "mismatch": [], "inconsistent": [],
                                       "prohibited_phrase": [], "accounting": []}
    declared: list[dict] = []
    open_defects = [(str(entry["value"]), re.compile(entry["context"]), entry)
                    for entry in config.get("open_defects", [])]
    counts = {"sections_with_a_claim": 0, "claim_values": 0, "verified": 0, "derived": 0,
              "exempt_tokens": 0, "exempt_claims": 0}
    field_values: dict[str, dict[str, list[int]]] = {}
    section_records: dict[str, set[str]] = {}
    section_claims: dict[str, list[dict]] = {}

    for block in blocks:
        section = block["section"]
        section_records.setdefault(section, set())
        section_records[section] |= block_records(block["text"], table_records, templates,
                                                  known)
        if block["kind"] == "heading":
            continue
        text = " ".join(block["text"].split())
        masked = mask_exempt_spans(text)
        for token, decimals, offset in numbers_of(masked):
            if is_exempt(token, text, offset, exemptions):
                counts["exempt_tokens"] += 1
                continue
            if any(token.replace(",", "") == value and context.search(text)
                   for value, context, _ in literal):
                counts["exempt_claims"] += 1
                continue
            entry = next((item for value, context, item in open_defects
                          if token.replace(",", "") == value and context.search(text)), None)
            if entry is not None:
                counts["open_defects"] = counts.get("open_defects", 0) + 1
                declared.append({"line": block["line"], "report_value": entry["value"],
                                 "record_value": entry["record_value"],
                                 "audit_item": entry["audit_item"], "reason": entry["reason"],
                                 "excerpt": excerpt(text)})
                continue
            section_claims.setdefault(section, []).append({
                "line": block["line"], "kind": block["kind"], "value": token.replace(",", ""),
                "decimals": decimals, "text": text})

    every_receipt = {record for records in section_records.values() for record in records
                     if not record.startswith("record:")}
    recap = [re.compile(pattern) for pattern in config.get("recap_sections", [])]
    bindings = [(re.compile(entry["pattern"]), set(entry["records"]))
                for entry in config.get("section_records", [])]
    conversions = config.get("unit_conversions", [])
    strict_scale = bool(config.get("scale_word_required", False))
    automatic_objects = bool(config.get("automatic_object_derivations", True))
    for section, claims in section_claims.items():
        counts["sections_with_a_claim"] += 1
        named = {record for record in section_records.get(section, set())
                 if not record.startswith("record:")}
        for pattern, records in bindings:
            if pattern.search(section):
                named |= records
                counts["bound_claims"] = counts.get("bound_claims", 0) + len(claims)
        if any(pattern.search(section) for pattern in recap):
            named |= every_receipt
            counts["recap_claims"] = counts.get("recap_claims", 0) + len(claims)
        counts["claim_values"] += len(claims)
        if not named:
            uncited = [claim for claim in claims if "@pmid" not in claim["text"]]
            counts["attributed_to_a_citation"] = (
                counts.get("attributed_to_a_citation", 0) + len(claims) - len(uncited))
            if uncited:
                findings["no_receipt"].append({
                    "section": section, "line": uncited[0]["line"],
                    "claim_values": len(uncited),
                    "values": sorted({claim["value"] for claim in uncited})[:12],
                    "excerpt": excerpt(uncited[0]["text"])})
            continue
        for claim in claims:
            value = claim["value"]
            decimals = claim["decimals"]
            matched = match_leaf(value, decimals, named, receipts,
                                 not strict_scale or writes_a_scale_word(claim["text"]))
            if matched:
                counts["verified"] += 1
                if matched[0]:
                    key = f"{matched[0]}#{matched[1]}"
                    field_values.setdefault(key, {}).setdefault(value, []).append(claim["line"])
                else:
                    counts["ambiguous_field"] = counts.get("ambiguous_field", 0) + 1
                continue
            inside = (derive_within_objects(value, decimals, named, receipts)
                      if automatic_objects else None)
            if inside:
                counts["derived"] += 1
                continue
            converted = resolve_unit_conversion(value, decimals, named, receipts, conversions)
            if converted:
                counts["converted"] = counts.get("converted", 0) + 1
                counts["derived"] += 1
                continue
            candidates = derivations.get(value, [])
            if candidates:
                resolved = False
                last: tuple[dict, float | None] | None = None
                for derived in candidates:
                    computed = resolve_derivation(derived, receipts, raw)
                    last = (derived, computed)
                    if (computed is not None
                            and abs(computed - float(value)) < claim_tolerance(value, decimals)):
                        counts["derived"] += 1
                        resolved = True
                        break
                if resolved:
                    continue
                derived, computed = last
                findings["mismatch"].append({
                    "line": claim["line"], "value": value, "kind": "derivation",
                    "report_value": value,
                    "receipt_value": "absent field" if computed is None
                                     else f"{computed:.{max(decimals, 4)}f}",
                    "detail": describe_derivation(derived),
                    "excerpt": excerpt(claim["text"])})
                continue
            findings["mismatch"].append({
                "line": claim["line"], "value": value, "kind": "no_matching_field",
                "section": section, "excerpt": excerpt(claim["text"]),
                "report_value": value,
                "receipt_value": nearest_value(value, named, receipts),
                "nearest": nearest_leaves(value, decimals, named, receipts)})

    for key, values in sorted(field_values.items()):
        distinct = distinct_numbers(values)
        if len(distinct) > 1:
            findings["inconsistent"].append({
                "quantity": key, "values": sorted(values),
                "sites": [{"value": value, "lines": sorted(set(lines))}
                          for value, lines in sorted(values.items())],
                "source": "receipt field"})

    findings["inconsistent"].extend(declared_quantity_findings(config, qmd))
    findings["prohibited_phrase"].extend(prohibited_findings(config, qmd))
    findings["accounting"] = accounting_findings(counts)
    counts["accounting"] = accounting(counts)
    total = sum(len(items) for items in findings.values())
    return {"counts": counts, "findings": findings, "total_findings": total,
            "open_defects": declared}


def accounting(counts: dict[str, int]) -> dict[str, object]:
    """Return the disjoint accounting of the claim values that this run read.

    Every claim value of a report reaches exactly one of three dispositions, which are a
    receipt field, a declared derivation and a citation in a section that names no record of
    this work. A token that an exemption rule removes and a claim that an exemption entry
    removes never become a claim value, so neither one is a component of the total. Both
    replicate reviews of 2026-09-19 read the four printed counts as four components of one
    total and found that they sum to 1,947 against a total of 1,935. The mapping below states
    which counts are components and which counts are outside the total.

    Args:
        counts: The count mapping of one run.

    Returns:
        A mapping with the three components, their sum, the total, the equality verdict and
        the two counts that stand outside the total.
    """
    components = {
        "verified_against_a_receipt_field": counts.get("verified", 0),
        "verified_against_a_declared_derivation": counts.get("derived", 0),
        "attributed_to_a_citation": counts.get("attributed_to_a_citation", 0),
    }
    component_sum = sum(components.values())
    total = counts.get("claim_values", 0)
    return {
        "components": components,
        "component_sum": component_sum,
        "claim_values": total,
        "components_sum_to_the_total": component_sum == total,
        "outside_the_total": {
            "exempt_tokens": counts.get("exempt_tokens", 0),
            "exempt_claims_with_a_written_reason": counts.get("exempt_claims", 0),
            "declared_open_defects": counts.get("open_defects", 0),
        },
    }


def accounting_findings(counts: dict[str, int]) -> list[dict]:
    """Return a finding when the three dispositions do not sum to the claim value total.

    Args:
        counts: The count mapping of one run.

    Returns:
        A list that holds one finding when the sum differs from the total, and no entry
        otherwise.
    """
    record = accounting(counts)
    if record["components_sum_to_the_total"]:
        return []
    return [{
        "components": record["components"],
        "component_sum": record["component_sum"],
        "claim_values": record["claim_values"],
        "reason": "The three dispositions of a claim value do not sum to the claim value total."}]


def distinct_numbers(values: dict[str, list[int]]) -> list[float]:
    """Give the distinct numeric values of one receipt field as the report prints them.

    Two printed forms of one number, for example 50 and 50.0, are one value. The comparison uses the
    tolerance of the coarser of the two printed precisions.

    Args:
        values: A map of printed value to the lines that print it.

    Returns:
        The distinct numeric values.
    """
    distinct: list[tuple[float, int]] = []
    for text in sorted(values):
        try:
            number = float(text)
        except ValueError:
            continue
        decimals = len(text.split(".")[1]) if "." in text else 0
        if not any(abs(number - seen) < max(tolerance(min(decimals, 6)), tolerance(seen_places))
                   for seen, seen_places in distinct):
            distinct.append((number, min(decimals, 6)))
    return [number for number, _ in distinct]


def claim_tolerance(value: str, decimals: int) -> float:
    """Give the comparison tolerance for one claim value as the report prints it.

    A value printed in scientific notation carries its precision in the mantissa, so the tolerance
    of 8.8e-06 is half of the last mantissa place scaled by the exponent and not half of a tenth.

    Args:
        value: The claim value as the report prints it.
        decimals: The number of decimal places of the mantissa.

    Returns:
        Half of the last printed place of the value.
    """
    base = tolerance(decimals)
    lowered = value.lower()
    if "e" not in lowered:
        return base
    try:
        exponent = int(lowered.split("e", 1)[1])
    except ValueError:
        return base
    return base * (10.0 ** exponent)


def tolerance(decimals: int) -> float:
    """Give the comparison tolerance for a value printed at one precision.

    Args:
        decimals: The number of decimal places that the report prints.

    Returns:
        Half of the last printed place.
    """
    return 0.5 * 10.0 ** (-decimals) + 1e-12


def full_token(text: str, offset: int, token: str) -> str:
    """Give the token with its thousands separators removed.

    Args:
        text: The block text.
        offset: The character offset of the token.
        token: The token string.

    Returns:
        The token with no comma.
    """
    del text, offset
    return token.replace(",", "")


def excerpt(text: str, width: int = 150) -> str:
    """Give a one-line excerpt of a block.

    Args:
        text: The block text.
        width: The character budget.

    Returns:
        The first characters of the block on one line.
    """
    flat = " ".join(text.split())
    return flat[:width] + ("..." if len(flat) > width else "")


def cited_records(text: str, templates: list[dict], known: set[str]) -> set[str]:
    """Give every record identifier that one piece of text cites by path or by template.

    Args:
        text: The block text or the caption text.
        templates: The record templates of the configuration. A template expands one caption
            placeholder, for example ``<slug>``, into the records that the placeholder stands for.
        known: The identifiers of the records that the checker loaded.

    Returns:
        The set of record identifiers. A receipt keeps its stem, a declared build record under
        ``reports/`` keeps its path with the extension removed, and a record that the checker did
        not load carries the prefix ``record:`` because the value check cannot read it.
    """
    named = set(RECEIPT_REFERENCE.findall(text))
    named |= set(TABLE_DATA_REFERENCE.findall(text))
    for path in BUILD_RECORD_REFERENCE.findall(text):
        key = path.rsplit(".", 1)[0]
        named.add(key if key in known else f"record:{path}")
    for record in RECORD_REFERENCE.findall(text):
        named.add(f"record:{record}")
    for template in templates:
        if re.search(template["pattern"], text):
            named |= set(template["records"])
    return named


def table_sources(blocks: list[dict], templates: list[dict],
                  known: set[str]) -> dict[str, set[str]]:
    """Map every table label to the records that its caption names.

    Args:
        blocks: The blocks of the report.
        templates: The record templates of the configuration.
        known: The identifiers of the records that the checker loaded.

    Returns:
        A map of table label to the set of record identifiers of its caption.
    """
    sources: dict[str, set[str]] = {}
    for block in blocks:
        label = TABLE_LABEL.search(block["text"])
        if not label:
            continue
        sources[f"{label.group(1)}-{label.group(2)}"] = cited_records(
            block["text"], templates, known)
    return sources


def block_records(text: str, table_records: dict[str, set[str]],
                  templates: list[dict], known: set[str]) -> set[str]:
    """Give every record that one block names, directly or through a table reference.

    Args:
        text: The block text.
        table_records: The map of table label to its records.
        templates: The record templates of the configuration.
        known: The identifiers of the records that the checker loaded.

    Returns:
        The set of record identifiers that the block can be checked against. A block that names a
        benchmark record under ``docs/benchmarks/`` and no receipt returns the record name, which
        the coverage check accepts and the value check cannot use.
    """
    named = cited_records(text, templates, known)
    for match in re.finditer(r"@((?:tbl|fig)-[A-Za-z0-9-]+)", text):
        named |= table_records.get(match.group(1), set())
    return named


def path_numbers(leaves: dict[str, object]) -> set[str]:
    """Give every number that appears as a key of one flattened receipt.

    A receipt stores a histogram as an object whose keys are the values that the run observed, for
    example the score of a spike. Such a key is a recorded value of the run, so the checker accepts
    it.

    Args:
        leaves: The flattened receipt.

    Returns:
        The set of numeric path segments.
    """
    found: set[str] = set()
    for path in leaves:
        for segment in path.replace("[", ".").replace("]", "").split("."):
            if segment and NUMBER.fullmatch(segment):
                found.add(segment.replace(",", ""))
    return found


def resolve_unit_conversion(value: str, decimals: int, named: set[str],
                            receipts: dict[str, dict[str, object]],
                            conversions: list[dict]) -> str | None:
    """Verify one claim value against a receipt field in another unit.

    A receipt records a duration in seconds and the report prints minutes, so the printed value is
    the recorded value divided by a declared constant. Each conversion names the leaf paths it
    applies to, so the rule stays narrow.

    Args:
        value: The claim value as the report prints it.
        decimals: The number of decimal places that the report prints.
        named: The record identifiers that the section names.
        receipts: The flattened records.
        conversions: The ``unit_conversions`` entries of the configuration.

    Returns:
        The record, the field and the conversion that produced the value, or None.
    """
    if not conversions:
        return None
    try:
        target = float(value)
    except ValueError:
        return None
    for entry in conversions:
        pattern = re.compile(entry["pattern"])
        divisor = float(entry["divide_by"])
        for record in sorted(named):
            for path, leaf in (receipts.get(record) or {}).items():
                if isinstance(leaf, bool) or not isinstance(leaf, (int, float)):
                    continue
                if not pattern.search(path):
                    continue
                if abs(float(leaf) / divisor - target) < tolerance(decimals):
                    return f"{record}#{path}={leaf} divided by {entry['divide_by']}"
    return None


def rounds_to(leaf: float, value: str, decimals: int) -> bool:
    """Say whether one numeric leaf prints as one claim value at the printed precision.

    The string forms of ``claim_forms`` miss a leaf that sits exactly half way, for example a
    recovery rate of 0.2435 that the report prints as 0.244, because the rounding of a binary float
    at that point goes the other way. The numeric comparison closes that gap.

    Args:
        leaf: The leaf value.
        value: The claim value as the report prints it.
        decimals: The number of decimal places that the report prints.

    Returns:
        True when the leaf lies inside half of the last printed place of the claim value.
    """
    try:
        target = float(value)
    except ValueError:
        return False
    if math.isnan(leaf) or math.isinf(leaf):
        return False
    return abs(leaf - target) <= claim_tolerance(value, decimals)


def match_leaf(value: str, decimals: int, named: set[str],
               receipts: dict[str, dict[str, object]],
               scale_words: bool = True) -> tuple[str, str] | None:
    """Find a receipt leaf that prints as one claim value.

    A numeric leaf is preferred over a string leaf that holds the value inside a sentence. When
    exactly one numeric leaf of the named receipts matches, the return value is that leaf and the
    consistency check of the caller can key on it. When several match, the value is verified and the
    caller receives the ambiguous marker so that it keys on nothing.

    Args:
        value: The claim value as the report prints it, with no thousands separator.
        decimals: The number of decimal places that the report prints.
        named: The receipt identifiers that the section names.
        receipts: The flattened receipts.
        scale_words: Accept a leaf printed in thousands, millions or billions.

    Returns:
        The receipt identifier and the JSON path of the matching leaf, the pair ``("", "")`` when
        the match is ambiguous or comes from a string leaf, or None when nothing matches.
    """
    numeric: list[tuple[str, str]] = []
    textual: list[tuple[str, str]] = []
    for receipt in sorted(named):
        leaves = receipts.get(receipt)
        if leaves is None:
            continue
        for path, leaf in leaves.items():
            if isinstance(leaf, (int, float)) and not isinstance(leaf, bool):
                if (value in claim_forms(leaf, decimals, scale_words)
                        or rounds_to(float(leaf), value, decimals)):
                    numeric.append((receipt, path))
                continue
            if value in claim_forms(leaf, decimals, scale_words):
                textual.append((receipt, path))
    if len(numeric) == 1:
        return numeric[0]
    if numeric or textual:
        return "", ""
    for receipt in sorted(named):
        leaves = receipts.get(receipt)
        if leaves is not None and value in path_numbers(leaves):
            return "", ""
    return None


def describe_derivation(spec: dict) -> str:
    """Give a one-line description of a declared derivation.

    Args:
        spec: The derivation entry.

    Returns:
        The operation and the receipt fields or rows that it reads.
    """
    if "fields" in spec:
        return f"{spec['operation']} of " + ", ".join(spec["fields"])
    rows = [spec.get("rows", ""), spec.get("null_rows", "")]
    return (f"{spec['operation']} at {spec.get('field', 'best_rank')} "
            f"at or below {spec.get('max_rank')} over "
            + ", ".join(item for item in rows if item))


def nearest_value(value: str, named: set[str],
                  receipts: dict[str, dict[str, object]]) -> str:
    """Give the value of the receipt leaf closest to one unmatched claim value.

    Args:
        value: The claim value as the report prints it.
        named: The record identifiers that the section names.
        receipts: The flattened receipts.

    Returns:
        The closest recorded value with its receipt and its field, or a message when the named
        records carry no numeric leaf.
    """
    nearest = nearest_leaves(value, 0, named, receipts, limit=1)
    return nearest[0] if nearest else "the named records carry no numeric field"


def nearest_leaves(value: str, decimals: int, named: set[str],
                   receipts: dict[str, dict[str, object]], limit: int = 3) -> list[str]:
    """Give the receipt leaves closest to one unmatched claim value.

    Args:
        value: The claim value.
        decimals: The number of decimal places that the report prints.
        named: The receipt identifiers that the block names.
        receipts: The flattened receipts.
        limit: The number of candidates to return.

    Returns:
        A list of ``receipt#path=value`` strings, closest first.
    """
    try:
        target = float(value)
    except ValueError:
        return []
    del decimals
    candidates: list[tuple[float, str]] = []
    for receipt in sorted(named):
        for path, leaf in (receipts.get(receipt) or {}).items():
            if isinstance(leaf, bool) or not isinstance(leaf, (int, float)):
                continue
            candidates.append((abs(float(leaf) - target), f"{receipt}#{path}={leaf}"))
    candidates.sort(key=lambda item: item[0])
    return [text for _, text in candidates[:limit]]


def declared_quantity_findings(config: dict, qmd: Path) -> list[dict]:
    """Check every declared quantity for one value across the report and the records.

    Args:
        config: The checker configuration.
        qmd: The build source path.

    Returns:
        A list of findings, one per quantity that carries two values.
    """
    findings: list[dict] = []
    for entry in config.get("quantities", []):
        pattern = re.compile(entry["pattern"])
        seen: dict[str, list[str]] = {}
        for relative in [str(qmd)] + entry.get("also_in", []):
            path = Path(relative) if Path(relative).is_absolute() else ROOT / relative
            if not path.exists():
                findings.append({"quantity": entry["name"], "values": [],
                                 "sites": [{"line": 0, "value": f"missing file {relative}"}],
                                 "source": "declared quantity"})
                continue
            for index, line in enumerate(path.read_text().splitlines(), start=1):
                for match in pattern.finditer(line):
                    seen.setdefault(normalize(match.group("value")), []).append(
                        f"{relative}:{index}")
        if not seen:
            findings.append({"quantity": entry["name"], "values": [],
                             "sites": [{"line": 0, "value": "no occurrence found"}],
                             "source": "declared quantity"})
            continue
        expected = normalize(str(entry["value"]))
        wrong = {value: sites for value, sites in seen.items() if value != expected}
        if wrong:
            findings.append({
                "quantity": entry["name"], "values": sorted(seen),
                "sites": [{"line": 0, "value": f"{value} at {', '.join(sites[:4])}"}
                          for value, sites in sorted(wrong.items())],
                "source": f"declared quantity, expected {expected}"})
    return findings


def normalize(value: str) -> str:
    """Normalize a captured quantity for comparison.

    Args:
        value: The captured text.

    Returns:
        The lower case text with no thousands separator and no surrounding space.
    """
    words = {"one": "1", "two": "2", "three": "3", "four": "4", "five": "5", "six": "6",
             "seven": "7", "eight": "8", "nine": "9", "ten": "10"}
    text = value.strip().lower().replace(",", "")
    return words.get(text, text)


def prohibited_findings(config: dict, qmd: Path) -> list[dict]:
    """Find every prohibited phrase in the report source.

    Args:
        config: The checker configuration.
        qmd: The build source path.

    Returns:
        A list of findings, one per occurrence.
    """
    findings: list[dict] = []
    text = qmd.read_text().splitlines()
    for entry in config.get("prohibited_phrases", []):
        pattern = re.compile(entry["pattern"], re.IGNORECASE)
        for index, line in enumerate(text, start=1):
            for match in pattern.finditer(line):
                findings.append({"line": index, "phrase": match.group(0),
                                 "rule": entry["reason"],
                                 "excerpt": excerpt(line)})
    return findings


def print_report(result: dict, quiet: bool) -> None:
    """Print the checker result.

    Args:
        result: The result mapping.
        quiet: Print the counts alone when no finding is open.
    """
    counts = result["counts"]
    findings = result["findings"]
    if not quiet or result["total_findings"]:
        for name, items in findings.items():
            for item in items:
                print(f"  {name.upper():<18} {json.dumps(item, sort_keys=True)}")
    print(f"claims checked {counts['claim_values']}, verified against a receipt field "
          f"{counts['verified']}, verified against a declared derivation {counts['derived']}, "
          f"attributed to a citation in a section that names no receipt "
          f"{counts.get('attributed_to_a_citation', 0)}, "
          f"sections with a claim {counts['sections_with_a_claim']}, "
          f"exempt tokens {counts['exempt_tokens']}, "
          f"exempt claims with a written reason {counts['exempt_claims']}")
    record = counts.get("accounting", {})
    if record:
        parts = record["components"]
        print(f"accounting: {parts['verified_against_a_receipt_field']} against a receipt field "
              f"+ {parts['verified_against_a_declared_derivation']} against a declared derivation "
              f"+ {parts['attributed_to_a_citation']} attributed to a citation "
              f"= {record['component_sum']} of {record['claim_values']} claim values; "
              f"{record['outside_the_total']['exempt_claims_with_a_written_reason']} exempt claims "
              f"and {record['outside_the_total']['exempt_tokens']} exempt tokens stand outside "
              f"that total")
    print(f"findings: no receipt {len(findings['no_receipt'])}, "
          f"value differs from the receipt {len(findings['mismatch'])}, "
          f"one quantity with two values {len(findings['inconsistent'])}, "
          f"prohibited phrase {len(findings['prohibited_phrase'])}, "
          f"accounting {len(findings.get('accounting', []))}")
    for item in result.get("open_defects", []):
        print(f"  OPEN DEFECT        line {item['line']}, the report prints "
              f"{item['report_value']} and the record holds {item['record_value']}, "
              f"item {item['audit_item']} owns the repair")


def check_one(qmd: Path, config_path: Path, receipt_dir: Path, quiet: bool,
              records_root: Path | None = None) -> dict | None:
    """Check one report and print its result.

    Args:
        qmd: The build source path.
        config_path: The checker configuration path.
        receipt_dir: The receipt directory.
        quiet: Print the counts alone when no finding is open.
        records_root: The directory that a declared record path resolves against.

    Returns:
        The result mapping, or None when the source or the configuration is absent.
    """
    if not qmd.exists():
        print(f"No report source at {qmd}.", file=sys.stderr)
        return None
    if not config_path.exists():
        print(f"No checker configuration at {config_path}.", file=sys.stderr)
        return None
    config = load_config(config_path)
    raw: dict[str, object] = {}
    receipts = load_receipts(receipt_dir, raw)
    receipts.update(load_declared_records(config, records_root))
    receipts = filter_leaf_paths(receipts, config.get("ignored_leaf_paths", []))
    result = check_report(qmd, config, receipts, raw)
    result["report"] = str(qmd)
    result["receipts_read"] = len(receipts)
    print(f"report {qmd.relative_to(ROOT) if qmd.is_relative_to(ROOT) else qmd}")
    print_report(result, quiet)
    return result


def main() -> int:
    """Check one report or several and report the findings.

    Returns:
        Zero when every report carries no finding, one when a report carries one, and two when a
        source or a configuration is absent.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("qmd", type=Path, nargs="+", help="the build source of each report")
    parser.add_argument("--config", type=Path, default=None,
                        help="the checker configuration, default <qmd dir>/claims.yaml")
    parser.add_argument("--receipts", type=Path, default=ROOT / "tracking" / "receipts")
    parser.add_argument("--records-root", type=Path, default=None,
                        help="the directory that a declared record path resolves against")
    parser.add_argument("--json", type=Path, default=None, help="write the result to this file")
    parser.add_argument("--quiet", action="store_true", help="print the counts alone on a pass")
    args = parser.parse_args()
    if args.config and len(args.qmd) > 1:
        print("--config names one configuration and does not apply to several reports.",
              file=sys.stderr)
        return 2
    results: list[dict] = []
    for entry in args.qmd:
        qmd = entry if entry.is_absolute() else ROOT / entry
        config_path = args.config or qmd.parent / "claims.yaml"
        result = check_one(qmd, config_path, args.receipts, args.quiet,
                           args.records_root)
        if result is None:
            return 2
        results.append(result)
    payload = results[0] if len(results) == 1 else {
        "reports": results,
        "total_findings": sum(item["total_findings"] for item in results)}
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        with args.json.open("w") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
    return 1 if payload["total_findings"] else 0


if __name__ == "__main__":
    sys.exit(main())
