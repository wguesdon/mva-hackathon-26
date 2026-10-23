#!/usr/bin/env python3
"""Measure the agreement between two readings of the Track 2 evidence records. Data class A.

The Track 2 edge provenance reading had one reader. A single reader classification is an unmeasured
instrument, and both the submission audit and the second jury round raised it. This script measures
the instrument. It computes raw agreement, Cohen's kappa with a large sample confidence interval and
the confusion matrix across the classification categories, and it names every item where the two
readings disagree.

The script holds three groups of functions.

* ``cohen_kappa`` and ``kappa_interval`` give the point estimate and the interval. The variance is
  the large sample variance of Fleiss, Cohen and Everitt (1969), which is the variance of the
  estimator around the true value and not the variance under the null of independence.
* ``kappa_power`` gives the power of a two sided test of one kappa value against another for a
  binary judgement, from the variance of Bloch and Kraemer (1989) that Donner and Eliasziw (1992)
  use for sample size. The design document quotes it before the reading runs.
* ``confusion`` and ``disagreements`` give the table and the item list.

The script reads only public accessions, PubMed identifiers and published sentences. No patient data
enters it.

Usage:
    uv run python scripts/track2_reading_agreement.py power
    uv run python scripts/track2_reading_agreement.py compare
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROVENANCE = ROOT / "reports" / "track2" / "provenance"
RECEIPTS = ROOT / "tracking" / "receipts"

FIRST_EDGES = ROOT / "tracking" / "track2_edge_provenance_2026_09_15.tsv"
FIRST_SENTENCES = PROVENANCE / "seed_edge_surfaces_2026_09_15.tsv"
SECOND_EDGES = PROVENANCE / "second_reading_edges_2026_09_18.tsv"
SECOND_EDGES_RUN1 = PROVENANCE / "second_reading_edges_run1_2026_09_18.tsv"
SECOND_SENTENCES = PROVENANCE / "second_reading_sentences_2026_09_18.tsv"
ADJUDICATION = PROVENANCE / "reading_adjudication_2026_09_18.tsv"
RECORD_CHECK = PROVENANCE / "reading_record_check_2026_09_18.tsv"
RECEIPT_PATH = RECEIPTS / "t2_reading_agreement_v1.json"

EDGE_CLASSES = ("SUPPORTED", "WRONG DIRECTION", "NOT AN INDICATION", "OFF TOPIC", "UNCLEAR")
SENTENCE_CLASSES = ("yes", "no")

# This configured subset of public edge identifiers defines one agreement summary.
# The identifiers carry no completed reading verdict.
CANDIDATE_EDGES = ("A1", "A2", "A3", "A4", "A8", "A9")

# The design document fixes these values before the reading runs. The power target of the binary
# judgement is a kappa of 0.80 against a null of 0.60, which is the floor that Landis and Koch call
# substantial agreement.
NULL_KAPPA = 0.60
ALTERNATIVE_KAPPA = 0.80
ALPHA = 0.05


def normal_cdf(value: float) -> float:
    """Give the standard normal distribution function at one point.

    Args:
        value: The point on the standard normal scale.

    Returns:
        The probability that a standard normal variate is at most the point.
    """
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


def normal_quantile(probability: float) -> float:
    """Give the standard normal quantile by bisection.

    Args:
        probability: A probability strictly between zero and one.

    Returns:
        The point whose standard normal distribution function equals the probability.

    Raises:
        ValueError: If the probability is not strictly between zero and one.
    """
    if not 0.0 < probability < 1.0:
        raise ValueError("The probability must be strictly between zero and one.")
    low, high = -40.0, 40.0
    for _ in range(200):
        middle = 0.5 * (low + high)
        if normal_cdf(middle) < probability:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def contingency(pairs: list[tuple[str, str]], classes: tuple[str, ...]) -> list[list[int]]:
    """Build the square table of the two readings.

    Args:
        pairs: One tuple of the first verdict and the second verdict for each item.
        classes: The class labels in the order of the rows and the columns.

    Returns:
        The counts, with the first reading on the rows and the second reading on the columns.
    """
    index = {name: position for position, name in enumerate(classes)}
    table = [[0 for _ in classes] for _ in classes]
    for first, second in pairs:
        table[index[first]][index[second]] += 1
    return table


def cohen_kappa(table: list[list[int]]) -> dict:
    """Compute Cohen's kappa and the two agreement proportions of one table.

    Args:
        table: The square contingency table of the two readings.

    Returns:
        The item count, the observed agreement, the expected agreement and kappa. Kappa is None when
        the expected agreement is one, which happens when both readings use one class only.
    """
    total = sum(sum(row) for row in table)
    if total == 0:
        return {"n": 0, "observed_agreement": None, "expected_agreement": None, "kappa": None}
    size = len(table)
    row_sums = [sum(row) for row in table]
    column_sums = [sum(table[row][column] for row in range(size)) for column in range(size)]
    observed = sum(table[position][position] for position in range(size)) / total
    expected = sum(row_sums[position] * column_sums[position] for position in range(size)) / (total * total)
    kappa = None if abs(1.0 - expected) < 1e-12 else (observed - expected) / (1.0 - expected)
    return {"n": total, "observed_agreement": observed, "expected_agreement": expected, "kappa": kappa}


def kappa_standard_error(table: list[list[int]]) -> float | None:
    """Give the large sample standard error of Cohen's kappa.

    The formula is the one of Fleiss, Cohen and Everitt (1969) for the variance of the estimator
    around the true value. It is the formula that a confidence interval needs. The variance under
    the null of independence is a different quantity and this function does not compute it.

    Args:
        table: The square contingency table of the two readings.

    Returns:
        The standard error, or None when the expected agreement is one or the table is empty.
    """
    result = cohen_kappa(table)
    total, kappa = result["n"], result["kappa"]
    if total == 0 or kappa is None:
        return None
    size = len(table)
    row_sums = [sum(row) / total for row in table]
    column_sums = [sum(table[row][column] for row in range(size)) / total for column in range(size)]
    observed = result["observed_agreement"]
    expected = result["expected_agreement"]
    if abs(1.0 - expected) < 1e-12:
        return None
    term_a = sum(
        table[position][position] / total * (1.0 - (row_sums[position] + column_sums[position]) * (1.0 - kappa)) ** 2
        for position in range(size)
    )
    term_b = 0.0
    for row in range(size):
        for column in range(size):
            if row == column:
                continue
            term_b += table[row][column] / total * (column_sums[row] + row_sums[column]) ** 2
    term_b *= (1.0 - kappa) ** 2
    term_c = (kappa - expected * (1.0 - kappa)) ** 2
    variance = (term_a + term_b - term_c) / (total * (1.0 - expected) ** 2)
    return math.sqrt(variance) if variance > 0 else 0.0


def kappa_interval(table: list[list[int]], alpha: float = ALPHA) -> dict:
    """Give Cohen's kappa with a two sided confidence interval.

    Args:
        table: The square contingency table of the two readings.
        alpha: The two sided error rate of the interval.

    Returns:
        The kappa result, the standard error and the two interval limits. A limit is clipped to the
        range of kappa, which is minus one to one.
    """
    result = cohen_kappa(table)
    error = kappa_standard_error(table)
    if result["kappa"] is None or error is None:
        return {**result, "standard_error": error, "ci_low": None, "ci_high": None, "alpha": alpha}
    half = normal_quantile(1.0 - alpha / 2.0) * error
    return {**result, "standard_error": error,
            "ci_low": max(-1.0, result["kappa"] - half), "ci_high": min(1.0, result["kappa"] + half),
            "alpha": alpha}


def binary_kappa_variance(kappa: float, prevalence: float, n: int) -> float:
    """Give the large sample variance of kappa for a binary judgement at a known prevalence.

    The expression is the one of Bloch and Kraemer (1989) that Donner and Eliasziw (1992) use for
    the sample size of a kappa study. It is a function of the true kappa, the marginal prevalence of
    the positive class and the item count, so it is available before any reading runs.

    Args:
        kappa: The true value of kappa.
        prevalence: The marginal proportion of the positive class.
        n: The item count.

    Returns:
        The variance of the estimator at that kappa.

    Raises:
        ValueError: If the prevalence is not strictly between zero and one, or n is not positive.
    """
    if not 0.0 < prevalence < 1.0:
        raise ValueError("The prevalence must be strictly between zero and one.")
    if n <= 0:
        raise ValueError("The item count must be positive.")
    spread = 2.0 * prevalence * (1.0 - prevalence)
    return ((1.0 - kappa) / n) * ((1.0 - kappa) * (1.0 - 2.0 * kappa) + kappa * (2.0 - kappa) / spread)


def kappa_power(n: int, prevalence: float, null_kappa: float = NULL_KAPPA,
                alternative_kappa: float = ALTERNATIVE_KAPPA, alpha: float = ALPHA) -> dict:
    """Give the power of a two sided test of one kappa value against another.

    Args:
        n: The item count of the reading.
        prevalence: The marginal proportion of the positive class.
        null_kappa: The kappa value of the null hypothesis.
        alternative_kappa: The kappa value of the alternative hypothesis.
        alpha: The two sided error rate.

    Returns:
        The two standard errors, the power and the half width of the interval at the alternative.
    """
    error_null = math.sqrt(binary_kappa_variance(null_kappa, prevalence, n))
    error_alternative = math.sqrt(binary_kappa_variance(alternative_kappa, prevalence, n))
    critical = normal_quantile(1.0 - alpha / 2.0)
    power = normal_cdf(((alternative_kappa - null_kappa) - critical * error_null) / error_alternative)
    return {"n": n, "prevalence": prevalence, "null_kappa": null_kappa,
            "alternative_kappa": alternative_kappa, "alpha": alpha,
            "standard_error_at_null": error_null, "standard_error_at_alternative": error_alternative,
            "power": power, "ci_half_width_at_alternative": critical * error_alternative}


def beta_quantile(probability: float, shape_a: float, shape_b: float) -> float:
    """Give the quantile of a beta distribution by bisection on its distribution function.

    Args:
        probability: The probability, strictly between zero and one.
        shape_a: The first shape parameter.
        shape_b: The second shape parameter.

    Returns:
        The quantile of the beta distribution.
    """
    low, high = 0.0, 1.0
    for _ in range(200):
        middle = 0.5 * (low + high)
        if regularized_beta(middle, shape_a, shape_b) < probability:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def regularized_beta(value: float, shape_a: float, shape_b: float, steps: int = 20000) -> float:
    """Give the regularized incomplete beta function by numerical integration.

    Args:
        value: The upper limit of the integral, between zero and one.
        shape_a: The first shape parameter.
        shape_b: The second shape parameter.
        steps: The number of trapezoid steps.

    Returns:
        The value of the regularized incomplete beta function.
    """
    if value <= 0.0:
        return 0.0
    if value >= 1.0:
        return 1.0
    log_norm = math.lgamma(shape_a + shape_b) - math.lgamma(shape_a) - math.lgamma(shape_b)
    width = value / steps
    total = 0.0
    for step in range(steps + 1):
        point = min(max(step * width, 1e-12), 1.0 - 1e-12)
        density = math.exp(log_norm + (shape_a - 1.0) * math.log(point)
                           + (shape_b - 1.0) * math.log(1.0 - point))
        total += density * (0.5 if step in (0, steps) else 1.0)
    return min(1.0, total * width)


def clopper_pearson(successes: int, trials: int, alpha: float = ALPHA) -> dict:
    """Give the exact two sided binomial interval of one proportion.

    The interval is needed where the agreement is perfect. The large sample standard error of kappa
    is zero at perfect agreement, so the normal interval on kappa is degenerate and states nothing.

    Args:
        successes: The number of agreements.
        trials: The number of items.
        alpha: The two sided error rate.

    Returns:
        The proportion and the two interval limits.
    """
    if trials == 0:
        return {"proportion": None, "ci_low": None, "ci_high": None, "alpha": alpha}
    low = 0.0 if successes == 0 else beta_quantile(alpha / 2.0, successes, trials - successes + 1)
    high = 1.0 if successes == trials else beta_quantile(1.0 - alpha / 2.0, successes + 1,
                                                         trials - successes)
    return {"proportion": successes / trials, "ci_low": low, "ci_high": high, "alpha": alpha}


def collapse(label: str, positive: str) -> str:
    """Map a class label onto a binary label for a collapsed comparison.

    Args:
        label: The class label of a reading.
        positive: The class label that becomes the positive class.

    Returns:
        The string "yes" when the label is the positive class, otherwise "no".
    """
    return "yes" if label == positive else "no"


def confusion(pairs: list[tuple[str, str]], classes: tuple[str, ...]) -> list[dict]:
    """Give the confusion matrix of two readings as a list of rows.

    Args:
        pairs: One tuple of the first verdict and the second verdict for each item.
        classes: The class labels in the order of the rows and the columns.

    Returns:
        One mapping per first reading class, with the count of each second reading class.
    """
    table = contingency(pairs, classes)
    rows = []
    for position, name in enumerate(classes):
        row = {"first_reading": name, "total": sum(table[position])}
        for column, other in enumerate(classes):
            row[other] = table[position][column]
        rows.append(row)
    return rows


def read_tsv(path: Path) -> list[dict]:
    """Read one tab separated file.

    Args:
        path: The file path.

    Returns:
        One mapping per data row.
    """
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(path: Path) -> str:
    """Give the first sixteen hexadecimal characters of the SHA-256 of a file.

    Args:
        path: The file path.

    Returns:
        The digest prefix, or the string "absent" when the file does not exist.
    """
    if not path.exists():
        return "absent"
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def sentence_key(row: dict) -> str:
    """Give the key of one seed sentence row.

    The key is the relation, the partner node and the index of the item inside that pair. It is the
    key that `scripts/track2_mention_audit.py` uses for the derived rule comparison.

    Args:
        row: A row of the seed surface table or of the second reading table.

    Returns:
        The key string.
    """
    return f"{row['relation']}|{row['partner']}|{row['item_index']}"


def load_sentence_pairs() -> tuple[list[tuple[str, str]], list[dict]]:
    """Pair the first and the second reading of the seed sentences.

    Returns:
        The verdict pairs and the rows where the two readings disagree.

    Raises:
        SystemExit: If an item of one reading has no partner in the other.
    """
    first = {sentence_key(row): row for row in read_tsv(FIRST_SENTENCES)}
    second = {sentence_key(row): row for row in read_tsv(SECOND_SENTENCES)}
    missing = set(first) ^ set(second)
    if missing:
        raise SystemExit(f"{len(missing)} seed sentence item(s) appear in one reading only.")
    pairs, differing = [], []
    for key in sorted(first):
        one, two = first[key]["names_syndrome"], second[key]["names_syndrome"]
        pairs.append((one, two))
        if one != two:
            differing.append({"key": key, "pmid": first[key]["pmid"],
                              "seed_surface": first[key]["seed_surface"],
                              "first_reading": one, "second_reading": two,
                              "second_reason": second[key].get("reason", ""),
                              "sentence": first[key]["sentence"]})
    return pairs, differing


def load_edge_pairs() -> tuple[list[tuple[str, str]], list[dict]]:
    """Pair the first and the second reading of the edges.

    Returns:
        The verdict pairs and the rows where the two readings disagree.

    Raises:
        SystemExit: If an edge of one reading has no partner in the other.
    """
    first = {row["edge_id"]: row for row in read_tsv(FIRST_EDGES)}
    second = {row["edge_id"]: row for row in read_tsv(SECOND_EDGES)}
    missing = set(first) ^ set(second)
    if missing:
        raise SystemExit(f"{len(missing)} edge(s) appear in one reading only.")
    pairs, differing = [], []
    for key in sorted(first):
        one, two = first[key]["classification"], second[key]["classification"]
        pairs.append((one, two))
        if one != two:
            differing.append({"edge_id": key, "edge": first[key]["edge"],
                              "relation": first[key]["relation"],
                              "evidence_kind": first[key]["evidence_kind"],
                              "sentences_read": first[key]["sentences_read"],
                              "first_reading": one, "second_reading": two,
                              "second_reason": second[key].get("reason", "")})
    return pairs, differing


def prevalence_of(pairs: list[tuple[str, str]], positive: str) -> float:
    """Give the marginal proportion of the positive class in the first reading.

    Args:
        pairs: One tuple of the first verdict and the second verdict for each item.
        positive: The class label that counts as positive.

    Returns:
        The proportion of items whose first reading verdict is the positive class.
    """
    return sum(1 for first, _ in pairs if first == positive) / len(pairs)


def block(pairs: list[tuple[str, str]], classes: tuple[str, ...]) -> dict:
    """Build the whole agreement block of one comparison.

    Args:
        pairs: One tuple of the first verdict and the second verdict for each item.
        classes: The class labels of the comparison.

    Returns:
        The kappa interval, the agreement counts and the confusion matrix.
    """
    table = contingency(pairs, classes)
    result = kappa_interval(table)
    agree = sum(1 for first, second in pairs if first == second)
    exact = clopper_pearson(agree, len(pairs))
    return {**result, "items_that_agree": agree, "items_that_differ": len(pairs) - agree,
            "agreement_exact_ci_low": exact["ci_low"], "agreement_exact_ci_high": exact["ci_high"],
            "kappa_interval_is_degenerate": result["standard_error"] == 0.0,
            "classes": list(classes), "confusion": confusion(pairs, classes),
            "first_reading_counts": dict(Counter(first for first, _ in pairs)),
            "second_reading_counts": dict(Counter(second for _, second in pairs))}


def power_report() -> dict:
    """Give the power of the design for both readings.

    Returns:
        One power block per reading, under the prevalence that the first reading recorded.
    """
    sentences = read_tsv(FIRST_SENTENCES)
    edges = read_tsv(FIRST_EDGES)
    sentence_prevalence = sum(1 for row in sentences if row["names_syndrome"] == "yes") / len(sentences)
    off_topic = sum(1 for row in edges if row["classification"] == "OFF TOPIC") / len(edges)
    return {
        "seed_sentence_reading": kappa_power(len(sentences), sentence_prevalence),
        "edge_reading_off_topic_collapse": kappa_power(len(edges), off_topic,
                                                       alternative_kappa=0.90),
        "note": ("Each reading is a census of its record, so the item count is the population size "
                 "and no larger study is available inside this record."),
    }


def run1_block() -> dict:
    """Compare the first run of the edge reading with the first reading and with run 2.

    Amendment 1 of the design records why run 1 is a sensitivity analysis and not the reported
    value. Run 1 named each entity by its accession alone and run 2 names it by the graph label.

    Returns:
        The two agreement blocks and the edges whose class the label changed, or a note when the
        run 1 file is absent.
    """
    if not SECOND_EDGES_RUN1.exists():
        return {"file_present": False}
    first = {row["edge_id"]: row["classification"] for row in read_tsv(FIRST_EDGES)}
    run1 = {row["edge_id"]: row["classification"] for row in read_tsv(SECOND_EDGES_RUN1)}
    run2 = {row["edge_id"]: row["classification"] for row in read_tsv(SECOND_EDGES)}
    keys = sorted(first)
    changed = [{"edge_id": key, "run1": run1[key], "run2": run2[key], "first_reading": first[key]}
               for key in keys if run1[key] != run2[key]]
    return {
        "file_present": True,
        "path": str(SECOND_EDGES_RUN1.relative_to(ROOT)),
        "sha256_prefix": digest(SECOND_EDGES_RUN1),
        "first_reading_against_run1": block([(first[key], run1[key]) for key in keys], EDGE_CLASSES),
        "run1_against_run2": block([(run1[key], run2[key]) for key in keys], EDGE_CLASSES),
        "edges_whose_class_the_entity_label_changed": changed,
        "count_of_edges_whose_class_the_entity_label_changed": len(changed),
    }


def seed_composition() -> dict:
    """Describe how the first reading of the seed sentences was produced.

    The `names_syndrome` column of the first reading is the output of the version 1 rule on most of
    its rows and a reader override on the rest. The composition changes how a perfect agreement is
    read, so the receipt records it.

    Returns:
        The sentence count, the count that the version 1 rule decided and the override count.
    """
    rows = read_tsv(FIRST_SENTENCES)
    overrides = sum(1 for row in rows if row["reader_override"].strip())
    return {"sentences": len(rows),
            "sentences_decided_by_the_version_1_rule": len(rows) - overrides,
            "sentences_with_a_reader_override": overrides,
            "override_reasons": dict(Counter(row["reader_override"] for row in rows
                                             if row["reader_override"].strip()))}


def adjudication_block() -> dict:
    """Summarize the adjudication of the disagreements.

    Returns:
        The counts of the adjudicated verdicts, or a note when the file is absent.
    """
    if not ADJUDICATION.exists():
        return {"file_present": False}
    rows = read_tsv(ADJUDICATION)
    picks = Counter(row["adjudicated_reading"] for row in rows)
    return {"file_present": True, "items": len(rows),
            "adjudicator_picked_the_first_reading": picks.get("first", 0),
            "adjudicator_picked_the_second_reading": picks.get("second", 0),
            "adjudicator_picked_neither": picks.get("neither", 0),
            "by_record": dict(Counter(row["record"] for row in rows)),
            "rows": rows}


def adjudicated_class_of_each_edge() -> dict[str, str]:
    """Give the adjudicated class of every edge.

    An edge where the two readings agree keeps that class. An edge where they differ takes the class
    of the model adjudicator, or the class of the record check when the adjudicator answered neither.

    Returns:
        One class per edge identifier.

    Raises:
        SystemExit: If an item is neither adjudicated nor resolved on the record.
    """
    first = {row["edge_id"]: row["classification"] for row in read_tsv(FIRST_EDGES)}
    second = {row["edge_id"]: row["classification"] for row in read_tsv(SECOND_EDGES)}
    adjudicated = {row["key"]: row for row in read_tsv(ADJUDICATION)
                   if row["record"] == "edges"} if ADJUDICATION.exists() else {}
    checked = {row["key"]: row for row in read_tsv(RECORD_CHECK)} if RECORD_CHECK.exists() else {}
    out = {}
    for edge_id, verdict in first.items():
        if verdict == second[edge_id]:
            out[edge_id] = verdict
            continue
        row = adjudicated.get(edge_id)
        if row is not None and row["adjudicated_class"]:
            out[edge_id] = row["adjudicated_class"]
            continue
        check = checked.get(edge_id)
        if check is None:
            raise SystemExit(f"{edge_id}: the item is neither adjudicated nor resolved.")
        out[edge_id] = check["verdict"]
    return out


def sensitivity_block() -> dict:
    """Give the counts of the first reading and of the adjudicated reading side by side.

    Returns:
        The class counts of both, over all edges, over the text-mined edges and over the curated
        edges, and the edges whose class the adjudication changed.
    """
    rows = read_tsv(FIRST_EDGES)
    first = {row["edge_id"]: row["classification"] for row in rows}
    kind = {row["edge_id"]: row["evidence_kind"] for row in rows}
    adjudicated = adjudicated_class_of_each_edge()
    mined = [edge_id for edge_id in first if kind[edge_id] == "pubtator3_sentence"]
    curated = [edge_id for edge_id in first if kind[edge_id] != "pubtator3_sentence"]
    changed = [{"edge_id": edge_id, "first_reading": first[edge_id],
                "adjudicated": adjudicated[edge_id]}
               for edge_id in sorted(first) if first[edge_id] != adjudicated[edge_id]]
    return {
        "edges": len(first),
        "text_mined_edges": len(mined),
        "curated_edges": len(curated),
        "first_reading_all": {name: sum(1 for value in first.values() if value == name)
                              for name in EDGE_CLASSES},
        "adjudicated_all": {name: sum(1 for value in adjudicated.values() if value == name)
                            for name in EDGE_CLASSES},
        "first_reading_text_mined": {name: sum(1 for edge_id in mined if first[edge_id] == name)
                                     for name in EDGE_CLASSES},
        "adjudicated_text_mined": {name: sum(1 for edge_id in mined if adjudicated[edge_id] == name)
                                   for name in EDGE_CLASSES},
        "first_reading_curated": {name: sum(1 for edge_id in curated if first[edge_id] == name)
                                  for name in EDGE_CLASSES},
        "adjudicated_curated": {name: sum(1 for edge_id in curated if adjudicated[edge_id] == name)
                                for name in EDGE_CLASSES},
        "edges_whose_class_the_adjudication_changed": changed,
        "count_of_edges_whose_class_the_adjudication_changed": len(changed),
        "candidate_edges": CANDIDATE_EDGES,
        "candidate_edges_count": len(CANDIDATE_EDGES),
        "candidate_edges_supported_in_the_first_reading":
            sum(1 for edge_id in CANDIDATE_EDGES if first[edge_id] == "SUPPORTED"),
        "candidate_edges_supported_after_adjudication":
            sum(1 for edge_id in CANDIDATE_EDGES if adjudicated[edge_id] == "SUPPORTED"),
        "candidate_edges_that_read_the_wrong_way_after_adjudication":
            sum(1 for edge_id in CANDIDATE_EDGES if adjudicated[edge_id] == "WRONG DIRECTION"),
    }


def record_check_block() -> dict:
    """Summarize the items that the record check resolved.

    Returns:
        The rows of the record check, or a note when the file is absent.
    """
    if not RECORD_CHECK.exists():
        return {"file_present": False}
    rows = read_tsv(RECORD_CHECK)
    return {"file_present": True, "items": len(rows),
            "premises_verified": sum(int(row["premises_verified"]) for row in rows),
            "rows": rows}


def main_power(_: argparse.Namespace) -> int:
    """Print the power of the design.

    Args:
        _: The parsed command line, unused.

    Returns:
        Zero.
    """
    print(json.dumps(power_report(), indent=2))
    return 0


def main_compare(args: argparse.Namespace) -> int:
    """Compare the two readings and write the receipt.

    Args:
        args: The parsed command line.

    Returns:
        Zero when both readings are present.
    """
    sentence_pairs, sentence_differ = load_sentence_pairs()
    edge_pairs, edge_differ = load_edge_pairs()
    off_topic_pairs = [(collapse(one, "OFF TOPIC"), collapse(two, "OFF TOPIC"))
                       for one, two in edge_pairs]
    supported_pairs = [(collapse(one, "SUPPORTED"), collapse(two, "SUPPORTED"))
                       for one, two in edge_pairs]
    receipt = {
        "receipt_id": "t2_reading_agreement_v1",
        "created_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "track": "t2_rank",
        "control": "agreement_of_two_independent_readings",
        "queue_item": "R3",
        "data_class_of_run": "A",
        "data_class_of_receipt": "A",
        "receipt_scope": ("Public entity accessions, PubMed identifiers and published sentences "
                          "only. No patient data enters this run."),
        "design": "reports/track2/provenance/second_reading_design.md",
        "design_commit": args.design_commit,
        "script": "scripts/track2_reading_agreement.py",
        "statistic": ("Cohen's kappa with the large sample variance of Fleiss, Cohen and Everitt "
                      "(1969). The interval is a normal interval on the kappa scale."),
        "confidence_level_percent": 100 * (1 - ALPHA),
        "first_reading_of_the_seed_sentences": seed_composition(),
        "inputs": {
            "first_reading_edges": {"path": str(FIRST_EDGES.relative_to(ROOT)),
                                    "sha256_prefix": digest(FIRST_EDGES)},
            "first_reading_sentences": {"path": str(FIRST_SENTENCES.relative_to(ROOT)),
                                        "sha256_prefix": digest(FIRST_SENTENCES)},
            "second_reading_edges": {"path": str(SECOND_EDGES.relative_to(ROOT)),
                                     "sha256_prefix": digest(SECOND_EDGES)},
            "second_reading_sentences": {"path": str(SECOND_SENTENCES.relative_to(ROOT)),
                                         "sha256_prefix": digest(SECOND_SENTENCES)},
        },
        "power_of_the_design": power_report(),
        "seed_sentence_reading": block(sentence_pairs, SENTENCE_CLASSES),
        "edge_reading_five_classes": block(edge_pairs, EDGE_CLASSES),
        "edge_reading_off_topic_collapse": block(off_topic_pairs, SENTENCE_CLASSES),
        "edge_reading_supported_collapse": block(supported_pairs, SENTENCE_CLASSES),
        "disagreements": {"seed_sentences": sentence_differ, "edges": edge_differ},
        "sensitivity_run1_of_the_edge_reading": run1_block(),
        "adjudication": adjudication_block(),
        "record_check": record_check_block(),
        "sensitivity_adjudicated_counts": sensitivity_block(),
    }
    RECEIPT_PATH.write_text(json.dumps(receipt, indent=1) + "\n")
    print(f"Wrote {RECEIPT_PATH.relative_to(ROOT)}")
    for name in ("seed_sentence_reading", "edge_reading_five_classes",
                 "edge_reading_off_topic_collapse", "edge_reading_supported_collapse"):
        item = receipt[name]
        print(f"  {name}: n={item['n']} agreement={item['observed_agreement']:.4f} "
              f"kappa={item['kappa']:.4f} ({item['ci_low']:.4f} to {item['ci_high']:.4f})")
    return 0


def main() -> int:
    """Dispatch the command line.

    Returns:
        The exit code of the chosen subcommand.
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("power", help="Print the power of the design.")
    compare = sub.add_parser("compare", help="Compare the two readings and write the receipt.")
    compare.add_argument("--design-commit", default="", help="The commit that records the design.")
    args = parser.parse_args()
    return {"power": main_power, "compare": main_compare}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
