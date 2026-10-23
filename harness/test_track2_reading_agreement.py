"""Test the agreement statistics of the Track 2 second reading against published values.

The script checks Cohen's kappa and its large sample standard error against the worked example of
Fleiss, Cohen and Everitt (1969), the exact binomial interval against two textbook values, the power
function against a hand computation, and the behaviour of the code at perfect agreement, where the
kappa interval is degenerate and the exact interval on raw agreement carries the precision.

Run it directly, or under `make harness-tests`.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import track2_reading_agreement as agreement  # noqa: E402  the scripts path is set above

# Fleiss, Cohen and Everitt (1969) and Fleiss (1981) page 218. Three categories and 100 items. The
# published kappa is 0.429 and the published 95 percent interval is 0.28 to 0.58.
FLEISS_TABLE = [[53, 5, 2], [11, 14, 5], [1, 6, 3]]


def check(condition: bool, message: str) -> None:
    """Stop the run when one condition does not hold.

    Args:
        condition: The condition to test.
        message: The message that names the failed check.

    Raises:
        SystemExit: If the condition does not hold.
    """
    if not condition:
        raise SystemExit(f"FAIL: {message}")


def test_kappa_matches_the_published_example() -> None:
    """Check kappa and its standard error against the published worked example."""
    result = agreement.kappa_interval(FLEISS_TABLE)
    check(abs(result["kappa"] - 0.4286) < 0.0005, f"kappa is {result['kappa']}")
    check(abs(result["observed_agreement"] - 0.70) < 1e-9,
          f"observed agreement is {result['observed_agreement']}")
    check(abs(result["standard_error"] - 0.0760) < 0.0005,
          f"the standard error is {result['standard_error']}")
    check(abs(result["ci_low"] - 0.28) < 0.01 and abs(result["ci_high"] - 0.58) < 0.01,
          f"the interval is {result['ci_low']} to {result['ci_high']}")


def test_kappa_is_zero_under_chance_agreement() -> None:
    """Check that a table with no association gives a kappa of zero."""
    result = agreement.cohen_kappa([[25, 25], [25, 25]])
    check(abs(result["kappa"]) < 1e-9, f"kappa is {result['kappa']}")


def test_perfect_agreement_gives_a_degenerate_interval() -> None:
    """Check that perfect agreement gives kappa one and a zero standard error."""
    result = agreement.kappa_interval([[10, 0], [0, 10]])
    check(abs(result["kappa"] - 1.0) < 1e-12, f"kappa is {result['kappa']}")
    check(result["standard_error"] == 0.0, f"the standard error is {result['standard_error']}")
    check(result["ci_low"] == 1.0 and result["ci_high"] == 1.0, "the interval is not degenerate")


def test_exact_binomial_interval_matches_textbook_values() -> None:
    """Check the exact binomial interval against two published values."""
    none_of_ten = agreement.clopper_pearson(0, 10)
    check(none_of_ten["ci_low"] == 0.0, f"the lower limit is {none_of_ten['ci_low']}")
    check(abs(none_of_ten["ci_high"] - 0.3085) < 0.0005,
          f"the upper limit is {none_of_ten['ci_high']}")
    two_of_ten = agreement.clopper_pearson(2, 10)
    check(abs(two_of_ten["ci_low"] - 0.0252) < 0.0005, f"the lower limit is {two_of_ten['ci_low']}")
    check(abs(two_of_ten["ci_high"] - 0.5561) < 0.0005,
          f"the upper limit is {two_of_ten['ci_high']}")
    all_of_many = agreement.clopper_pearson(230, 230)
    check(abs(all_of_many["ci_low"] - 0.9841) < 0.0005,
          f"the lower limit is {all_of_many['ci_low']}")
    check(all_of_many["ci_high"] == 1.0, f"the upper limit is {all_of_many['ci_high']}")


def test_power_matches_the_design_table() -> None:
    """Check the power function against the two values of the design document."""
    sentences = agreement.kappa_power(230, 49 / 230)
    check(abs(sentences["power"] - 0.9317) < 0.0005, f"the power is {sentences['power']}")
    edges = agreement.kappa_power(43, 12 / 43, alternative_kappa=0.90)
    check(abs(edges["power"] - 0.6672) < 0.0005, f"the power is {edges['power']}")


def test_confusion_matrix_counts_every_item() -> None:
    """Check that the confusion matrix keeps every item and the class order."""
    pairs = [("yes", "yes"), ("yes", "no"), ("no", "no"), ("no", "no")]
    rows = agreement.confusion(pairs, ("yes", "no"))
    check(sum(row["total"] for row in rows) == len(pairs), "the matrix loses an item")
    check(rows[0]["yes"] == 1 and rows[0]["no"] == 1, f"the first row is {rows[0]}")
    check(rows[1]["no"] == 2, f"the second row is {rows[1]}")


def test_collapse_maps_one_class_to_the_positive_label() -> None:
    """Check the collapse of five classes onto a binary judgement."""
    check(agreement.collapse("OFF TOPIC", "OFF TOPIC") == "yes", "the positive class is wrong")
    check(agreement.collapse("SUPPORTED", "OFF TOPIC") == "no", "the negative class is wrong")


def main() -> int:
    """Run every test of this file.

    Returns:
        Zero when every test passes.
    """
    test_kappa_matches_the_published_example()
    test_kappa_is_zero_under_chance_agreement()
    test_perfect_agreement_gives_a_degenerate_interval()
    test_exact_binomial_interval_matches_textbook_values()
    test_power_matches_the_design_table()
    test_confusion_matrix_counts_every_item()
    test_collapse_maps_one_class_to_the_positive_label()
    print("OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
