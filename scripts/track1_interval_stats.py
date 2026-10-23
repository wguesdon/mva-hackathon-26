#!/usr/bin/env python3
"""Confidence intervals for the Track 1 benchmark records. Python standard library only. Class A.

Two kinds of quantity appear in the Track 1 records. A proportion such as pair top-1 is the mean of a
Bernoulli trial, and the Wilson score interval applies to it. A mean such as the mean rank points or
the mean F-max is not a Bernoulli mean, so this module gives it a percentile bootstrap interval and
the record names which interval belongs to which quantity.
"""

from __future__ import annotations

import math
import random
import statistics

BOOTSTRAP_RESAMPLES = 10000
BOOTSTRAP_SEED = 20260917
Z_95 = 1.959963984540054


def wilson(successes: int, trials: int, z: float = Z_95) -> dict:
    """Return the Wilson score interval of a proportion.

    Args:
        successes: The number of successes.
        trials: The number of trials.
        z: The standard normal quantile of the interval.

    Returns:
        The point estimate, the lower bound, the upper bound and the two counts.
    """
    if trials == 0:
        return {"estimate": None, "low": None, "high": None, "successes": 0, "trials": 0,
                "method": "Wilson score interval at 95 percent"}
    proportion = successes / trials
    denominator = 1.0 + z * z / trials
    centre = (proportion + z * z / (2 * trials)) / denominator
    spread = z * math.sqrt(proportion * (1 - proportion) / trials
                           + z * z / (4 * trials * trials)) / denominator
    return {"estimate": round(proportion, 4), "low": round(max(0.0, centre - spread), 4),
            "high": round(min(1.0, centre + spread), 4), "successes": successes, "trials": trials,
            "method": "Wilson score interval at 95 percent"}


def bootstrap_mean(values: list[float], resamples: int = BOOTSTRAP_RESAMPLES,
                   seed: int = BOOTSTRAP_SEED, digits: int = 3) -> dict:
    """Return the percentile bootstrap interval of a mean.

    Args:
        values: The observations.
        resamples: The number of bootstrap resamples.
        seed: The seed of the resampling.
        digits: The number of decimal digits to keep.

    Returns:
        The mean, the lower bound, the upper bound and the method.
    """
    if not values:
        return {"estimate": None, "low": None, "high": None, "n": 0,
                "method": "percentile bootstrap at 95 percent"}
    rng = random.Random(seed)
    size = len(values)
    means = []
    for _ in range(resamples):
        sample = [values[rng.randrange(size)] for _ in range(size)]
        means.append(statistics.mean(sample))
    means.sort()
    low = means[int(0.025 * resamples)]
    high = means[min(resamples - 1, int(0.975 * resamples))]
    return {"estimate": round(statistics.mean(values), digits), "low": round(low, digits),
            "high": round(high, digits), "n": size,
            "method": f"percentile bootstrap at 95 percent, {resamples} resamples, seed {seed}"}


def arm_intervals(results: list[dict]) -> dict:
    """Return every reported quantity of one arm with its interval.

    Args:
        results: The per-scenario results of one arm, with `rank_points`, `f_max` and
            `full_match_rank`.

    Returns:
        The four quantities with their intervals.
    """
    return {
        "mean_rank_points": bootstrap_mean([float(r["rank_points"]) for r in results], digits=2),
        "mean_f_max": bootstrap_mean([float(r["f_max"]) for r in results], digits=4),
        "pair_top1": wilson(sum(1 for r in results if r["full_match_rank"] == 1), len(results)),
        "pair_top10": wilson(sum(1 for r in results if r["full_match_rank"] is not None),
                             len(results)),
    }
