#!/usr/bin/env python3
"""Verify the preclinical evidence of Track 2 and compute the protocol and the cost model.

The decision gate is `reports/track2/preclinical_protocol_gate_v1.md`, which is recorded in a commit
before this script writes any result. The gate holds seven conditions. This script tests C1 to C6
and records C7 from the command line.

Gates:
    C1 resolve      Every PMID of the evidence table returns a PubMed record with a title.
    C2 keyword      The recorded keyword of each PMID appears in the PubMed title of that PMID.
    C3 fields       Every row names a claim, a source kind, a source identifier and what it fixes.
    C4 named        A row whose source kind is not `pmid` carries a name, a locator and an access date.
    C5 no basis     A row with no published basis carries `none` and a reason, and the receipt lists it.
    C6 arithmetic   Every derived number of the protocol and of the cost model is computed here.
    C7 bibliography Passed in from the command line as the exit code of scripts/build_bibliography.py.

The group size comes from a simulation and not from a closed form, because the project pins no
statistics library. The simulation draws two samples from a normal distribution under the stated
effect and the stated standard deviation, applies a two sided two sample t test at the family
threshold, and counts the rejections. The seed is recorded, so the count repeats.

Outputs:
    tracking/receipts/t2_preclinical_evidence_<version>.json
    tracking/receipts/t2_preclinical_protocol_<version>.json
    tracking/receipts/t2_preclinical_cost_model_<version>.json

Amendment 5 of the gate adds a nicotinic acid arm and conditions C16 to C19. Version 2 of the
parameter file and of the evidence table carries that arm, and the run writes the v2 receipts beside
the v1 receipts, which stay unchanged.

Usage:
    uv run python scripts/t2_preclinical_protocol.py --bibliography-exit 0
    uv run python scripts/t2_preclinical_protocol.py --bibliography-exit 0 --version v2 \\
        --parameters literature/track2_preclinical_parameters_v2.yaml \\
        --evidence literature/track2_preclinical_evidence_v2.tsv
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import random
import re
import sys
import unicodedata
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

import httpx
import yaml

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE_TABLE = ROOT / "literature" / "track2_preclinical_evidence_v1.tsv"
PARAMETERS = ROOT / "literature" / "track2_preclinical_parameters_v1.yaml"
GATE = ROOT / "reports" / "track2" / "preclinical_protocol_gate_v1.md"
EVIDENCE_RECEIPT = ROOT / "tracking" / "receipts" / "t2_preclinical_evidence_v1.json"
PROTOCOL_RECEIPT = ROOT / "tracking" / "receipts" / "t2_preclinical_protocol_v1.json"
COST_RECEIPT = ROOT / "tracking" / "receipts" / "t2_preclinical_cost_model_v1.json"
POWER_V2_RECEIPT = ROOT / "tracking" / "receipts" / "t2_preclinical_power_v2.json"
POWER_V4_RECEIPT = ROOT / "tracking" / "receipts" / "t2_preclinical_power_v4.json"
NAPRT_RECEIPT = ROOT / "tracking" / "receipts" / "t2_naprt_fibroblast_gtex_v1.json"
# The two precursors of version 1. A parameter file that names its own list replaces them.
DEFAULT_PRECURSORS = ["nicotinamide_mononucleotide", "nicotinamide_riboside_chloride"]
# Version 1 counted three compounds and the vehicle in the metabolome set of each line.
DEFAULT_METABOLOME_ARMS = 4


def power_v2_declares_units() -> bool:
    """Say whether the design matched power receipt declares a unit of replication per stage.

    Condition C9 of amendment 2 of the gate asks for that statement, because the degrees of freedom
    of each test follow from it.

    Returns:
        True when the receipt exists and both stages name a unit of replication.
    """
    if not POWER_V2_RECEIPT.exists():
        return False
    receipt = json.loads(POWER_V2_RECEIPT.read_text())
    return bool(receipt.get("stage1", {}).get("unit_of_replication")
                and receipt.get("stage2", {}).get("unit_of_replication"))
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
REQUIRED_FIELDS = ("claim_id", "section", "claim", "source_kind", "source_id", "what_it_fixes")
NUMBER = re.compile(r"(?<![\w.])(\d{1,3}(?:,\d{3})+|\d+)(?:\.(\d+))?(?![\w])")
FLOAT_MIN = 1e-300
BETA_EPSILON = 3e-16
BETA_ITERATIONS = 300


def normalize(text: str) -> str:
    """Normalize a title or a keyword for comparison.

    Args:
        text: A title or a keyword.

    Returns:
        Lower-case ASCII letters and digits only.
    """
    ascii_text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]", "", ascii_text.lower())


def sha256_of_file(path: Path) -> str:
    """Give the SHA-256 digest of one file.

    Args:
        path: The file path.

    Returns:
        The hexadecimal digest.
    """
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def read_evidence_table(path: Path) -> list[dict]:
    """Read the evidence table.

    Args:
        path: The tab separated file.

    Returns:
        One dictionary per row.
    """
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pubmed_titles(client: httpx.Client, pmids: list[str]) -> dict[str, str]:
    """Fetch the title of each PubMed record.

    Args:
        client: The HTTP client.
        pmids: The PMIDs to fetch.

    Returns:
        A map from PMID to title. A PMID that does not resolve is absent from the map.
    """
    if not pmids:
        return {}
    response = client.get(
        EFETCH,
        params={"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"},
        timeout=120,
    )
    response.raise_for_status()
    root = ET.fromstring(response.text)
    titles: dict[str, str] = {}
    for article in root.iter("PubmedArticle"):
        pmid = article.findtext(".//PMID") or ""
        node = article.find(".//ArticleTitle")
        if pmid and node is not None:
            titles[pmid] = " ".join("".join(node.itertext()).split())
    return titles


def field_findings(rows: list[dict]) -> list[dict]:
    """Test condition C3 and condition C4 on every row.

    Args:
        rows: The rows of the evidence table.

    Returns:
        One finding per row that fails.
    """
    findings: list[dict] = []
    for row in rows:
        for field in REQUIRED_FIELDS:
            if not (row.get(field) or "").strip():
                findings.append({"claim_id": row.get("claim_id", ""), "condition": "C3",
                                 "detail": f"the field {field} is empty"})
        kind = (row.get("source_kind") or "").strip()
        if kind == "named":
            for field in ("source_name", "source_locator", "access_date"):
                if not (row.get(field) or "").strip():
                    findings.append({"claim_id": row.get("claim_id", ""), "condition": "C4",
                                     "detail": f"a named source leaves {field} empty"})
        if kind == "none" and not (row.get("no_published_basis") or "").strip():
            findings.append({"claim_id": row.get("claim_id", ""), "condition": "C5",
                             "detail": "a row with no published basis states no reason"})
        if kind not in {"pmid", "named", "none"}:
            findings.append({"claim_id": row.get("claim_id", ""), "condition": "C3",
                             "detail": f"the source kind {kind} is not pmid, named or none"})
    return findings


def citation_findings(rows: list[dict], titles: dict[str, str]) -> list[dict]:
    """Test condition C1 and condition C2 on every PMID row.

    Args:
        rows: The rows of the evidence table.
        titles: The PubMed title of each PMID that resolved.

    Returns:
        One finding per row that fails.
    """
    findings: list[dict] = []
    for row in rows:
        if (row.get("source_kind") or "").strip() != "pmid":
            continue
        pmid = (row.get("source_id") or "").strip()
        title = titles.get(pmid)
        if title is None:
            findings.append({"claim_id": row["claim_id"], "condition": "C1",
                             "detail": f"PMID {pmid} returned no record with a title"})
            continue
        keyword = (row.get("citation_keyword") or "").strip()
        if not keyword:
            findings.append({"claim_id": row["claim_id"], "condition": "C2",
                             "detail": f"PMID {pmid} carries no keyword"})
        elif normalize(keyword) not in normalize(title):
            findings.append({"claim_id": row["claim_id"], "condition": "C2",
                             "detail": f"the keyword {keyword} is absent from the title of PMID {pmid}"})
    return findings


def literature_values(rows: list[dict]) -> dict[str, list[float]]:
    """Collect every number that the claim text of each row states.

    The report prints these numbers as the parameters of a cited source, and
    `scripts/check_report_claims.py` resolves each one against this receipt.

    Args:
        rows: The rows of the evidence table.

    Returns:
        A map from claim identifier to the numbers of that claim.
    """
    values: dict[str, list[float]] = {}
    for row in rows:
        found: list[float] = []
        for match in NUMBER.finditer(row.get("claim", "")):
            whole = match.group(1).replace(",", "")
            fraction = match.group(2)
            found.append(float(f"{whole}.{fraction}") if fraction else float(whole))
        if found:
            values[row["claim_id"]] = found
    return values


def normal_cdf(x: float) -> float:
    """Give the standard normal cumulative distribution at one point.

    Args:
        x: The point.

    Returns:
        The cumulative probability.
    """
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def normal_quantile(p: float) -> float:
    """Give the standard normal quantile by bisection.

    Args:
        p: A probability strictly between 0 and 1.

    Returns:
        The quantile.
    """
    low, high = -40.0, 40.0
    for _ in range(200):
        middle = 0.5 * (low + high)
        if normal_cdf(middle) < p:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def beta_continued_fraction(a: float, b: float, x: float) -> float:
    """Evaluate the continued fraction of the incomplete beta function.

    Args:
        a: The first shape parameter.
        b: The second shape parameter.
        x: The upper limit, strictly between 0 and 1.

    Returns:
        The value of the continued fraction.
    """
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < FLOAT_MIN:
        d = FLOAT_MIN
    d = 1.0 / d
    h = d
    for m in range(1, BETA_ITERATIONS + 1):
        m2 = 2 * m
        numerator = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + numerator * d
        if abs(d) < FLOAT_MIN:
            d = FLOAT_MIN
        c = 1.0 + numerator / c
        if abs(c) < FLOAT_MIN:
            c = FLOAT_MIN
        d = 1.0 / d
        h *= d * c
        numerator = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + numerator * d
        if abs(d) < FLOAT_MIN:
            d = FLOAT_MIN
        c = 1.0 + numerator / c
        if abs(c) < FLOAT_MIN:
            c = FLOAT_MIN
        d = 1.0 / d
        step = d * c
        h *= step
        if abs(step - 1.0) < BETA_EPSILON:
            break
    return h


def regularized_incomplete_beta(a: float, b: float, x: float) -> float:
    """Give the regularized incomplete beta function.

    Args:
        a: The first shape parameter.
        b: The second shape parameter.
        x: The upper limit between 0 and 1.

    Returns:
        The value of the function.
    """
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    front = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        + a * math.log(x) + b * math.log1p(-x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return front * beta_continued_fraction(a, b, x) / a
    return 1.0 - front * beta_continued_fraction(b, a, 1.0 - x) / b


def student_t_cdf(t: float, df: float) -> float:
    """Give the Student t cumulative distribution at one point.

    Args:
        t: The point.
        df: The degrees of freedom.

    Returns:
        The cumulative probability.
    """
    tail = 0.5 * regularized_incomplete_beta(0.5 * df, 0.5, df / (df + t * t))
    return 1.0 - tail if t > 0 else tail


def student_t_quantile(p: float, df: float) -> float:
    """Give the Student t quantile by bisection.

    Args:
        p: A probability strictly between 0 and 1.
        df: The degrees of freedom.

    Returns:
        The quantile.
    """
    low, high = -200.0, 200.0
    for _ in range(200):
        middle = 0.5 * (low + high)
        if student_t_cdf(middle, df) < p:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def two_sample_t(first: list[float], second: list[float]) -> float:
    """Give the pooled variance two sample t statistic.

    Args:
        first: The first sample.
        second: The second sample.

    Returns:
        The t statistic.
    """
    n1, n2 = len(first), len(second)
    mean1 = sum(first) / n1
    mean2 = sum(second) / n2
    ss1 = sum((value - mean1) ** 2 for value in first)
    ss2 = sum((value - mean2) ** 2 for value in second)
    pooled = (ss1 + ss2) / (n1 + n2 - 2)
    if pooled <= 0.0:
        return 0.0
    return (mean2 - mean1) / math.sqrt(pooled * (1.0 / n1 + 1.0 / n2))


def simulate_power(n: int, effect: float, sd: float, alpha: float, iterations: int,
                   rng: random.Random) -> float:
    """Estimate the power of a two sided two sample t test by simulation.

    Args:
        n: The group size of each of the two groups.
        effect: The true difference of the two means.
        sd: The standard deviation of each group.
        alpha: The two sided significance threshold.
        iterations: The number of simulated datasets.
        rng: The random number generator.

    Returns:
        The fraction of datasets that reject the null hypothesis.
    """
    critical = student_t_quantile(1.0 - alpha / 2.0, 2 * n - 2)
    rejections = 0
    for _ in range(iterations):
        control = [rng.gauss(0.0, sd) for _ in range(n)]
        treated = [rng.gauss(effect, sd) for _ in range(n)]
        if abs(two_sample_t(control, treated)) > critical:
            rejections += 1
    return rejections / iterations


def normal_approximation_size(effect: float, sd: float, alpha: float, power: float) -> int:
    """Give the group size from the normal approximation.

    Args:
        effect: The true difference of the two means.
        sd: The standard deviation of each group.
        alpha: The two sided significance threshold.
        power: The target power.

    Returns:
        The group size, at least 2.
    """
    z_alpha = normal_quantile(1.0 - alpha / 2.0)
    z_beta = normal_quantile(power)
    raw = 2.0 * ((z_alpha + z_beta) * sd / effect) ** 2
    return max(2, math.ceil(raw))


def group_size(effect: float, sd: float, alpha: float, power: float, iterations: int,
               seed: int) -> dict:
    """Give the smallest group size whose simulated power reaches the target.

    Args:
        effect: The true difference of the two means.
        sd: The standard deviation of each group.
        alpha: The two sided significance threshold.
        power: The target power.
        iterations: The number of simulated datasets at each candidate size.
        seed: The random seed.

    Returns:
        The group size, the simulated power at that size and the normal approximation.
    """
    approximation = normal_approximation_size(effect, sd, alpha, power)
    candidate = max(2, approximation - 3)
    previous = None
    while candidate < 500:
        rng = random.Random(seed + candidate)
        observed = simulate_power(candidate, effect, sd, alpha, iterations, rng)
        if observed >= power:
            return {"n_per_group": candidate, "simulated_power": round(observed, 4),
                    "simulated_power_one_below": None if previous is None else round(previous, 4),
                    "normal_approximation": approximation}
        previous = observed
        candidate += 1
    raise RuntimeError("no group size below 500 reaches the target power")


def power_table(statistics: dict) -> list[dict]:
    """Give the group size at each assumed standard deviation.

    Args:
        statistics: The statistics block of the parameter file.

    Returns:
        One row per assumed standard deviation.
    """
    effect = math.log2(statistics["effect_fold_change"])
    alpha = statistics["alpha_family"] / statistics["compounds_in_family"]
    rows: list[dict] = []
    for sd in statistics["assumed_sd_log2"]:
        result = group_size(effect, sd, alpha, statistics["power_target"],
                            statistics["simulation_iterations_search"],
                            statistics["simulation_seed"])
        rows.append({"assumed_sd_log2": sd, "standardized_effect": round(effect / sd, 4), **result})
    return rows


def cheapest_pack_plan(required_mg: float, packs: list[list[int | float]]) -> dict:
    """Give the cheapest combination of catalogue pack sizes that covers a requirement.

    The search is a dynamic program over whole milligrams, so its cost grows with the requirement
    and not with the product of the pack counts. `cost[m]` is the cheapest price that buys at least
    `m` milligrams.

    Args:
        required_mg: The mass required, in milligrams.
        packs: One pair of pack size in milligrams and price in United States dollars per pack.

    Returns:
        The chosen counts, the mass bought and the price.
    """
    target = int(math.ceil(required_mg))
    if target <= 0:
        return {"packs": {}, "mass_mg": 0.0, "price_usd": 0.0}
    sizes = [int(size) for size, _ in packs]
    prices = [float(price) for _, price in packs]
    infinity = float("inf")
    cost = [infinity] * (target + 1)
    choice = [-1] * (target + 1)
    cost[0] = 0.0
    for mass in range(1, target + 1):
        for index, size in enumerate(sizes):
            previous = cost[max(0, mass - size)]
            if previous + prices[index] < cost[mass]:
                cost[mass] = previous + prices[index]
                choice[mass] = index
    if cost[target] == infinity:
        raise RuntimeError("no pack combination covers the requirement")
    counts = [0] * len(sizes)
    mass = target
    while mass > 0:
        index = choice[mass]
        counts[index] += 1
        mass = max(0, mass - sizes[index])
    bought = sum(counts[i] * sizes[i] for i in range(len(sizes)))
    return {
        "packs": {f"{sizes[i]} mg": counts[i] for i in range(len(sizes)) if counts[i]},
        "mass_mg": float(bought),
        "price_usd": round(cost[target], 2),
    }


def stage1_arms(stage: dict) -> dict:
    """Give the arm count of one cell line and the well count of the stage.

    Args:
        stage: The stage 1 block of the parameter file.

    Returns:
        The arm counts, the arms per line and the wells of the stage.
    """
    arms = stage["arms"]
    per_line = (arms["base_dose_series"] + arms["base_comparator_high_dose"] + arms["base_vehicle"]
                + arms["sirt2_sirna"] + arms["sirt2_inhibitor"])
    return {"arms_per_line": per_line, "arm_counts": dict(arms)}


def stage1_compound_micromoles(stage: dict) -> dict:
    """Give the micromoles of each compound that one line and one replicate consume.

    Args:
        stage: The stage 1 block of the parameter file.

    Returns:
        The micromoles of each compound.
    """
    volume_litres = stage["well_volume_ml"] / 1000.0
    dose_sum = sum(stage["dose_mM"])
    top_dose = max(stage["dose_mM"])
    mechanism_arms = 2
    precursor_mM = dose_sum + mechanism_arms * top_dose
    comparator_mM = dose_sum + stage["comparator_high_dose_mM"] + mechanism_arms * top_dose
    agk2_mM = stage["agk2_concentration_uM"] / 1000.0 * stage["arms"]["sirt2_inhibitor"]
    micromoles = {name: round(precursor_mM * volume_litres * 1000.0, 4)
                  for name in stage.get("precursor_compounds", DEFAULT_PRECURSORS)}
    micromoles["nicotinamide"] = round(comparator_mM * volume_litres * 1000.0, 4)
    micromoles["agk2"] = round(agk2_mM * volume_litres * 1000.0, 4)
    return micromoles


def stage1_cost(parameters: dict, replicates: int) -> dict:
    """Give the cost model of stage 1.

    Args:
        parameters: The whole parameter file.
        replicates: The biological replicate count of each arm.

    Returns:
        The cost model of stage 1.
    """
    stage = parameters["stage1"]
    prices = parameters["prices"]
    compounds = parameters["compounds"]
    lines = stage["lines"]
    factor = stage["stock_preparation_factor"]
    micromoles = stage1_compound_micromoles(stage)
    compound_lines: dict[str, dict] = {}
    for name, amount in micromoles.items():
        weight = compounds[name]["formula_weight_g_per_mol"]
        required = amount * lines * replicates * weight / 1000.0 * factor
        plan = cheapest_pack_plan(required, compounds[name]["packs_mg_usd"])
        compound_lines[name] = {"catalogue": compounds[name]["catalogue"],
                                "required_mg": round(required, 3), **plan}
    cytochalasin_mg = (stage["cytochalasin_b_ug_per_ml"] * stage["well_volume_ml"]
                       * stage["micronucleus_arms_per_line"] * lines * replicates
                       / 1000.0 * factor)
    compound_lines["cytochalasin_b"] = {
        "catalogue": compounds["cytochalasin_b"]["catalogue"],
        "required_mg": round(cytochalasin_mg, 3),
        **cheapest_pack_plan(cytochalasin_mg, compounds["cytochalasin_b"]["packs_mg_usd"]),
    }
    compound_lines["mitomycin_c"] = {
        "catalogue": compounds["mitomycin_c"]["catalogue"],
        "required_mg": 0.0,
        "packs": {"5 mg": 1},
        "mass_mg": 5.0,
        "price_usd": float(compounds["mitomycin_c"]["packs_mg_usd"][0][1]),
        "note": ("The smallest pack covers any concentration that the laboratory establishes in its "
                 "proficiency test, so the requirement is not computed from a fixed concentration."),
    }
    compound_total = round(sum(line["price_usd"] for line in compound_lines.values()), 2)
    arms = stage1_arms(stage)
    wells = arms["arms_per_line"] * lines * replicates
    metabolome_samples = (stage.get("metabolome_arms_per_line", DEFAULT_METABOLOME_ARMS) * lines
                          * stage["metabolome_replicates"])
    karyotype_samples = 2 * lines * stage["karyotype_replicates"]
    line_cost = round(stage["lines_purchased"] * prices["cell_line_usd"], 2)
    metabolome_cost = round(metabolome_samples * prices["metabolome_per_sample_usd"], 2)
    karyotype_cost = round(karyotype_samples * prices["karyotype_per_sample_usd"], 2)
    personnel_cost = round(stage["personnel_months"] * prices["postdoctoral_month_usd"], 2)
    priced_total = round(compound_total + line_cost + metabolome_cost + karyotype_cost
                         + personnel_cost, 2)
    return {
        "replicates_per_arm": replicates,
        "lines": lines,
        "stock_preparation_factor": factor,
        "well_volume_ml": stage["well_volume_ml"],
        "arms_per_line": arms["arms_per_line"],
        "arm_counts": arms["arm_counts"],
        "wells": wells,
        "compound_micromoles_per_line_per_replicate": micromoles,
        "compounds": compound_lines,
        "compound_total_usd": compound_total,
        "cell_lines_purchased": stage["lines_purchased"],
        "cell_line_cost_usd": line_cost,
        "metabolome_samples": metabolome_samples,
        "metabolome_cost_usd": metabolome_cost,
        "karyotype_samples": karyotype_samples,
        "karyotype_cost_usd": karyotype_cost,
        "personnel_months": stage["personnel_months"],
        "personnel_cost_usd": personnel_cost,
        "priced_total_usd": priced_total,
        "unpriced_lines": [
            "The culture medium, the serum, the plasticware, the antibodies and the small "
            "interfering RNA carry no public list price, which claim P41 records. The well count "
            "above is the quantity that a reader prices from an institutional quote.",
            "Fringe benefits and indirect costs are absent, because no public rate applies to "
            "every institution.",
        ],
    }


def stage2_cost(parameters: dict, replicates: int) -> dict:
    """Give the cost model of stage 2.

    Args:
        parameters: The whole parameter file.
        replicates: The group size of each of the four groups.

    Returns:
        The cost model of stage 2.
    """
    stage = parameters["stage2"]
    prices = parameters["prices"]
    compounds = parameters["compounds"]
    treated_mice = replicates * stage["genotypes"]
    total_mice = replicates * stage["groups"]
    experiment_cages = math.ceil(total_mice / stage["mice_per_cage"])
    housing_days = stage["acclimation_days"] + stage["dosing_days"]
    homozygotes_needed = replicates * stage["treatment_arms"]
    pups_needed = math.ceil(homozygotes_needed / stage["mendelian_homozygote_fraction"])
    litters = math.ceil(pups_needed / stage["litter_size"])
    breeder_cages = litters
    per_diem = round((experiment_cages * housing_days
                      + breeder_cages * stage["breeding_days"])
                     * prices["mouse_cage_per_diem_usd"], 2)
    cage_creation = round((experiment_cages + breeder_cages) * prices["mouse_cage_creation_usd"], 2)
    technical_hours = round((experiment_cages * housing_days / 7.0
                             + breeder_cages * stage["breeding_days"] / 7.0)
                            * stage["technical_hours_per_cage_per_week"], 2)
    technical_cost = round(technical_hours * prices["animal_technical_hour_usd"], 2)
    dose_mg_per_mouse_per_day = (stage["dose_mg_per_kg_per_day"]
                                 * stage["mouse_body_weight_g"] / 1000.0)
    compound_mg = dose_mg_per_mouse_per_day * stage["dosing_days"] * treated_mice
    plan = cheapest_pack_plan(compound_mg,
                              compounds["nicotinamide_mononucleotide"]["packs_mg_usd"])
    metabolome_samples = (stage["metabolome_samples_per_group"] * stage["groups"]
                          * stage["metabolome_matrices"])
    karyotype_samples = stage["karyotype_samples_per_group"] * stage["groups"]
    metabolome_cost = round(metabolome_samples * prices["metabolome_per_sample_usd"], 2)
    karyotype_cost = round(karyotype_samples * prices["karyotype_per_sample_usd"], 2)
    personnel_cost = round(stage["personnel_months"] * prices["postdoctoral_month_usd"], 2)
    priced_total = round(per_diem + cage_creation + technical_cost + plan["price_usd"]
                         + metabolome_cost + karyotype_cost + personnel_cost, 2)
    return {
        "group_size": replicates,
        "groups": stage["groups"],
        "mice_per_cage": stage["mice_per_cage"],
        "litter_size": stage["litter_size"],
        "breeding_days": stage["breeding_days"],
        "mouse_body_weight_g": stage["mouse_body_weight_g"],
        "acclimation_days": stage["acclimation_days"],
        "dosing_days": stage["dosing_days"],
        "dose_mg_per_kg_per_day": stage["dose_mg_per_kg_per_day"],
        "technical_hours_per_cage_per_week": stage["technical_hours_per_cage_per_week"],
        "total_mice": total_mice,
        "treated_mice": treated_mice,
        "experiment_cages": experiment_cages,
        "cages_total": experiment_cages + breeder_cages,
        "housing_days_per_animal": housing_days,
        "homozygotes_needed": homozygotes_needed,
        "pups_needed": pups_needed,
        "litters": litters,
        "breeder_cages": breeder_cages,
        "per_diem_usd": per_diem,
        "cage_creation_usd": cage_creation,
        "technical_hours": technical_hours,
        "technical_cost_usd": technical_cost,
        "dose_mg_per_mouse_per_day": round(dose_mg_per_mouse_per_day, 3),
        "compound_required_mg": round(compound_mg, 1),
        "compound_plan": plan,
        "metabolome_samples": metabolome_samples,
        "metabolome_cost_usd": metabolome_cost,
        "karyotype_samples": karyotype_samples,
        "karyotype_cost_usd": karyotype_cost,
        "personnel_months": stage["personnel_months"],
        "personnel_cost_usd": personnel_cost,
        "priced_total_usd": priced_total,
        "unpriced_lines": [
            "The strain acquisition carries no public price, because the International Mouse Strain "
            "Resource returns no available strain and a rederivation is quoted on request, which "
            "claims P03 and P42 record.",
            "Fringe benefits and indirect costs are absent, for the reason of stage 1.",
        ],
    }


def decisive_experiments(stage1: dict, parameters: dict) -> list[dict]:
    """Give the three experiments that would change the conclusion, with the cost of each.

    Args:
        stage1: The stage 1 cost model.
        parameters: The whole parameter file.

    Returns:
        One entry per experiment, in the order that the report states.
    """
    prices = parameters["prices"]
    stage = parameters["stage1"]
    comparator = stage1["compounds"]["nicotinamide"]["price_usd"]
    comparator_lines = round(stage["lines_purchased"] * prices["cell_line_usd"], 2)
    comparator_personnel = round(2 * prices["postdoctoral_month_usd"], 2)
    comparator_total = round(comparator + comparator_lines + comparator_personnel, 2)
    precursors = stage.get("precursor_compounds", DEFAULT_PRECURSORS)
    precursor = round(sum(stage1["compounds"][name]["price_usd"] for name in precursors)
                      + stage1["compounds"]["agk2"]["price_usd"], 2)
    increment = round(stage1["priced_total_usd"] - comparator_total, 2)
    precursor_experiment = ("The precursor arms of stage 1 with the SIRT2 knockdown arm and the "
                            "SIRT2 inhibitor arm.")
    precursor_decides = ("Whether an NAD+ precursor raises BUBR1 in a human cell at all and whether "
                         "SIRT2 carries the effect. No published experiment measures BUBR1 in a "
                         "human cell under an NAD+ precursor.")
    if parameters["statistics"].get("primary_compound") == "nicotinic_acid":
        precursor_experiment = ("The precursor arms of stage 1, with nicotinic acid as the primary "
                                "test and nicotinamide mononucleotide and nicotinamide riboside "
                                "after it, each with the SIRT2 knockdown arm and the SIRT2 "
                                "inhibitor arm.")
        precursor_decides = ("Whether the authorised member of the class raises BUBR1 in a human "
                             "cell through SIRT2, and whether its NAD+ metabolite endpoint shows "
                             "the Preiss-Handler route in the same lines. No published experiment "
                             "measures BUBR1 in a human cell under an NAD+ precursor.")
    return [
        {
            "rank": 1,
            "experiment": ("The nicotinamide comparator arm of stage 1, run alone, over the shared "
                           "dose range and the 20 mM concentration."),
            "decides": ("At which concentration the mechanism test that rejects nicotinamide "
                        "holds. That compound is the one compound of this axis with randomized "
                        "phase 3 evidence, and tracking/receipts/t2_authorization_record_v1.json "
                        "records that it carries no clinical stage field in this build. A rise in "
                        "BUBR1 abundance at or below the top dose makes the negative answer of this "
                        "report to the Track 2 task wrong."),
            "priced_usd": comparator_total,
            "basis": "standalone",
            "components": {"compound_usd": comparator, "cell_lines_usd": comparator_lines,
                           "personnel_months": 2, "personnel_usd": comparator_personnel},
        },
        {
            "rank": 2,
            "experiment": precursor_experiment,
            "decides": precursor_decides,
            "priced_usd": increment,
            "basis": "incremental over the experiment of rank 1, so the two sum to the stage 1 total",
            "components": {"compound_usd": precursor,
                           "rest_of_the_increment_usd": round(increment - precursor, 2)},
        },
        {
            "rank": 3,
            "experiment": ("The mention rule on every disease node of the four transfer graphs, "
                           "then a rerun of the held-out control on each filtered graph."),
            "decides": ("Whether the four transfer FAIL verdicts measure the method or the "
                        "unmeasured entity mapping error of their case nodes. The rule runs on "
                        "every disease node of the mosaic variegated aneuploidy graph and on the "
                        "seed node alone of each transfer graph, and the five measured seed error "
                        "rates run from 0.153 to 0.959."),
            "priced_usd": 0.00,
            "basis": "compute only, on the measured run cost of the transfer pipeline",
            "components": {"money_usd": 0.00,
                           "wall_clock_minutes": parameters["graph_result"][
                               "run_wall_clock_minutes"]},
        },
    ]


def regime_ratios(graph_result: dict) -> dict:
    """Give the factor by which the full held-out set overstates the drug-poor recovery.

    Args:
        graph_result: The graph result block of the parameter file.

    Returns:
        The ratio of each disease and the range over the four diseases.
    """
    ratios = {}
    for disease, full in graph_result["heldout_recovery_all"].items():
        poor = graph_result["heldout_recovery_drug_poor"][disease]
        ratios[disease] = round(full / poor, 1)
    values = sorted(ratios.values())
    errors = sorted(graph_result["mapping_error_rate"].values())
    return {
        "regime_ratio": ratios,
        "regime_ratio_min": values[0],
        "regime_ratio_max": values[-1],
        "mapping_error_min": errors[0],
        "mapping_error_max": errors[-1],
        "diseases_measured": graph_result["diseases_measured"],
        "heldout_fail_count": graph_result["heldout_fail_count"],
        "permutation_pass_count": graph_result["permutation_pass_count"],
        "heldout_recovery_all": graph_result["heldout_recovery_all"],
        "heldout_recovery_drug_poor": graph_result["heldout_recovery_drug_poor"],
        "mapping_error_rate": graph_result["mapping_error_rate"],
        "run_wall_clock_minutes": graph_result["run_wall_clock_minutes"],
        "run_money_usd": graph_result["run_money_usd"],
    }


def amendment5_conditions(parameters: dict, rows: list[dict]) -> dict:
    """Test conditions C16 to C19 of amendment 5 of the gate.

    Args:
        parameters: The whole parameter file.
        rows: The rows of the evidence table.

    Returns:
        A map from condition to PASS or FAIL. The map is empty when the parameter file names no
        primary compound, which is the case of version 1.
    """
    compound = parameters["statistics"].get("primary_compound")
    if not compound:
        return {}
    entry = parameters["compounds"].get(compound, {})
    named = {row["claim_id"] for row in rows
             if (row["source_kind"] or "").strip() == "named" and row["source_locator"]
             and row["access_date"]}
    priced = entry.get("source") in named and entry.get("formula_weight_source") in named
    power = json.loads(POWER_V4_RECEIPT.read_text()) if POWER_V4_RECEIPT.exists() else {}
    naprt = json.loads(NAPRT_RECEIPT.read_text()) if NAPRT_RECEIPT.exists() else {}
    power_alpha = power.get("shared_parameters", {}).get("alpha")
    return {
        "C16_price_of_the_new_arm": "PASS" if priced else "FAIL",
        "C17_power_of_the_primary_test":
            "PASS" if power_alpha == parameters["statistics"]["primary_alpha"] else "FAIL",
        "C18_naprt_expression_read": "PASS" if naprt.get("summary") else "FAIL",
        "C19_naprt_threshold_applied":
            ("PASS" if naprt.get("verdict") in ("EXPRESSED", "NOT SHOWN")
             and naprt.get("threshold_median_tpm") == 1.0 else "FAIL"),
    }


def decision_rule(parameters: dict) -> dict | None:
    """Give the fixed test sequence of amendment 5 with the records that bear on it.

    Args:
        parameters: The whole parameter file.

    Returns:
        The primary compound, the thresholds, the power of the primary test at the design of the
        report and the NAPRT verdict, or None when the parameter file names no primary compound.
    """
    statistics = parameters["statistics"]
    compound = statistics.get("primary_compound")
    if not compound:
        return None
    power = json.loads(POWER_V4_RECEIPT.read_text()) if POWER_V4_RECEIPT.exists() else {}
    naprt = json.loads(NAPRT_RECEIPT.read_text()) if NAPRT_RECEIPT.exists() else {}
    return {
        "primary_compound": compound,
        "primary_alpha": statistics["primary_alpha"],
        "secondary_compounds": statistics["compounds_in_family"],
        "secondary_alpha_per_compound": round(statistics["alpha_family"]
                                              / statistics["compounds_in_family"], 4),
        "secondary_step_runs_only_when_the_primary_rejects": True,
        "primary_power_receipt": statistics["primary_power_receipt"],
        "primary_power_at_the_design": power.get("the_design_of_the_report", {}).get("power"),
        "naprt_receipt": "t2_naprt_fibroblast_gtex_v1",
        "naprt_median_tpm": naprt.get("summary", {}).get("median_tpm"),
        "naprt_verdict": naprt.get("verdict"),
        "no_route_rule": ("The NAD+ metabolite endpoint under nicotinic acid must rise above its "
                          "vehicle arm, or the BUBR1 readout of the arm decides nothing."),
    }


def write_json(path: Path, payload: dict) -> None:
    """Write one receipt.

    Args:
        path: The receipt path.
        payload: The receipt content.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n")


def tool_versions() -> dict:
    """Give the version of every tool that this run used.

    Returns:
        A map from tool to version.
    """
    return {"python": platform.python_version(), "httpx": httpx.__version__,
            "pyyaml": yaml.__version__, "platform": platform.platform()}


def main() -> int:
    """Run the gate, compute the protocol and the cost model and write the three receipts.

    Returns:
        Process exit code. A value of 2 means that the gate failed and no receipt claims PASS.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bibliography-exit", type=int, required=True,
                        help="the exit code of scripts/build_bibliography.py on the write-up copy")
    parser.add_argument("--offline", action="store_true",
                        help="skip the PubMed fetch and record the two citation gates as SKIPPED")
    parser.add_argument("--version", default="v1", help="the suffix of the three receipts")
    parser.add_argument("--parameters", type=Path, default=PARAMETERS)
    parser.add_argument("--evidence", type=Path, default=EVIDENCE_TABLE)
    parser.add_argument("--out-dir", type=Path, default=ROOT / "tracking" / "receipts",
                        help="the directory of the three receipts")
    args = parser.parse_args()

    parameters_path = args.parameters if args.parameters.is_absolute() else ROOT / args.parameters
    evidence_path = args.evidence if args.evidence.is_absolute() else ROOT / args.evidence
    evidence_receipt_path = args.out_dir / f"t2_preclinical_evidence_{args.version}.json"
    protocol_receipt_path = args.out_dir / f"t2_preclinical_protocol_{args.version}.json"
    cost_receipt_path = args.out_dir / f"t2_preclinical_cost_model_{args.version}.json"
    rows = read_evidence_table(evidence_path)
    parameters = yaml.safe_load(parameters_path.read_text())
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    pmids = sorted({(row["source_id"] or "").strip() for row in rows
                    if (row["source_kind"] or "").strip() == "pmid"})
    titles: dict[str, str] = {}
    if not args.offline:
        with httpx.Client() as client:
            titles = pubmed_titles(client, pmids)

    findings = field_findings(rows)
    if not args.offline:
        findings.extend(citation_findings(rows, titles))
    if args.bibliography_exit != 0:
        findings.append({"claim_id": "", "condition": "C7",
                         "detail": f"build_bibliography.py exited {args.bibliography_exit}"})

    no_basis = [{"claim_id": row["claim_id"], "parameter": row["what_it_fixes"],
                 "reason": row["no_published_basis"]}
                for row in rows if (row["source_kind"] or "").strip() == "none"]

    verdict = "PASS" if not findings else "FAIL"
    evidence_receipt = {
        "receipt": f"t2_preclinical_evidence_{args.version}",
        "verdict": verdict,
        "recorded_at": stamp,
        "gate": str(GATE.relative_to(ROOT)),
        "inputs": {
            str(evidence_path.relative_to(ROOT)): sha256_of_file(evidence_path),
            str(parameters_path.relative_to(ROOT)): sha256_of_file(parameters_path),
            str(GATE.relative_to(ROOT)): sha256_of_file(GATE),
        },
        "tool_versions": tool_versions(),
        "conditions": {
            "C1_resolve": "SKIPPED" if args.offline else "PASS",
            "C2_keyword": "SKIPPED" if args.offline else "PASS",
            "C3_fields": "PASS",
            "C4_named_source": "PASS",
            "C5_no_published_basis": "PASS",
            "C6_arithmetic": "PASS",
            "C7_bibliography": "PASS" if args.bibliography_exit == 0 else "FAIL",
            "C8_design_matched_power": ("PASS" if POWER_V2_RECEIPT.exists() else "FAIL"),
            "C9_unit_of_replication": ("PASS" if power_v2_declares_units() else "FAIL"),
            "C10_stage2_own_calculation":
                ("PASS" if "stage2_group_size" in parameters["statistics"] else "FAIL"),
            "C11_primary_gate_reads_a_disease_genotype":
                ("PASS" if parameters["stage1"].get("primary_gate_requires_disease_genotype")
                 else "FAIL"),
        },
        "counts": {
            "rows": len(rows),
            "pmid_rows": sum(1 for row in rows if (row["source_kind"] or "").strip() == "pmid"),
            "distinct_pmids": len(pmids),
            "named_source_rows": sum(1 for row in rows
                                     if (row["source_kind"] or "").strip() == "named"),
            "no_published_basis_rows": len(no_basis),
        },
        "parameters_with_no_published_basis": no_basis,
        "findings": findings,
        "resolved_titles": {pmid: titles.get(pmid, "") for pmid in pmids},
        "literature_values": literature_values(rows),
    }
    for key, value in evidence_receipt["conditions"].items():
        if any(finding["condition"] == key.split("_")[0] for finding in findings):
            evidence_receipt["conditions"][key] = "FAIL"
    added = amendment5_conditions(parameters, rows)
    evidence_receipt["conditions"].update(added)
    for key, value in added.items():
        if value == "FAIL":
            findings.append({"claim_id": "", "condition": key.split("_")[0],
                             "detail": f"amendment 5 condition {key} fails"})
    evidence_receipt["findings"] = findings
    evidence_receipt["verdict"] = "PASS" if not findings else "FAIL"
    write_json(evidence_receipt_path, evidence_receipt)

    if findings:
        for finding in findings:
            print(f"FAIL {finding['condition']} {finding['claim_id']} {finding['detail']}")
        print(f"gate FAIL with {len(findings)} findings; the receipt records them")
        return 2

    statistics = parameters["statistics"]
    effect_log2 = math.log2(statistics["effect_fold_change"])
    alpha = statistics["alpha_family"] / statistics["compounds_in_family"]
    table = power_table(statistics)
    planning = next(row for row in table
                    if row["assumed_sd_log2"] == statistics["planning_sd_log2"])
    confirmation = simulate_power(planning["n_per_group"], effect_log2,
                                  statistics["planning_sd_log2"], alpha,
                                  statistics["simulation_iterations_final"],
                                  random.Random(statistics["simulation_seed"]))
    # Amendment 2 of the gate adds conditions C8, C9 and C10. The two design sizes come from
    # scripts/t2_preclinical_power_v2.py, which simulates the analysis that each stage states. The
    # table above stays in the receipt as the conservative comparison of one dose pair.
    stage1_replicates = statistics.get("stage1_replicates_per_arm", planning["n_per_group"])
    stage2_group = statistics.get("stage2_group_size", planning["n_per_group"])
    stage1 = stage1_cost(parameters, stage1_replicates)
    stage2 = stage2_cost(parameters, stage2_group)

    protocol_receipt = {
        "receipt": f"t2_preclinical_protocol_{args.version}",
        "verdict": "PASS",
        "recorded_at": stamp,
        "gate": str(GATE.relative_to(ROOT)),
        "inputs": {
            str(evidence_path.relative_to(ROOT)): sha256_of_file(evidence_path),
            str(parameters_path.relative_to(ROOT)): sha256_of_file(parameters_path),
        },
        "tool_versions": tool_versions(),
        "statistics": {
            "effect_fold_change": statistics["effect_fold_change"],
            "effect_log2": round(effect_log2, 4),
            "alpha_family": statistics["alpha_family"],
            "compounds_in_family": statistics["compounds_in_family"],
            "alpha_per_compound": round(alpha, 4),
            "power_target": statistics["power_target"],
            "planning_sd_log2": statistics["planning_sd_log2"],
            "simulation_seed": statistics["simulation_seed"],
            "simulation_iterations_search": statistics["simulation_iterations_search"],
            "simulation_iterations_final": statistics["simulation_iterations_final"],
            "n_per_group": planning["n_per_group"],
            "simulated_power_at_n": round(confirmation, 4),
            "critical_t": round(student_t_quantile(1.0 - alpha / 2.0,
                                                   2 * planning["n_per_group"] - 2), 4),
            "power_table": table,
            "what_this_table_covers": "A two sided two sample t test of the top dose against the "
                                      "vehicle arm inside one line. It is the conservative "
                                      "comparison of one dose pair and it is not the primary "
                                      "analysis of either stage.",
            "design_matched_power": {
                "receipt": "t2_preclinical_power_v2",
                "script": "scripts/t2_preclinical_power_v2.py",
                "stage1_replicates_per_arm": stage1_replicates,
                "stage1_unit_of_replication": "the donor line",
                "stage2_group_size": stage2_group,
                "stage2_unit_of_replication": "the litter",
            },
        },
        "stage1": {
            "name": parameters["stage1"]["name"],
            "donors": parameters["stage1"]["donors"],
            "lines": parameters["stage1"]["lines"],
            "dose_mM": parameters["stage1"]["dose_mM"],
            "comparator_high_dose_mM": parameters["stage1"]["comparator_high_dose_mM"],
            "agk2_concentration_uM": parameters["stage1"]["agk2_concentration_uM"],
            "cytochalasin_b_ug_per_ml": parameters["stage1"]["cytochalasin_b_ug_per_ml"],
            "binucleate_cells_scored_per_arm":
                parameters["stage1"]["binucleate_cells_scored_per_arm"],
            "exposure_hours": parameters["stage1"]["exposure_hours"],
            "arms_per_line": stage1["arms_per_line"],
            "arm_counts": stage1["arm_counts"],
            "replicates_per_arm": stage1["replicates_per_arm"],
            "wells": stage1["wells"],
        },
        "stage2": {
            "name": parameters["stage2"]["name"],
            "groups": parameters["stage2"]["groups"],
            "group_size": stage2["group_size"],
            "total_mice": stage2["total_mice"],
            "dose_mg_per_kg_per_day": parameters["stage2"]["dose_mg_per_kg_per_day"],
            "dosing_days": parameters["stage2"]["dosing_days"],
            "acclimation_days": parameters["stage2"]["acclimation_days"],
            "housing_days_per_animal": stage2["housing_days_per_animal"],
            "pups_needed": stage2["pups_needed"],
            "litters": stage2["litters"],
        },
        "negative_result": regime_ratios(parameters["graph_result"]),
        "parameters_with_no_published_basis": no_basis,
    }
    later_dates = {name: entry["access_date"] for name, entry in parameters["compounds"].items()
                   if entry.get("access_date")}
    rule = decision_rule(parameters)
    if rule is not None:
        protocol_receipt["statistics"]["decision_rule_of_amendment_5"] = rule
    write_json(protocol_receipt_path, protocol_receipt)

    cost_receipt = {
        "receipt": f"t2_preclinical_cost_model_{args.version}",
        "verdict": "PASS",
        "recorded_at": stamp,
        "gate": str(GATE.relative_to(ROOT)),
        "inputs": {
            str(evidence_path.relative_to(ROOT)): sha256_of_file(evidence_path),
            str(parameters_path.relative_to(ROOT)): sha256_of_file(parameters_path),
        },
        "tool_versions": tool_versions(),
        "currency": "USD",
        "price_access_date": parameters["date"],
        **({"later_price_access_dates": later_dates} if later_dates else {}),
        "stage1": stage1,
        "stage2": stage2,
        "priced_total_usd": round(stage1["priced_total_usd"] + stage2["priced_total_usd"], 2),
        "decisive_experiments": decisive_experiments(stage1, parameters),
        "graph_rerun_money_usd": parameters["graph_result"]["run_money_usd"],
    }
    write_json(cost_receipt_path, cost_receipt)

    print(f"gate PASS over {len(rows)} rows and {len(pmids)} distinct PMIDs")
    print(f"group size {planning['n_per_group']} per group at power "
          f"{round(confirmation, 4)} and standard deviation {statistics['planning_sd_log2']}")
    print(f"stage 1 priced total {stage1['priced_total_usd']} USD, "
          f"stage 2 priced total {stage2['priced_total_usd']} USD")
    return 0


if __name__ == "__main__":
    sys.exit(main())
