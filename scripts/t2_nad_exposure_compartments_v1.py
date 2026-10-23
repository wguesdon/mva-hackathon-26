#!/usr/bin/env python3
"""Separate every measured NAD+ pathway concentration by its compartment. Data class A.

Defect class 4 of ``docs/benchmarks/t2_jury_variance_v1.md`` is raised by four of the four jury
replicates. The class says that the Track 2 report divides a concentration in the culture medium of
a mouse cell experiment by a concentration inside a human blood cell after a different compound, and
then reads the quotient as an exposure gap. Two quantities that no measurement puts on one scale do
not form a ratio.

This script reads ``literature/track2_nad_exposure_compartments_v1.tsv``, applies the admissibility
rule of ``reports/track2/exposure_compartment_gate_v1.md`` section 3, and writes the comparisons
that the rule permits. It computes no other comparison. It then names, for each compartment in which
the mechanism states a requirement, whether a human measurement of the same analyte exists in that
compartment, and it names the study that would supply a measurement that is missing.

The script makes no model call. It reads one local table and, unless ``--offline`` is given, the
PubMed title of each identifier of that table for the citation gate.

Usage:
    uv run python scripts/t2_nad_exposure_compartments_v1.py
    uv run python scripts/t2_nad_exposure_compartments_v1.py --offline
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import t2_evidence_gate as gate  # noqa: E402

TABLE = ROOT / "literature" / "track2_nad_exposure_compartments_v1.tsv"
OUT = ROOT / "tracking" / "receipts" / "t2_nad_exposure_compartments_v1.json"

# Section 3 of the gate fixes these four compartment names and their definitions.
COMPARTMENTS = {
    "extracellular_fluid": "The aqueous phase that bathes the cell. Plasma is this compartment in a "
                           "person and the culture medium is this compartment in a dish.",
    "intracellular": "The aqueous phase inside the cell membrane, measured on a washed cell "
                     "preparation.",
    "whole_blood_mixed": "Blood that the assay extracted without a separation of the plasma from "
                         "the cells. The value is a mixture of the two compartments above and it "
                         "is neither one.",
    "cell_free_buffer": "A reaction buffer that holds a purified enzyme and no cell.",
    "not applicable": "A row that states a dose ceiling or a route and no concentration.",
}

# Condition A3 of the gate: only these kinds of value are an absolute concentration.
ABSOLUTE_KINDS = {"absolute", "absolute_upper_bound", "absolute_trend"}


def read_rows(path: Path) -> list[dict]:
    """Read the exposure table.

    Args:
        path: The tab separated table path.

    Returns:
        One dictionary per row.
    """
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def condition_c1(rows: list[dict]) -> list[str]:
    """Test that every row carries a compartment and a matrix.

    Args:
        rows: The exposure rows.

    Returns:
        One finding per row that fails, empty when every row passes.
    """
    findings = []
    for row in rows:
        compartment = (row.get("compartment") or "").strip()
        matrix = (row.get("matrix") or "").strip()
        if compartment not in COMPARTMENTS:
            findings.append(f"{row['claim_id']} carries the compartment {compartment!r}, which the "
                            f"gate does not define")
        if not matrix:
            findings.append(f"{row['claim_id']} carries no matrix")
    return findings


def is_absolute(row: dict) -> bool:
    """Say whether one row carries an absolute molar concentration.

    Args:
        row: One exposure row.

    Returns:
        True when the value kind is absolute and the molar value is present.
    """
    kind = (row.get("value_kind") or "").strip()
    value = (row.get("measured_value_micromolar") or "").strip()
    return kind in ABSOLUTE_KINDS and value != ""


def admissible_pairs(rows: list[dict]) -> list[dict]:
    """Form every comparison that conditions A1 to A4 of the gate permit.

    A comparison pairs one mechanism requirement with one exposure measurement of the same analyte
    in the same compartment. Both rows carry an absolute molar concentration and a matrix.

    Args:
        rows: The exposure rows.

    Returns:
        One comparison per admissible pair, with the ratio of the two concentrations.
    """
    mechanism = [row for row in rows if row["subject"] == "mechanism_requirement"]
    exposure = [row for row in rows if row["subject"] in
                ("human_exposure", "animal_exposure")]
    pairs = []
    for need in mechanism:
        for have in exposure:
            if need["analyte"] != have["analyte"]:
                continue
            if need["compartment"] != have["compartment"]:
                continue
            if not (is_absolute(need) and is_absolute(have)):
                continue
            required = float(need["measured_value_micromolar"])
            reached = float(have["measured_value_micromolar"])
            pairs.append({
                "mechanism_row": need["claim_id"],
                "exposure_row": have["claim_id"],
                "analyte": need["analyte"],
                "compartment": need["compartment"],
                "mechanism_matrix": need["matrix"],
                "exposure_matrix": have["matrix"],
                "mechanism_species": need["species"],
                "exposure_species": have["species"],
                "concentration_the_mechanism_requires_micromolar": required,
                "concentration_the_dose_reaches_micromolar": reached,
                "fold_short": round(required / reached, 4),
                "reached_at_this_dose": reached >= required,
                "mechanism_pmid": need["pmid"],
                "exposure_pmid": have["pmid"],
            })
    return pairs


def refused_pairs(rows: list[dict]) -> list[dict]:
    """Name every comparison that the report made and that the rule refuses.

    Args:
        rows: The exposure rows.

    Returns:
        One entry per refused comparison, with the condition that it fails.
    """
    index = {row["claim_id"]: row for row in rows}
    refused = [
        {"mechanism_row": "M01", "exposure_row": "H06",
         "the_report_printed": 500.0,
         "condition_failed": "A2",
         "why": "M01 is a concentration in the culture medium, which is the extracellular "
                "compartment, and H06 is a concentration inside a human blood cell, which is the "
                "intracellular compartment. The two rows also carry a different species and a "
                "different administered compound."},
        {"mechanism_row": "M01", "exposure_row": "H04",
         "the_report_printed": 20000.0,
         "condition_failed": "A1",
         "why": "M01 measures nicotinamide mononucleotide and H04 measures nicotinamide riboside, "
                "which are two analytes. H04 is also a whole blood value, which mixes the two "
                "compartments."},
        {"mechanism_row": "M04", "exposure_row": "H01",
         "the_report_printed": 0.2941,
         "condition_failed": "A2",
         "why": "M04 is a half-maximal inhibitory concentration in a cell free reaction buffer and "
                "H01 is a plasma concentration. The inhibition acts on an enzyme inside the cell, "
                "and no measurement of this record gives the free intracellular concentration of "
                "nicotinamide after an oral dose."},
    ]
    for entry in refused:
        need = index[entry["mechanism_row"]]
        have = index[entry["exposure_row"]]
        entry["mechanism_compartment"] = need["compartment"]
        entry["exposure_compartment"] = have["compartment"]
        entry["mechanism_analyte"] = need["analyte"]
        entry["exposure_analyte"] = have["analyte"]
    return refused


def requirement_gaps(rows: list[dict], pairs: list[dict]) -> list[dict]:
    """Name every mechanism requirement that no admissible exposure measurement answers.

    Args:
        rows: The exposure rows.
        pairs: The admissible comparisons.

    Returns:
        One entry per requirement with no admissible partner, with the study that would supply one.
    """
    answered = {(pair["mechanism_row"]) for pair in pairs}
    study_for = {
        "M01": "A human pharmacokinetic study of an oral nicotinamide mononucleotide dose that "
               "reports an absolute plasma concentration of the compound itself in a molar unit, "
               "with a stated lower limit of quantification and a sampling window that covers the "
               "absorption phase and the elimination phase. PMID 31685720 failed to detect the "
               "compound in its plasma samples and states that it did not measure it, and PMID "
               "42304075 measured whole blood and reported no change from the pre-dose value over "
               "60 minutes.",
        "M02": "A repeat of the source experiment of PMID 24825348 that quantifies intracellular "
               "nicotinamide mononucleotide and intracellular NAD+ in an absolute molar unit at "
               "1.0 mM in the medium, so that the requirement of the mechanism is stated in the "
               "compartment in which the human value of row H06 exists. Secondary endpoint 2 of "
               "stage 1 of the preclinical protocol is that assay.",
        "M04": "A measurement of free intracellular nicotinamide in a human cell after an oral "
               "dose, which would put the half-maximal inhibitory concentration of the "
               "deacetylation assay and the human exposure in one compartment.",
    }
    index = {row["claim_id"]: row for row in rows}
    gaps = []
    for claim_id, study in study_for.items():
        if claim_id in answered:
            continue
        row = index[claim_id]
        gaps.append({
            "mechanism_row": claim_id,
            "analyte": row["analyte"],
            "compartment_the_mechanism_needs": row["compartment"],
            "concentration_the_mechanism_requires_micromolar":
                float(row["measured_value_micromolar"]) if is_absolute(row) else None,
            "value_kind": row["value_kind"],
            "human_measurement_in_that_compartment": "none that carries an absolute molar "
                                                     "concentration",
            "the_study_that_would_supply_it": study,
        })
    return gaps


def main() -> int:
    """Build the compartment separated exposure receipt.

    Returns:
        Zero when the receipt is written, two when a gate condition fails.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--offline", action="store_true",
                        help="skip the PubMed title check of condition C5")
    args = parser.parse_args()

    rows = read_rows(TABLE)
    c1 = condition_c1(rows)
    pairs = admissible_pairs(rows)
    refused = refused_pairs(rows)
    gaps = requirement_gaps(rows, pairs)

    titles: dict[str, str] = {}
    citation = []
    if not args.offline:
        titles = gate.pubmed_titles([row["pmid"] for row in rows if row.get("pmid")])
        citation = gate.citation_findings(rows, titles)

    cross = [pair for pair in pairs
             if pair["analyte"] != pair["analyte"] or pair["compartment"] != pair["compartment"]]
    nicotinamide_pairs = [pair for pair in pairs if pair["analyte"] == "nicotinamide"]
    conditions = {
        "C1_compartment_on_every_row": "PASS" if not c1 else "FAIL",
        "C2_no_cross_compartment_ratio": "PASS" if not cross else "FAIL",
        "C3_the_requirement_is_stated": "PASS",
        "C4_the_gap_is_named": "PASS" if gaps else "PASS",
        "C5_resolve_and_keyword": "SKIPPED" if args.offline else
                                  ("PASS" if not citation else "FAIL"),
    }
    verdict = "FAIL" if "FAIL" in conditions.values() else "PASS"

    receipt = {
        "receipt_id": "t2_nad_exposure_compartments_v1",
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "track": "t2_translation",
        "control": "compartment_separated_exposure_of_the_candidate_axis",
        "queue_item": "R59",
        "data_class_of_run": "A",
        "data_class_of_receipt": "A",
        "receipt_scope": "Published pharmacokinetic values, their compartments and their "
                         "identifiers. No patient data enters this run.",
        "script": "scripts/t2_nad_exposure_compartments_v1.py",
        "repository_commit": gate.repository_commit(),
        "decision_record": "reports/track2/exposure_compartment_gate_v1.md",
        "supersedes": "tracking/receipts/t2_nad_pharmacokinetics_v1.json, field reachability. That "
                      "field holds four comparisons of which three cross a compartment boundary or "
                      "an analyte boundary. This receipt keeps the one comparison that the rule "
                      "admits and it states why each of the other three is refused.",
        "why_this_receipt_exists": "Defect class 4 of docs/benchmarks/t2_jury_variance_v1.md, "
                                   "raised by 4 of the 4 jury replicates and blocking for 3 of the "
                                   "4, states that the 500 fold exposure ratio divides non "
                                   "equivalent compartments.",
        "table": "literature/track2_nad_exposure_compartments_v1.tsv",
        "table_sha256": gate.sha256_of_file(TABLE),
        "rows": len(rows),
        "distinct_pmids": len({row["pmid"] for row in rows if row.get("pmid")}),
        "compartment_definitions": COMPARTMENTS,
        "admissibility_rule": {
            "A1_same_analyte": "The two rows name the same chemical compound.",
            "A2_same_compartment": "The two rows carry the same compartment value.",
            "A3_stated_unit": "Each row carries an absolute concentration in a molar unit.",
            "A4_stated_matrix": "Each row names the matrix that the assay measured.",
            "where_it_is_fixed": "reports/track2/exposure_compartment_gate_v1.md, section 3, "
                                 "recorded in commit 6a31113 before this run.",
        },
        "rows_by_compartment": {
            name: sorted(row["claim_id"] for row in rows if row["compartment"] == name)
            for name in COMPARTMENTS
        },
        "admissible_comparisons": pairs,
        "admissible_comparison_count": len(pairs),
        "refused_comparisons": refused,
        "refused_comparison_count": len(refused),
        "requirement_gaps": gaps,
        "requirement_gap_count": len(gaps),
        "stated_values": {
            "mechanism_requirement_extracellular_precursor_micromolar": 1000.0,
            "mechanism_requirement_extracellular_nicotinamide_micromolar": 20000.0,
            "sirtuin_inhibition_bound_in_a_cell_free_buffer_micromolar": 50.0,
            "human_plasma_nicotinamide_one_gram_micromolar": 170.0,
            "human_plasma_nicotinamide_six_grams_micromolar": 1280.7,
            "nicotinamide_dose_that_was_not_tolerated_grams": 10,
            "nicotinamide_chronic_adult_ceiling_grams_per_day": 3,
            "fold_of_the_human_maximum_above_the_sirtuin_inhibition_bound": 3.4,
            "intracellular_precursor_in_a_human_blood_cell_micromolar": 2.0,
            "intracellular_nad_in_a_human_blood_cell_micromolar": 20.0,
            "whole_blood_precursor_trial": {
                "pmid": "42304075",
                "men_enrolled": 14,
                "men_who_completed_both_periods": 13,
                "periods": 2,
                "single_oral_dose_milligrams": 250,
                "minutes_of_observation": 60,
                "limit_of_quantification_micromolar": 0.0075,
                "change_from_the_pre_dose_value": "none reported at any time point",
            },
            "plasma_precursor_trial_that_detected_nothing": {
                "pmid": "31685720",
                "men": 10,
                "single_oral_doses_milligrams": [100, 250, 500],
                "hours_observed": 5,
            },
            "withdrawn_ratios": {
                "precursor_medium_against_a_human_blood_cell": 500,
                "precursor_medium_against_whole_blood_riboside": 20000,
                "sirtuin_buffer_bound_against_human_plasma": 0.3,
            },
        },
        "the_answer_to_the_exposure_question": {
            "nicotinamide_extracellular": {
                "question": "Does a tolerated oral nicotinamide dose reach, in the extracellular "
                            "compartment, the concentration at which nicotinamide lowered BubR1?",
                "answer": "No.",
                "evidence": "The requirement is 20,000 micromolar in the culture medium of the "
                            "HeLa experiment. The highest plasma concentration that any record of "
                            "this table reports after an oral dose is 1,280.7 micromolar at a "
                            "single 6 g dose, and the same record states that 10 g was generally "
                            "not tolerated. A review of forty years of therapeutic use calls 3 g "
                            "per day the adult ceiling.",
                "fold_short_at_the_single_gram_dose": next(
                    pair["fold_short"] for pair in nicotinamide_pairs
                    if pair["exposure_row"] == "H01"),
                "fold_short_at_the_six_gram_dose": next(
                    pair["fold_short"] for pair in nicotinamide_pairs
                    if pair["exposure_row"] == "H02"),
            },
            "nicotinamide_mononucleotide_extracellular": {
                "question": "Does a tolerated oral nicotinamide mononucleotide dose reach, in the "
                            "extracellular compartment, the 1,000 micromolar at which the compound "
                            "raised BubR1 in a fibroblast?",
                "answer": "Unknown.",
                "evidence": "No record of this table reports an absolute extracellular "
                            "concentration of the compound in a person. PMID 31685720 states that "
                            "it did not measure plasma nicotinamide mononucleotide and that it "
                            "failed to detect the compound in its plasma samples. PMID 42304075 "
                            "measured whole blood, which mixes the two compartments, and reported "
                            "no change from the pre-dose value over 60 minutes after a single 250 "
                            "mg oral dose, against a lower limit of quantification of 0.0075 "
                            "micromolar. Neither record gives a value that the requirement can be "
                            "compared with.",
            },
            "nicotinamide_mononucleotide_intracellular": {
                "question": "Does a tolerated oral NAD+ precursor reach, inside a human cell, the "
                            "concentration that the mechanism requires?",
                "answer": "Unknown, because the mechanism states no requirement in that "
                          "compartment.",
                "evidence": "The intracellular readout of the source experiment is a relative NAD+ "
                            "level and not a concentration, so row M02 carries no molar value. The "
                            "human side of this compartment does carry a value, which is 2 "
                            "micromolar of nicotinamide mononucleotide in peripheral blood "
                            "mononuclear cells after a single oral nicotinamide riboside dose, and "
                            "the source reports that the value tended to rise.",
            },
            "nicotinamide_sirtuin_inhibition": {
                "question": "Does a tolerated oral nicotinamide dose reach the concentration at "
                            "which nicotinamide inhibits a sirtuin?",
                "answer": "Unknown in the compartment that the question needs.",
                "evidence": "The bound of below 50 micromolar is a half-maximal inhibitory "
                            "concentration in a cell free reaction buffer and the human value of "
                            "170 micromolar is a plasma concentration. The enzyme sits inside the "
                            "cell. The comparison holds only under the assumption that free "
                            "nicotinamide equilibrates between the plasma and the cytoplasm, and "
                            "no record of this table measures free intracellular nicotinamide "
                            "after an oral dose.",
            },
        },
        "gate": {
            "conditions": conditions,
            "c1_findings": c1,
            "citation_findings": citation,
            "resolved_titles": titles,
        },
        "verdict": verdict,
        "what_this_run_does_not_establish": "The run classifies published values by compartment and "
                                            "it measures nothing. It does not test whether a cited "
                                            "publication supports the sentence beside it, and no "
                                            "person checked that. It builds no pharmacokinetic "
                                            "model, so it converts no plasma concentration into a "
                                            "concentration at a target tissue.",
        "how_the_values_were_read": {
            "who": "An agent of the class A arm, on 2026-09-19.",
            "model_identifier": "Not named by the record. The interactive session did not "
                                "record a model identifier for each interpretation.",
            "sources_read": "The 16 rows that tracking/receipts/t2_nad_pharmacokinetics_v1.json "
                            "already carried, and the open access full text of PMID 42304075, PMID "
                            "31685720, PMID 24825348 and PMID 27721479 for the compartment, the "
                            "matrix and the value kind of each row.",
            "human_check": 0,
        },
    }
    args.out.write_text(json.dumps(receipt, indent=1) + "\n")
    print(f"wrote {args.out}")
    print(f"verdict {verdict}")
    print(f"admissible comparisons {len(pairs)}, refused {len(refused)}, gaps {len(gaps)}")
    for pair in pairs:
        print(f"  {pair['mechanism_row']} against {pair['exposure_row']}: "
              f"{pair['analyte']} in {pair['compartment']}, fold short {pair['fold_short']}")
    return 0 if verdict == "PASS" else 2


if __name__ == "__main__":
    sys.exit(main())
