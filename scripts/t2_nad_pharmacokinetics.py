#!/usr/bin/env python3
"""Record the measured exposure of nicotinamide and of two NAD+ precursors, and test reachability.

Defect 10 of `reports/jury/jury_review_track2_2026-09-18T134500Z.md` states that the Track 2 report
labels an in vitro concentration a clinical exposure while it carries no pharmacokinetic record that
maps an oral dose onto a concentration. Rank 9 of section 9 of `reports/audit/repair_round_status.md`
names the same gap.

This script writes the record. It holds one row per measured exposure, each with the compound, the
species, the dose, the route, the measured concentration, the half life and the PMID. It then
divides the concentration that the proposed mechanism was measured at by the highest concentration
that a tolerated dose reaches, and it states whether the mechanism concentration is reachable.

The arithmetic is a division of two recorded values and it crosses no threshold that a gate fixes.
The reachability verdict is the sign of that ratio against one, which is a property of the two
measurements and not a decision of this work.

Usage:
    uv run python scripts/t2_nad_pharmacokinetics.py
    uv run python scripts/t2_nad_pharmacokinetics.py --offline
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import t2_evidence_gate as gate  # noqa: E402

# The molar mass that converts a mass concentration of the 1993 and 1995 radiotherapy studies into
# the molar concentration that every other record of this table carries. PubChem CID 936.
NICOTINAMIDE_G_PER_MOL = 122.12
NICOTINAMIDE_PUBCHEM_CID = 936

# Every number that the report states from this record is a named leaf here, so that
# `scripts/check_report_claims.py` resolves it against one field and not against a number inside a
# sentence of the table.
STATED_VALUES = {
    "nicotinamide_one_gram_oral_human": {
        "pmid": "38093631",
        "volunteers": 6,
        "dose_grams": 1,
        "maximum_plasma_concentration_micromolar": 170,
        "maximum_plasma_concentration_first_quartile": 141,
        "maximum_plasma_concentration_third_quartile": 184,
        "time_to_maximum_hours": 1,
        "terminal_half_life_hours": 2.1,
        "terminal_half_life_first_quartile": 1.8,
        "terminal_half_life_third_quartile": 2.6,
        "area_under_the_curve_micromolar_hours": 510,
        "area_under_the_curve_first_quartile": 423,
        "area_under_the_curve_third_quartile": 725,
    },
    "nicotinamide_six_grams_oral_human": {
        "pmid": "8532910",
        "patients": 12,
        "dose_grams": 6,
        "maximum_plasma_concentration_micrograms_per_millilitre": 156.4,
        "standard_error_micrograms_per_millilitre": 33.6,
        "maximum_plasma_concentration_micromolar": 1280.7,
        "dose_that_was_not_tolerated_grams": 10,
    },
    "nicotinamide_mouse_intraperitoneal": {
        "pmid": "8356223",
        "matched_dose_milligrams_per_kilogram": 171,
        "human_peak_micrograms_per_millilitre": 160,
        "human_peak_micromolar": 1310.2,
        "dose_range_low_milligrams_per_kilogram": 100,
        "dose_range_high_milligrams_per_kilogram": 1000,
    },
    "nicotinamide_tolerated_ceiling_human": {
        "pmid": "11126400",
        "adult_daily_dose_grams": 3,
    },
    "nicotinamide_riboside_steady_state_human": {
        "pmid": "29211728",
        "volunteers": 8,
        "dose_milligrams": 1000,
        "doses_per_day": 2,
        "whole_blood_riboside_baseline_micromolar": 0.023,
        "whole_blood_riboside_trough_micromolar": 0.029,
        "whole_blood_riboside_peak_micromolar": 0.05,
        "whole_blood_riboside_average_steady_state_micromolar": 0.04,
        "riboside_half_life_hours": 2.7,
        "riboside_elimination_rate_constant_per_hour": 0.26,
        "whole_blood_nad_baseline_micromolar": 27,
        "whole_blood_nad_average_steady_state_micromolar": 50,
        "whole_blood_nad_mean_rise_micromolar": 26.7,
    },
    "nicotinamide_riboside_single_dose_human": {
        "pmid": "27721479",
        "subjects": 12,
        "doses_milligrams": [100, 300, 1000],
        "peak_cell_mononucleotide_micromolar": 2,
        "peak_cell_nad_micromolar": 20,
        "nad_rise_at_100_milligrams_micromolar": 4,
        "nad_rise_at_the_higher_doses_micromolar": 6.5,
    },
    "nicotinamide_riboside_eight_weeks_human": {
        "pmid": "31278280",
        "doses_milligrams_per_day": [100, 300, 1000],
        "nad_rise_percent": [22, 51, 142],
        "weeks": 8,
    },
    "nicotinamide_mononucleotide_human": {
        "single_dose_pmid": "31685720",
        "single_doses_milligrams": [100, 250, 500],
        "men": 10,
        "hours_observed": 5,
        "repeated_dose_pmid": "36002548",
        "repeated_dose_milligrams_per_day": 1250,
        "repeated_dose_weeks": 4,
        "adults": 31,
    },
    "nicotinamide_mononucleotide_mouse": {
        "intraperitoneal_pmid": "24825348",
        "intraperitoneal_dose_milligrams_per_kilogram_per_day": 500,
        "intraperitoneal_days": 7,
        "ages_months": [3, 30],
        "oral_pmid": "28068222",
        "oral_doses_milligrams_per_kilogram_per_day": [100, 300],
        "oral_months": 12,
    },
}

ROWS = [
    {
        "claim_id": "K01",
        "subject": "mechanism_concentration",
        "compound": "nicotinamide mononucleotide",
        "species": "mouse",
        "system": "wild-type mouse embryonic fibroblasts, and Sirt2 knockout fibroblasts as the "
                  "mechanism control",
        "dose": "1.0 mM in the culture medium for 24 hours",
        "route": "culture medium",
        "measured_quantity": "BubR1 abundance by immunoblot",
        "measured_value_micromolar": 1000.0,
        "half_life_hours": "",
        "pmid": "24825348",
        "citation_keyword": "BubR1",
        "what_the_row_fixes": "The concentration at which the proposed mechanism raised the primary "
                              "endpoint of stage 1 in any system",
    },
    {
        "claim_id": "K02",
        "subject": "mechanism_concentration",
        "compound": "nicotinamide",
        "species": "human",
        "system": "HeLa cells",
        "dose": "20 mM in the culture medium for 16 hours",
        "route": "culture medium",
        "measured_quantity": "BubR1 abundance by immunoblot against tubulin",
        "measured_value_micromolar": 20000.0,
        "half_life_hours": "",
        "pmid": "24825348",
        "citation_keyword": "BubR1",
        "what_the_row_fixes": "The concentration at which nicotinamide lowered the primary endpoint "
                              "of stage 1, which is the measurement that rejects the compound",
    },
    {
        "claim_id": "K03",
        "subject": "mechanism_concentration",
        "compound": "nicotinamide",
        "species": "in vitro",
        "system": "yeast Sir2 and human SIRT1 deacetylation assay",
        "dose": "the half-maximal inhibitory concentration of the assay",
        "route": "reaction buffer",
        "measured_quantity": "noncompetitive inhibition of deacetylation",
        "measured_value_micromolar": 50.0,
        "half_life_hours": "",
        "pmid": "12297502",
        "citation_keyword": "nicotinamide",
        "what_the_row_fixes": "The upper bound of the concentration at which nicotinamide inhibits a "
                              "sirtuin, which the paper reports as below 50 micromolar",
    },
    {
        "claim_id": "K04",
        "subject": "human_exposure",
        "compound": "nicotinamide",
        "species": "human",
        "system": "6 nonpregnant healthy female volunteers",
        "dose": "1 g, single",
        "route": "oral",
        "measured_quantity": "plasma maximum concentration, Table 2, median with the first and the "
                             "third quartile",
        "measured_value_micromolar": 170.0,
        "half_life_hours": 2.1,
        "pmid": "38093631",
        "citation_keyword": "Nicotinamide",
        "what_the_row_fixes": "The one record that gives the plasma concentration, the time to the "
                              "maximum and the terminal half life of an oral nicotinamide dose in "
                              "the molar unit of the protocol",
    },
    {
        "claim_id": "K05",
        "subject": "human_exposure",
        "compound": "nicotinamide",
        "species": "human",
        "system": "12 patients under radiation therapy and hyperthermia",
        "dose": "3, 4, 5, 6 and 10 g, single, escalating",
        "route": "oral",
        "measured_quantity": "mean maximum recorded plasma concentration at 6 g, with one standard "
                             "error of the mean, converted from 156.4 micrograms per millilitre",
        "measured_value_micromolar": 1280.7,
        "half_life_hours": "",
        "pmid": "8532910",
        "citation_keyword": "Nicotinamide pharmacokinetics",
        "what_the_row_fixes": "The plasma concentration of the highest oral dose that the record "
                              "calls well tolerated, and the dose at which tolerance ends",
    },
    {
        "claim_id": "K06",
        "subject": "human_exposure",
        "compound": "nicotinamide",
        "species": "mouse",
        "system": "female CDF1 mice with a C3H mammary carcinoma",
        "dose": "100 to 1000 mg per kg, and 171 mg per kg for the matched peak",
        "route": "intraperitoneal",
        "measured_quantity": "the intraperitoneal dose whose peak plasma concentration equals the "
                             "peak of 6 g of oral nicotinamide in a person, which the study gives "
                             "as 160 micrograms per millilitre",
        "measured_value_micromolar": 1310.2,
        "half_life_hours": "",
        "pmid": "8356223",
        "citation_keyword": "Nicotinamide pharmacokinetics",
        "what_the_row_fixes": "The one record that puts a mouse dose and a human dose of "
                              "nicotinamide on the same plasma concentration",
    },
    {
        "claim_id": "K07",
        "subject": "tolerated_dose",
        "compound": "nicotinamide",
        "species": "human",
        "system": "a review of forty years of therapeutic use",
        "dose": "3 g per day",
        "route": "oral",
        "measured_quantity": "the adult daily dose above which the review calls nicotinamide a drug "
                             "with toxic potential, with reversible hepatotoxicity at very high "
                             "doses",
        "measured_value_micromolar": "",
        "half_life_hours": "",
        "pmid": "11126400",
        "citation_keyword": "high-dose nicotinamide",
        "what_the_row_fixes": "The ceiling of a chronic oral nicotinamide dose in an adult",
    },
    {
        "claim_id": "K08",
        "subject": "human_exposure",
        "compound": "nicotinamide riboside",
        "species": "human",
        "system": "8 healthy volunteers at steady state on day 9",
        "dose": "1000 mg twice daily",
        "route": "oral",
        "measured_quantity": "whole blood nicotinamide riboside at the 3 hour peak, mean with the "
                             "standard deviation, against a baseline of 0.023 micromolar",
        "measured_value_micromolar": 0.05,
        "half_life_hours": 2.7,
        "pmid": "29211728",
        "citation_keyword": "nicotinamide riboside",
        "what_the_row_fixes": "The highest measured blood concentration of nicotinamide riboside "
                              "itself at a dose that produced no adverse event",
    },
    {
        "claim_id": "K09",
        "subject": "human_exposure",
        "compound": "nicotinamide riboside",
        "species": "human",
        "system": "the same 8 healthy volunteers",
        "dose": "1000 mg twice daily",
        "route": "oral",
        "measured_quantity": "whole blood NAD+ at steady state, mean with the standard deviation, "
                             "against a baseline of 27 micromolar",
        "measured_value_micromolar": 50.0,
        "half_life_hours": "",
        "pmid": "29211728",
        "citation_keyword": "nicotinamide riboside",
        "what_the_row_fixes": "That the dose raises the pool that the mechanism needs, and by how "
                              "much",
    },
    {
        "claim_id": "K10",
        "subject": "human_exposure",
        "compound": "nicotinamide mononucleotide",
        "species": "human",
        "system": "12 healthy subjects, peripheral blood mononuclear cells",
        "dose": "100, 300 and 1000 mg of oral nicotinamide riboside, single",
        "route": "oral",
        "measured_quantity": "the averaged peak concentration of nicotinamide mononucleotide inside "
                             "the cell at the 300 mg and the 1000 mg dose",
        "measured_value_micromolar": 2.0,
        "half_life_hours": "",
        "pmid": "27721479",
        "citation_keyword": "Nicotinamide riboside",
        "what_the_row_fixes": "The highest measured intracellular concentration of the compound of "
                              "the mechanism in a person",
    },
    {
        "claim_id": "K11",
        "subject": "human_exposure",
        "compound": "nicotinamide riboside",
        "species": "human",
        "system": "overweight but otherwise healthy men and women, 8 weeks",
        "dose": "100, 300 and 1000 mg per day",
        "route": "oral",
        "measured_quantity": "the rise of whole blood NAD+ at the three doses, which is 22, 51 and "
                             "142 percent, with no difference in adverse events against placebo",
        "measured_value_micromolar": "",
        "half_life_hours": "",
        "pmid": "31278280",
        "citation_keyword": "Nicotinamide Riboside Chloride",
        "what_the_row_fixes": "The dose response of the NAD+ pool and the tolerated chronic dose of "
                              "nicotinamide riboside",
    },
    {
        "claim_id": "K12",
        "subject": "route",
        "compound": "nicotinamide riboside and nicotinamide mononucleotide",
        "species": "mouse",
        "system": "isotope tracer flux measurement across tissues",
        "dose": "the tracer dose of the study",
        "route": "intravenous against oral",
        "measured_quantity": "the form in which each compound reaches a tissue",
        "measured_value_micromolar": "",
        "half_life_hours": "",
        "pmid": "29685734",
        "citation_keyword": "NAD Synthesis-Breakdown Fluxes",
        "what_the_row_fixes": "That the intravenous route delivers the intact molecule to several "
                              "tissues and that the same compound given orally is metabolised to "
                              "nicotinamide in the liver first",
    },
    {
        "claim_id": "K13",
        "subject": "mouse_exposure",
        "compound": "nicotinamide mononucleotide",
        "species": "mouse",
        "system": "C57BL/6 mice aged 3 months and 30 months",
        "dose": "500 mg per kg per day for 7 days",
        "route": "intraperitoneal",
        "measured_quantity": "NAD+ in heart and testis, and the testis BubR1 immunoblot of the aged "
                             "animals against the young animals",
        "measured_value_micromolar": "",
        "half_life_hours": "",
        "pmid": "24825348",
        "citation_keyword": "BubR1",
        "what_the_row_fixes": "The one exposure in any animal that moved the primary endpoint of "
                              "stage 2, and the route that delivered it",
    },
    {
        "claim_id": "K14",
        "subject": "mouse_exposure",
        "compound": "nicotinamide mononucleotide",
        "species": "mouse",
        "system": "C57BL/6N mice, 12 months of administration",
        "dose": "100 and 300 mg per kg per day",
        "route": "oral, in the drinking water",
        "measured_quantity": "tissue NAD+ and a set of age-associated measurements, with no "
                             "observed toxicity",
        "measured_value_micromolar": "",
        "half_life_hours": "",
        "pmid": "28068222",
        "citation_keyword": "Nicotinamide Mononucleotide",
        "what_the_row_fixes": "The oral dose and the route that stage 2 of the protocol copies. No "
                              "measurement of this study reads BubR1",
    },
    {
        "claim_id": "K15",
        "subject": "human_exposure",
        "compound": "nicotinamide mononucleotide",
        "species": "human",
        "system": "10 healthy Japanese men",
        "dose": "100, 250 and 500 mg, single",
        "route": "oral",
        "measured_quantity": "plasma N-methyl-2-pyridone-5-carboxamide and "
                             "N-methyl-4-pyridone-5-carboxamide over 5 hours, which rose with the "
                             "dose, with no change in heart rate, blood pressure, oxygen saturation "
                             "or body temperature",
        "measured_value_micromolar": "",
        "half_life_hours": "",
        "pmid": "31685720",
        "citation_keyword": "nicotinamide mononucleotide",
        "what_the_row_fixes": "That an oral dose of the compound leaves the body as a nicotinamide "
                              "metabolite, and that the study measured no concentration of the "
                              "compound itself",
    },
    {
        "claim_id": "K16",
        "subject": "tolerated_dose",
        "compound": "nicotinamide mononucleotide",
        "species": "human",
        "system": "31 healthy adults aged 20 to 65, 4 weeks",
        "dose": "1250 mg once daily",
        "route": "oral",
        "measured_quantity": "anthropometry, haematology, biochemistry, urine and body composition, "
                             "with no severe adverse event",
        "measured_value_micromolar": "",
        "half_life_hours": "",
        "pmid": "36002548",
        "citation_keyword": "nicotinamide mononucleotide",
        "what_the_row_fixes": "The highest daily oral dose of the compound that a trial has reported "
                              "as tolerated",
    },
]


def reachability(rows: list[dict]) -> dict:
    """Divide each mechanism concentration by the highest concentration a tolerated dose reaches.

    Args:
        rows: The evidence rows.

    Returns:
        One entry per comparison, with the two concentrations, the ratio and the verdict.
    """
    by_id = {row["claim_id"]: row for row in rows}
    comparisons = [
        {
            "comparison": "nicotinamide_mononucleotide_mechanism_against_human_intracellular",
            "mechanism_row": "K01",
            "exposure_row": "K10",
            "question": "Does an oral NAD+ precursor reach, inside a human cell, the concentration "
                        "at which nicotinamide mononucleotide raised BubR1?",
        },
        {
            "comparison": "nicotinamide_mononucleotide_mechanism_against_human_blood_riboside",
            "mechanism_row": "K01",
            "exposure_row": "K08",
            "question": "Does oral nicotinamide riboside at 1000 mg twice daily reach that "
                        "concentration in blood as the riboside itself?",
        },
        {
            "comparison": "nicotinamide_rejection_against_human_plasma",
            "mechanism_row": "K02",
            "exposure_row": "K04",
            "question": "Does an oral nicotinamide dose reach the concentration at which BubR1 fell "
                        "in cells?",
        },
        {
            "comparison": "nicotinamide_sirtuin_inhibition_against_human_plasma",
            "mechanism_row": "K03",
            "exposure_row": "K04",
            "question": "Does an oral nicotinamide dose reach the concentration at which "
                        "nicotinamide inhibits a sirtuin?",
        },
    ]
    out = []
    for item in comparisons:
        needed = by_id[item["mechanism_row"]]["measured_value_micromolar"]
        reached = by_id[item["exposure_row"]]["measured_value_micromolar"]
        ratio = round(needed / reached, 4)
        out.append({**item,
                    "concentration_the_mechanism_was_measured_at_micromolar": needed,
                    "concentration_a_tolerated_dose_reaches_micromolar": reached,
                    "fold_short": ratio,
                    "reachable_at_a_tolerated_dose": ratio <= 1.0})
    ladder = {
        "protocol_comparator_steps_millimolar": [0.01, 0.0316, 0.1, 0.316, 1.0, 20.0],
        "human_maximum_plasma_concentration_millimolar": 0.17,
        "step_that_brackets_the_human_maximum": 0.1,
        "fold_of_the_top_non_pharmacological_step_above_the_human_maximum": round(1000 / 170, 4),
        "fold_of_the_rejection_concentration_above_the_human_maximum": round(20000 / 170, 4),
        "fold_of_the_human_maximum_above_the_sirtuin_inhibition_bound": round(170 / 50, 4),
    }
    return {"comparisons": out, "concentration_ladder": ladder,
            "verdict": "NOT_REACHABLE" if any(not item["reachable_at_a_tolerated_dose"]
                                              for item in out[:2]) else "REACHABLE"}


def main() -> int:
    """Write the pharmacokinetic table and its receipt.

    Returns:
        Zero when the gate passes, and 2 when a citation fails.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path,
                        default=ROOT / "literature" / "track2_nad_pharmacokinetics_v1.tsv")
    parser.add_argument("--receipt", type=Path,
                        default=ROOT / "tracking" / "receipts" /
                        "t2_nad_pharmacokinetics_v1.json")
    parser.add_argument("--offline", action="store_true",
                        help="Write the table and skip the citation gate.")
    args = parser.parse_args()

    digest = gate.write_table(args.table, ROWS)
    titles = {} if args.offline else gate.pubmed_titles([row["pmid"] for row in ROWS])
    findings = [] if args.offline else gate.citation_findings(ROWS, titles)
    result = reachability(ROWS)

    receipt = {
        "receipt_id": args.receipt.stem,
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "track": "t2_rank",
        "control": "measured_exposure_of_nicotinamide_and_of_two_nad_precursors",
        "queue_item": "R48",
        "data_class_of_run": "A",
        "data_class_of_receipt": "A",
        "receipt_scope": "Published pharmacokinetic values and their identifiers. No patient data "
                         "enters this run.",
        "script": "scripts/t2_nad_pharmacokinetics.py",
        "repository_commit": gate.repository_commit(),
        "decision_gate": "None. The run states no verdict that a threshold decides. It divides two "
                         "recorded concentrations and reports the ratio, and the reachability field "
                         "is the sign of that ratio against one.",
        "why_this_receipt_exists": "Defect 10 of reports/jury/jury_review_track2_2026-09-18T134500Z"
                                   ".md states that the report labels an unlinked in vitro "
                                   "concentration a clinical exposure, that no cited "
                                   "pharmacokinetic record maps the oral human dose onto a "
                                   "concentration, and that the report supplies no target plasma or "
                                   "tissue concentration.",
        "table": str(args.table.relative_to(ROOT)),
        "table_sha256": digest,
        "rows": len(ROWS),
        "distinct_pmids": len({row["pmid"] for row in ROWS}),
        "nicotinamide_grams_per_mole": NICOTINAMIDE_G_PER_MOL,
        "nicotinamide_pubchem_cid": NICOTINAMIDE_PUBCHEM_CID,
        "evidence": ROWS,
        "stated_values": STATED_VALUES,
        "reachability": result,
        "how_the_values_were_read": {
            "who": "An agent of the class A arm, on 2026-09-18.",
            "model_identifier": "Not named by the record. The interactive session did not "
                                "record a model identifier for each interpretation.",
            "sources_read": "The PubMed abstract of every record of the table, through "
                            "scripts/t2_pubmed_abstracts.py, and the open access full text of six "
                            "records, through scripts/t2_europepmc_fulltext.py and two fetches of "
                            "the PubMed Central page of a record that the Europe PMC service does "
                            "not serve as full text, which are PMID 38093631 and PMID 23242215.",
            "human_check": 0,
            "what_no_check_covers": "No person read a value, an abstract or a full text. The "
                                    "citation gate resolves every PMID and finds the recorded "
                                    "keyword in its title, and it does not test whether the "
                                    "publication supports the statement beside it.",
        },
        "gate": {
            "C1_resolve": "SKIPPED" if args.offline else ("PASS" if not findings else "FAIL"),
            "C2_keyword": "SKIPPED" if args.offline else ("PASS" if not findings else "FAIL"),
            "findings": findings,
            "resolved_titles": titles,
        },
        "verdict": "PASS" if not findings else "FAIL",
        "what_this_run_does_not_establish": "The gate tests that each PMID resolves and that the "
                                            "recorded keyword is in its title. It does not test "
                                            "whether the cited publication supports the sentence "
                                            "beside it, and no person and no second model checked "
                                            "that. The table also holds no measured free "
                                            "nicotinamide concentration in a tissue under either "
                                            "precursor, because no publication that this work read "
                                            "reports one.",
    }
    args.receipt.write_text(json.dumps(receipt, indent=1) + "\n")
    print(f"wrote {args.table}")
    print(f"wrote {args.receipt}")
    print(json.dumps(result, indent=1))
    print(f"gate findings: {len(findings)}")
    for finding in findings:
        print(f"  {finding}")
    return 0 if not findings else 2


if __name__ == "__main__":
    sys.exit(main())
