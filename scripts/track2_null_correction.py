"""Recompute every Track 2 control probability under the corrected null model.

The gate is section 15 of `docs/benchmarks/t2_null_model_v1.md`. Two defects of the recorded
probabilities are repaired here.

The first defect is the resolution. A comparison against B permuted graphs cannot resolve a
probability below `1 / (B + 1)`. The Poisson binomial tail of the recorded runs reaches a smaller
value, and that value is a property of the model and not a measurement of the permutation
distribution.

The second defect is the over-dispersion. The permuted recovery count is more variable than the
independence assumption of the Poisson binomial allows, so the recorded tail is smaller than the
true probability of the same count. The correction is a beta-binomial whose first two moments are
the first two moments of the permuted count, and `scripts/track2_controls.py` holds it.

The script reranks nothing and it launches no pipeline. It reads the `by_top_k` block of each
held-out receipt, or the stored case rows of a receipt written before that block existed, and it
reads the survivor count of each permutation control receipt. It writes one receipt with one row per
control, per graph, per stratum and per size of the top list, and it prints the table that the report
carries.

The data class of every input and of the receipt is A. No patient file is opened.

Usage:
    uv run python scripts/track2_null_correction.py
    uv run python scripts/track2_null_correction.py --receipt-id t2_null_correction_v1
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_controls as controls  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RECEIPT_DIR = ROOT / "tracking" / "receipts"
DEFAULT_RECEIPT_ID = "t2_null_correction_v1"
# The endpoint that the decision record fixed before any layer result existed. Every graph and every
# stratum reports it, and the mosaic variegated aneuploidy co-mention arm also reports the grid of
# Table 4.
REPORT_ENDPOINT = 10
GRID_TOP_K = (1, 3, 10, 25)

# One entry per held-out receipt. The graph, the gene layer and the scope of the null identify the
# control that the row belongs to.
HELD_OUT_SPECS = (
    {"graph": "Mosaic variegated aneuploidy", "arm": "Co-mention", "null_scope": "all_edges",
     "receipt": "t2_heldout_real_v3.json", "grid": True},
    {"graph": "Mosaic variegated aneuploidy", "arm": "Functional", "null_scope": "all_edges",
     "receipt": "t2_heldout_functional_v3.json", "grid": False},
    {"graph": "Mosaic variegated aneuploidy", "arm": "Union", "null_scope": "all_edges",
     "receipt": "t2_heldout_union_v3.json", "grid": False},
    {"graph": "Mosaic variegated aneuploidy", "arm": "Co-mention",
     "null_scope": "indication_edges", "receipt": "t2_heldout_indication_null_v2.json",
     "grid": False},
    {"graph": "Nijmegen breakage syndrome", "arm": "Co-mention", "null_scope": "all_edges",
     "receipt": "t2_transfer_heldout_nijmegen_breakage_syndrome_v2.json", "grid": False},
    {"graph": "Duchenne muscular dystrophy", "arm": "Co-mention", "null_scope": "all_edges",
     "receipt": "t2_transfer_heldout_duchenne_muscular_dystrophy_v2.json", "grid": False},
    {"graph": "Cystic fibrosis", "arm": "Co-mention", "null_scope": "all_edges",
     "receipt": "t2_transfer_heldout_cystic_fibrosis_v1.json", "grid": False},
)

# One entry per degree-preserving permutation control. The receipt path differs between the mosaic
# variegated aneuploidy run and the three transfer runs.
PERMUTATION_SPECS = (
    {"graph": "Mosaic variegated aneuploidy", "receipt": "t2_permutation_real_v2.json",
     "field": "result"},
    {"graph": "Cystic fibrosis", "receipt": "t2_transfer_cystic_fibrosis_v1.json",
     "field": "permutation_control"},
    {"graph": "Nijmegen breakage syndrome",
     "receipt": "t2_transfer_nijmegen_breakage_syndrome_v1.json", "field": "permutation_control"},
    {"graph": "Duchenne muscular dystrophy",
     "receipt": "t2_transfer_duchenne_muscular_dystrophy_v1.json", "field": "permutation_control"},
)

STRATUM_LABELS = {"full": "every eligible disease", "restricted": "drug-poor, 9 or fewer drugs",
                  "drug_rich": "drug-rich, 10 or more drugs",
                  "mesh_supplementary": "MeSH supplementary records"}


def stratum_label(name: str) -> str:
    """Return the report label of one stratum key.

    Args:
        name: The stratum key of the receipt.

    Returns:
        The label, or the key itself when no label is declared.
    """
    return STRATUM_LABELS.get(name, name)


def format_probability(value: float | None) -> str:
    """Return one probability in the form that the report prints.

    Args:
        value: The probability, or `None` when the stratum reports none.

    Returns:
        A string with three decimal places at or above 0.1 and three significant digits below it,
        or "not reported".
    """
    if value is None:
        return "not reported"
    if value >= 0.1:
        return f"{value:.3f}"
    return f"{value:.3g}"


def reported_probability(decision: float | None, floor: float | None) -> str:
    """Return the decision probability as a point value or as the bound of its own design.

    A comparison against B permuted graphs cannot resolve a probability below `1 / (B + 1)`, so a
    decision probability at that floor states a bound and not a point value.

    Args:
        decision: The decision probability of the row.
        floor: The smallest probability that the permutation count of the row reaches.

    Returns:
        The value that the report prints.
    """
    if decision is None:
        return "not reported"
    if floor is not None and decision <= floor:
        return f"{floor:.3f} or less"
    return format_probability(decision)


def held_out_row(spec: dict, record: dict, top_k: int, name: str, summary: dict) -> dict:
    """Build one row of the corrected table from one stratum summary.

    Args:
        spec: The entry of `HELD_OUT_SPECS` that names the graph and the arm.
        record: The resolution record of the receipt.
        top_k: The size of the top list of the row.
        name: The stratum key.
        summary: The resolution block of that stratum.

    Returns:
        The row.
    """
    floor = summary["smallest_reachable_permutation_probability"]
    decision = summary.get("decision_probability")
    quantile = summary["permutation_quantile_probability"]
    corrected_upper = summary.get("dispersion_corrected_probability_at_null_upper_limit")
    conservative = (None if quantile is None or corrected_upper is None
                    else max(quantile, corrected_upper))
    return {"control": "held-out indication recovery", "graph": spec["graph"],
            "arm": spec["arm"], "null_scope": spec["null_scope"], "stratum": name,
            "stratum_label": stratum_label(name), "top_k": top_k,
            "receipt_id": record["receipt_id"], "cases": summary["cases"],
            "permutations": record["permutations"], "resolution_floor": floor,
            "recovery_rate": summary["recovery_rate"],
            "null_rate_smoothed": summary["null_rate_smoothed"],
            "null_count_dispersion_factor": summary.get("null_count_dispersion_factor"),
            "dispersion_factor_applied": summary.get("dispersion_factor_applied"),
            "probability_uncorrected": summary["chance_probability"],
            "probability_uncorrected_at_null_upper_limit":
                summary["chance_probability_at_null_upper_limit"],
            "probability_corrected": summary.get("dispersion_corrected_probability"),
            "probability_corrected_at_null_upper_limit":
                summary.get("dispersion_corrected_probability_at_null_upper_limit"),
            "permutation_quantile_probability": quantile,
            "decision_probability": decision,
            "decision_probability_reported": reported_probability(decision, floor),
            "decision_probability_at_null_upper_limit": conservative,
            "decision_probability_at_null_upper_limit_reported":
                reported_probability(conservative, floor),
            "probability_resolution": summary["probability_resolution"],
            "verdict_as_recorded": summary.get("verdict_as_recorded"),
            "verdict": summary["verdict"]}


def permutation_row(spec: dict, receipt_dir: Path) -> dict:
    """Build one row of the corrected table from one permutation control receipt.

    The permutation control reports an empirical p for each candidate of the real top ten and it
    assumes nothing, so the correction of the held-out tail does not apply to it. The row states its
    own resolution floor, which is `1 / (permutations + 1)`.

    Args:
        spec: The entry of `PERMUTATION_SPECS`.
        receipt_dir: The receipt directory.

    Returns:
        The row.

    Raises:
        KeyError: When the receipt does not carry the named field.
    """
    record = json.loads((receipt_dir / spec["receipt"]).read_text(encoding="utf-8"))
    block = record[spec["field"]]
    permutations = block["permutations"]
    floor = 1.0 / (permutations + 1)
    return {"control": "degree-preserving permutation", "graph": spec["graph"],
            "arm": "Co-mention", "null_scope": "all_edges", "stratum": "top_ten_candidates",
            "stratum_label": "the 10 candidates of the real top ten", "top_k": block["top_k"],
            "receipt_id": record["receipt_id"], "cases": block["top_k"],
            "permutations": permutations, "resolution_floor": float(f"{floor:.3g}"),
            "candidates_total": block["candidates_total"], "survivors": block["survivors"],
            "decision_probability": None,
            "decision_probability_reported": f"{floor:.3f} or less per candidate",
            "probability_resolution": "resolved",
            "verdict_as_recorded": block["verdict"], "verdict": block["verdict"]}


def collect_held_out_rows(receipt_dir: Path) -> list[dict]:
    """Read every held-out receipt and return every stratum of every size of the top list.

    Each row carries `in_report_table`. The report table holds the pre-specified top ten endpoint of
    every graph and every arm, and it leaves out the MeSH supplementary stratum, which carries no
    indication edge on any graph and returns INSUFFICIENT_CASES everywhere. Every other row enters
    the census of the dispersion and of the verdicts.

    Args:
        receipt_dir: The receipt directory.

    Returns:
        One row per graph, arm, stratum and size of the top list.
    """
    rows: list[dict] = []
    for spec in HELD_OUT_SPECS:
        path = receipt_dir / spec["receipt"]
        if not path.exists():
            continue
        record = controls.receipt_resolution(path, GRID_TOP_K, controls.DEFAULT_CONFIDENCE,
                                             ("drug_rich", "mesh_supplementary"))
        for top_k in GRID_TOP_K:
            for name, summary in record["by_top_k"].get(str(top_k), {}).items():
                row = held_out_row(spec, record, top_k, name, summary)
                row["in_report_table"] = (top_k == REPORT_ENDPOINT
                                          and name != "mesh_supplementary")
                rows.append(row)
    return rows


def verdict_census(rows: list[dict]) -> dict:
    """Count the verdicts of every held-out comparison of this work.

    Args:
        rows: Every row of the corrected table.

    Returns:
        The number of comparisons of each verdict, and the rows that read PASS.
    """
    held = [row for row in rows if row["control"] == "held-out indication recovery"]
    counts: dict[str, int] = {}
    for row in held:
        counts[row["verdict"]] = counts.get(row["verdict"], 0) + 1
    passes = [{"graph": row["graph"], "arm": row["arm"], "null_scope": row["null_scope"],
               "stratum": row["stratum"], "top_k": row["top_k"],
               "recovery_rate": row["recovery_rate"], "null_rate_smoothed": row["null_rate_smoothed"],
               "decision_probability_reported": row["decision_probability_reported"]}
              for row in held if row["verdict"] == "PASS"]
    return {"comparisons": len(held), "by_verdict": counts, "pass_rows": passes}


def changed_rows(rows: list[dict]) -> list[dict]:
    """Return every row whose decision probability differs from the recorded Poisson binomial value.

    Args:
        rows: The rows of the corrected table.

    Returns:
        One entry per changed row, with the old value and the new value.
    """
    changes = []
    for row in rows:
        if row["control"] != "held-out indication recovery":
            continue
        old = row["probability_uncorrected"]
        new = row["decision_probability"]
        if old is None or new is None or old == new:
            continue
        changes.append({"graph": row["graph"], "arm": row["arm"],
                        "null_scope": row["null_scope"], "stratum": row["stratum"],
                        "top_k": row["top_k"], "receipt_id": row["receipt_id"],
                        "probability_uncorrected": old,
                        "decision_probability": new,
                        "decision_probability_reported": row["decision_probability_reported"],
                        "direction": "higher" if new > old else "lower"})
    return changes


def moved_verdicts(rows: list[dict]) -> list[dict]:
    """Return every row whose verdict differs from the verdict of its source receipt.

    Args:
        rows: The rows of the corrected table.

    Returns:
        One entry per moved verdict.
    """
    return [{"graph": row["graph"], "arm": row["arm"], "stratum": row["stratum"],
             "top_k": row["top_k"], "verdict_as_recorded": row["verdict_as_recorded"],
             "verdict": row["verdict"]}
            for row in rows
            if row.get("verdict_as_recorded") and row["verdict_as_recorded"] != row["verdict"]]


def dispersion_summary(rows: list[dict], permutations: int | None = None,
                       null_scope: str | None = None) -> dict:
    """Describe the dispersion of the permuted count over every row that measures it.

    Args:
        rows: The rows of the corrected table.
        permutations: Read only the rows of this permutation count, or every row when `None`.
        null_scope: Read only the rows of this null scope, or every scope when `None`.

    Returns:
        The number of rows, the smallest and the largest factor, the number above 1 and the number
        that the correction changes.
    """
    factors = [(row["null_count_dispersion_factor"], row) for row in rows
               if row.get("null_count_dispersion_factor") is not None
               and (permutations is None or row.get("permutations") == permutations)
               and (null_scope is None or row.get("null_scope") == null_scope)]
    if not factors:
        return {"rows": 0}
    values = [value for value, _ in factors]
    low = min(factors, key=lambda item: item[0])
    high = max(factors, key=lambda item: item[0])
    deviations = [math.log(value) / math.sqrt(2.0 / (row["permutations"] - 1))
                  for value, row in factors if value > 0 and row["permutations"] > 1]
    return {"rows": len(values), "smallest": low[0], "largest": high[0],
            "more_than_1_96_sd_above_one": sum(1 for value in deviations if value > 1.96),
            "more_than_1_96_sd_below_one": sum(1 for value in deviations if value < -1.96),
            "smallest_row": f"{low[1]['graph']}, {low[1]['arm']}, {low[1]['stratum']}, "
                            f"top {low[1]['top_k']}",
            "largest_row": f"{high[1]['graph']}, {high[1]['arm']}, {high[1]['stratum']}, "
                           f"top {high[1]['top_k']}",
            "above_one": sum(1 for value in values if value > 1.0),
            "corrected": sum(1 for _, row in factors
                             if (row.get("dispersion_factor_applied") or 0) > 1.0)}


def table_rows(rows: list[dict]) -> list[list[str]]:
    """Render the rows of the report table as strings.

    Args:
        rows: The rows of the corrected table.

    Returns:
        One list of cells per row.
    """
    rendered = []
    for row in rows:
        recovery = row.get("recovery_rate")
        null = row.get("null_rate_smoothed")
        if row["control"] == "degree-preserving permutation":
            recovery_cell = f"{row['survivors']} of {row['top_k']} survive"
            null_cell = "not applicable"
        else:
            recovery_cell = "no rate" if recovery is None else f"{recovery:.3f}"
            null_cell = "none" if null is None else f"{null:.3f}"
        rendered.append([
            row["control"], row["graph"], row["arm"],
            "indication only" if row["null_scope"] == "indication_edges" else "all edges",
            row["stratum_label"], str(row["top_k"]), str(row["cases"]),
            str(row["permutations"]),
            f"{row['resolution_floor']:.3f}" if row["resolution_floor"] else "none",
            recovery_cell, null_cell,
            row["decision_probability_reported"], row["verdict"]])
    return rendered


def markdown_table(rows: list[dict]) -> str:
    """Return the report table as Markdown.

    Args:
        rows: The rows of the corrected table.

    Returns:
        The table.
    """
    header = ["Control", "Graph", "Gene layer", "Null scope", "Stratum", "k", "Cases",
              "Permutations", "Floor", "Recovery", "Null", "p", "Verdict"]
    lines = ["| " + " | ".join(header) + " |",
             "|" + "|".join(["---"] * len(header)) + "|"]
    for cells in table_rows([row for row in rows if row.get("in_report_table", True)]):
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    """Build the command line parser.

    Returns:
        The parser.
    """
    parser = argparse.ArgumentParser(description="Correct every Track 2 control probability.")
    parser.add_argument("--receipt-dir", type=Path, default=RECEIPT_DIR,
                        help="The directory that holds the control receipts to read.")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="The directory of the written receipt. It defaults to --receipt-dir.")
    parser.add_argument("--receipt-id", default=DEFAULT_RECEIPT_ID,
                        help="The receipt identifier and file stem.")
    parser.add_argument("--table", type=Path, default=None,
                        help="Write the Markdown table to this file as well.")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Recompute every probability and write the receipt.

    Args:
        argv: The command line arguments.

    Returns:
        0 on success, 2 when no held-out receipt is present.
    """
    args = build_parser().parse_args(argv)
    rows = collect_held_out_rows(args.receipt_dir)
    if not rows:
        print("No held-out receipt was found.", file=sys.stderr)
        return 2
    rows.extend(permutation_row(spec, args.receipt_dir) for spec in PERMUTATION_SPECS)
    record = {
        "receipt_id": args.receipt_id,
        "note": ("Every probability of every Track 2 control under the corrected null model of "
                 "section 15 of docs/benchmarks/t2_null_model_v1.md. The decision probability is "
                 "the larger of the permutation quantile and the beta-binomial tail whose first "
                 "two moments are the first two moments of the permuted recovery count."),
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "track": "t2_rank", "control": "corrected_control_probabilities",
        "data_class_of_run": "A", "data_class_of_receipt": "A",
        "receipt_scope": "Public graph identifiers only. No patient data enters this run.",
        "scripts": {"analysis": "scripts/track2_null_correction.py",
                    "correction": "scripts/track2_controls.py"},
        "decision_record": "docs/benchmarks/t2_null_model_v1.md, section 15",
        "parameters": {"alpha": controls.DEFAULT_ALPHA,
                       "min_recovery": controls.DEFAULT_MIN_RECOVERY,
                       "min_cases": controls.DEFAULT_MIN_CASES,
                       "report_endpoint": REPORT_ENDPOINT, "grid_top_k": list(GRID_TOP_K),
                       "minimum_dispersion_factor": controls.MIN_DISPERSION_FACTOR,
                       "confidence": controls.DEFAULT_CONFIDENCE},
        "rows": rows,
        "changed_probabilities": changed_rows(rows),
        "moved_verdicts": moved_verdicts(rows),
        "verdicts": verdict_census(rows),
        "counts": {"held_out_comparisons": sum(
                       1 for row in rows
                       if row["control"] == "held-out indication recovery"),
                   "permutation_controls": sum(
                       1 for row in rows
                       if row["control"] == "degree-preserving permutation"),
                   "rows_of_the_report_table": sum(
                       1 for row in rows if row.get("in_report_table", True)),
                   "changed_probabilities": len(changed_rows(rows)),
                   "moved_verdicts": len(moved_verdicts(rows))},
        "dispersion_of_the_reported_rows": dispersion_summary(
            [row for row in rows if row.get("in_report_table")]),
        "dispersion_of_the_all_edge_grid_at_200_permutations": dispersion_summary(
            rows, 200, "all_edges"),
        "dispersion_of_the_indication_only_grid_at_200_permutations": dispersion_summary(
            rows, 200, "indication_edges"),
        "table_markdown": markdown_table(rows)}
    path = controls.write_receipt(record, args.receipt_id, args.out_dir or args.receipt_dir)
    print(record["table_markdown"])
    print()
    print(f"rows: {len(rows)}, changed probabilities: {len(record['changed_probabilities'])}, "
          f"moved verdicts: {len(record['moved_verdicts'])}")
    print(f"verdicts: {json.dumps(record['verdicts']['by_verdict'])} of "
          f"{record['verdicts']['comparisons']} held-out comparisons")
    print(f"dispersion of the all-edge grid at 200 permutations: "
          f"{json.dumps(record['dispersion_of_the_all_edge_grid_at_200_permutations'])}")
    print(f"dispersion of the indication-only grid at 200 permutations: "
          f"{json.dumps(record['dispersion_of_the_indication_only_grid_at_200_permutations'])}")
    if args.table:
        args.table.write_text(record["table_markdown"] + "\n", encoding="utf-8")
    print(f"receipt: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
