"""Fit the cost model of a Track 2 run and write the coefficients that the runner reads.

Section 6 of `docs/decisions/0009_track2_one_entry_point_and_a_prospective_disease.md` fixes the
form of the model before the fit. Two relations are fitted with ordinary least squares.

The first relation predicts the held-out control, which is the stage that dominates the cost of a
large graph. The response is `log10` of the seconds of one ranking inside the control and the
predictor is `log10` of the number of edges of the graph that the control ranks. One ranking is the
unit of work of the control, because the control ranks the graph once for each disease that carries
an indication edge and once more for each of those diseases on each permutation of the null.

The second relation predicts the network half of the cost. The response is the number of external
calls of a run and the predictor is the number of nodes of its graph, fitted through the origin,
because a run with no node sends no call.

The module reads the measured runs from the receipts and the run logs alone. It invents no value.

Usage:
    uv run python scripts/track2_cost_model_fit.py --version 1 \\
        --run cystic_fibrosis --run nijmegen_breakage_syndrome \\
        --run duchenne_muscular_dystrophy
"""

from __future__ import annotations

import argparse
import datetime as dt
import math
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402

RECEIPT_DIR = common.ROOT / "tracking" / "receipts"
TRANSFER_DIR = common.ROOT / "reports" / "track2" / "transfer"
MODEL_DIR = common.ROOT / "literature"
SECONDS_PER_MINUTE = 60.0

# The run that is not a transfer run. Its cost point has its own receipt, because the mosaic
# variegated aneuploidy held-out control of the report ran at a different permutation count and its
# statistical fields are not the fields of a cost measurement.
MVA_COST_POINT = "t2_cost_point_mva_heldout_v1"


def ols(xs: list[float], ys: list[float]) -> dict:
    """Fit an ordinary least squares line.

    Args:
        xs: The predictor values.
        ys: The response values.

    Returns:
        The slope, the intercept, the coefficient of determination and the number of points.

    Raises:
        ValueError: When fewer than two points are given, or when every predictor is equal.
    """
    if len(xs) < 2 or len(xs) != len(ys):
        raise ValueError("an ordinary least squares fit needs at least two paired points.")
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    sxx = sum((x - mean_x) ** 2 for x in xs)
    if sxx == 0.0:
        raise ValueError("every predictor value is equal, so the slope is not defined.")
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = mean_y - slope * mean_x
    fitted = [intercept + slope * x for x in xs]
    ss_residual = sum((y - f) ** 2 for y, f in zip(ys, fitted))
    ss_total = sum((y - mean_y) ** 2 for y in ys)
    return {
        "slope": round(slope, 4),
        "intercept": round(intercept, 4),
        "r_squared": round(1.0 - ss_residual / ss_total, 4) if ss_total else None,
        "points": len(xs),
    }


def origin_slope(xs: list[float], ys: list[float]) -> dict:
    """Fit a line through the origin.

    Args:
        xs: The predictor values.
        ys: The response values.

    Returns:
        The slope, the coefficient of determination against the mean of the response, and the
        number of points.

    Raises:
        ValueError: When no point is given, or when every predictor is zero.
    """
    if not xs or len(xs) != len(ys):
        raise ValueError("a fit through the origin needs at least one paired point.")
    denominator = sum(x * x for x in xs)
    if denominator == 0.0:
        raise ValueError("every predictor value is zero, so the slope is not defined.")
    slope = sum(x * y for x, y in zip(xs, ys)) / denominator
    mean_y = sum(ys) / len(ys)
    ss_residual = sum((y - slope * x) ** 2 for x, y in zip(xs, ys))
    ss_total = sum((y - mean_y) ** 2 for y in ys)
    return {
        "slope": round(slope, 4),
        "r_squared": round(1.0 - ss_residual / ss_total, 4) if ss_total else None,
        "points": len(xs),
    }


def network_point(slug: str) -> dict:
    """Read the network counters and the whole run cost of one disease.

    Args:
        slug: The disease slug.

    Returns:
        The label, the call count, the byte count, the run wall clock and the peak memory, or an
        empty mapping when neither run receipt exists.
    """
    for name in (f"t2_transfer_{slug}_v1.json", f"t2_run_disease_{slug}_v1.json"):
        receipt = common.read_json(RECEIPT_DIR / name)
        if not receipt:
            continue
        cost = receipt["cost"]
        log = common.read_json(TRANSFER_DIR / slug / "run_log.json") or {}
        label = receipt["disease"].get("label", slug)
        return {
            "label": label,
            "api_calls_total": cost["api_calls_total"],
            "api_bytes_total": cost["api_bytes_total"],
            "run_wall_seconds": cost["wall_seconds_total"],
            "peak_memory_mib": cost.get("peak_memory_mib") or peak_memory_of(log.get("stages", {})),
            "source": f"tracking/receipts/{name}",
        }
    return {}


def control_point(slug: str) -> dict:
    """Read the held-out control measurement of one disease.

    A same day cost measurement is preferred over the control of the original run, because every
    such measurement ran on one machine under one load. The record names which one it used.

    Args:
        slug: The disease slug.

    Returns:
        The graph size, the control work and the wall clock of the control.

    Raises:
        FileNotFoundError: When no measurement of the control exists.
    """
    cost_receipt = common.read_json(RECEIPT_DIR / f"t2_cost_point_{slug}_heldout_v1.json")
    if cost_receipt:
        measured = cost_receipt["measurement"]
        return {
            "graph_edges": measured["graph_edges"], "graph_nodes": measured["graph_nodes"],
            "ranked_edges": measured["ranked_edges"], "diseases": measured["diseases"],
            "permutations": measured["permutations"],
            "heldout_wall_seconds": measured["heldout_wall_seconds"],
            "heldout_peak_memory_mib": measured["peak_memory_mib"],
            "measured_on": cost_receipt["created_utc"][:10],
            "measurement": "a same day cost measurement",
            "source": f"tracking/receipts/t2_cost_point_{slug}_heldout_v1.json",
        }
    heldout = common.read_json(RECEIPT_DIR / f"t2_transfer_heldout_{slug}_v1.json")
    log = common.read_json(TRANSFER_DIR / slug / "run_log.json") or {}
    stage = log.get("stages", {}).get("heldout", {})
    if not heldout or not stage:
        raise FileNotFoundError(
            f"no held-out measurement for {slug}. Run `uv run python scripts/track2_cost_point.py "
            f"--graph-dir reports/track2/transfer/{slug}/graph --slug {slug}` first.")
    manifest = common.read_json(TRANSFER_DIR / slug / "graph" / "manifest.json") or {}
    return {
        "graph_edges": manifest.get("edge_count"), "graph_nodes": manifest.get("node_count"),
        "ranked_edges": heldout["graph"]["edges"],
        "diseases": heldout["result"]["full"]["cases"],
        "permutations": heldout["result"]["full"]["permutations"],
        "heldout_wall_seconds": stage["wall_seconds"],
        "heldout_peak_memory_mib": stage.get("peak_memory_mib"),
        "measured_on": heldout["created_utc"][:10],
        "measurement": "the control stage of the original run",
        "source": f"tracking/receipts/t2_transfer_heldout_{slug}_v1.json",
    }


def disease_point(slug: str) -> dict:
    """Assemble the cost point of one disease.

    Args:
        slug: The disease slug.

    Returns:
        The measured point.

    Raises:
        FileNotFoundError: When no held-out measurement of the disease exists.
    """
    control = control_point(slug)
    network = network_point(slug)
    point = build_point(slug, network.get("label", slug), control["graph_edges"],
                        control["graph_nodes"], control["diseases"], control["permutations"],
                        control["heldout_wall_seconds"], network.get("api_calls_total"),
                        network.get("api_bytes_total"), network.get("run_wall_seconds"),
                        control["source"])
    point["ranked_edges"] = control["ranked_edges"]
    point["heldout_peak_memory_mib"] = control["heldout_peak_memory_mib"]
    point["run_peak_memory_mib"] = network.get("peak_memory_mib")
    point["measured_on"] = control["measured_on"]
    point["measurement"] = control["measurement"]
    point["network_source"] = network.get("source")
    return point


def peak_memory_of(stages: dict) -> float | None:
    """Return the largest stage peak memory of a run.

    Args:
        stages: The stage records of a run log.

    Returns:
        The peak in mebibytes, or None when no stage recorded one.
    """
    peaks = [record["peak_memory_mib"] for record in stages.values()
             if record.get("peak_memory_mib") is not None]
    return max(peaks) if peaks else None


def mva_point() -> dict:
    """Read the cost point of the mosaic variegated aneuploidy run.

    The mosaic variegated aneuploidy run predates the run log, so its network counters are not
    recorded and only its held-out control carries a measured wall clock.

    Returns:
        The measured point.

    Raises:
        FileNotFoundError: When the cost point receipt is missing.
    """
    path = RECEIPT_DIR / f"{MVA_COST_POINT}.json"
    receipt = common.read_json(path)
    if not receipt:
        raise FileNotFoundError(
            f"{path} is missing. Run `uv run python scripts/track2_cost_point.py` first.")
    measured = receipt["measurement"]
    point = build_point("mosaic_variegated_aneuploidy", "Mosaic variegated aneuploidy",
                        measured["graph_edges"], measured["graph_nodes"], measured["diseases"],
                        measured["permutations"], measured["heldout_wall_seconds"],
                        None, None, None, str(path.relative_to(common.ROOT)))
    point["ranked_edges"] = measured.get("ranked_edges")
    point["heldout_peak_memory_mib"] = measured.get("peak_memory_mib")
    point["run_peak_memory_mib"] = None
    point["measured_on"] = receipt["created_utc"][:10]
    point["measurement"] = "a same day cost measurement"
    point["network_source"] = None
    return point


def build_point(slug: str, label: str, edges: int, nodes: int, diseases: int, permutations: int,
                heldout_seconds: float, api_calls: int | None, api_bytes: int | None,
                run_seconds: float | None, source: str) -> dict:
    """Assemble one measured cost point.

    Args:
        slug: The disease slug.
        label: The disease name.
        edges: The number of edges of the ranking graph.
        nodes: The number of nodes of the ranking graph.
        diseases: The number of diseases that carry an indication edge.
        permutations: The permutation count of the held-out control.
        heldout_seconds: The wall clock of the held-out control.
        api_calls: The number of external calls of the whole run, when the run recorded one.
        api_bytes: The number of response bytes of the whole run, when the run recorded one.
        run_seconds: The wall clock of the whole run, when the run recorded one.
        source: The record that holds the values.

    Returns:
        The point.
    """
    rankings = diseases * (permutations + 1)
    return {
        "slug": slug,
        "label": label,
        "graph_edges": edges,
        "graph_nodes": nodes,
        "diseases_with_an_indication": diseases,
        "heldout_permutations": permutations,
        "heldout_rankings": rankings,
        "heldout_wall_seconds": round(heldout_seconds, 1),
        "heldout_wall_minutes": round(heldout_seconds / SECONDS_PER_MINUTE, 1),
        "seconds_per_ranking": round(heldout_seconds / rankings, 4) if rankings else None,
        "api_calls_total": api_calls,
        "api_bytes_total": api_bytes,
        "run_wall_minutes": round(run_seconds / SECONDS_PER_MINUTE, 1) if run_seconds else None,
        "source": source,
    }


def fit_points(points: list[dict]) -> dict:
    """Fit the two relations of the cost model.

    Args:
        points: The measured cost points.

    Returns:
        The two fits, with the residual of each point and the range of validity.
    """
    xs = [math.log10(point["graph_edges"]) for point in points]
    ys = [math.log10(point["seconds_per_ranking"]) for point in points]
    ranking = ols(xs, ys)
    residuals = []
    for point, x, y in zip(points, xs, ys):
        predicted = ranking["intercept"] + ranking["slope"] * x
        residuals.append({
            "slug": point["slug"],
            "observed_seconds_per_ranking": point["seconds_per_ranking"],
            "predicted_seconds_per_ranking": round(10.0 ** predicted, 4),
            "ratio_predicted_over_observed": round(10.0 ** (predicted - y), 3),
        })
    ranking["residuals"] = residuals
    ranking["largest_absolute_ratio"] = round(
        max(max(row["ratio_predicted_over_observed"], 1.0 / row["ratio_predicted_over_observed"])
            for row in residuals), 3)
    ranking["valid_edges_from"] = min(point["graph_edges"] for point in points)
    ranking["valid_edges_to"] = max(point["graph_edges"] for point in points)

    with_calls = [point for point in points if point.get("api_calls_total")]
    calls = origin_slope([float(point["graph_nodes"]) for point in with_calls],
                         [float(point["api_calls_total"]) for point in with_calls])
    calls["runs"] = [point["slug"] for point in with_calls]
    return {"ranking_seconds": ranking, "calls_per_node": calls}


def predict_seconds_per_ranking(model: dict, edges: int) -> float:
    """Return the predicted seconds of one ranking for a graph of a given size.

    Args:
        model: The `ranking_seconds` block of the cost model.
        edges: The number of edges of the graph.

    Returns:
        The predicted seconds.
    """
    return 10.0 ** (model["intercept"] + model["slope"] * math.log10(max(edges, 1)))


def write_model(version: int, fits: dict, points: list[dict]) -> Path:
    """Write the pinned coefficient file that the runner reads.

    Args:
        version: The model version.
        fits: The fitted relations.
        points: The measured points behind them.

    Returns:
        The written path.
    """
    path = MODEL_DIR / f"track2_cost_model_v{version}.yaml"
    payload = {
        "version": version,
        "created": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "decision_record": ("docs/decisions/"
                            "0009_track2_one_entry_point_and_a_prospective_disease.md"),
        "fitted_on": [point["slug"] for point in points],
        "ranking_seconds": fits["ranking_seconds"],
        "calls_per_node": fits["calls_per_node"],
        "points": points,
    }
    header = ("# The Track 2 cost model. scripts/track2_cost_model_fit.py writes this file and\n"
              "# scripts/track2_run_disease.py reads it to project the cost of a run before the\n"
              "# run spends the compute. Every value comes from a measured run.\n")
    path.write_text(header + yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return path


def main() -> int:
    """Fit the cost model and write the coefficient file and the receipt.

    Returns:
        0 on success, 1 when a measured run is missing.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--version", type=int, required=True, help="The model version to write.")
    parser.add_argument("--run", action="append", default=[],
                        help="A transfer run slug. Repeat the option for several runs.")
    parser.add_argument("--with-mva", action="store_true",
                        help="Add the mosaic variegated aneuploidy cost point.")
    args = parser.parse_args()

    points: list[dict] = []
    try:
        if args.with_mva:
            points.append(mva_point())
        for slug in args.run:
            points.append(disease_point(slug))
    except FileNotFoundError as error:
        print(f"FAIL: {error}")
        return 1
    points.sort(key=lambda point: point["graph_edges"])

    fits = fit_points(points)
    model_path = write_model(args.version, fits, points)
    receipt = {
        "receipt_id": f"t2_cost_model_v{args.version}",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "track": "t2_rank",
        "control": "cost_model_of_a_run_for_a_new_disease",
        "queue_item": "R17",
        "decision_record": ("docs/decisions/"
                            "0009_track2_one_entry_point_and_a_prospective_disease.md"),
        "data_class_of_run": "A",
        "data_class_of_receipt": "A",
        "receipt_scope": ("Graph counts, wall clock times and call counters of public graph runs "
                          "only. No patient data enters this fit."),
        "script": "scripts/track2_cost_model_fit.py",
        "model_file": str(model_path.relative_to(common.ROOT)),
        "model_file_sha256": common.sha256_of(model_path),
        "points": points,
        "ranking_seconds": fits["ranking_seconds"],
        "calls_per_node": fits["calls_per_node"],
    }
    path = common.write_json(RECEIPT_DIR / f"t2_cost_model_v{args.version}.json", receipt)
    fit = fits["ranking_seconds"]
    print(f"ranking seconds: log10(s) = {fit['intercept']} + {fit['slope']} log10(edges), "
          f"R squared {fit['r_squared']} on {fit['points']} runs")
    print(f"  valid from {fit['valid_edges_from']} to {fit['valid_edges_to']} edges, "
          f"largest absolute ratio {fit['largest_absolute_ratio']}")
    for row in fit["residuals"]:
        print(f"    {row['slug']}: observed {row['observed_seconds_per_ranking']} s, predicted "
              f"{row['predicted_seconds_per_ranking']} s, ratio "
              f"{row['ratio_predicted_over_observed']}")
    calls = fits["calls_per_node"]
    print(f"calls per node: {calls['slope']}, R squared {calls['r_squared']} on "
          f"{calls['points']} runs")
    print(f"model: {model_path}")
    print(f"receipt: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
