"""Run the whole Track 2 method on one rare disease, from a disease identifier to two controls.

This is the single entry point of the Track 2 method. It takes a disease identifier, resolves every
other identifier that the pipeline needs, checks its prerequisites, builds the graph, ranks the
candidates, runs the held-out control and the permutation control, and writes one receipt with the
result and the cost of the run. Every parameter carries the default that
`docs/decisions/0009_track2_one_entry_point_and_a_prospective_disease.md` fixes, so a reader who
gives only a disease identifier runs the recorded configuration.

The run has four phases.

| Phase | What it does | What stops it |
|---|---|---|
| Preflight | Reads the tool versions, the pinned release files, the four API hosts, the free disk space and the write permission of each output directory. | Exit status 2 and the name of the missing item. |
| Resolve | Turns the disease identifier into the causal gene, the MeSH identifier, the PubTator3 accession and the graph node identifier. | Exit status 3 and the identifier that a person must supply. |
| Build and rank | Runs thirteen stages, from the PubTator3 relations of the seed gene to the ranking. | Exit status 1 and the failing stage. |
| Control | Projects the cost of the held-out control, then runs it and the permutation control. | The controls report FAIL as a result and the run continues. |

The projection between the third phase and the fourth is the guard that lets a reader stop before
the expensive stage. `--estimate-only` stops there, and `--max-heldout-minutes` stops there when the
projection exceeds a budget that the reader sets.

Usage:
    uv run python scripts/track2_run_disease.py --disease OMIM:268300
    uv run python scripts/track2_run_disease.py --disease "Bloom syndrome" --estimate-only
    uv run python scripts/track2_run_disease.py --check-only
"""

from __future__ import annotations

import argparse
import datetime as dt
import shutil
import subprocess
import sys
import time
from pathlib import Path

import httpx
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402
import track2_cost_model_fit as cost_model  # noqa: E402
import track2_disease_resolver as resolver  # noqa: E402
import track2_source_versions as versions  # noqa: E402
import track2_transfer_pipeline as pipeline  # noqa: E402

ROOT = common.ROOT
RESOURCES = ROOT / "literature" / "track2_runner_resources.yaml"
DEFAULT_MODEL = ROOT / "literature" / "track2_cost_model_v2.yaml"
DEFAULT_BASE = ROOT / "reports" / "track2" / "transfer"
RECEIPT_DIR = ROOT / "tracking" / "receipts"

# The stages that build the graph and the ranking, and the stages that run the two controls. The
# runner prints the cost projection between the two groups.
BUILD_STAGES = tuple(stage for stage in pipeline.STAGE_ORDER if stage not in ("heldout",
                                                                             "permutation"))
CONTROL_STAGES = ("heldout", "permutation")

EXIT_OK = 0
EXIT_STAGE_FAILED = 1
EXIT_PREFLIGHT_FAILED = 2
EXIT_RESOLUTION_FAILED = 3
EXIT_OVER_BUDGET = 4

GIB = 1024 ** 3
SECONDS_PER_MINUTE = 60.0

# NCBI answers a burst with HTTP 429, so a probe that follows an earlier run retries before it
# reports the host as unreachable.
PROBE_ATTEMPTS = 3
PROBE_RETRY_PAUSE_S = 4.0

# The wait before the runner repeats a stage that a transient network refusal ended.
STAGE_RETRY_PAUSE_S = 60.0
DEFAULT_STAGE_RETRIES = 3


def repository_commit() -> str:
    """Return the commit of the repository that runs the pipeline.

    Returns:
        The commit hash, with the suffix `-dirty` when the working tree carries a change, or an
        empty string when the directory is not a git repository.
    """
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                              text=True, check=False, timeout=30)
        status = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, capture_output=True,
                                text=True, check=False, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return ""
    if head.returncode != 0:
        return ""
    return head.stdout.strip() + ("-dirty" if status.stdout.strip() else "")


def load_resources(path: Path = RESOURCES) -> dict:
    """Read the resource manifest.

    Args:
        path: The manifest file.

    Returns:
        The manifest.

    Raises:
        FileNotFoundError: When the manifest is missing.
    """
    if not path.is_file():
        raise FileNotFoundError(f"{path} is missing, so the runner cannot check its resources.")
    return yaml.safe_load(path.read_text()) or {}


def check_tool(name: str, command: list[str]) -> dict:
    """Check that one command line tool answers.

    Args:
        name: The tool name for the report.
        command: The command that prints its version.

    Returns:
        The check record, with `ok`, the version text and the repair instruction.
    """
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return {"check": name, "ok": False, "detail": "the command did not run",
                "fix": f"install {name} and put it on PATH."}
    if result.returncode != 0:
        return {"check": name, "ok": False, "detail": result.stderr.strip()[:200],
                "fix": f"install {name} and put it on PATH."}
    return {"check": name, "ok": True, "detail": result.stdout.strip().splitlines()[0]}


def check_python(minimum: str) -> dict:
    """Check the Python version of the interpreter that runs the pipeline.

    Args:
        minimum: The smallest accepted version, for example `3.11`.

    Returns:
        The check record.
    """
    wanted = tuple(int(part) for part in minimum.split("."))
    found = sys.version_info[:len(wanted)]
    return {
        "check": "python",
        "ok": found >= wanted,
        "detail": ".".join(str(part) for part in sys.version_info[:3]),
        "fix": f"use Python {minimum} or later. `uv python install {minimum}` installs one.",
    }


def check_release_file(entry: dict) -> dict:
    """Check that a pinned release file is present and unchanged.

    Args:
        entry: One record of the `release_files` block of the manifest.

    Returns:
        The check record.
    """
    path = ROOT / entry["path"]
    fix = (f"download release {entry['release']} with "
           f"`curl -L -o {entry['path']} {entry['download']}`.")
    if not path.is_file():
        return {"check": entry["id"], "ok": False, "detail": f"{entry['path']} is missing",
                "fix": fix}
    digest = common.sha256_of(path)
    if digest != entry["sha256"]:
        return {"check": entry["id"], "ok": False,
                "detail": (f"{entry['path']} has the SHA-256 {digest[:16]} and release "
                           f"{entry['release']} has {entry['sha256'][:16]}"),
                "fix": fix}
    return {"check": entry["id"], "ok": True, "detail": f"release {entry['release']}, {digest[:16]}"}


def probe_once(client: httpx.Client, entry: dict) -> tuple[bool, str]:
    """Send one probe to an API host.

    Args:
        client: The HTTP client.
        entry: One record of the `endpoints` block of the manifest.

    Returns:
        Whether the host answered, and the detail for the report.
    """
    try:
        if entry.get("probe_method") == "POST":
            response = client.post(entry["probe"], content=entry.get("probe_body", ""),
                                   headers={"Content-Type": "application/json"})
        else:
            response = client.get(entry["probe"])
    except httpx.HTTPError as error:
        return False, type(error).__name__
    return response.status_code < 400, f"HTTP {response.status_code}"


def check_endpoint(client: httpx.Client, entry: dict, attempts: int = PROBE_ATTEMPTS) -> dict:
    """Check that one API host answers.

    NCBI answers a burst of calls with HTTP 429, and a preflight probe that follows an earlier run
    can meet that status while the host is healthy. The check therefore probes again after a pause
    before it reports a failure.

    Args:
        client: The HTTP client.
        entry: One record of the `endpoints` block of the manifest.
        attempts: The number of probes before the check fails.

    Returns:
        The check record.
    """
    fix = (f"the host {entry['host']} supplies {entry['what_it_supplies'].lower()} "
           "Check the network connection and any proxy, then run the check again.")
    detail = ""
    for attempt in range(attempts):
        ok, detail = probe_once(client, entry)
        if ok:
            return {"check": entry["id"], "ok": True, "detail": detail}
        if attempt + 1 < attempts:
            time.sleep(PROBE_RETRY_PAUSE_S * (attempt + 1))
    return {"check": entry["id"], "ok": False, "detail": f"{detail} on {attempts} probes",
            "fix": fix}


def check_disk(base: Path, minimum_gib: int) -> dict:
    """Check the free disk space of the output directory.

    Args:
        base: The output base directory.
        minimum_gib: The smallest accepted free space in gibibytes.

    Returns:
        The check record.
    """
    target = base
    while not target.exists() and target != target.parent:
        target = target.parent
    free_gib = shutil.disk_usage(target).free / GIB
    return {
        "check": "free_disk",
        "ok": free_gib >= minimum_gib,
        "detail": f"{free_gib:.1f} GiB free at {target}",
        "fix": (f"the raw responses of one run occupy between 40 MB and 245 MB and the graph "
                f"tables occupy up to 100 MB. Free {minimum_gib} GiB under {base}."),
    }


def check_writable(base: Path) -> dict:
    """Check that the runner can write under the output base directory.

    Args:
        base: The output base directory.

    Returns:
        The check record.
    """
    try:
        base.mkdir(parents=True, exist_ok=True)
        probe = base / ".write_probe"
        probe.write_text("probe", encoding="utf-8")
        probe.unlink()
    except OSError as error:
        return {"check": "writable_output", "ok": False, "detail": str(error)[:200],
                "fix": f"give the account write permission on {base}."}
    return {"check": "writable_output", "ok": True, "detail": str(base)}


def check_cost_model(path: Path) -> dict:
    """Check that the cost model file is present.

    Args:
        path: The cost model file.

    Returns:
        The check record.
    """
    if not path.is_file():
        return {"check": "cost_model", "ok": False, "detail": f"{path} is missing",
                "fix": ("run `uv run python scripts/track2_cost_model_fit.py --version 1 --run "
                        "cystic_fibrosis --run nijmegen_breakage_syndrome --run "
                        "duchenne_muscular_dystrophy` to fit it from the recorded runs.")}
    return {"check": "cost_model", "ok": True, "detail": str(path.relative_to(ROOT))}


def preflight(base: Path, model_path: Path, resources: dict) -> list[dict]:
    """Run every prerequisite check.

    Args:
        base: The output base directory.
        model_path: The cost model file.
        resources: The resource manifest.

    Returns:
        One record per check, in the order that the runner reports them.
    """
    settings = resources.get("preflight", {})
    checks = [
        check_tool("uv", ["uv", "--version"]),
        check_python(str(settings.get("min_python", "3.11"))),
        check_cost_model(model_path),
    ]
    checks.extend(check_release_file(entry) for entry in resources.get("release_files", []))
    checks.append(check_writable(base))
    checks.append(check_disk(base, int(settings.get("min_free_disk_gib", 3))))
    timeout = float(settings.get("probe_timeout_s", 30))
    with httpx.Client(timeout=timeout, follow_redirects=True,
                      headers={"User-Agent": common.USER_AGENT}) as client:
        checks.extend(check_endpoint(client, entry) for entry in resources.get("endpoints", []))
    return checks


def report_checks(checks: list[dict]) -> bool:
    """Print the result of every check and return whether all of them passed.

    Args:
        checks: The check records.

    Returns:
        True when every check passed.
    """
    failed = [check for check in checks if not check["ok"]]
    for check in checks:
        mark = "ok  " if check["ok"] else "FAIL"
        print(f"  [{mark}] {check['check']}: {check['detail']}")
    for check in failed:
        print(f"FAIL: {check['check']}. {check.get('fix', '')}")
    return not failed


def projection(model: dict, manifest: dict, permutations: int) -> dict:
    """Project what the held-out control of a built graph will cost.

    Args:
        model: The cost model file content.
        manifest: The manifest of the built graph.
        permutations: The permutation count of the control.

    Returns:
        The projected rankings, seconds and minutes, with the range over which the model holds.
    """
    block = model["ranking_seconds"]
    edges = int(manifest["edge_count"])
    diseases = int(manifest["diseases_with_an_indication"])
    rankings = diseases * (permutations + 1)
    per_ranking = cost_model.predict_seconds_per_ranking(block, edges)
    seconds = per_ranking * rankings
    return {
        "graph_edges": edges,
        "graph_nodes": int(manifest["node_count"]),
        "diseases_with_an_indication": diseases,
        "permutations": permutations,
        "rankings": rankings,
        "model_version": model["version"],
        "model_fitted_on": model["fitted_on"],
        "predicted_seconds_per_ranking": round(per_ranking, 4),
        "predicted_wall_seconds": round(seconds, 1),
        "predicted_wall_minutes": round(seconds / SECONDS_PER_MINUTE, 1),
        "model_valid_edges_from": block["valid_edges_from"],
        "model_valid_edges_to": block["valid_edges_to"],
        "edges_inside_the_measured_range": (block["valid_edges_from"] <= edges
                                            <= block["valid_edges_to"]),
        "largest_absolute_ratio_of_the_fit": block["largest_absolute_ratio"],
    }


def print_projection(values: dict) -> None:
    """Print the cost projection of the held-out control.

    Args:
        values: The record that `projection` returned.
    """
    print("cost projection of the held-out control")
    print(f"  graph: {values['graph_edges']} edges, {values['graph_nodes']} nodes, "
          f"{values['diseases_with_an_indication']} diseases with an indication")
    print(f"  work: {values['rankings']} rankings at {values['permutations']} permutations")
    print(f"  model v{values['model_version']}: "
          f"{values['predicted_seconds_per_ranking']} seconds per ranking")
    print(f"  projected wall clock: {values['predicted_wall_minutes']} minutes")
    if not values["edges_inside_the_measured_range"]:
        print(f"  WARNING: the model was fitted between {values['model_valid_edges_from']} and "
              f"{values['model_valid_edges_to']} edges and this graph is outside that range, so "
              "the projection is an extrapolation.")
    print(f"  the largest absolute error of the fit on its own runs is a factor of "
          f"{values['largest_absolute_ratio_of_the_fit']}")


def source_versions() -> list[dict]:
    """Read the release identifier of every source that the run read.

    Returns:
        One record per source, or an empty list when a source did not answer.
    """
    try:
        return versions.collect()
    except (httpx.HTTPError, ValueError, KeyError):
        return []


def first_failed_stage(log: dict, stages: list[str]) -> str:
    """Return the first stage of a list that ended with a non zero status.

    Args:
        log: The run log of the disease.
        stages: The stages to inspect.

    Returns:
        The stage name, or an empty string when every stage returned zero.
    """
    for name in stages:
        record = log.get("stages", {}).get(name)
        if record and record.get("exit_status") not in (0, None):
            return name
    return ""


def run_with_retries(row: dict, base: Path, stages: list[str], force: bool,
                     retries: int) -> tuple[dict, str]:
    """Run a group of stages and repeat a stage that a transient refusal ended.

    Every fetch stage writes its cache as it goes and skips the work that the cache already holds,
    so a repeat resumes rather than starts again. The first run of this entry point met a run of
    HTTP 429 answers from NCBI that exhausted the attempt budget of one call, and the stage stopped
    with 178 calls spent and nothing kept. The retry and the cache checkpoint together close that.

    Args:
        row: The resolved disease row.
        base: The output base directory.
        stages: The stages to run, in order.
        force: True to run a stage whose output already exists.
        retries: The number of extra attempts of a failed stage.

    Returns:
        The run log and the name of the stage that still fails, or an empty string.
    """
    log = pipeline.run_disease(row, base, stages, force)
    failed = first_failed_stage(log, stages)
    attempt = 0
    while failed and attempt < retries:
        attempt += 1
        print(f"retry {attempt} of {retries}: the stage {failed} returned a non zero status. "
              f"Waiting {STAGE_RETRY_PAUSE_S} seconds, then running the remaining stages again.")
        time.sleep(STAGE_RETRY_PAUSE_S)
        remaining = stages[stages.index(failed):]
        log = pipeline.run_disease(row, base, remaining, force=True)
        failed = first_failed_stage(log, remaining)
    return log, failed


def run_receipt(row: dict, log: dict, projected: dict | None, checks: list[dict],
                started: str) -> dict:
    """Assemble the receipt of one run of the entry point.

    Args:
        row: The resolved disease row.
        log: The run log of the pipeline.
        projected: The cost projection, or None when the run stopped before it.
        checks: The preflight check records.
        started: The start time of the run.

    Returns:
        The receipt.
    """
    stages = log.get("stages", {})
    by_host: dict[str, dict[str, int]] = {}
    for record in stages.values():
        for host, counter in record.get("api_calls", {}).items():
            entry = by_host.setdefault(host, {"calls": 0, "bytes": 0, "failed": 0})
            entry["calls"] += counter.get("calls", 0)
            entry["bytes"] += counter.get("bytes", 0)
            entry["failed"] += counter.get("failed", 0)
    peaks = [record["peak_memory_mib"] for record in stages.values()
             if record.get("peak_memory_mib") is not None]
    seconds = sum(record.get("wall_seconds", 0.0) for record in stages.values())
    return {
        "receipt_id": f"t2_run_disease_{row['slug']}_v1",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "track": "t2_rank",
        "control": "one_command_run_of_the_track2_method_on_one_disease",
        "queue_item": "R17",
        "decision_record": ("docs/decisions/"
                            "0009_track2_one_entry_point_and_a_prospective_disease.md"),
        "data_class_of_run": "A",
        "data_class_of_receipt": "A",
        "receipt_scope": ("Public accessions, public counts and public identifiers only. No "
                          "patient data enters this run."),
        "script": "scripts/track2_run_disease.py",
        "repository_commit": repository_commit(),
        "started_utc": started,
        "disease": row,
        "preflight": checks,
        "cost_projection": projected,
        "cost": {
            "wall_seconds_total": round(seconds, 1),
            "wall_minutes_total": round(seconds / SECONDS_PER_MINUTE, 1),
            "peak_memory_mib": max(peaks) if peaks else None,
            "peak_memory_stage": max(stages.items(),
                                     key=lambda item: item[1].get("peak_memory_mib") or 0.0)[0]
                                 if peaks else None,
            "by_host": by_host,
            "api_calls_total": sum(entry["calls"] for entry in by_host.values()),
            "api_bytes_total": sum(entry["bytes"] for entry in by_host.values()),
            "money_usd": 0.0,
            "money_note": ("PubTator3, the NCBI E-utilities, Open Targets and ChEMBL serve every "
                           "endpoint of this runner with no charge and with no API key. The run "
                           "is local, on one CPU, with no cloud instance and no cloud storage, so "
                           "the run has no cloud spend."),
        },
        "per_stage": [
            {"stage": name, "wall_seconds": record.get("wall_seconds"),
             "peak_memory_mib": record.get("peak_memory_mib"),
             "exit_status": record.get("exit_status"),
             "api_calls": record.get("api_call_total"), "api_bytes": record.get("api_bytes_total")}
            for name, record in stages.items()
        ],
    }


def main() -> int:
    """Run the Track 2 method on one disease from the command line.

    Returns:
        0 on success, 1 when a build stage failed, 2 when a prerequisite is missing, 3 when the
        disease identifier did not resolve and 4 when the projected cost exceeds the budget.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--disease", help="An OMIM identifier, a MeSH identifier or a disease "
                                          "name, for example OMIM:268300.")
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE,
                        help="The base output directory.")
    parser.add_argument("--hpo-dir", type=Path, default=resolver.DEFAULT_HPO_DIR,
                        help="The directory of the pinned Human Phenotype Ontology release.")
    parser.add_argument("--cost-model", type=Path, default=DEFAULT_MODEL,
                        help="The cost model that projects the held-out control.")
    parser.add_argument("--check-only", action="store_true",
                        help="Run the preflight check and stop.")
    parser.add_argument("--estimate-only", action="store_true",
                        help="Build the graph, print the cost projection and stop.")
    parser.add_argument("--max-heldout-minutes", type=float, default=0.0,
                        help="Stop before the held-out control when the projection exceeds this "
                             "budget. 0 turns the budget off.")
    parser.add_argument("--stage-retries", type=int, default=DEFAULT_STAGE_RETRIES,
                        help="The number of extra attempts of a build stage that a transient "
                             "network refusal ended.")
    parser.add_argument("--force", action="store_true", help="Run a stage that already ran.")
    args = parser.parse_args()

    started = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    started_clock = time.time()
    try:
        resources = load_resources()
    except FileNotFoundError as error:
        print(f"FAIL: {error}")
        return EXIT_PREFLIGHT_FAILED

    print("preflight")
    checks = preflight(args.base, args.cost_model, resources)
    if not report_checks(checks):
        return EXIT_PREFLIGHT_FAILED
    if args.check_only:
        print("preflight: every check passed.")
        return EXIT_OK
    if not args.disease:
        print("FAIL: --disease is required unless --check-only is given.")
        return EXIT_PREFLIGHT_FAILED

    try:
        row = resolver.resolve(args.disease, args.hpo_dir)
    except resolver.ResolutionError as error:
        print(f"FAIL: {error}")
        return EXIT_RESOLUTION_FAILED
    row_path = resolver.write_row(row)
    print(f"disease: {row['label']} ({row['omim']}), genes "
          f"{', '.join(gene['symbol'] for gene in row['genes'])}, node {row['node']}")
    print(f"  resolved row: {row_path.relative_to(ROOT)}")

    log, failed = run_with_retries(row, args.base, list(BUILD_STAGES), args.force,
                                   args.stage_retries)
    if failed:
        print(f"FAIL: the stage {failed} returned a non zero status on "
              f"{args.stage_retries + 1} attempts. The run stopped. Every completed stage is on "
              "disk, so the same command resumes the run.")
        common.write_json(RECEIPT_DIR / f"t2_run_disease_{row['slug']}_v1.json",
                          run_receipt(row, log, None, checks, started))
        return EXIT_STAGE_FAILED

    paths = pipeline.disease_paths(args.base, row)
    manifest = common.read_json(paths["manifest"]) or {}
    model = yaml.safe_load(args.cost_model.read_text())
    projected = projection(model, manifest, pipeline.HELDOUT_PERMUTATIONS)
    print_projection(projected)

    if args.estimate_only:
        common.write_json(RECEIPT_DIR / f"t2_run_disease_{row['slug']}_v1.json",
                          run_receipt(row, log, projected, checks, started))
        print("estimate only: the two controls did not run.")
        return EXIT_OK
    if args.max_heldout_minutes and projected["predicted_wall_minutes"] > args.max_heldout_minutes:
        print(f"FAIL: the projected {projected['predicted_wall_minutes']} minutes exceed the "
              f"budget of {args.max_heldout_minutes} minutes. Raise --max-heldout-minutes or run "
              "the control on a machine with more compute.")
        common.write_json(RECEIPT_DIR / f"t2_run_disease_{row['slug']}_v1.json",
                          run_receipt(row, log, projected, checks, started))
        return EXIT_OVER_BUDGET

    log = pipeline.run_disease(row, args.base, list(CONTROL_STAGES), args.force)
    receipt = run_receipt(row, log, projected, checks, started)
    receipt["source_versions"] = source_versions()
    path = common.write_json(RECEIPT_DIR / f"t2_run_disease_{row['slug']}_v1.json", receipt)
    cost = receipt["cost"]
    print(f"run: {round((time.time() - started_clock) / SECONDS_PER_MINUTE, 1)} minutes of wall "
          f"clock in this process, {cost['wall_minutes_total']} minutes over the stages")
    print(f"  peak memory: {cost['peak_memory_mib']} MiB in the stage {cost['peak_memory_stage']}")
    print(f"  external calls: {cost['api_calls_total']}, response bytes "
          f"{cost['api_bytes_total']}, money {cost['money_usd']} USD")
    print(f"receipt: {path}")
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
