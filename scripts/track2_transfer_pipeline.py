"""Run the whole Track 2 pipeline for one rare disease and record what the run cost.

The Track 2 report measured one disease. This driver runs the same stages for any disease that
`literature/track2_transfer_diseases.yaml` names, from the first PubTator3 call to the two
controls, and it records the wall clock time, the number of API calls per host and the number of
response bytes of every stage. The run needs no hand step. The only input is the disease row, and
that row holds only public accessions.

Every stage is a separate process, so a stage that already produced its output is skipped on a
second run. `--force` runs a stage again. `--only` runs one stage.

The stages run in this order.

| Stage | What it does |
|---|---|
| `relations` | The PubTator3 relations of the seed gene and of the neighbours that it reveals. |
| `entities` | The source database identifier of every disease and chemical accession. |
| `counts` | The publication count of every entity, for the specificity weight. |
| `consequences` | The disease states that the seed gene produces, and what acts on them. |
| `open_targets` | The compounds and the Reactome pathways of every gene. |
| `chembl` | The indications of every compound. |
| `evidence` | The PMIDs and the evidence sentence of each entity pair. |
| `graph_unfiltered` | The graph with the mention rule off. |
| `surface_forms` | The mention rule of the disease node, derived from its MeSH entry terms. |
| `mention_fetch` | Every sentence of every literature edge of the disease node. |
| `mention_apply` | The verdict of the rule on every sentence and on every edge. |
| `graph_filtered` | The graph with the mention rule on. |
| `rank` | The ranking from the disease node. |
| `heldout` | The held-out indication recovery control with its degree-matched null. |
| `permutation` | The degree-preserving permutation control. |

Usage:
    uv run python scripts/track2_transfer_pipeline.py --disease cystic_fibrosis
    uv run python scripts/track2_transfer_pipeline.py --disease all --only rank
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402
import track2_peak_memory as peak_memory  # noqa: E402

ROOT = common.ROOT
DEFAULT_CONFIG = ROOT / "literature" / "track2_transfer_diseases.yaml"
DEFAULT_BASE = ROOT / "reports" / "track2" / "transfer"
SURFACE_FORM_DIR = ROOT / "literature" / "track2_surface_forms"
RECEIPT_DIR = ROOT / "tracking" / "receipts"

STAGE_ORDER = ("relations", "entities", "counts", "consequences", "open_targets", "chembl",
               "evidence", "graph_unfiltered", "surface_forms", "mention_fetch", "mention_apply",
               "graph_filtered", "rank", "heldout", "permutation")

# `scripts/track2_controls.py` returns 1 when a control reads FAIL. A failed control is a result of
# the run and not a failure of the stage, so the driver records the verdict and continues.
CONTROL_STAGES = ("heldout", "permutation")

# The parameters that ADR 0007 fixed before any result existed.
NEIGHBOUR_MIN_PUBLICATIONS = 5
MAX_NEIGHBOURS = 40
TIER2_CHEMICAL_MIN = 20
EVIDENCE_MAX_PAIRS = 2600
EVIDENCE_MAX_PAGES = 3
CONSEQUENCE_MIN_SEEDS = 1
CONSEQUENCE_MIN_NPMI = 0.07
CONSEQUENCE_LIMIT = 12
HELDOUT_TOP_K = 10
HELDOUT_PERMUTATIONS = 5
PERMUTATION_COUNT = 200
RANDOM_SEED = 20260915
MAX_DISEASE_DRUGS = 9
MAX_DISEASE_GENES = 10


def load_diseases(path: Path) -> list[dict]:
    """Read the disease rows.

    Args:
        path: The configuration file.

    Returns:
        The disease records in file order.
    """
    data = yaml.safe_load(path.read_text()) or {}
    return list(data.get("diseases", []))


def select_diseases(rows: list[dict], name: str) -> list[dict]:
    """Choose the disease rows that a command line name selects.

    Args:
        rows: Every disease record.
        name: A slug, or `all`.

    Returns:
        The selected records.

    Raises:
        ValueError: When no record carries the slug.
    """
    if name == "all":
        return rows
    chosen = [row for row in rows if row["slug"] == name]
    if not chosen:
        raise ValueError(f"no disease with the slug {name}.")
    return chosen


def disease_paths(base: Path, disease: dict) -> dict[str, Path]:
    """Build the output paths of one disease.

    Args:
        base: The base output directory.
        disease: The disease record.

    Returns:
        The directories and the files that the stages read and write.
    """
    home = base / disease["slug"]
    raw = home / "raw"
    graph = home / "graph"
    return {
        "home": home, "raw": raw, "graph": graph,
        "provenance": home / "provenance", "ranking": home / "ranking",
        "relations_index": raw / "pubtator_relations_index.json",
        "entities": raw / "pubtator_entities.json",
        "counts": raw / "pubtator_entity_counts.json",
        "consequences": raw / "pubtator_consequences.json",
        "open_targets": raw / "open_targets_targets.json",
        "chembl": raw / "chembl_indications.json",
        "evidence": raw / "pubtator_evidence.json",
        "mention_cache": raw / "edge_provenance_evidence.json",
        "nodes_tsv": graph / "nodes.tsv",
        "edges_tsv": graph / "edges.tsv",
        "drug_evidence_tsv": graph / "drug_evidence.tsv",
        "manifest": graph / "manifest.json",
        "rule": SURFACE_FORM_DIR / f"{disease['slug']}_v1.yaml",
        "sentences": home / "provenance" / "mention_sentences.tsv",
        "verdicts": home / "provenance" / "mention_edge_verdicts.tsv",
        "mention_report": home / "provenance" / "mention_report.json",
        "unfiltered_manifest": home / "graph_manifest_unfiltered.json",
        "run_log": home / "run_log.json",
    }


def stage_command(stage: str, disease: dict, paths: dict[str, Path]) -> list[str]:
    """Build the command of one stage.

    Args:
        stage: The stage name.
        disease: The disease record.
        paths: The paths of `disease_paths`.

    Returns:
        The command as an argument list.

    Raises:
        ValueError: When the stage name is unknown.
    """
    genes = [entry["symbol"] for entry in disease["genes"]]
    raw = str(paths["raw"])
    graph = str(paths["graph"])
    node = disease["node"]
    if stage == "relations":
        return ["scripts/fetch_pubtator_relations.py", "--genes", *genes, "--raw-dir", raw,
                "--neighbour-min-publications", str(NEIGHBOUR_MIN_PUBLICATIONS),
                "--max-neighbours", str(MAX_NEIGHBOURS)]
    if stage == "entities":
        return ["scripts/fetch_pubtator_entities.py", "--raw-dir", raw,
                "--tier2-chemical-min-publications", str(TIER2_CHEMICAL_MIN)]
    if stage == "counts":
        return ["scripts/fetch_pubtator_entity_counts.py", "--raw-dir", raw,
                "--tier2-chemical-min-publications", str(TIER2_CHEMICAL_MIN)]
    if stage == "consequences":
        return ["scripts/fetch_consequence_relations.py", "--raw-dir", raw,
                "--tier2-chemical-min-publications", str(TIER2_CHEMICAL_MIN),
                "--min-seeds", str(CONSEQUENCE_MIN_SEEDS),
                "--min-npmi", str(CONSEQUENCE_MIN_NPMI), "--limit", str(CONSEQUENCE_LIMIT)]
    if stage == "open_targets":
        return ["scripts/fetch_open_targets_drugs.py", "--raw-dir", raw,
                "--tier2-chemical-min-publications", str(TIER2_CHEMICAL_MIN)]
    if stage == "chembl":
        return ["scripts/fetch_chembl_indications.py", "--raw-dir", raw]
    if stage == "evidence":
        return ["scripts/fetch_pubtator_evidence.py", "--raw-dir", raw,
                "--tier2-chemical-min-publications", str(TIER2_CHEMICAL_MIN),
                "--max-pairs", str(EVIDENCE_MAX_PAIRS), "--max-pages", str(EVIDENCE_MAX_PAGES)]
    if stage == "graph_unfiltered":
        return ["scripts/build_track2_graph.py", "--raw-dir", raw, "--graph-dir", graph,
                "--tier2-chemical-min-publications", str(TIER2_CHEMICAL_MIN),
                "--seed-mention-filter", "off", "--seed-node", node]
    if stage == "surface_forms":
        return ["scripts/track2_disease_surface_forms.py", "--mesh-id", disease["mesh_id"],
                "--node", node, "--accession", disease["accession"],
                "--label", disease["label"], "--out", str(paths["rule"])]
    if stage == "mention_fetch":
        return ["scripts/track2_mention_audit.py", "fetch", "--seed", node,
                "--edges", str(paths["edges_tsv"]), "--entities", str(paths["entities"]),
                "--cache", str(paths["mention_cache"])]
    if stage == "mention_apply":
        return ["scripts/track2_mention_audit.py", "apply", "--rule", str(paths["rule"]),
                "--cache", str(paths["mention_cache"]), "--entities", str(paths["entities"]),
                "--sentences-out", str(paths["sentences"]),
                "--verdicts-out", str(paths["verdicts"]),
                "--report-out", str(paths["mention_report"]), "--update-rule"]
    if stage == "graph_filtered":
        return ["scripts/build_track2_graph.py", "--raw-dir", raw, "--graph-dir", graph,
                "--tier2-chemical-min-publications", str(TIER2_CHEMICAL_MIN),
                "--seed-mention-filter", "strict", "--seed-node", node,
                "--seed-rule", str(paths["rule"]),
                "--seed-evidence", str(paths["mention_cache"]),
                "--seed-verdicts-out", str(paths["verdicts"])]
    if stage == "rank":
        return ["scripts/track2_rank.py", "--seed", node, "--nodes", str(paths["nodes_tsv"]),
                "--edges", str(paths["edges_tsv"]),
                "--drug-evidence", str(paths["drug_evidence_tsv"]),
                "--out-dir", str(paths["ranking"])]
    if stage == "heldout":
        return ["scripts/track2_controls.py", "heldout", "--nodes", str(paths["nodes_tsv"]),
                "--edges", str(paths["edges_tsv"]),
                "--drug-evidence", str(paths["drug_evidence_tsv"]),
                "--top-k", str(HELDOUT_TOP_K), "-n", str(HELDOUT_PERMUTATIONS),
                "--random-seed", str(RANDOM_SEED), "--regime", "drug_poor",
                "--max-disease-drugs", str(MAX_DISEASE_DRUGS),
                "--max-disease-genes", str(MAX_DISEASE_GENES),
                "--receipt-id", f"t2_transfer_heldout_{disease['slug']}_v1",
                "--note", f"Transfer run of ADR 0007 for {disease['label']}."]
    if stage == "permutation":
        return ["scripts/track2_controls.py", "permutation", "--seed", node,
                "--nodes", str(paths["nodes_tsv"]), "--edges", str(paths["edges_tsv"]),
                "--drug-evidence", str(paths["drug_evidence_tsv"]),
                "-n", str(PERMUTATION_COUNT), "--top-k", str(HELDOUT_TOP_K),
                "--random-seed", str(RANDOM_SEED),
                "--receipt-id", f"t2_transfer_permutation_{disease['slug']}_v1",
                "--note", f"Transfer run of ADR 0007 for {disease['label']}."]
    raise ValueError(f"unknown stage {stage}.")


def stage_output(stage: str, disease: dict, paths: dict[str, Path]) -> Path:
    """Return the file whose presence means that a stage already ran.

    Args:
        stage: The stage name.
        disease: The disease record.
        paths: The paths of `disease_paths`.

    Returns:
        The output path.
    """
    by_stage = {
        "relations": paths["relations_index"], "entities": paths["entities"],
        "counts": paths["counts"], "consequences": paths["consequences"],
        "open_targets": paths["open_targets"], "chembl": paths["chembl"],
        "evidence": paths["evidence"], "graph_unfiltered": paths["unfiltered_manifest"],
        "surface_forms": paths["rule"], "mention_fetch": paths["mention_cache"],
        "mention_apply": paths["mention_report"], "graph_filtered": paths["manifest"],
        "rank": paths["ranking"],
        "heldout": RECEIPT_DIR / f"t2_transfer_heldout_{disease['slug']}_v1.json",
        "permutation": RECEIPT_DIR / f"t2_transfer_permutation_{disease['slug']}_v1.json",
    }
    return by_stage[stage]


def run_stage(stage: str, disease: dict, paths: dict[str, Path], log_dir: Path) -> dict:
    """Run one stage in its own process and record what it cost.

    Args:
        stage: The stage name.
        disease: The disease record.
        paths: The paths of `disease_paths`.
        log_dir: The directory of the API call counter files.

    Returns:
        The stage record, with the wall clock time, the peak memory, the exit status and the API
        counters.
    """
    command = stage_command(stage, disease, paths)
    counter_path = log_dir / f"{disease['slug']}_{stage}_calls.json"
    if counter_path.is_file():
        counter_path.unlink()
    environment = dict(os.environ)
    environment[common.CALL_LOG_VARIABLE] = str(counter_path)
    started = time.time()
    process = subprocess.Popen(["uv", "run", "python", *command], cwd=ROOT, env=environment,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    sampler = peak_memory.start_sampler(process.pid)
    stdout, stderr = process.communicate()
    peak_mib = peak_memory.stop_sampler(sampler)
    seconds = time.time() - started
    counters = common.read_json(counter_path) or {}
    record = {
        "stage": stage,
        "command": "uv run python " + " ".join(command),
        "started_utc": dt.datetime.fromtimestamp(started, dt.timezone.utc)
                         .isoformat(timespec="seconds"),
        "wall_seconds": round(seconds, 1),
        "peak_memory_mib": peak_mib,
        "exit_status": process.returncode,
        "api_calls": counters,
        "api_call_total": sum(entry["calls"] for entry in counters.values()),
        "api_bytes_total": sum(entry["bytes"] for entry in counters.values()),
        "stdout_tail": stdout.strip().splitlines()[-12:],
        "stderr_tail": stderr.strip().splitlines()[-6:],
    }
    return record


def merge_log(path: Path, log: dict, disease: dict) -> dict:
    """Write the run log after a stage, keeping the stages that another process wrote.

    Two stages of one disease can run in two processes, for example the two controls on a large
    graph. Each process holds its own copy of the log, so the copy is merged with the file on disk
    before the write and no stage record is lost.

    Args:
        path: The run log file.
        log: The log of this process.
        disease: The disease record.

    Returns:
        The merged log.
    """
    on_disk = common.read_json(path) or {}
    stages = dict(on_disk.get("stages", {}))
    stages.update(log.get("stages", {}))
    merged = dict(on_disk)
    merged["disease"] = disease
    merged["stages"] = stages
    merged["updated_utc"] = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    common.write_json(path, merged)
    return merged


def snapshot_unfiltered(paths: dict[str, Path]) -> None:
    """Keep the manifest of the build that runs with the mention rule off.

    Args:
        paths: The paths of `disease_paths`.
    """
    manifest = common.read_json(paths["manifest"])
    if manifest is not None:
        common.write_json(paths["unfiltered_manifest"], manifest)


def run_disease(disease: dict, base: Path, stages: list[str], force: bool) -> dict:
    """Run the selected stages of one disease.

    Args:
        disease: The disease record.
        base: The base output directory.
        stages: The stages to run, in order.
        force: True to run a stage whose output already exists.

    Returns:
        The run log of the disease.
    """
    paths = disease_paths(base, disease)
    for key in ("home", "raw", "graph", "provenance", "ranking"):
        paths[key].mkdir(parents=True, exist_ok=True)
    log_dir = paths["home"] / "call_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log = common.read_json(paths["run_log"]) or {"disease": disease, "stages": {}}
    log["disease"] = disease
    for stage in stages:
        output = stage_output(stage, disease, paths)
        if not force and output.exists() and stage in log.get("stages", {}):
            print(f"[{disease['slug']}] {stage}: already done")
            continue
        print(f"[{disease['slug']}] {stage}: running")
        record = run_stage(stage, disease, paths, log_dir)
        log.setdefault("stages", {})[stage] = record
        if stage == "graph_unfiltered" and record["exit_status"] == 0:
            snapshot_unfiltered(paths)
        log = merge_log(paths["run_log"], log, disease)
        status = "ok" if record["exit_status"] == 0 else f"FAIL {record['exit_status']}"
        print(f"[{disease['slug']}] {stage}: {status} in {record['wall_seconds']} s, "
              f"{record['api_call_total']} API calls, {record['api_bytes_total']} bytes, "
              f"peak {record['peak_memory_mib']} MiB")
        for line in record["stdout_tail"]:
            print(f"    {line}")
        if record["exit_status"] != 0:
            for line in record["stderr_tail"]:
                print(f"    ! {line}")
            if stage not in CONTROL_STAGES:
                break
    return log


def main() -> int:
    """Run the transfer pipeline from the command line.

    Returns:
        0 when every stage of every selected disease returned 0, 1 otherwise.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG,
                        help="The disease configuration file.")
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE,
                        help="The base output directory.")
    parser.add_argument("--disease", default="all", help="A disease slug, or `all`.")
    parser.add_argument("--only", action="append", default=[],
                        help="Run only this stage. Repeat the option for several stages.")
    parser.add_argument("--from-stage", default="", help="Start at this stage.")
    parser.add_argument("--force", action="store_true", help="Run a stage that already ran.")
    args = parser.parse_args()

    stages = list(args.only) if args.only else list(STAGE_ORDER)
    if args.from_stage:
        stages = list(STAGE_ORDER[STAGE_ORDER.index(args.from_stage):])
    failures = 0
    for disease in select_diseases(load_diseases(args.config), args.disease):
        log = run_disease(disease, args.base, stages, args.force)
        failures += sum(1 for name, record in log.get("stages", {}).items()
                        if name not in CONTROL_STAGES
                        and record.get("exit_status") not in (0, None))
    print(f"stages with a non zero exit status: {failures}")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
