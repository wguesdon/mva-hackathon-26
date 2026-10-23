"""Fetch the approved and trial indications of every compound in the Track 2 graph.

The graph needs an edge that joins a disease to a compound that treats it. Without that edge the
held-out indication control of `scripts/track2_controls.py` has no case to remove, so the ranking
cannot be validated.

The source is the ChEMBL `drug_indication` endpoint. A record maps a molecule to an indication and
carries both a MeSH identifier and an EFO identifier, plus the highest clinical phase reached for
that indication. The MeSH identifier matters here. The disease nodes of the literature layer carry
MeSH identifiers too, so a ChEMBL record joins straight onto a disease that already holds gene
edges. An indication whose disease is not already in the graph is dropped, because a compound that
reaches a disease by no other route teaches the control nothing.

The endpoint gives no PMID. It gives clinical trial identifiers, EMA records and ATC codes. The
edge therefore carries an empty PMID list and counts the distinct references as its evidence.

Usage:
    uv run python scripts/fetch_chembl_indications.py
    uv run python scripts/fetch_chembl_indications.py --refresh
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402

CHEMBL_URL = "https://www.ebi.ac.uk/chembl/api/data/drug_indication.json"
CACHE_NAME = "chembl_indications.json"
BATCH_SIZE = 20
PAGE_LIMIT = 1000


def graph_drug_ids(targets: dict) -> list[str]:
    """Collect the ChEMBL molecule identifiers that the Open Targets layer introduced.

    Args:
        targets: The Open Targets cache, keyed by gene symbol.

    Returns:
        The molecule identifiers in sorted order.
    """
    found = {drug.get("drug_id", "") for record in targets.values()
             for drug in record.get("drugs", [])}
    return sorted(identifier for identifier in found if identifier)


def graph_disease_mesh(entities: dict) -> set[str]:
    """Collect the MeSH identifiers of the disease nodes of the literature layer.

    Args:
        entities: The resolved PubTator3 entity cache, keyed by accession.

    Returns:
        The MeSH identifiers, for example `{"D001943"}`.
    """
    found: set[str] = set()
    for accession, entry in entities.items():
        biotype = entry.get("biotype") or common.accession_kind(accession)
        if biotype == "disease" and entry.get("database") == "ncbi_mesh" and entry.get("identifier"):
            found.add(entry["identifier"])
    return found


def reference_ids(record: dict) -> set[str]:
    """Collect the distinct reference identifiers of one indication record.

    A reference identifier field holds a comma separated list of clinical trial identifiers, so
    each identifier is split out and counted once.

    Args:
        record: One `drug_indications` record.

    Returns:
        The distinct reference identifiers, which may be empty.
    """
    found: set[str] = set()
    for reference in record.get("indication_refs") or []:
        for identifier in str(reference.get("ref_id", "") or "").split(","):
            if identifier.strip():
                found.add(identifier.strip())
    return found


def reference_count(record: dict) -> int:
    """Count the distinct references that support one indication record.

    Args:
        record: One `drug_indications` record.

    Returns:
        The number of distinct reference identifiers, and at least 1.
    """
    return max(len(reference_ids(record)), 1)


def keep_indications(records: list[dict], disease_mesh: set[str]) -> list[dict]:
    """Keep the indication records whose disease already belongs to the graph.

    Args:
        records: The raw `drug_indications` records.
        disease_mesh: The MeSH identifiers of the graph disease nodes.

    Returns:
        One record per molecule and disease pair, with the highest phase and the reference count.
    """
    merged: dict[tuple[str, str], dict] = {}
    references: dict[tuple[str, str], set[str]] = {}
    for record in records:
        mesh_id = record.get("mesh_id", "") or ""
        molecule = record.get("molecule_chembl_id", "") or ""
        if not molecule or mesh_id not in disease_mesh:
            continue
        key = (molecule, mesh_id)
        phase = float(record.get("max_phase_for_ind", 0) or 0)
        references.setdefault(key, set()).update(reference_ids(record))
        entry = merged.get(key)
        if entry is None:
            merged[key] = {"molecule_chembl_id": molecule, "mesh_id": mesh_id,
                           "mesh_heading": record.get("mesh_heading", "") or "",
                           "efo_id": record.get("efo_id", "") or "",
                           "efo_term": record.get("efo_term", "") or "",
                           "max_phase_for_ind": phase, "reference_count": 0}
            continue
        entry["max_phase_for_ind"] = max(entry["max_phase_for_ind"], phase)
    for key, entry in merged.items():
        entry["reference_count"] = max(len(references.get(key, set())), 1)
    return [merged[key] for key in sorted(merged)]


def fetch_batch(client, molecules: list[str]) -> list[dict]:
    """Retrieve every indication record of one batch of molecules.

    Args:
        client: The HTTP client.
        molecules: The ChEMBL molecule identifiers of the batch.

    Returns:
        The raw `drug_indications` records of the batch.
    """
    records: list[dict] = []
    offset = 0
    while True:
        payload = common.get_json(client, CHEMBL_URL,
                                  {"molecule_chembl_id__in": ",".join(molecules),
                                   "limit": PAGE_LIMIT, "offset": offset})
        batch = payload.get("drug_indications", []) or []
        records.extend(batch)
        meta = payload.get("page_meta", {}) or {}
        if not meta.get("next") or not batch:
            break
        offset += PAGE_LIMIT
        time.sleep(common.DEFAULT_PAUSE_S)
    return records


def main() -> int:
    """Fetch the indications of every graph compound and cache the joined records.

    Returns:
        0 on success, 1 when a required cache file is missing.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--raw-dir", type=Path, default=common.RAW_DIR,
                        help="The raw cache directory.")
    parser.add_argument("--batch-size", type=int, default=BATCH_SIZE,
                        help="The number of molecules per request.")
    parser.add_argument("--refresh", action="store_true",
                        help="Fetch every molecule again instead of reading the cache.")
    args = parser.parse_args()

    drug_cache = common.read_json(args.raw_dir / "open_targets_targets.json") or {}
    entity_cache = common.read_json(args.raw_dir / "pubtator_entities.json") or {}
    if not drug_cache or not entity_cache:
        print("FAIL: run the PubTator3 and the Open Targets fetch scripts first.")
        return 1
    molecules = graph_drug_ids(drug_cache.get("targets", {}))
    disease_mesh = graph_disease_mesh(entity_cache.get("entities", {}))
    print(f"  {len(molecules)} compounds in the graph, "
          f"{len(disease_mesh)} disease nodes carry a MeSH identifier")

    path = args.raw_dir / CACHE_NAME
    cached = common.read_json(path) or {}
    raw_by_molecule = {} if args.refresh else dict(cached.get("raw_by_molecule", {}))
    todo = [m for m in molecules if m not in raw_by_molecule]
    print(f"  {len(todo)} compounds to query now")

    with common.make_client() as client:
        batches = [todo[start:start + args.batch_size]
                   for start in range(0, len(todo), args.batch_size)]
        for number, batch in enumerate(batches, start=1):
            records = fetch_batch(client, batch)
            for molecule in batch:
                raw_by_molecule[molecule] = []
            for record in records:
                molecule = record.get("molecule_chembl_id", "")
                if molecule in raw_by_molecule:
                    raw_by_molecule[molecule].append(record)
            if number % 5 == 0:
                print(f"    {number}/{len(batches)} batches fetched", flush=True)
            time.sleep(common.DEFAULT_PAUSE_S)

    flat = [record for records in raw_by_molecule.values() for record in records]
    indications = keep_indications(flat, disease_mesh)
    diseases = {entry["mesh_id"] for entry in indications}
    compounds = {entry["molecule_chembl_id"] for entry in indications}
    payload = {"fetched_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
               "endpoint": CHEMBL_URL, "compounds_requested": len(molecules),
               "raw_indication_records": len(flat),
               "indications_on_graph_diseases": len(indications),
               "diseases_with_an_indication": len(diseases),
               "compounds_with_an_indication": len(compounds),
               "indications": indications, "raw_by_molecule": raw_by_molecule}
    common.write_json(path, payload)
    print(f"{len(flat)} raw indication records. {len(indications)} join a graph disease, "
          f"covering {len(diseases)} diseases and {len(compounds)} compounds. Cache: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
