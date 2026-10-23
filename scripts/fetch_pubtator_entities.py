"""Resolve each PubTator3 disease and chemical accession to its source database identifier.

A relation record names an entity by an accession such as `@DISEASE_Aneuploidy`. The graph needs
a stable identifier such as `MESH:D000782`. The autocomplete endpoint returns that identifier.
The script keeps a hit only when the returned `_id` equals the accession, so a near match never
becomes a wrong identifier.

A gene needs no call. The graph identifies a gene by its symbol, which the accession already holds.

Usage:
    uv run python scripts/fetch_pubtator_entities.py
    uv run python scripts/fetch_pubtator_entities.py --refresh
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402

AUTOCOMPLETE_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/entity/autocomplete/"
CACHE_NAME = "pubtator_entities.json"
RESOLVED_KINDS = ("disease", "chemical")

# The number of accessions between two writes of the cache. A stage that stops on a network refusal
# keeps the work of every completed block, so a second run resumes.
CHECKPOINT_EVERY = 50


def pick_hit(accession: str, hits: list[dict]) -> dict:
    """Choose the autocomplete hit that is the entity itself.

    Args:
        accession: The accession that was queried, for example `@CHEMICAL_Nocodazole`.
        hits: The autocomplete results.

    Returns:
        The matching hit, or an empty mapping when no hit has the same accession.
    """
    for hit in hits or []:
        if hit.get("_id") == accession:
            return hit
    return {}


def resolve_accession(client, accession: str) -> dict:
    """Call the autocomplete endpoint for one accession.

    Args:
        client: The HTTP client.
        accession: The accession to resolve.

    Returns:
        A record with the database identifier, the database name and the display name. The
        identifier is empty when the endpoint returned no exact match.
    """
    label = common.accession_label(accession)
    payload = common.get_json(client, AUTOCOMPLETE_URL, {"query": label, "limit": 10})
    hits = payload if isinstance(payload, list) else []
    hit = pick_hit(accession, hits)
    return {"accession": accession, "identifier": hit.get("db_id", ""),
            "database": hit.get("db", ""), "name": hit.get("name", label),
            "biotype": hit.get("biotype", common.accession_kind(accession))}


def write_cache(path: Path, wanted: list[str], entries: dict) -> int:
    """Write the resolved accessions to the cache file.

    The write also runs while the loop is still going, so a stage that stops on a network refusal
    keeps the accessions that it already resolved and a second run resumes from them.

    Args:
        path: The cache file.
        wanted: Every accession that the graph needs.
        entries: The records resolved so far.

    Returns:
        The number of wanted accessions that carry a database identifier.
    """
    hits = sum(1 for accession in wanted if entries.get(accession, {}).get("identifier"))
    common.write_json(path, {
        "fetched_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "endpoint": AUTOCOMPLETE_URL, "requested": len(wanted), "resolved": hits,
        "entities": entries,
    })
    return hits


def accessions_to_resolve(nodes: set[str]) -> list[str]:
    """Select the accessions that need a database identifier.

    Args:
        nodes: Every node accession of the graph.

    Returns:
        The disease and chemical accessions, in sorted order.
    """
    wanted = [a for a in nodes if common.accession_kind(a) in RESOLVED_KINDS]
    return sorted(wanted)


def main() -> int:
    """Resolve every disease and chemical accession of the graph and cache the result.

    Returns:
        0 on success, 1 when the relation cache is missing.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--raw-dir", type=Path, default=common.RAW_DIR,
                        help="The raw cache directory.")
    parser.add_argument("--tier2-chemical-min-publications", type=int, default=20,
                        help="The smallest publication count of a second tier chemical.")
    parser.add_argument("--refresh", action="store_true",
                        help="Resolve every accession again instead of reading the cache.")
    args = parser.parse_args()

    index = common.read_json(args.raw_dir / "pubtator_relations_index.json")
    if not index:
        print("FAIL: run scripts/fetch_pubtator_relations.py first.")
        return 1
    seeds = index["seed_genes"]
    neighbours = index["derived_neighbours"]
    records = common.load_relation_cache(args.raw_dir, seeds + neighbours, "pubtator_relations")
    nodes, _ = common.select_graph(records, seeds, neighbours,
                                   args.tier2_chemical_min_publications)
    nodes.update(common.consequence_nodes(args.raw_dir,
                                          args.tier2_chemical_min_publications))
    wanted = accessions_to_resolve(nodes)

    path = args.raw_dir / CACHE_NAME
    cached = common.read_json(path) or {}
    entries = {} if args.refresh else dict(cached.get("entities", {}))
    todo = [a for a in wanted if a not in entries]
    print(f"  {len(wanted)} disease and chemical accessions, {len(todo)} to resolve now")

    try:
        with common.make_client() as client:
            for number, accession in enumerate(todo, start=1):
                entries[accession] = resolve_accession(client, accession)
                if number % CHECKPOINT_EVERY == 0:
                    write_cache(path, wanted, entries)
                    print(f"    {number}/{len(todo)} done")
                time.sleep(common.DEFAULT_PAUSE_S)
    finally:
        hits = write_cache(path, wanted, entries)

    print(f"{hits} of {len(wanted)} accessions carry a database identifier. Cache: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
