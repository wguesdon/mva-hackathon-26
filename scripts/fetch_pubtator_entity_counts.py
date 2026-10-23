"""Fetch the total publication count of each PubTator3 entity in the graph.

A raw co-mention count measures how much an entity is studied, not how specifically two entities
relate. AKT1 returns 30,502 relations because AKT1 appears in 355,552 publications, not because
AKT1 has a close tie to BUB1B. A specificity weight needs the total count of each entity, so that
the observed co-mention count can be compared against the count that the two totals predict.

The PubTator3 search endpoint gives that total. A query with a single accession as the whole text
returns a `count` field, which is the number of publications that mention the entity.

The corpus size is the other term of the weight. PubTator3 publishes no size for its own index, so
the script reads the PubMed record count from the NCBI EInfo endpoint and records it with the date.
Section 7 of `reports/track2/graph_build_notes.md` states what that approximation costs.

Usage:
    uv run python scripts/fetch_pubtator_entity_counts.py
    uv run python scripts/fetch_pubtator_entity_counts.py --refresh
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from itertools import repeat
from pathlib import Path

import httpx

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402

SEARCH_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/"
EINFO_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi"
CACHE_NAME = "pubtator_entity_counts.json"
SAVE_EVERY = 50

# The single entity query is heavier than a pair query, and NCBI throttles a burst of them.
# Two workers with the standard pause stay under the limit.
DEFAULT_WORKERS = 2


def parse_count(payload: dict) -> int:
    """Read the publication count of a single entity search.

    Args:
        payload: The decoded search response.

    Returns:
        The count, or 0 when the response carries none.
    """
    if not isinstance(payload, dict):
        return 0
    return int(payload.get("count", 0) or 0)


def fetch_count(client, accession: str) -> int:
    """Retrieve the number of publications that mention one entity.

    Args:
        client: The HTTP client.
        accession: The accession, for example `@GENE_BUB1B`.

    Returns:
        The publication count.
    """
    payload = common.get_json(client, SEARCH_URL, {"text": accession})
    time.sleep(common.DEFAULT_PAUSE_S)
    return parse_count(payload)


def parse_corpus_size(payload: dict) -> tuple[int, str]:
    """Read the PubMed record count and its update date from an EInfo response.

    Args:
        payload: The decoded EInfo response.

    Returns:
        The record count and the update date. The count is 0 when the response is not usable.
    """
    entries = ((payload or {}).get("einforesult", {}) or {}).get("dbinfo", []) or []
    if not entries:
        return 0, ""
    return int(entries[0].get("count", 0) or 0), entries[0].get("lastupdate", "")


def fetch_corpus_size(client) -> tuple[int, str]:
    """Retrieve the number of records in PubMed.

    Args:
        client: The HTTP client.

    Returns:
        The record count and the update date.
    """
    payload = common.get_json(client, EINFO_URL, {"db": "pubmed", "retmode": "json"})
    return parse_corpus_size(payload)


def wanted_accessions(nodes: set[str]) -> list[str]:
    """Select the accessions that need a total publication count.

    Args:
        nodes: Every node accession of the graph.

    Returns:
        The gene, disease and chemical accessions, in sorted order.
    """
    kinds = {"gene", "disease", "chemical"}
    return sorted(a for a in nodes if common.accession_kind(a) in kinds)


def main() -> int:
    """Fetch the total count of every graph entity and the corpus size, then cache both.

    Returns:
        0 on success, 1 when the relation cache is missing or the corpus size failed.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--raw-dir", type=Path, default=common.RAW_DIR,
                        help="The raw cache directory.")
    parser.add_argument("--tier2-chemical-min-publications", type=int, default=20,
                        help="The smallest publication count of a second tier chemical.")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                        help="The number of parallel requests to the search endpoint.")
    parser.add_argument("--refresh", action="store_true",
                        help="Fetch every entity again instead of reading the cache.")
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
    wanted = wanted_accessions(nodes)

    path = args.raw_dir / CACHE_NAME
    cached = common.read_json(path) or {}
    counts = {} if args.refresh else dict(cached.get("counts", {}))
    todo = [a for a in wanted if a not in counts]
    print(f"  {len(wanted)} entities in the graph, {len(todo)} to count now")

    with common.make_client() as client:
        corpus_size, corpus_date = fetch_corpus_size(client)
        if corpus_size <= 0:
            print("FAIL: the NCBI EInfo endpoint returned no PubMed record count.")
            return 1
        print(f"  corpus size {corpus_size} records, updated {corpus_date}")
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            results = pool.map(fetch_count, repeat(client), todo)
            try:
                for number, (accession, count) in enumerate(zip(todo, results), start=1):
                    counts[accession] = count
                    if number % SAVE_EVERY == 0:
                        common.write_json(path, {"counts": counts})
                        print(f"    {number}/{len(todo)} entities counted", flush=True)
            except httpx.HTTPError as error:
                common.write_json(path, {"counts": counts})
                print(f"  the endpoint stopped answering after {len(counts)} entities: {error}")
                print("  the cache holds the counts so far. Run the script again to continue.")

    non_zero = sum(1 for a in wanted if counts.get(a, 0) > 0)
    payload = {"fetched_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
               "endpoint": SEARCH_URL, "corpus_endpoint": EINFO_URL,
               "corpus_size": corpus_size, "corpus_date": corpus_date,
               "entities_requested": len(wanted), "entities_with_a_count": non_zero,
               "counts": counts}
    common.write_json(path, payload)
    print(f"{non_zero} of {len(wanted)} entities carry a publication count. Cache: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
