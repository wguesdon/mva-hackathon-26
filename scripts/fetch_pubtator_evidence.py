"""Fetch the PMIDs and the evidence sentence of each entity pair in the literature layer.

The relations endpoint gives a publication count but no PMID. The search endpoint gives both.
A query of the form `relations:<type>|<accession_1>|<accession_2>` returns one result per
publication. Each result carries the PMID and a `text_hl` field. That field holds the sentence
that the two entities share, with the PubTator3 annotation markup around each mention. The
markup is removed here, so the graph keeps the sentence as the authors wrote it.

The search endpoint pages ten results at a time. The script stops at `--max-pages` pages per
pair, so a pair with many publications keeps its strongest evidence without a long run.

A fetch of every pair is not practical, so the order matters. A rare disease has a small
literature, and a plain publication count would rank every MVA edge below the cancer edges of the
hub genes. The script therefore fetches first each pair that touches a seed gene or a disease that
every seed gene shares. See `priority_accessions`.

Usage:
    uv run python scripts/fetch_pubtator_evidence.py
    uv run python scripts/fetch_pubtator_evidence.py --max-pairs 2000 --max-pages 3 --refresh
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from itertools import repeat
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402

SEARCH_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/"
CACHE_NAME = "pubtator_evidence.json"
SAVE_EVERY = 100

# NCBI asks for no more than three requests a second without an API key. Four workers on an
# endpoint that answers in about two seconds stays under that rate.
DEFAULT_WORKERS = 4


def pair_id(record: dict) -> str:
    """Build the cache key of one relation record.

    Args:
        record: A relation record with a type, a source and a target.

    Returns:
        The key, for example `associate|@GENE_BUB1B|@DISEASE_Aneuploidy`.
    """
    return f"{record.get('type', '')}|{record.get('source', '')}|{record.get('target', '')}"


def parse_results(payload: dict) -> list[dict]:
    """Read the PMID and the plain sentence of every result on one search page.

    Args:
        payload: The decoded search response.

    Returns:
        One record per result, with the PMID, the sentence and the title.
    """
    out: list[dict] = []
    for result in payload.get("results", []) or []:
        pmid = str(result.get("pmid", "") or "")
        if not pmid:
            continue
        out.append({"pmid": pmid, "sentence": common.plain_sentence(result.get("text_hl", "")),
                    "title": result.get("title", ""), "year": str(result.get("date", ""))[:4]})
    return out


def merge_identifiers(payload: dict, found: dict[str, str]) -> None:
    """Collect the database identifier of each marked entity on one search page.

    Args:
        payload: The decoded search response.
        found: The mapping of accession to identifier to update in place.
    """
    for result in payload.get("results", []) or []:
        for accession, identifier in common.highlight_identifiers(result.get("text_hl", "")).items():
            found.setdefault(accession, identifier)


def fetch_pair(client, record: dict, max_pages: int) -> dict:
    """Retrieve the publications of one entity pair.

    Args:
        client: The HTTP client.
        record: The relation record.
        max_pages: The largest number of result pages to read.

    Returns:
        A record with the reported total, the evidence list and the identifiers that the
        highlights revealed.
    """
    text = common.pair_query(record.get("source", ""), record.get("target", ""),
                             record.get("type", "ANY"))
    evidence: list[dict] = []
    identifiers: dict[str, str] = {}
    total = 0
    pages = 1
    page = 1
    while page <= min(max_pages, pages):
        payload = common.get_json(client, SEARCH_URL, {"text": text, "page": page})
        if not isinstance(payload, dict):
            break
        total = int(payload.get("count", 0) or 0)
        pages = int(payload.get("total_pages", 1) or 1)
        evidence.extend(parse_results(payload))
        merge_identifiers(payload, identifiers)
        page += 1
        time.sleep(common.DEFAULT_PAUSE_S)
    return {"query": text, "search_count": total, "pages_read": page - 1,
            "total_pages": pages, "evidence": evidence, "identifiers": identifiers}


def priority_accessions(records_by_gene: dict[str, list[dict]], seeds: list[str],
                        min_seed_links: int = 2) -> set[str]:
    """Choose the entities whose edges must carry provenance whatever their publication count.

    A plain publication count ranks a hub gene of cancer biology above every MVA edge, because
    the MVA literature is small. The disease under study would then reach the graph with no PMID,
    which is the one thing the design forbids. The priority set holds the seed genes and any
    disease that several seed genes share. The disease is derived from the relations, not named
    here, so the rule still works if the seed list changes.

    Args:
        records_by_gene: The relation records of each gene, keyed by symbol.
        seeds: The seed gene symbols.
        min_seed_links: The number of seed genes a disease must link to. The caller passes the
            seed count, so by default only a disease that every seed gene shares ranks first.

    Returns:
        The accessions whose edges are fetched first.
    """
    chosen = {f"@GENE_{symbol}" for symbol in seeds}
    seeds_of_disease: dict[str, set[str]] = {}
    for symbol in seeds:
        for record in records_by_gene.get(symbol, []):
            for side in ("source", "target"):
                accession = record.get(side, "")
                if accession.startswith("@DISEASE_"):
                    seeds_of_disease.setdefault(accession, set()).add(symbol)
    for accession, linked in seeds_of_disease.items():
        if len(linked) >= min_seed_links:
            chosen.add(accession)
    return chosen


def rank_records(edges: list[dict], priority: set[str] | None = None) -> list[dict]:
    """Order the relation records so that the pair that matters most comes first.

    A priority pair comes before every other pair. Inside each group the best supported pair
    comes first, then the key, so the order is stable between runs.

    Args:
        edges: The selected relation records.
        priority: The accessions that rank first, from `priority_accessions`.

    Returns:
        The records in fetch order.
    """
    chosen = priority or set()
    return sorted(edges, key=lambda r: (
        0 if (r.get("source", "") in chosen or r.get("target", "") in chosen) else 1,
        -int(r.get("publications", 0) or 0), pair_id(r)))


def main() -> int:
    """Fetch the evidence of the best supported entity pairs and cache the result.

    Returns:
        0 on success, 1 when the relation cache is missing.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--raw-dir", type=Path, default=common.RAW_DIR,
                        help="The raw cache directory.")
    parser.add_argument("--tier2-chemical-min-publications", type=int, default=20,
                        help="The smallest publication count of a second tier chemical.")
    parser.add_argument("--min-seed-links", type=int, default=0,
                        help="The number of seed genes a disease must share to rank first. "
                             "0 means every seed gene.")
    parser.add_argument("--max-pairs", type=int, default=1400,
                        help="The largest number of entity pairs to fetch evidence for.")
    parser.add_argument("--max-pages", type=int, default=3,
                        help="The largest number of search result pages per pair.")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                        help="The number of parallel requests to the search endpoint.")
    parser.add_argument("--refresh", action="store_true",
                        help="Fetch every pair again instead of reading the cache.")
    args = parser.parse_args()

    index = common.read_json(args.raw_dir / "pubtator_relations_index.json")
    if not index:
        print("FAIL: run scripts/fetch_pubtator_relations.py first.")
        return 1
    seeds = index["seed_genes"]
    neighbours = index["derived_neighbours"]
    records = common.load_relation_cache(args.raw_dir, seeds + neighbours, "pubtator_relations")
    _, edges = common.select_graph(records, seeds, neighbours,
                                   args.tier2_chemical_min_publications)
    consequence_edges = common.consequence_records(args.raw_dir,
                                                   args.tier2_chemical_min_publications)
    edges = edges + consequence_edges
    priority = priority_accessions(records, seeds, args.min_seed_links or len(seeds))
    priority.update(common.consequence_accessions(args.raw_dir))
    wanted = rank_records(edges, priority)[:args.max_pairs]

    path = args.raw_dir / CACHE_NAME
    cached = common.read_json(path) or {}
    pairs = {} if args.refresh else dict(cached.get("pairs", {}))
    todo = [r for r in wanted if pair_id(r) not in pairs]
    in_priority = sum(1 for r in wanted
                      if r.get("source", "") in priority or r.get("target", "") in priority)
    print(f"  {len(edges)} selected relations, {len(priority)} priority entities, "
          f"{len(wanted)} pairs requested ({in_priority} of them priority), "
          f"{len(todo)} to fetch now")

    with common.make_client() as client:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            results = pool.map(fetch_pair, repeat(client), todo, repeat(args.max_pages))
            for number, (record, result) in enumerate(zip(todo, results), start=1):
                pairs[pair_id(record)] = result
                if number % SAVE_EVERY == 0:
                    common.write_json(path, {"pairs": pairs})
                    print(f"    {number}/{len(todo)} pairs fetched", flush=True)

    with_pmids = sum(1 for key in pairs if pairs[key]["evidence"])
    with_sentence = sum(1 for key in pairs
                        for item in pairs[key]["evidence"] if item["sentence"])
    payload = {"fetched_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
               "endpoint": SEARCH_URL, "max_pages": args.max_pages,
               "pairs_requested": len(wanted), "priority_pairs": in_priority,
               "pairs_with_pmids": with_pmids,
               "evidence_items_with_sentence": with_sentence, "pairs": pairs}
    common.write_json(path, payload)
    print(f"{with_pmids} of {len(pairs)} pairs carry a PMID. "
          f"{with_sentence} evidence items carry a sentence. Cache: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
