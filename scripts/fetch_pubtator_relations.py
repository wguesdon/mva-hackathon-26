"""Fetch the PubTator3 relations of the MVA genes and of the neighbours that they reveal.

The script runs in two rounds. The first round fetches the relations of each seed gene. The
second round reads the gene neighbours out of those relations and fetches the relations of each
neighbour. The neighbour list is never hardcoded. It comes from the publication counts that the
first round returned.

The endpoint is `https://www.ncbi.nlm.nih.gov/research/pubtator3-api/relations?e1=@GENE_<SYMBOL>`.
It answers with a list of records. Each record holds a relation type, a source accession, a target
accession and a publication count.

Usage:
    uv run python scripts/fetch_pubtator_relations.py
    uv run python scripts/fetch_pubtator_relations.py --refresh --neighbour-min-publications 5
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402

RELATIONS_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/relations"
SEED_GENES = ("BUB1B", "CEP57", "TRIP13")
OUT_SUBDIR = "pubtator_relations"
INDEX_NAME = "pubtator_relations_index.json"


def relations_path(raw_dir: Path, symbol: str) -> Path:
    """Return the cache path of one gene.

    Args:
        raw_dir: The raw cache directory.
        symbol: The gene symbol.

    Returns:
        The JSON file that holds the raw response.
    """
    return raw_dir / OUT_SUBDIR / f"{symbol}.json"


def fetch_gene(client, raw_dir: Path, symbol: str, refresh: bool) -> list[dict]:
    """Fetch the relations of one gene, or read them from the cache.

    Args:
        client: The HTTP client.
        raw_dir: The raw cache directory.
        symbol: The gene symbol.
        refresh: True to ignore the cache and call the endpoint.

    Returns:
        The raw relation records.
    """
    path = relations_path(raw_dir, symbol)
    if not refresh:
        cached = common.read_json(path)
        if isinstance(cached, list):
            return cached
    payload = common.get_json(client, RELATIONS_URL, {"e1": f"@GENE_{symbol}"})
    records = payload if isinstance(payload, list) else []
    common.write_json(path, records)
    time.sleep(common.DEFAULT_PAUSE_S)
    return records


def gene_neighbour_counts(records_by_gene: dict[str, list[dict]],
                          seeds: tuple[str, ...] | list[str]) -> Counter:
    """Sum the publication counts of every gene that shares a relation with a seed gene.

    Args:
        records_by_gene: The relation records of each seed gene, keyed by symbol.
        seeds: The seed gene symbols.

    Returns:
        A counter of gene accession to the total publication count over the seed genes.
    """
    seed_accessions = {f"@GENE_{symbol}" for symbol in seeds}
    totals: Counter = Counter()
    for symbol in seeds:
        for record in records_by_gene.get(symbol, []):
            publications = int(record.get("publications", 0) or 0)
            for side in ("source", "target"):
                accession = record.get(side, "")
                if accession.startswith("@GENE_") and accession not in seed_accessions:
                    totals[accession] += publications
    return totals


def derive_neighbours(records_by_gene: dict[str, list[dict]], seeds: tuple[str, ...] | list[str],
                      min_publications: int, max_neighbours: int) -> list[str]:
    """Select the gene neighbours that the seed relations reveal.

    A neighbour qualifies when the publication counts that link it to the seed genes reach the
    threshold. The selection sorts by that total, then by symbol, so the result is stable.

    Args:
        records_by_gene: The relation records of each seed gene, keyed by symbol.
        seeds: The seed gene symbols.
        min_publications: The smallest total publication count that a neighbour must reach.
        max_neighbours: The largest number of neighbours to keep.

    Returns:
        The neighbour gene symbols in rank order.
    """
    totals = gene_neighbour_counts(records_by_gene, seeds)
    ranked = sorted(totals.items(), key=lambda item: (-item[1], item[0]))
    chosen = [accession for accession, total in ranked if total >= min_publications]
    symbols = [common.split_accession(accession)[1] for accession in chosen[:max_neighbours]]
    return [symbol for symbol in symbols if symbol]


def main() -> int:
    """Fetch the seed relations, derive the neighbours and fetch the neighbour relations.

    Returns:
        0 on success, 1 when a seed gene returned no relation.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--genes", nargs="+", default=list(SEED_GENES),
                        help="The seed gene symbols.")
    parser.add_argument("--neighbour-min-publications", type=int, default=5,
                        help="The smallest total publication count of a derived neighbour.")
    parser.add_argument("--max-neighbours", type=int, default=40,
                        help="The largest number of derived neighbours to fetch.")
    parser.add_argument("--raw-dir", type=Path, default=common.RAW_DIR,
                        help="The raw cache directory.")
    parser.add_argument("--refresh", action="store_true",
                        help="Ignore the cache and call the endpoint again.")
    args = parser.parse_args()

    records_by_gene: dict[str, list[dict]] = {}
    with common.make_client() as client:
        for symbol in args.genes:
            records_by_gene[symbol] = fetch_gene(client, args.raw_dir, symbol, args.refresh)
            print(f"  seed {symbol:<10} {len(records_by_gene[symbol]):>5} relations")
        empty = [symbol for symbol in args.genes if not records_by_gene[symbol]]
        if empty:
            print(f"FAIL: no relation returned for {', '.join(empty)}")
            return 1
        neighbours = derive_neighbours(records_by_gene, args.genes,
                                       args.neighbour_min_publications, args.max_neighbours)
        print(f"  derived {len(neighbours)} neighbour genes from the seed relations")
        for symbol in neighbours:
            records_by_gene[symbol] = fetch_gene(client, args.raw_dir, symbol, args.refresh)
            print(f"  neighbour {symbol:<10} {len(records_by_gene[symbol]):>5} relations")

    index = {
        "fetched_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "endpoint": RELATIONS_URL,
        "seed_genes": list(args.genes),
        "neighbour_min_publications": args.neighbour_min_publications,
        "max_neighbours": args.max_neighbours,
        "derived_neighbours": neighbours,
        "relation_counts": {symbol: len(records) for symbol, records in records_by_gene.items()},
    }
    path = common.write_json(args.raw_dir / INDEX_NAME, index)
    total = sum(len(records) for records in records_by_gene.values())
    print(f"{len(records_by_gene)} genes, {total} relation records. Index: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
