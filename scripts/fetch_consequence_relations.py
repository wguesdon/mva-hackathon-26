"""Select the downstream consequences of the MVA genes and fetch what acts on them.

No approved compound binds BUB1B, CEP57 or TRIP13, so a search that starts at the gene has no
answer. The answerable question starts at the effect. BUBR1 loss produces aneuploidy and
chromosomal instability, and the literature does record compounds that act on those states.

The consequence set is derived, not named. The script takes every disease entity that PubTator3
returns for the seed genes, keeps the entities that several seed genes share, and ranks them by
normalized pointwise mutual information against the seed genes. A cancer term such as
`@DISEASE_Neoplasms` appears beside every gene in the literature and scores near zero. A term that
describes the mechanism scores high. The threshold is a command line argument, so the rule is
visible and a reader can change it.

The script then calls the relations endpoint with each kept consequence as `e1`. That call returns
relation types that a gene query never returns, including `cause` and `treat`. Those are the edges
that let a ranker travel from the disease state to a compound.

Usage:
    uv run python scripts/fetch_consequence_relations.py
    uv run python scripts/fetch_consequence_relations.py --min-npmi 0.2 --refresh
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402

RELATIONS_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/relations"
OUT_SUBDIR = "pubtator_consequences"
CACHE_NAME = "pubtator_consequences.json"


def seed_disease_links(records_by_gene: dict[str, list[dict]],
                       seeds: list[str]) -> dict[str, dict[str, int]]:
    """Collect the publication count that joins each disease entity to each seed gene.

    Args:
        records_by_gene: The relation records of each gene, keyed by symbol.
        seeds: The seed gene symbols.

    Returns:
        For each disease accession, the publication count against each seed gene that links to it.
    """
    links: dict[str, dict[str, int]] = {}
    for symbol in seeds:
        accession = f"@GENE_{symbol}"
        for record in records_by_gene.get(symbol, []):
            source, target = record.get("source", ""), record.get("target", "")
            if accession not in (source, target):
                continue
            other = target if source == accession else source
            if not other.startswith("@DISEASE_"):
                continue
            counts = links.setdefault(other, {})
            counts[symbol] = counts.get(symbol, 0) + int(record.get("publications", 0) or 0)
    return links


def score_consequences(links: dict[str, dict[str, int]], counts: dict[str, int],
                       seeds: list[str], corpus_size: int, min_seeds: int) -> list[dict]:
    """Score each candidate disease entity by its specificity to the seed genes.

    The score is the highest normalized pointwise mutual information that the entity reaches
    against any single seed gene. The highest value is used rather than the mean, because CEP57
    and TRIP13 have a much smaller literature than BUB1B and a mean would hide a real tie.

    Args:
        links: The publication count of each disease against each seed gene.
        counts: The total publication count of each entity.
        seeds: The seed gene symbols.
        corpus_size: The number of records in the corpus.
        min_seeds: The number of seed genes a candidate must link to.

    Returns:
        One record per candidate, sorted by score, highest first.
    """
    scored: list[dict] = []
    for accession, per_seed in links.items():
        if len(per_seed) < min_seeds:
            continue
        best = None
        best_seed = ""
        for symbol, pair_count in per_seed.items():
            value = common.npmi(pair_count, counts.get(accession, 0),
                                counts.get(f"@GENE_{symbol}", 0), corpus_size)
            if value is not None and (best is None or value > best):
                best = value
                best_seed = symbol
        if best is None:
            continue
        scored.append({"accession": accession, "label": common.accession_label(accession),
                       "npmi": round(best, 4), "best_seed": best_seed,
                       "seed_count": len(per_seed),
                       "publications": sum(per_seed.values()),
                       "entity_publications": counts.get(accession, 0)})
    return sorted(scored, key=lambda r: (-r["npmi"], r["accession"]))


def select_consequences(scored: list[dict], min_npmi: float, min_publications: int,
                        limit: int) -> list[dict]:
    """Keep the candidates that reach both thresholds.

    The specificity threshold removes a term that appears beside every gene. The publication
    threshold removes a term that shares two papers with two seed genes by coincidence, which
    scores well only because the term itself is rare.

    Args:
        scored: The scored candidates.
        min_npmi: The smallest score that a consequence must reach.
        min_publications: The smallest co-mention count against the seed genes.
        limit: The largest number of consequences to keep.

    Returns:
        The kept records in score order.
    """
    return [record for record in scored
            if record["npmi"] >= min_npmi
            and record["publications"] >= min_publications][:limit]


def fetch_consequence(client, raw_dir: Path, accession: str, refresh: bool) -> list[dict]:
    """Fetch the relations of one consequence entity, or read them from the cache.

    Args:
        client: The HTTP client.
        raw_dir: The raw cache directory.
        accession: The consequence accession, for example `@DISEASE_Aneuploidy`.
        refresh: True to ignore the cache and call the endpoint.

    Returns:
        The raw relation records.
    """
    _, body = common.split_accession(accession)
    path = raw_dir / OUT_SUBDIR / f"{body}.json"
    if not refresh:
        cached = common.read_json(path)
        if isinstance(cached, list):
            return cached
    payload = common.get_json(client, RELATIONS_URL, {"e1": accession})
    records = payload if isinstance(payload, list) else []
    common.write_json(path, records)
    time.sleep(common.DEFAULT_PAUSE_S)
    return records


def main() -> int:
    """Select the consequence entities and fetch the relations of each one.

    Returns:
        0 on success, 1 when a required cache file is missing.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--raw-dir", type=Path, default=common.RAW_DIR,
                        help="The raw cache directory.")
    parser.add_argument("--tier2-chemical-min-publications", type=int, default=20,
                        help="The smallest publication count of a second tier chemical.")
    parser.add_argument("--min-seeds", type=int, default=2,
                        help="The number of seed genes a consequence must link to.")
    parser.add_argument("--min-npmi", type=float, default=0.07,
                        help="The smallest specificity score of a consequence.")
    parser.add_argument("--min-publications", type=int, default=5,
                        help="The smallest co-mention count against the seed genes.")
    parser.add_argument("--limit", type=int, default=12,
                        help="The largest number of consequence entities to keep.")
    parser.add_argument("--refresh", action="store_true",
                        help="Fetch every consequence again instead of reading the cache.")
    args = parser.parse_args()

    index = common.read_json(args.raw_dir / "pubtator_relations_index.json")
    count_cache = common.read_json(args.raw_dir / "pubtator_entity_counts.json") or {}
    if not index:
        print("FAIL: run scripts/fetch_pubtator_relations.py first.")
        return 1
    if not count_cache.get("counts"):
        print("FAIL: run scripts/fetch_pubtator_entity_counts.py first.")
        return 1
    seeds = index["seed_genes"]
    neighbours = index["derived_neighbours"]
    records = common.load_relation_cache(args.raw_dir, seeds + neighbours, "pubtator_relations")
    links = seed_disease_links(records, seeds)
    scored = score_consequences(links, count_cache["counts"], seeds,
                                count_cache["corpus_size"], args.min_seeds)
    chosen = select_consequences(scored, args.min_npmi, args.min_publications,
                                 args.limit)
    print(f"  {len(links)} disease entities link to a seed gene, "
          f"{len(scored)} link to {args.min_seeds} or more, {len(chosen)} pass the threshold")
    for record in chosen:
        print(f"    npmi {record['npmi']:+.3f}  seeds {record['seed_count']}  "
              f"{record['label']}")

    fetched: dict[str, int] = {}
    with common.make_client() as client:
        for record in chosen:
            rows = fetch_consequence(client, args.raw_dir, record["accession"], args.refresh)
            fetched[record["accession"]] = len(rows)
            print(f"    {record['label']:<40} {len(rows):>5} relations")

    payload = {"fetched_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
               "endpoint": RELATIONS_URL, "min_seeds": args.min_seeds,
               "min_npmi": args.min_npmi, "min_publications": args.min_publications,
               "limit": args.limit,
               "corpus_size": count_cache["corpus_size"],
               "candidates_scored": len(scored), "consequences": chosen,
               "relation_counts": fetched, "all_candidates": scored}
    path = common.write_json(args.raw_dir / CACHE_NAME, payload)
    total = sum(fetched.values())
    print(f"{len(chosen)} consequence entities, {total} relation records. Cache: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
