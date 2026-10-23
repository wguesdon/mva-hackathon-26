"""Fetch the compounds and the Reactome pathways of every gene in the Track 2 graph.

The source is the Open Targets Platform GraphQL API. Two calls run for each batch of genes.
The first call maps a gene symbol to an Ensembl gene identifier. The second call reads the
`drugAndClinicalCandidates` field and the `pathways` field of each target.

The `drugAndClinicalCandidates` field replaced the older `knownDrugs` field. A row gives the
ChEMBL identifier of the compound, the compound name, the mechanism of action, the action type
and the highest clinical stage that the compound reached against that target.

Usage:
    uv run python scripts/fetch_open_targets_drugs.py
    uv run python scripts/fetch_open_targets_drugs.py --refresh
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402

GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"
CACHE_NAME = "open_targets_targets.json"
BATCH_SIZE = 20

MAP_QUERY = """
query MapSymbols($terms: [String!]!) {
  mapIds(queryTerms: $terms, entityNames: ["target"]) {
    mappings { term hits { id name entity } }
  }
}
"""

TARGET_QUERY = """
query TargetEvidence($ids: [String!]!) {
  targets(ensemblIds: $ids) {
    id
    approvedSymbol
    pathways { pathwayId pathway topLevelTerm }
    drugAndClinicalCandidates {
      count
      rows {
        maxClinicalStage
        drug {
          id
          name
          drugType
          maximumClinicalStage
          mechanismsOfAction { rows { actionType mechanismOfAction targetName } }
        }
      }
    }
  }
}
"""


def batched(items: list[str], size: int) -> list[list[str]]:
    """Split a list into batches of a fixed size.

    Args:
        items: The values to split.
        size: The largest batch size.

    Returns:
        The batches in order.
    """
    return [items[start:start + size] for start in range(0, len(items), size)]


def parse_mappings(payload: dict) -> dict[str, str]:
    """Read the symbol to Ensembl identifier mapping from a `mapIds` response.

    A term keeps a hit only when the hit name equals the term, so a partner gene of a read
    through locus never takes the place of the gene that was asked for.

    Args:
        payload: The decoded GraphQL response.

    Returns:
        The Ensembl gene identifier of each symbol that resolved.
    """
    found: dict[str, str] = {}
    mappings = (payload.get("data") or {}).get("mapIds") or {}
    for mapping in mappings.get("mappings", []) or []:
        term = mapping.get("term", "")
        for hit in mapping.get("hits", []) or []:
            if hit.get("entity") == "target" and hit.get("name", "").upper() == term.upper():
                found[term] = hit.get("id", "")
                break
    return {term: identifier for term, identifier in found.items() if identifier}


def action_of(drug: dict) -> tuple[str, str]:
    """Read the action type and the mechanism text of a compound.

    Args:
        drug: The `drug` object of one candidate row.

    Returns:
        The action type in lower case and the mechanism text. Both are empty when the
        compound carries no mechanism record.
    """
    mechanisms = (drug.get("mechanismsOfAction") or {}).get("rows") or []
    if not mechanisms:
        return "", ""
    first = mechanisms[0]
    return (first.get("actionType") or "").lower(), first.get("mechanismOfAction") or ""


def parse_targets(payload: dict) -> list[dict]:
    """Flatten a `targets` response into one record per gene.

    Args:
        payload: The decoded GraphQL response.

    Returns:
        One record per target, with its pathways and its compound rows.
    """
    out: list[dict] = []
    targets = (payload.get("data") or {}).get("targets") or []
    for target in targets:
        candidates = target.get("drugAndClinicalCandidates") or {}
        drugs = []
        for row in candidates.get("rows", []) or []:
            drug = row.get("drug") or {}
            action, mechanism = action_of(drug)
            drugs.append({"drug_id": drug.get("id", ""), "drug_name": drug.get("name", ""),
                          "drug_type": drug.get("drugType", "") or "",
                          "action_type": action, "mechanism_of_action": mechanism,
                          "max_clinical_stage_on_target": row.get("maxClinicalStage", "") or "",
                          "max_clinical_stage_overall": drug.get("maximumClinicalStage", "") or ""})
        pathways = [{"pathway_id": p.get("pathwayId", ""), "pathway": p.get("pathway", ""),
                     "top_level_term": p.get("topLevelTerm", "")}
                    for p in target.get("pathways", []) or []]
        out.append({"ensembl_id": target.get("id", ""), "symbol": target.get("approvedSymbol", ""),
                    "drug_count": int(candidates.get("count", 0) or 0),
                    "drugs": drugs, "pathways": pathways})
    return out


def main() -> int:
    """Resolve every graph gene on Open Targets and cache its compounds and pathways.

    Returns:
        0 on success, 1 when the relation cache is missing or the API returned an error.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--raw-dir", type=Path, default=common.RAW_DIR,
                        help="The raw cache directory.")
    parser.add_argument("--tier2-chemical-min-publications", type=int, default=20,
                        help="The smallest publication count of a second tier chemical.")
    parser.add_argument("--refresh", action="store_true",
                        help="Query every gene again instead of reading the cache.")
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
    symbols = sorted({common.split_accession(a)[1] for a in nodes
                      if common.accession_kind(a) == "gene"})

    path = args.raw_dir / CACHE_NAME
    cached = common.read_json(path) or {}
    targets = {} if args.refresh else dict(cached.get("targets", {}))
    todo = [s for s in symbols if s not in targets]
    print(f"  {len(symbols)} genes in the graph, {len(todo)} to query now")

    mapped_total = 0
    with common.make_client() as client:
        for batch in batched(todo, BATCH_SIZE):
            mapping = parse_mappings(common.post_json(
                client, GRAPHQL_URL, {"query": MAP_QUERY, "variables": {"terms": batch}}))
            time.sleep(common.DEFAULT_PAUSE_S)
            for symbol in batch:
                if symbol not in mapping:
                    targets[symbol] = {"symbol": symbol, "ensembl_id": "", "drug_count": 0,
                                       "drugs": [], "pathways": [], "resolved": False}
            ids = sorted(set(mapping.values()))
            if not ids:
                continue
            payload = common.post_json(client, GRAPHQL_URL,
                                       {"query": TARGET_QUERY, "variables": {"ids": ids}})
            if payload.get("errors"):
                print(f"FAIL: Open Targets returned {payload['errors'][0].get('message', '')}")
                return 1
            by_symbol = {record["symbol"]: record for record in parse_targets(payload)}
            for symbol, ensembl_id in mapping.items():
                record = by_symbol.get(symbol)
                if record is None:
                    targets[symbol] = {"symbol": symbol, "ensembl_id": ensembl_id,
                                       "drug_count": 0, "drugs": [], "pathways": [],
                                       "resolved": False}
                    continue
                record["resolved"] = True
                targets[symbol] = record
                mapped_total += 1
            print(f"    {len(targets)}/{len(symbols)} genes cached")
            time.sleep(common.DEFAULT_PAUSE_S)

    with_drugs = sum(1 for s in symbols if targets.get(s, {}).get("drugs"))
    drug_ids = {d["drug_id"] for s in symbols for d in targets.get(s, {}).get("drugs", [])}
    payload = {"fetched_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
               "endpoint": GRAPHQL_URL, "genes_requested": len(symbols),
               "genes_with_a_compound": with_drugs, "distinct_compounds": len(drug_ids),
               "targets": targets}
    common.write_json(path, payload)
    print(f"{with_drugs} of {len(symbols)} genes have a compound. "
          f"{len(drug_ids)} distinct compounds. Cache: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
