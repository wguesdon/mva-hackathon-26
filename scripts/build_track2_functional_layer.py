"""Build an optional functional interaction layer for the Track 2 knowledge graph from STRING.

The default graph derives its gene neighbours and its gene to gene edges from literature
co-mention. This script builds a second gene to gene layer from STRING, which measures function
rather than citation. The layer is an option of the ranker (`--gene-layer`). It never replaces the
default graph, and the default graph tables under `reports/track2/graph/` stay unchanged.

The source is STRING v12.0, the current release, for Homo sapiens (taxon 9606). STRING 12.5 is a
preview release on the STRING website, so it is not used. STRING distributes all data under the
Creative Commons BY 4.0 licence (https://string-db.org/cgi/access, section Licensing).

Only two evidence channels enter the score, `experimental` and `database`. The `textmining`
channel is left out, because it counts co-mention in the literature, and that is the signal of the
default layer. The two channels are combined with the STRING rule: the prior 0.041 is removed from
each channel score, the channels combine as `1 - prod(1 - s)`, and the prior is added back once.
An edge enters the layer when that score reaches 700, the STRING high confidence threshold.

The script writes these files to `reports/track2/functional/`.

* `nodes.tsv`: the genes that the functional neighbour set adds, and their compounds and pathways.
* `edges.tsv`: the STRING gene to gene edges, plus the compound and pathway edges of the new genes.
* `drug_evidence.tsv`: the Open Targets drug records of the new genes.
* `manifest.json`: the release, the URLs, the SHA-256 of each download, the licence, the parameters,
  the identifier mapping rate and the neighbour sets.

The STRING downloads go to `tracking/cache/string_v12_0/`, which git ignores. The Open Targets
records of the new genes go to `reports/track2/raw/open_targets_functional_targets.json`, also
ignored, so the default Open Targets cache is never edited.

Usage:
    uv run python scripts/build_track2_functional_layer.py
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import sys
import time
from collections import defaultdict
from pathlib import Path

import httpx

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_track2_graph as build  # noqa: E402
import fetch_open_targets_drugs as open_targets  # noqa: E402
import track2_common as common  # noqa: E402

STRING_VERSION = "12.0"
STRING_TAXON = "9606"
STRING_BASE_URL = "https://stringdb-downloads.org/download"
STRING_LICENCE = "CC BY 4.0"
STRING_LICENCE_URL = "https://string-db.org/cgi/access"
STRING_FILES = {
    "links_detailed": f"protein.links.detailed.v{STRING_VERSION}/"
                      f"{STRING_TAXON}.protein.links.detailed.v{STRING_VERSION}.txt.gz",
    "info": f"protein.info.v{STRING_VERSION}/{STRING_TAXON}.protein.info.v{STRING_VERSION}.txt.gz",
    "aliases": f"protein.aliases.v{STRING_VERSION}/"
               f"{STRING_TAXON}.protein.aliases.v{STRING_VERSION}.txt.gz",
}
CACHE_DIR = common.ROOT / "tracking" / "cache" / "string_v12_0"
OUT_DIR = common.ROOT / "reports" / "track2" / "functional"
OT_CACHE_NAME = "open_targets_functional_targets.json"

# The parameters of the layer. They were fixed here before any held-out result of the layer existed,
# and they stay the default of every option, so every receipt of the report reproduces with no
# option. ADR 0014 sweeps `min_score` on a development split through `--min-score`, and the sweep
# writes each layer to its own directory.
CHANNELS = ("experimental", "database")
PRIOR = 0.041
MIN_SCORE = 700
MAX_NEIGHBOURS = 40

RELATION = "functional_interaction"
SOURCE = "string_v12_0_expdb"
NODE_SOURCE = "string_v12_0"
EDGE_HEADER = build.EDGE_HEADER + ["string_expdb_score"]

# The alias tiers of the symbol mapping, in order. The first tier holds the sources that name the
# current human gene symbol. The second tier holds the previous HGNC symbol, because PubTator3 can
# carry a symbol that HGNC has since renamed. A free text synonym source is not used, because it
# can name a protein family.
ALIAS_TIERS = (
    ("current_symbol", ("Ensembl_HGNC_symbol", "Ensembl_HGNC", "BioMart_HUGO",
                        "Ensembl_EntrezGene", "UniProt_GN_Name")),
    ("previous_symbol", ("Ensembl_HGNC_prev_symbol",)),
)


def combine_channels(scores: list[int], prior: float = PRIOR) -> int:
    """Combine STRING channel scores into one score with the STRING prior correction.

    Args:
        scores: The channel scores as integers from 0 to 1000.
        prior: The prior probability of an interaction that STRING removes from each channel.

    Returns:
        The combined score as an integer from 0 to 1000. A pair with no channel evidence returns 0.
    """
    if not any(scores):
        return 0
    miss = 1.0
    for score in scores:
        value = score / 1000.0
        corrected = max(0.0, (value - prior) / (1.0 - prior))
        miss *= 1.0 - corrected
    combined = 1.0 - miss
    combined = combined * (1.0 - prior) + prior
    return int(round(combined * 1000))


def download(client: httpx.Client, url: str, path: Path) -> None:
    """Stream one file to disk through a temporary name.

    Args:
        client: The HTTP client.
        url: The source URL.
        path: The destination file.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".part")
    with client.stream("GET", url) as response:
        response.raise_for_status()
        with partial.open("wb") as handle:
            for chunk in response.iter_bytes(1 << 20):
                handle.write(chunk)
    partial.replace(path)


def sha256_stream(path: Path) -> str:
    """Compute the SHA-256 digest of a large file in chunks.

    Args:
        path: The file.

    Returns:
        The hexadecimal digest.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        chunk = handle.read(1 << 20)
        while chunk:
            digest.update(chunk)
            chunk = handle.read(1 << 20)
    return digest.hexdigest()


def fetch_string_files(cache_dir: Path) -> dict[str, dict]:
    """Download the STRING files that are absent from the cache and record their checksums.

    Args:
        cache_dir: The cache directory.

    Returns:
        For each file key, the URL, the local path, the byte size and the SHA-256 digest.
    """
    records: dict[str, dict] = {}
    with httpx.Client(timeout=httpx.Timeout(600, connect=30), follow_redirects=True,
                      headers={"User-Agent": common.USER_AGENT}) as client:
        for key, relative in STRING_FILES.items():
            url = f"{STRING_BASE_URL}/{relative}"
            path = cache_dir / Path(relative).name
            if not path.is_file():
                print(f"  download {url}")
                download(client, url, path)
            records[key] = {"url": url, "path": str(path.relative_to(common.ROOT)),
                            "bytes": path.stat().st_size, "sha256": sha256_stream(path)}
    return records


def read_info(path: Path) -> dict[str, str]:
    """Read the preferred name of every STRING protein.

    Args:
        path: The gzip protein.info file.

    Returns:
        The preferred name of each STRING protein identifier.
    """
    names: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                names[parts[0]] = parts[1]
    return names


def read_symbol_aliases(path: Path, tiers: tuple = ALIAS_TIERS) -> dict[str, dict[str, set[str]]]:
    """Read the gene symbol aliases of every STRING protein, one table per alias tier.

    Args:
        path: The gzip protein.aliases file.
        tiers: The tier names and the alias sources of each tier.

    Returns:
        For each tier name, the STRING protein identifiers of each upper case alias.
    """
    tier_of_source = {source: name for name, sources in tiers for source in sources}
    tables: dict[str, dict[str, set[str]]] = {name: defaultdict(set) for name, _ in tiers}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            for source in parts[2].split(" "):
                if source in tier_of_source:
                    tables[tier_of_source[source]][parts[1].upper()].add(parts[0])
    return tables


def map_symbols(symbols: list[str], names: dict[str, str],
                aliases: dict[str, dict[str, set[str]]]
                ) -> tuple[dict[str, str], dict[str, str], list[str]]:
    """Map graph gene symbols to STRING protein identifiers.

    A symbol maps first through the preferred name. A symbol without a preferred name match maps
    through the alias tiers in order, and a tier maps a symbol only when the alias names exactly one
    protein. An ambiguous alias does not map, so no identifier is guessed.

    Args:
        symbols: The graph gene symbols.
        names: The preferred name of each STRING protein.
        aliases: For each tier name in order, the proteins of each upper case alias.

    Returns:
        The protein of each mapped symbol, the method of each mapped symbol, and the unmapped
        symbols in sorted order.
    """
    by_name: dict[str, list[str]] = defaultdict(list)
    for protein, name in names.items():
        by_name[name.upper()].append(protein)
    mapped: dict[str, str] = {}
    method: dict[str, str] = {}
    unmapped: list[str] = []
    for symbol in sorted(symbols):
        key = symbol.upper()
        if len(by_name.get(key, [])) == 1:
            mapped[symbol] = by_name[key][0]
            method[symbol] = "preferred_name"
            continue
        for tier, table in aliases.items():
            if len(table.get(key, set())) == 1:
                mapped[symbol] = next(iter(table[key]))
                method[symbol] = tier
                break
        else:
            unmapped.append(symbol)
    return mapped, method, unmapped


def read_functional_links(path: Path, channels: tuple[str, ...] = CHANNELS,
                          min_score: int = MIN_SCORE) -> dict[tuple[str, str], int]:
    """Read the STRING pairs whose combined score over the chosen channels reaches the threshold.

    The detailed file lists each pair in both orders. The result keeps one unordered pair.

    Args:
        path: The gzip protein.links.detailed file.
        channels: The channel columns that enter the score.
        min_score: The smallest combined score of a kept pair.

    Returns:
        The combined score of each kept unordered pair of STRING protein identifiers.
    """
    links: dict[tuple[str, str], int] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().split()
        columns = [header.index(channel) for channel in channels]
        for line in handle:
            parts = line.split()
            first, second = parts[0], parts[1]
            if first > second:
                continue
            score = combine_channels([int(parts[column]) for column in columns])
            if score >= min_score:
                links[(first, second)] = score
    return links


def derive_functional_neighbours(links: dict[tuple[str, str], int], seeds: dict[str, str],
                                 names: dict[str, str],
                                 max_neighbours: int = MAX_NEIGHBOURS) -> list[dict]:
    """Select the functional neighbours of the seed genes.

    The rule mirrors the co-mention rule of `scripts/fetch_pubtator_relations.py`. A gene scores the
    sum of its STRING scores over the seed genes, and the selection sorts by that sum and then by
    symbol, so the result is stable.

    Args:
        links: The kept STRING pairs and their scores.
        seeds: The STRING protein of each seed symbol.
        names: The preferred name of each STRING protein.
        max_neighbours: The largest number of neighbours to keep.

    Returns:
        One record per neighbour with the symbol, the protein, the summed score and the seeds.
    """
    seed_proteins = {protein: symbol for symbol, protein in seeds.items()}
    totals: dict[str, int] = defaultdict(int)
    partners: dict[str, set[str]] = defaultdict(set)
    for (first, second), score in links.items():
        for mine, other in ((first, second), (second, first)):
            if mine in seed_proteins and other not in seed_proteins:
                totals[other] += score
                partners[other].add(seed_proteins[mine])
    ranked = sorted(totals.items(), key=lambda item: (-item[1], names.get(item[0], item[0])))
    return [{"symbol": names.get(protein, protein), "string_id": protein, "score_sum": total,
             "seeds": sorted(partners[protein])}
            for protein, total in ranked[:max_neighbours]]


def functional_edges(links: dict[tuple[str, str], int], symbol_of: dict[str, str]) -> list[dict]:
    """Build the gene to gene edges of the layer between the genes of the universe.

    Args:
        links: The kept STRING pairs and their scores.
        symbol_of: The graph symbol of each STRING protein of the universe.

    Returns:
        The edge rows. The evidence count is 1, the convention of the other database edges of the
        graph, and the STRING score sits in its own column.
    """
    rows = []
    for (first, second), score in sorted(links.items()):
        if first not in symbol_of or second not in symbol_of:
            continue
        source_id, target_id = sorted((f"GENE:{symbol_of[first]}", f"GENE:{symbol_of[second]}"))
        if source_id == target_id:
            continue
        rows.append({"source_id": source_id, "target_id": target_id, "relation": RELATION,
                     "evidence_count": 1, "specificity_npmi": "", "pmids": "", "source": SOURCE,
                     "string_expdb_score": score})
    return rows


def fetch_new_targets(symbols: list[str], cache_path: Path) -> dict:
    """Fetch the Open Targets compounds and pathways of the genes that the layer adds.

    Args:
        symbols: The new gene symbols.
        cache_path: The separate cache file of these genes.

    Returns:
        The target records, keyed by symbol.
    """
    cached = common.read_json(cache_path) or {}
    targets = dict(cached.get("targets", {}))
    todo = [symbol for symbol in symbols if symbol not in targets]
    with common.make_client() as client:
        for batch in open_targets.batched(todo, open_targets.BATCH_SIZE):
            mapping = open_targets.parse_mappings(common.post_json(
                client, open_targets.GRAPHQL_URL,
                {"query": open_targets.MAP_QUERY, "variables": {"terms": batch}}))
            time.sleep(common.DEFAULT_PAUSE_S)
            for symbol in batch:
                targets.setdefault(symbol, {"symbol": symbol, "ensembl_id": "", "drug_count": 0,
                                            "drugs": [], "pathways": [], "resolved": False})
            ids = sorted(set(mapping.values()))
            if not ids:
                continue
            payload = common.post_json(client, open_targets.GRAPHQL_URL,
                                       {"query": open_targets.TARGET_QUERY,
                                        "variables": {"ids": ids}})
            if payload.get("errors"):
                raise RuntimeError(f"Open Targets returned {payload['errors'][0].get('message')}")
            by_symbol = {record["symbol"]: record for record in open_targets.parse_targets(payload)}
            for symbol in mapping:
                if symbol in by_symbol:
                    by_symbol[symbol]["resolved"] = True
                    targets[symbol] = by_symbol[symbol]
            time.sleep(common.DEFAULT_PAUSE_S)
    common.write_json(cache_path, {
        "fetched_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "endpoint": open_targets.GRAPHQL_URL, "targets": targets})
    return targets


def read_graph_genes(graph_dir: Path) -> tuple[list[str], dict[str, dict]]:
    """Read the gene symbols and every node row of the default graph.

    Args:
        graph_dir: The default graph directory.

    Returns:
        The gene symbols and the node rows keyed by node identifier.
    """
    nodes: dict[str, dict] = {}
    with (graph_dir / "nodes.tsv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            nodes[row["node_id"]] = dict(row)
    genes = sorted(node_id.split(":", 1)[1] for node_id, row in nodes.items()
                   if row["kind"] == "gene")
    return genes, nodes


def jaccard(first: set[str], second: set[str]) -> float:
    """Return the Jaccard index of two sets.

    Args:
        first: The first set.
        second: The second set.

    Returns:
        The size of the intersection over the size of the union, or 0 for two empty sets.
    """
    union = first | second
    return len(first & second) / len(union) if union else 0.0


def main() -> int:
    """Download STRING, map the graph genes, build the layer and write its tables.

    Returns:
        0 on success, 1 when an input is missing.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--graph-dir", type=Path, default=common.GRAPH_DIR,
                        help="The default graph directory.")
    parser.add_argument("--cache-dir", type=Path, default=CACHE_DIR,
                        help="The STRING download cache.")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR, help="The layer directory.")
    parser.add_argument("--raw-dir", type=Path, default=common.RAW_DIR,
                        help="The raw cache directory of the Open Targets records.")
    parser.add_argument("--min-score", type=int, default=MIN_SCORE,
                        help="The smallest combined STRING score of an edge of the layer.")
    parser.add_argument("--max-neighbours", type=int, default=MAX_NEIGHBOURS,
                        help="The largest number of functional neighbours of one seed gene.")
    args = parser.parse_args()

    manifest_graph = common.read_json(args.graph_dir / "manifest.json") or {}
    seeds = manifest_graph.get("seed_genes", [])
    comention_neighbours = manifest_graph.get("derived_neighbour_genes", [])
    if not seeds:
        print("FAIL: the default graph manifest names no seed gene.")
        return 1
    downloads = fetch_string_files(args.cache_dir)
    names = read_info(args.cache_dir / Path(STRING_FILES["info"]).name)
    aliases = read_symbol_aliases(args.cache_dir / Path(STRING_FILES["aliases"]).name)
    genes, graph_nodes = read_graph_genes(args.graph_dir)
    mapped, method, unmapped = map_symbols(genes, names, aliases)
    print(f"  mapped {len(mapped)} of {len(genes)} graph genes to STRING v{STRING_VERSION}")

    links = read_functional_links(args.cache_dir / Path(STRING_FILES["links_detailed"]).name,
                                  min_score=args.min_score)
    print(f"  {len(links)} human pairs reach {args.min_score} on {', '.join(CHANNELS)}")
    seed_map = {symbol: mapped[symbol] for symbol in seeds if symbol in mapped}
    neighbours = derive_functional_neighbours(links, seed_map, names,
                                              max_neighbours=args.max_neighbours)

    symbol_of = {protein: symbol for symbol, protein in mapped.items()}
    graph_symbols = set(genes)
    new_genes: list[str] = []
    for record in neighbours:
        if record["string_id"] in symbol_of:
            record["symbol"] = symbol_of[record["string_id"]]
            continue
        if record["symbol"] in graph_symbols:
            continue
        symbol_of[record["string_id"]] = record["symbol"]
        new_genes.append(record["symbol"])
    edges = functional_edges(links, symbol_of)

    targets = fetch_new_targets(sorted(new_genes), args.raw_dir / OT_CACHE_NAME)
    work_nodes = {node_id: dict(row) for node_id, row in graph_nodes.items()}
    for symbol in new_genes:
        work_nodes[f"GENE:{symbol}"] = {"node_id": f"GENE:{symbol}", "kind": "gene",
                                        "label": symbol, "source": NODE_SOURCE}
    drug_edges, phases = build.build_drug_layer(targets, work_nodes)
    pathway_edges = build.build_pathway_layer(targets, work_nodes)
    new_node_rows = sorted((row for node_id, row in work_nodes.items()
                            if node_id not in graph_nodes),
                           key=lambda row: (row["kind"], row["node_id"]))
    all_edges = edges + [dict(row, string_expdb_score="") for row in drug_edges + pathway_edges]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    nodes_path = common.write_tsv(args.out_dir / "nodes.tsv", build.NODE_HEADER, new_node_rows)
    edges_path = common.write_tsv(args.out_dir / "edges.tsv", EDGE_HEADER, all_edges)
    drugs_path = common.write_tsv(args.out_dir / "drug_evidence.tsv", build.DRUG_HEADER, phases)

    functional_set = {record["symbol"] for record in neighbours}
    comention_set = set(comention_neighbours)
    seed_edges = {seed: sorted({row["target_id"] if row["source_id"] == f"GENE:{seed}"
                                else row["source_id"] for row in edges
                                if f"GENE:{seed}" in (row["source_id"], row["target_id"])})
                  for seed in seeds}
    manifest = {
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "source": {"name": "STRING", "version": STRING_VERSION, "taxon": STRING_TAXON,
                   "release_status": "current release since 2023-07-26; 12.5 is a preview",
                   "licence": STRING_LICENCE, "licence_url": STRING_LICENCE_URL,
                   "citation": "Szklarczyk D et al. The STRING database in 2023. Nucleic Acids "
                               "Res 2023;51(D1):D638-D646. PMID 36370105",
                   "downloads": downloads},
        "parameters": {"channels": list(CHANNELS), "prior": PRIOR, "min_score": args.min_score,
                       "max_neighbours": args.max_neighbours, "relation": RELATION,
                       "edge_evidence_count": 1, "alias_tiers": {name: list(sources) for name, sources in ALIAS_TIERS},
                       "combination": "STRING rule: remove the prior from each channel, combine "
                                      "as 1 - prod(1 - s), add the prior back once"},
        "mapping": {"graph_genes": len(genes), "mapped": len(mapped),
                    "mapping_rate": round(len(mapped) / len(genes), 4) if genes else 0.0,
                    "by_method": {name: sum(1 for value in method.values() if value == name)
                                  for name in ["preferred_name"] + [n for n, _ in ALIAS_TIERS]},
                    "alias_mapped": {symbol: value for symbol, value in sorted(method.items())
                                     if value != "preferred_name"},
                    "unmapped": unmapped},
        "seed_genes": seeds,
        "neighbours": {
            "comention": comention_neighbours,
            "functional": neighbours,
            "functional_not_in_default_graph": sorted(new_genes),
            "overlap": sorted(functional_set & comention_set),
            "overlap_count": len(functional_set & comention_set),
            "jaccard": round(jaccard(functional_set, comention_set), 4),
        },
        "seed_functional_partners_in_layer": seed_edges,
        "counts": {"human_pairs_at_threshold": len(links),
                   "gene_gene_edges": len(edges), "new_gene_nodes": len(new_genes),
                   "new_nodes_total": len(new_node_rows), "new_drug_edges": len(drug_edges),
                   "new_pathway_edges": len(pathway_edges), "new_drug_rows": len(phases)},
        "files": {name: {"sha256": common.sha256_of(path), "bytes": path.stat().st_size}
                  for name, path in (("nodes.tsv", nodes_path), ("edges.tsv", edges_path),
                                     ("drug_evidence.tsv", drugs_path))},
    }
    path = common.write_json(args.out_dir / "manifest.json", manifest)
    print(f"  functional neighbours {len(neighbours)}, overlap with co-mention "
          f"{len(functional_set & comention_set)}, Jaccard {manifest['neighbours']['jaccard']}")
    print(f"  {len(edges)} gene to gene edges, {len(new_genes)} new genes, "
          f"{len(drug_edges)} new drug edges, {len(pathway_edges)} new pathway edges")
    print(f"Manifest: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
