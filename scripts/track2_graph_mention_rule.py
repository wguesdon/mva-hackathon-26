"""Run the derived mention rule on every disease node of a Track 2 graph.

`scripts/track2_seed_filter.py` decides whether one mention of a disease concept names the disease,
and `scripts/track2_mention_audit.py` runs that decision over one disease node. Both run on the seed
node of a build. Every other disease node of the same build keeps the entity mapping error that the
rule exists for, and the held-out control of `scripts/track2_controls.py` evaluates those nodes.

This module runs the rule over every disease node of a graph. ADR 0010 fixes its scope, its
parameters and its verdict rule.

Three steps run in order.

``rules`` reads the MeSH entry terms of every disease node and writes one rule file per node with
`scripts/track2_disease_surface_forms.py`. A node whose identifier carries no MeSH concept receives
no rule and keeps every edge.

``fetch`` retrieves every result page of every literature edge of every disease node from the
PubTator3 search endpoint, to at most 20 pages, and caches one file per node. A second run sends no
request for a node that the cache already holds.

``apply`` decides every cached sentence with the rule of its own node, writes a per-node table, an
edge verdict table and a sentence table, and writes the distribution of the share of literature
edges that no sentence supports.

`scripts/build_track2_graph.py --mention-filter-scope graph` calls `apply_graph_filter` of this
module and removes every literature edge of a disease node that the rule rejects.

Usage:
    uv run python scripts/track2_graph_mention_rule.py rules \\
        --graph-dir reports/track2/graph_unfiltered
    uv run python scripts/track2_graph_mention_rule.py fetch \\
        --graph-dir reports/track2/graph_unfiltered
    uv run python scripts/track2_graph_mention_rule.py apply \\
        --graph-dir reports/track2/graph_unfiltered
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import statistics
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import httpx
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402
import track2_disease_surface_forms as surface  # noqa: E402
import track2_mention_audit as audit  # noqa: E402
import track2_seed_filter as seed_filter  # noqa: E402

DEFAULT_GRAPH_DIR = common.ROOT / "reports" / "track2" / "graph_unfiltered"
DEFAULT_RULES_DIR = common.ROOT / "literature" / "track2_surface_forms" / "mva_graph"
DEFAULT_CACHE_DIR = common.RAW_DIR / "graph_mention_cache"
DEFAULT_ENTITIES = common.RAW_DIR / "pubtator_entities.json"
DEFAULT_OUT_DIR = common.ROOT / "reports" / "track2" / "graph_mentions"

DISEASE_KIND = "disease"
MESH_PREFIX = "DISEASE:MESH:"
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
LITERATURE_SOURCES = audit.LITERATURE_SOURCES
MAX_PAGES_ALL = audit.MAX_PAGES_ALL
DEFAULT_WORKERS = audit.DEFAULT_WORKERS

NODE_HEADER = ["node_id", "label", "mesh_id", "has_rule", "name_reachable", "literature_edges",
               "edges_with_sentence", "edges_without_sentence", "sentences", "sentences_accepted",
               "edges_kept", "edges_dropped", "share_unsupported_before", "share_unsupported_after",
               "abbreviation_only_drops", "restorable_drops"]
EDGE_HEADER = ["disease_node", "relation", "source_accession", "target_accession", "sentences",
               "accepted_sentences", "verdict", "abbreviation_only", "restorable",
               "accepted_forms", "rejected_forms"]
FORM_HEADER = ["disease_node", "form", "mentions", "accepted", "verdict"]

# The rule rejects a mention whose matched text is only an abbreviation of the concept and whose
# publication writes no long form. Every false reject of the reader comparison of the mosaic
# variegated aneuploidy node falls in that class, so the class is the measurable proxy of the false
# reject cost on a node that no reader read.
ABBREVIATION_REASON = seed_filter.REASON_UNCORROBORATED


def safe_stem(node_id: str) -> str:
    """Return a file stem for one node identifier.

    Args:
        node_id: The graph node identifier, for example `DISEASE:MESH:C536987`.

    Returns:
        The identifier with each colon read as an underscore, in lower case.
    """
    return node_id.replace(":", "_").lower()


def mesh_id_of(node_id: str) -> str:
    """Return the MeSH identifier of one disease node.

    Args:
        node_id: The graph node identifier.

    Returns:
        The MeSH identifier, or an empty string when the node carries no MeSH concept.
    """
    if not node_id.startswith(MESH_PREFIX):
        return ""
    return node_id[len(MESH_PREFIX):]


def disease_nodes(nodes_path: Path) -> dict[str, dict]:
    """Read every disease node of the graph.

    Args:
        nodes_path: The path of `nodes.tsv`.

    Returns:
        The node rows, keyed by node identifier, in table order.
    """
    out: dict[str, dict] = {}
    with nodes_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["kind"] == DISEASE_KIND:
                out[row["node_id"]] = row
    return out


def literature_edges_by_node(edges_path: Path, nodes: dict[str, dict]) -> dict[str, list[dict]]:
    """Group the literature edges of the graph by the disease node that they touch.

    Args:
        edges_path: The path of `edges.tsv`.
        nodes: The disease node rows, keyed by node identifier.

    Returns:
        One list of distinct source, target and relation records per disease node.
    """
    out: dict[str, list[dict]] = {node: [] for node in nodes}
    seen: set[tuple[str, str, str]] = set()
    with edges_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["source"] not in LITERATURE_SOURCES:
                continue
            key = (row["source_id"], row["target_id"], row["relation"])
            if key in seen:
                continue
            seen.add(key)
            for endpoint in (row["source_id"], row["target_id"]):
                if endpoint in out:
                    out[endpoint].append(row)
    return out


def rule_path_for(rules_dir: Path, node_id: str) -> Path:
    """Return the rule file path of one disease node.

    Args:
        rules_dir: The directory that holds one rule file per node.
        node_id: The graph node identifier.

    Returns:
        The rule file path.
    """
    return rules_dir / f"{safe_stem(node_id)}.yaml"


def cache_path_for(cache_dir: Path, node_id: str) -> Path:
    """Return the sentence cache path of one disease node.

    Args:
        cache_dir: The directory that holds one cache file per node.
        node_id: The graph node identifier.

    Returns:
        The cache file path.
    """
    return cache_dir / f"{safe_stem(node_id)}.json"


def accession_for_node(node_id: str, accessions: dict[str, str]) -> str:
    """Return the PubTator3 accession of one node identifier.

    Args:
        node_id: The node identifier.
        accessions: The mapping of `track2_mention_audit.accession_by_node`.

    Returns:
        The accession, or an empty string when the node does not resolve.
    """
    return audit.accession_for(node_id, accessions)


def search_mesh_uid(mesh_id: str) -> str:
    """Look up the NCBI unique identifier of a MeSH concept with the search endpoint.

    `track2_disease_surface_forms.mesh_uid` builds the unique identifier from the MeSH identifier,
    which is the rule for a record that NCBI created before 2015. A later record such as D000077192
    carries an unrelated unique identifier and the built value returns the body "cannot get document
    summary". The search endpoint answers the `MeSH Unique ID` field for both classes of record.

    Args:
        mesh_id: The MeSH identifier, for example `D000077192`.

    Returns:
        The unique identifier, or an empty string when the search returns no record.
    """
    with common.make_client() as client:
        payload = common.get_json(client, ESEARCH_URL,
                                  {"db": "mesh", "retmode": "json",
                                   "term": f"{mesh_id}[MeSH Unique ID]"})
    ids = ((payload or {}).get("esearchresult") or {}).get("idlist") or []
    return str(ids[0]) if ids else ""


def entry_terms_by_uid(uid: str) -> tuple[list[str], str]:
    """Read the entry terms of one MeSH record by its unique identifier.

    Args:
        uid: The NCBI unique identifier of the record.

    Returns:
        The entry terms in record order, and the retrieval time.

    Raises:
        ValueError: When the record holds no entry term.
    """
    with common.make_client() as client:
        payload = common.get_json(client, surface.ESUMMARY_URL,
                                  {"db": "mesh", "id": uid, "retmode": "json"})
    record = ((payload or {}).get("result") or {}).get(uid) or {}
    terms = [str(term).strip() for term in record.get("ds_meshterms", []) if str(term).strip()]
    if not terms:
        raise ValueError(f"MeSH record {uid} returned no entry term.")
    return terms, dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def entry_terms_for(mesh_id: str) -> tuple[list[str], str, str, str]:
    """Read the entry terms of one MeSH concept, with the search endpoint as a fallback.

    Args:
        mesh_id: The MeSH identifier.

    Returns:
        The entry terms, the retrieval time, the unique identifier that answered and the route
        that produced that identifier, which is `formula` or `esearch`.

    Raises:
        ValueError: When neither route returns an entry term.
    """
    try:
        terms, fetched_at = surface.fetch_entry_terms(mesh_id)
        return terms, fetched_at, surface.mesh_uid(mesh_id), "formula"
    except ValueError:
        uid = search_mesh_uid(mesh_id)
        if not uid:
            raise
        terms, fetched_at = entry_terms_by_uid(uid)
        return terms, fetched_at, uid, "esearch"


def write_rule_file(node_id: str, label: str, accession: str, rules_dir: Path) -> tuple[Path, int]:
    """Fetch the MeSH entry terms of one disease node and write its rule file.

    Args:
        node_id: The graph node identifier.
        label: The node label.
        accession: The PubTator3 accession of the node.
        rules_dir: The output directory.

    Returns:
        The rule file path and the number of entry terms.

    Raises:
        ValueError: When the MeSH record holds no entry term.
    """
    mesh_id = mesh_id_of(node_id)
    terms, fetched_at, uid, route = entry_terms_for(mesh_id)
    rule = surface.build_rule(mesh_id, node_id, accession, label or terms[0], terms, fetched_at)
    rule["source"]["mesh_uid"] = uid
    rule["source"]["mesh_uid_route"] = route
    path = rule_path_for(rules_dir, node_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(rule, sort_keys=False, allow_unicode=True, width=98))
    return path, len(terms)


def run_rules(args: argparse.Namespace) -> int:
    """Write one rule file for every disease node that carries a MeSH concept.

    Args:
        args: The parsed command line.

    Returns:
        0 on success.
    """
    nodes = disease_nodes(args.graph_dir / "nodes.tsv")
    accessions = audit.accession_by_node(args.entities)
    written, skipped, failed, no_mesh = 0, 0, [], []
    for node_id, row in sorted(nodes.items()):
        if not mesh_id_of(node_id):
            no_mesh.append(node_id)
            continue
        path = rule_path_for(args.rules_dir, node_id)
        if path.is_file() and not args.refresh:
            skipped += 1
            continue
        try:
            written_path, terms = write_rule_file(node_id, row.get("label", ""),
                                                  accession_for_node(node_id, accessions),
                                                  args.rules_dir)
        except (httpx.HTTPError, ValueError) as error:
            failed.append(f"{node_id} ({type(error).__name__})")
            continue
        written += 1
        if written % 25 == 0:
            print(f"  {written} rule files written, last {written_path.name} with {terms} terms")
        time.sleep(common.DEFAULT_PAUSE_S)
    print(f"disease nodes {len(nodes)}, rule files written {written}, already present {skipped}, "
          f"no MeSH concept {len(no_mesh)}, failed {len(failed)}")
    if no_mesh:
        print(f"  no MeSH concept: {', '.join(no_mesh)}")
    if failed:
        print(f"  failed: {', '.join(failed)}")
    return 0


def fetch_node(spec: tuple[str, str, str, str]) -> tuple[str, dict | None, str]:
    """Retrieve every page of one entity pair with its own client.

    Args:
        spec: The cache key, the relation and the two accessions.

    Returns:
        The cache key, the pages, and the error name. The pages are None when the endpoint refused
        the query, and the error name is empty on success.
    """
    return audit.fetch_one(spec)


def pair_specs(edges: list[dict], node_id: str,
               accessions: dict[str, str]) -> tuple[list[tuple[str, str, str, str]], list[str]]:
    """Build the retrieval specification of every literature edge of one disease node.

    Args:
        edges: The literature edge rows of the node.
        node_id: The disease node identifier.
        accessions: The mapping of `track2_mention_audit.accession_by_node`.

    Returns:
        The specifications and the edge keys whose endpoints have no accession.
    """
    specs: list[tuple[str, str, str, str]] = []
    unresolved: list[str] = []
    seen: set[str] = set()
    for edge in edges:
        acc_1 = accession_for_node(edge["source_id"], accessions)
        acc_2 = accession_for_node(edge["target_id"], accessions)
        if not acc_1 or not acc_2:
            unresolved.append(f"{edge['source_id']}|{edge['target_id']}|{edge['relation']}")
            continue
        key = f"{edge['relation']}|{acc_1}|{acc_2}|{MAX_PAGES_ALL}"
        if key in seen:
            continue
        seen.add(key)
        specs.append((key, edge["relation"], acc_1, acc_2))
    if not specs and not unresolved:
        return specs, unresolved
    return specs, unresolved


def run_fetch(args: argparse.Namespace) -> int:
    """Fill the sentence cache of every literature edge of every disease node.

    Args:
        args: The parsed command line.

    Returns:
        0 on success.
    """
    nodes = disease_nodes(args.graph_dir / "nodes.tsv")
    by_node = literature_edges_by_node(args.graph_dir / "edges.tsv", nodes)
    accessions = audit.accession_by_node(args.entities)
    order = sorted(nodes, key=lambda node: (-len(by_node[node]), node))
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    done, pairs_fetched, refused_total, unresolved_total = 0, 0, 0, 0
    started = time.monotonic()
    for node_id in order:
        path = cache_path_for(args.cache_dir, node_id)
        cache = common.read_json(path) or {}
        pairs = cache.get("pairs", {})
        refused = list(cache.get("refused_pairs", []))
        fetched_at = (cache.get("fetched_at")
                      or dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"))
        specs, unresolved = pair_specs(by_node[node_id], node_id, accessions)
        todo = [spec for spec in specs if spec[0] not in pairs and spec[0] not in refused]
        done += 1
        if not todo:
            unresolved_total += len(unresolved)
            refused_total += len(refused)
            pairs_fetched += len(pairs)
            continue
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            for key, pages, error in pool.map(fetch_node, todo):
                if pages is None:
                    refused.append(key)
                else:
                    pairs[key] = pages
        common.write_json(path, {"fetched_at": fetched_at, "endpoint": audit.SEARCH_URL,
                                 "node": node_id, "unresolved_edges": unresolved,
                                 "refused_pairs": refused, "pairs": pairs})
        pairs_fetched += len(pairs)
        refused_total += len(refused)
        unresolved_total += len(unresolved)
        elapsed = time.monotonic() - started
        print(f"  [{done}/{len(order)}] {node_id} {len(pairs)} pairs, {len(refused)} refused, "
              f"{len(unresolved)} unresolved, {elapsed:.0f}s elapsed", flush=True)
    print(f"nodes {len(order)}, cached pairs {pairs_fetched}, refused {refused_total}, "
          f"edges with no accession {unresolved_total}")
    return 0


def load_rules(rules_dir: Path, nodes: dict[str, dict]) -> dict[str, seed_filter.SurfaceRule]:
    """Read the rule file of every disease node that has one.

    Args:
        rules_dir: The directory that holds one rule file per node.
        nodes: The disease node rows, keyed by node identifier.

    Returns:
        The rule of each node that has a rule file.
    """
    rules: dict[str, seed_filter.SurfaceRule] = {}
    for node_id in nodes:
        path = rule_path_for(rules_dir, node_id)
        if not path.is_file():
            continue
        try:
            rules[node_id] = seed_filter.load_rule(path)
        except ValueError:
            # Amendment 1 of ADR 0010 names this class. Every entry term of the concept is one
            # short content token, so the distinctiveness rule produces no token set and no accept
            # pattern, and no matched text can name the concept. The node then carries no usable
            # rule. It joins the nodes with no rule at all, so its edges stay in the graph and the
            # run counts it, rather than losing every edge to a rule that can accept nothing.
            continue
    return rules


def node_pair_verdicts(cache: dict, rule: seed_filter.SurfaceRule) -> dict[tuple, dict]:
    """Decide every cached pair of one disease node.

    Args:
        cache: The cache record of the node.
        rule: The rule of the node.

    Returns:
        The records of `track2_seed_filter.pair_verdicts`.
    """
    return seed_filter.pair_verdicts(cache.get("pairs", {}) or {}, rule)


def name_reachable(rule: seed_filter.SurfaceRule) -> bool:
    """Say whether any name of the concept can pass the distinctiveness rule.

    A mention is accepted when the content tokens of its matched text and the content tokens of one
    entry term contain one another and the shared part is distinctive. A concept whose every entry
    term reduces to one short content token fails that test on its own preferred term, so no long
    form of it is ever accepted and no abbreviation of it is ever corroborated by a long form.
    Amendment 1 of ADR 0010 adds this count.

    Args:
        rule: The rule of the node.

    Returns:
        True when at least one entry term token set passes the distinctiveness rule.
    """
    if rule.accept_patterns:
        return True
    return any(seed_filter.distinctive(set(entry), rule) for entry in rule.token_sets)


def near_miss(form: str, rule: seed_filter.SurfaceRule) -> bool:
    """Say whether a rejected matched text is a name of the concept that the threshold refused.

    The matched text contains an entry term or sits inside one, and the shared content tokens fail
    the distinctiveness rule on their own. "hypoxia" against the entry term `Hypoxia` is the case.

    Args:
        form: The matched text of the mention.
        rule: The rule of the node.

    Returns:
        True when a containment holds and the shared part is not distinctive.
    """
    tokens = set(surface.content_tokens(form))
    if not tokens or rule.accept_patterns:
        return False
    for entry in rule.token_sets:
        shared = tokens & entry
        if shared and (tokens <= entry or entry <= tokens) \
                and not seed_filter.distinctive(shared, rule):
            return True
    return False


def restorable(record: dict, rule: seed_filter.SurfaceRule) -> bool:
    """Say whether a dropped edge carries a sentence that a reader could restore.

    Two classes of rejection carry the false reject cost of the rule. The first is a mention whose
    matched text is an abbreviation of the concept with no long form in its publication. The second
    is a near miss, which is a name of the concept that the distinctiveness rule refused. An edge
    with neither class holds no mention that the rule connected to the concept at all.

    Args:
        record: One record of `track2_seed_filter.pair_verdicts`.
        rule: The rule of the node.

    Returns:
        True when the edge is dropped and one of its sentences falls in either class.
    """
    if record["verdict"] != "drop":
        return False
    for item in record["items"]:
        if item["reason"] == ABBREVIATION_REASON:
            return True
        if any(near_miss(form, rule) for form in item["rejected_forms"]):
            return True
    return False


def abbreviation_only(record: dict) -> bool:
    """Say whether an edge was dropped only because its abbreviations carry no long form.

    Args:
        record: One record of `track2_seed_filter.pair_verdicts`.

    Returns:
        True when the edge is dropped and at least one of its sentences was rejected with the
        uncorroborated abbreviation reason.
    """
    if record["verdict"] != "drop":
        return False
    return any(item["reason"] == ABBREVIATION_REASON for item in record["items"])


def node_summary(node_id: str, row: dict, edges: list[dict], rule: seed_filter.SurfaceRule | None,
                 verdicts: dict[tuple, dict]) -> dict:
    """Summarize the rule on one disease node.

    Args:
        node_id: The disease node identifier.
        row: The node row of `nodes.tsv`.
        edges: The literature edge rows of the node.
        rule: The rule of the node, or None when the node has no MeSH concept.
        verdicts: The records of `node_pair_verdicts`.

    Returns:
        One table row.
    """
    items = [item for record in verdicts.values() for item in record["items"]]
    kept = sum(1 for record in verdicts.values() if record["verdict"] == "keep")
    with_sentence = len(verdicts)
    total = len(edges)
    return {
        "node_id": node_id,
        "label": row.get("label", ""),
        "mesh_id": mesh_id_of(node_id),
        "has_rule": "yes" if rule is not None else "no",
        "literature_edges": total,
        "edges_with_sentence": with_sentence,
        "edges_without_sentence": total - with_sentence,
        "sentences": len(items),
        "sentences_accepted": sum(1 for item in items if item["accepted"]),
        "edges_kept": kept,
        "edges_dropped": with_sentence - kept,
        "share_unsupported_before": (round(1 - kept / with_sentence, 4) if with_sentence else None),
        "share_unsupported_after": (0.0 if with_sentence else None),
        "abbreviation_only_drops": sum(1 for record in verdicts.values()
                                       if abbreviation_only(record)),
        "restorable_drops": (sum(1 for record in verdicts.values() if restorable(record, rule))
                             if rule is not None else 0),
        "name_reachable": ("yes" if rule is not None and name_reachable(rule)
                           else "no" if rule is not None else ""),
    }


def edge_rows_of(node_id: str, verdicts: dict[tuple, dict],
                 rule: seed_filter.SurfaceRule | None) -> list[dict]:
    """Build one verdict row per literature edge of one disease node.

    Args:
        node_id: The disease node identifier.
        verdicts: The records of `node_pair_verdicts`.
        rule: The rule of the node, or None when the node has no MeSH concept.

    Returns:
        The rows, sorted by verdict and accepted sentence count.
    """
    rows = [{"disease_node": node_id, "relation": record["relation"],
             "source_accession": record["source_accession"],
             "target_accession": record["target_accession"],
             "sentences": record["sentences"],
             "accepted_sentences": record["accepted_sentences"],
             "verdict": record["verdict"],
             "abbreviation_only": "yes" if abbreviation_only(record) else "no",
             "restorable": ("yes" if rule is not None and restorable(record, rule) else "no"),
             "accepted_forms": " | ".join(record["accepted_forms"]),
             "rejected_forms": " | ".join(record["rejected_forms"])}
            for record in verdicts.values()]
    return sorted(rows, key=lambda row: (row["verdict"], -row["accepted_sentences"],
                                         row["relation"], row["source_accession"],
                                         row["target_accession"]))


def form_rows_of(node_id: str, verdicts: dict[tuple, dict],
                 rule: seed_filter.SurfaceRule) -> list[dict]:
    """Build one row per observed surface form of one disease node.

    Args:
        node_id: The disease node identifier.
        verdicts: The records of `node_pair_verdicts`.
        rule: The rule of the node.

    Returns:
        The rows, in descending order of mention count.
    """
    return [{"disease_node": node_id, "form": record["form"], "mentions": record["mentions"],
             "accepted": record["accepted"], "verdict": record["verdict"]}
            for record in audit.form_records(verdicts, rule) if record["form"]]


def near_miss_counts(form_rows: list[dict],
                     rules: dict[str, seed_filter.SurfaceRule]) -> dict:
    """Count the rejected mentions that the distinctiveness rule refused although a name matched.

    Amendment 1 of ADR 0010 adds this count. A near miss is a matched text that contains an entry
    term of the concept or sits inside one, and whose shared content tokens fail the distinctiveness
    rule on their own.

    Args:
        form_rows: The rows of `form_rows_of`.
        rules: The rule of each node that has one.

    Returns:
        The total rejected mentions of the class, the number of distinct forms, and the fifteen
        forms with the largest count.
    """
    totals: Counter[str] = Counter()
    for row in form_rows:
        rule = rules.get(row["disease_node"])
        missed = int(row["mentions"]) - int(row["accepted"])
        if rule is None or missed <= 0 or not near_miss(row["form"], rule):
            continue
        totals[f'{row["form"]} ({row["disease_node"]})'] += missed
    return {"mentions": sum(totals.values()), "distinct_node_and_form_pairs": len(totals),
            "largest": [{"form_and_node": key, "rejected_mentions": value}
                        for key, value in totals.most_common(15)]}


def most_rejected(form_rows: list[dict], count: int) -> list[dict]:
    """Return the surface forms that the rule rejected most often across the whole graph.

    Args:
        form_rows: The rows of `form_rows_of`.
        count: The number of forms to return.

    Returns:
        One record per form, with the rejected mention count and the nodes that carry it.
    """
    rejected: Counter[str] = Counter()
    nodes: dict[str, set[str]] = {}
    for row in form_rows:
        missed = int(row["mentions"]) - int(row["accepted"])
        if missed <= 0:
            continue
        rejected[row["form"]] += missed
        nodes.setdefault(row["form"], set()).add(row["disease_node"])
    return [{"form": form, "rejected_mentions": total, "disease_nodes": len(nodes[form])}
            for form, total in rejected.most_common(count)]


def retrieval_block(cache_dir: Path) -> dict:
    """Summarize the retrieval from the cache files that it wrote.

    Every cache file carries the start time of the run in its `fetched_at` field, and the
    modification time of the last file is the end of the run, so the two give the wall clock with no
    record outside the cache itself.

    Args:
        cache_dir: The directory that holds one cache file per node.

    Returns:
        The node count, the pair count, the result count, the refused pairs, the byte size and the
        wall clock of the retrieval.
    """
    files = sorted(cache_dir.glob("*.json")) if cache_dir.is_dir() else []
    if not files:
        return {"nodes": 0}
    pairs, results, refused, unresolved, size = 0, 0, 0, 0, 0
    started: list[str] = []
    last = 0.0
    for path in files:
        record = common.read_json(path) or {}
        cached = record.get("pairs", {}) or {}
        pairs += len(cached)
        results += sum(len(pair.get("items", [])) for pair in cached.values())
        refused += len(record.get("refused_pairs", []) or [])
        unresolved += len(record.get("unresolved_edges", []) or [])
        size += path.stat().st_size
        last = max(last, path.stat().st_mtime)
        if record.get("fetched_at"):
            started.append(str(record["fetched_at"]))
    start = min(started) if started else ""
    finished = dt.datetime.fromtimestamp(last, dt.timezone.utc).isoformat(timespec="seconds")
    minutes = None
    if start:
        begin = dt.datetime.fromisoformat(start)
        minutes = round((dt.datetime.fromisoformat(finished) - begin).total_seconds() / 60.0, 1)
    return {"nodes": len(files), "entity_pairs": pairs, "search_results": results,
            "entity_pairs_refused": refused, "edges_with_no_accession": unresolved,
            "cache_bytes": size, "started_utc": start, "finished_utc": finished,
            "wall_clock_minutes": minutes, "page_cap": MAX_PAGES_ALL,
            "workers": DEFAULT_WORKERS, "money_cost_usd": 0.0}


def percentile_of(values: list[float], target: float) -> float | None:
    """Return the percentile of one rate inside a distribution of rates.

    Args:
        values: The rates.
        target: The rate to place.

    Returns:
        The share of the rates that fall below the target, in percent.
    """
    if not values:
        return None
    return round(100.0 * sum(1 for value in values if value < target) / len(values), 1)


def quantile_block(values: list[float]) -> dict:
    """Return the distribution of a list of rates.

    Args:
        values: The rates.

    Returns:
        The count, the minimum, the three quartiles, the maximum and the mean.
    """
    if not values:
        return {"n": 0}
    ordered = sorted(values)
    block = {"n": len(ordered), "min": round(ordered[0], 4), "max": round(ordered[-1], 4),
             "median": round(statistics.median(ordered), 4),
             "mean": round(statistics.fmean(ordered), 4)}
    if len(ordered) >= 4:
        quarters = statistics.quantiles(ordered, n=4, method="inclusive")
        block["q1"] = round(quarters[0], 4)
        block["q3"] = round(quarters[2], 4)
    if len(ordered) >= 10:
        tenths = statistics.quantiles(ordered, n=10, method="inclusive")
        block["p90"] = round(tenths[8], 4)
    return block


def extremes(rows: list[dict], key: str, count: int, reverse: bool) -> list[dict]:
    """Return the rows at one end of a distribution.

    Args:
        rows: The node summary rows.
        key: The column to order on.
        count: The number of rows to return.
        reverse: True for the largest values.

    Returns:
        The rows, each reduced to the fields that a record needs.
    """
    usable = [row for row in rows if row[key] is not None]
    ordered = sorted(usable, key=lambda row: (row[key], -row["literature_edges"], row["node_id"]),
                     reverse=reverse)
    return [{"node_id": row["node_id"], "label": row["label"], key: row[key],
             "literature_edges": row["literature_edges"],
             "edges_with_sentence": row["edges_with_sentence"],
             "edges_kept": row["edges_kept"]}
            for row in ordered[:count]]


def collect(graph_dir: Path, rules_dir: Path, cache_dir: Path,
            entities_path: Path) -> tuple[dict, dict, dict, dict]:
    """Read every node, its rule, its cache and its verdicts.

    Args:
        graph_dir: The graph directory.
        rules_dir: The rule directory.
        cache_dir: The sentence cache directory.
        entities_path: The resolved entity cache.

    Returns:
        The disease node rows, the literature edges of each node, the rule of each node and the
        verdicts of each node.
    """
    nodes = disease_nodes(graph_dir / "nodes.tsv")
    by_node = literature_edges_by_node(graph_dir / "edges.tsv", nodes)
    rules = load_rules(rules_dir, nodes)
    verdicts: dict[str, dict] = {}
    for node_id in nodes:
        rule = rules.get(node_id)
        if rule is None:
            verdicts[node_id] = {}
            continue
        cache = common.read_json(cache_path_for(cache_dir, node_id)) or {}
        verdicts[node_id] = node_pair_verdicts(cache, rule)
    return nodes, by_node, rules, verdicts


def dropped_edge_keys_graph(verdicts: dict[str, dict], entities: dict) -> set[tuple[str, str, str]]:
    """List the graph edge keys that the rule removes over every disease node.

    Args:
        verdicts: The verdicts of each node.
        entities: The resolved entity cache, keyed by accession.

    Returns:
        The source, target and relation of every removed edge, in both orders.
    """
    keys: set[tuple[str, str, str]] = set()
    for node_verdicts in verdicts.values():
        keys |= seed_filter.dropped_edge_keys(node_verdicts, entities)
    return keys


def apply_graph_filter(edges: list[dict], rules_dir: Path, cache_dir: Path, entities: dict,
                       graph_dir: Path) -> tuple[list[dict], dict]:
    """Remove every literature edge of a disease node that the rule of that node rejects.

    Args:
        edges: The edge rows of the whole graph.
        rules_dir: The directory that holds one rule file per node.
        cache_dir: The sentence cache directory.
        entities: The resolved entity cache, keyed by accession.
        graph_dir: The graph directory that holds the node table of the unfiltered build.

    Returns:
        The kept edge rows and the record of the filter for the manifest.

    Raises:
        FileNotFoundError: When no rule file and no cache file exists.
    """
    nodes = disease_nodes(graph_dir / "nodes.tsv")
    rules = load_rules(rules_dir, nodes)
    if not rules:
        raise FileNotFoundError(f"no rule file under {rules_dir}; run "
                                "scripts/track2_graph_mention_rule.py rules first")
    verdicts: dict[str, dict] = {}
    for node_id in nodes:
        rule = rules.get(node_id)
        cache = common.read_json(cache_path_for(cache_dir, node_id)) or {} if rule else {}
        verdicts[node_id] = node_pair_verdicts(cache, rule) if rule else {}
    removed = dropped_edge_keys_graph(verdicts, entities)
    kept = [row for row in edges
            if not (row["source"] in LITERATURE_SOURCES
                    and (row["source_id"], row["target_id"], row["relation"]) in removed)]
    all_records = [record for node in verdicts.values() for record in node.values()]
    items = [item for record in all_records for item in record["items"]]
    return kept, {
        "mode": "graph",
        "rules_dir": seed_filter.repo_path(rules_dir),
        "cache_dir": seed_filter.repo_path(cache_dir),
        "disease_nodes": len(nodes),
        "disease_nodes_with_a_rule": len(rules),
        "literature_edges_decided": len(all_records),
        "literature_edges_kept": sum(1 for record in all_records if record["verdict"] == "keep"),
        "literature_edges_dropped": sum(1 for record in all_records
                                        if record["verdict"] == "drop"),
        "sentences": len(items),
        "sentences_accepted": sum(1 for item in items if item["accepted"]),
        "edge_rows_removed": len(edges) - len(kept),
        "rejection_reasons": dict(Counter(item["reason"] for item in items
                                          if not item["accepted"]).most_common()),
    }


def run_apply(args: argparse.Namespace) -> int:
    """Decide every cached sentence of every disease node and write the tables.

    Args:
        args: The parsed command line.

    Returns:
        0 on success, 1 when no rule file exists.
    """
    nodes, by_node, rules, verdicts = collect(args.graph_dir, args.rules_dir, args.cache_dir,
                                              args.entities)
    if not rules:
        print(f"FAIL: no rule file under {args.rules_dir}.")
        return 1
    node_rows, edge_rows, form_rows, sentence_rows = [], [], [], []
    entities = (common.read_json(args.entities) or {}).get("entities", {})
    for node_id, row in nodes.items():
        rule = rules.get(node_id)
        node_verdicts = verdicts[node_id]
        node_rows.append(node_summary(node_id, row, by_node[node_id], rule, node_verdicts))
        edge_rows.extend(edge_rows_of(node_id, node_verdicts, rule))
        if rule is None:
            continue
        form_rows.extend(form_rows_of(node_id, node_verdicts, rule))
        for sentence in audit.sentence_rows(node_verdicts, entities, rule):
            sentence["disease_node"] = node_id
            sentence_rows.append(sentence)
    node_rows.sort(key=lambda row: (-row["literature_edges"], row["node_id"]))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    node_path = common.write_tsv(args.out_dir / "node_summary.tsv", NODE_HEADER, node_rows)
    edge_path = common.write_tsv(args.out_dir / "edge_verdicts.tsv", EDGE_HEADER, edge_rows)
    form_path = common.write_tsv(args.out_dir / "surface_forms.tsv", FORM_HEADER, form_rows)
    sentence_path = common.write_tsv(args.sentences_out,
                                     ["disease_node"] + audit.SENTENCE_HEADER, sentence_rows)
    report = summarize(node_rows, form_rows, sentence_rows, rules, args)
    report["near_miss_rejections"] = near_miss_counts(form_rows, rules)
    report["tables"] = {name: {"path": seed_filter.repo_path(path),
                               "sha256_prefix": common.sha256_of(path)[:16],
                               "rows": count}
                        for name, path, count in (("node_summary", node_path, len(node_rows)),
                                                  ("edge_verdicts", edge_path, len(edge_rows)),
                                                  ("surface_forms", form_path, len(form_rows)),
                                                  ("sentences", sentence_path,
                                                   len(sentence_rows)))}
    common.write_json(args.report_out, report)
    print(f"disease nodes {report['disease_nodes']}, with a rule "
          f"{report['disease_nodes_with_a_rule']}")
    print(f"literature edges {report['literature_edges']}, decided "
          f"{report['literature_edges_decided']}, kept {report['literature_edges_kept']}, "
          f"dropped {report['literature_edges_dropped']}")
    print(f"sentences {report['sentences']}, accepted {report['sentences_accepted']}")
    print(f"share of edges that no sentence supports: {report['share_unsupported']}")
    print(f"nodes that lose every edge: {report['nodes_that_lose_every_edge']}")
    print(f"report: {args.report_out}")
    return 0


def summarize(node_rows: list[dict], form_rows: list[dict], sentence_rows: list[dict],
              rules: dict[str, seed_filter.SurfaceRule], args: argparse.Namespace) -> dict:
    """Build the report record of one graph-wide run.

    Args:
        node_rows: The rows of `node_summary`.
        form_rows: The rows of `form_rows_of`.
        sentence_rows: The sentence rows.
        rules: The rule of each node that has one.
        args: The parsed command line.

    Returns:
        The report record.
    """
    decided = [row for row in node_rows if row["edges_with_sentence"]]
    shares = [row["share_unsupported_before"] for row in decided]
    weighted_drop = sum(row["edges_dropped"] for row in node_rows)
    weighted_decided = sum(row["edges_with_sentence"] for row in node_rows)
    lost_all = [row for row in decided if row["edges_kept"] == 0]
    return {
        "run_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "graph_dir": seed_filter.repo_path(args.graph_dir),
        "rules_dir": seed_filter.repo_path(args.rules_dir),
        "cache_dir": seed_filter.repo_path(args.cache_dir),
        "page_cap": MAX_PAGES_ALL,
        "disease_nodes": len(node_rows),
        "disease_nodes_with_a_rule": len(rules),
        "disease_nodes_without_a_rule": [row["node_id"] for row in node_rows
                                         if row["has_rule"] == "no"],
        "literature_edges": sum(row["literature_edges"] for row in node_rows),
        "literature_edges_decided": weighted_decided,
        "literature_edges_without_a_sentence": sum(row["edges_without_sentence"]
                                                   for row in node_rows),
        "literature_edges_kept": sum(row["edges_kept"] for row in node_rows),
        "literature_edges_dropped": weighted_drop,
        "sentences": len(sentence_rows),
        "sentences_accepted": sum(1 for row in sentence_rows if row["rule_accepts"] == "yes"),
        "rejection_reasons": dict(Counter(row["reason"] for row in sentence_rows
                                          if row["rule_accepts"] == "no").most_common()),
        "share_unsupported": quantile_block(shares),
        "nodes_that_lose_no_edge": sum(1 for value in shares if value == 0.0),
        "nodes_at_or_below_0_1": sum(1 for value in shares if value <= 0.1),
        "nodes_at_or_above_0_5": sum(1 for value in shares if value >= 0.5),
        "percentile_of_the_four_earlier_rates": {
            str(rate): percentile_of(shares, rate) for rate in (0.153, 0.346, 0.414, 0.959)},
        "retrieval": retrieval_block(args.cache_dir),
        "share_unsupported_pooled": (round(weighted_drop / weighted_decided, 4)
                                     if weighted_decided else None),
        "nodes_with_the_most_literature_edges": [
            {"node_id": row["node_id"], "label": row["label"],
             "literature_edges": row["literature_edges"],
             "edges_with_sentence": row["edges_with_sentence"],
             "edges_kept": row["edges_kept"],
             "share_unsupported_before": row["share_unsupported_before"]}
            for row in sorted(decided, key=lambda r: (-r["literature_edges"], r["node_id"]))[:12]],
        "nodes_at_the_low_extreme": extremes(decided, "share_unsupported_before", 10, False),
        "nodes_at_the_high_extreme": extremes(decided, "share_unsupported_before", 10, True),
        "nodes_that_lose_every_edge": len(lost_all),
        "nodes_that_lose_every_edge_list": [{"node_id": row["node_id"], "label": row["label"],
                                             "literature_edges": row["literature_edges"]}
                                            for row in sorted(
                                                lost_all,
                                                key=lambda r: (-r["literature_edges"],
                                                               r["node_id"]))[:25]],
        "abbreviation_only_drops": sum(row["abbreviation_only_drops"] for row in node_rows),
        "restorable_drops": sum(row["restorable_drops"] for row in node_rows),
        "sentences_rejected_as_uncorroborated_abbreviation": sum(
            1 for row in sentence_rows if row["reason"] == ABBREVIATION_REASON),
        "nodes_that_no_name_reaches": [
            {"node_id": row["node_id"], "label": row["label"],
             "literature_edges": row["literature_edges"], "edges_kept": row["edges_kept"]}
            for row in node_rows if row["name_reachable"] == "no"],
        "nodes_that_no_name_reaches_totals": {
            "nodes": sum(1 for row in node_rows if row["name_reachable"] == "no"),
            "literature_edges": sum(row["literature_edges"] for row in node_rows
                                    if row["name_reachable"] == "no"),
            "edges_kept": sum(row["edges_kept"] for row in node_rows
                              if row["name_reachable"] == "no")},
        "distinct_surface_forms": len({row["form"] for row in form_rows}),
        "distinct_node_and_form_pairs": len(form_rows),
        "most_rejected_forms": most_rejected(form_rows, 15),
    }


def build_parser() -> argparse.ArgumentParser:
    """Build the command line parser.

    Returns:
        The parser with one subcommand per step.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    steps = parser.add_subparsers(dest="step", required=True)
    for name, help_text in (("rules", "Write one rule file per disease node."),
                            ("fetch", "Cache every sentence of every disease node."),
                            ("apply", "Decide every cached sentence of every disease node.")):
        step = steps.add_parser(name, help=help_text)
        step.add_argument("--graph-dir", type=Path, default=DEFAULT_GRAPH_DIR,
                          help="The graph directory of the unfiltered build.")
        step.add_argument("--rules-dir", type=Path, default=DEFAULT_RULES_DIR,
                          help="The directory that holds one rule file per node.")
        step.add_argument("--entities", type=Path, default=DEFAULT_ENTITIES,
                          help="The resolved entity cache.")
        if name == "rules":
            step.add_argument("--refresh", action="store_true",
                              help="Write a rule file again when one is already present.")
        if name in ("fetch", "apply"):
            step.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR,
                              help="The sentence cache directory.")
        if name == "fetch":
            step.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                              help="The number of concurrent requests to the search endpoint.")
        if name == "apply":
            step.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                              help="The output directory of the tables.")
            step.add_argument("--report-out", type=Path,
                              default=DEFAULT_OUT_DIR / "graph_mention_report.json",
                              help="The report file to write.")
            step.add_argument("--sentences-out", type=Path,
                              default=common.RAW_DIR / "graph_mention_sentences.tsv",
                              help="The sentence table to write. The default sits with the raw "
                                   "cache, because the table reaches 40 MB and the repository "
                                   "commits the edge verdict table and the surface form table "
                                   "instead.")
    return parser


def main() -> int:
    """Run one step of the graph-wide mention rule.

    Returns:
        The exit status of the step.
    """
    args = build_parser().parse_args()
    if args.step == "rules":
        return run_rules(args)
    if args.step == "fetch":
        return run_fetch(args)
    return run_apply(args)


if __name__ == "__main__":
    sys.exit(main())
