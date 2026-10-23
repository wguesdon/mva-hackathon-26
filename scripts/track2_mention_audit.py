"""Read every sentence of every literature edge of a disease node and apply the mention rule.

The Track 2 graph joins a disease node to a gene or to a compound when PubTator3 reports that the
two entities co-occur. The recogniser maps a stretch of text onto a concept, and it maps a general
phrase onto a rare disease concept often enough to change a ranking. `scripts/track2_seed_filter.py`
decides whether one mention names the disease, from a rule file that
`scripts/track2_disease_surface_forms.py` derives from the MeSH entry terms of the concept.

This module runs that decision over a whole disease node, in two steps.

``fetch`` reads the literature edges of the node out of `edges.tsv` and retrieves every result page
of each edge from the PubTator3 search endpoint. The pages are cached, so a second run sends no
request.

``apply`` decides every cached sentence, writes a sentence table, writes an edge verdict table,
writes the observed surface forms with a count and a verdict back into the rule file, and, when a
table of reader verdicts is given, measures the rule against the reader.

Usage:
    uv run python scripts/track2_mention_audit.py fetch --seed DISEASE:MESH:D003550 \\
        --raw-dir reports/track2/cystic_fibrosis/raw --edges reports/track2/cystic_fibrosis/graph/edges.tsv
    uv run python scripts/track2_mention_audit.py apply --rule literature/track2_surface_forms/d003550_v1.yaml
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
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
import track2_seed_filter as seed_filter  # noqa: E402

SEARCH_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/"
LITERATURE_SOURCES = ("pubtator3", "pubtator3_consequence")
MAX_PAGES_ALL = 20
SAVE_EVERY = 50

# NCBI asks for no more than three requests a second without an API key. The search endpoint
# answers in about two seconds, so four workers stay under that rate.
DEFAULT_WORKERS = 4
SENTENCE_HEADER = ["relation", "partner", "item_index", "pmid", "surface_form", "rule_accepts",
                   "reason", "sentence"]
NO_FORM = "(no marked mention)"


def literature_edges(edges_path: Path, seed: str) -> list[dict]:
    """Read the literature edges of one disease node from the graph table.

    Args:
        edges_path: The path of `edges.tsv`.
        seed: The disease node identifier.

    Returns:
        One record per distinct source, target and relation, in table order.
    """
    seen: set[tuple[str, str, str]] = set()
    rows: list[dict] = []
    with edges_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["source"] not in LITERATURE_SOURCES:
                continue
            if seed not in (row["source_id"], row["target_id"]):
                continue
            key = (row["source_id"], row["target_id"], row["relation"])
            if key in seen:
                continue
            seen.add(key)
            rows.append(row)
    return rows


def accession_by_node(entities_path: Path) -> dict[str, str]:
    """Map each graph node identifier back to its PubTator3 accession.

    Args:
        entities_path: The resolved entity cache.

    Returns:
        The accession of each node identifier.
    """
    entities = (common.read_json(entities_path) or {}).get("entities", {})
    out: dict[str, str] = {}
    for accession, entry in entities.items():
        node = common.node_id_for_accession(accession, entry.get("identifier", ""),
                                            entry.get("database", ""))
        if node:
            out.setdefault(node, accession)
    return out


def accession_for(node_id: str, accessions: dict[str, str]) -> str:
    """Return the PubTator3 accession of one node identifier.

    Args:
        node_id: The node identifier.
        accessions: The mapping of `accession_by_node`.

    Returns:
        The accession, or an empty string when the node does not resolve.
    """
    if node_id.startswith("GENE:"):
        return f"@GENE_{node_id.split(':', 1)[1]}"
    return accessions.get(node_id, "")


def fetch_pages(client, relation: str, acc_1: str, acc_2: str, max_pages: int) -> dict:
    """Retrieve every result page of one entity pair from the PubTator3 search endpoint.

    Args:
        client: The HTTP client.
        relation: The relation name.
        acc_1: The first accession.
        acc_2: The second accession.
        max_pages: The largest number of pages to read.

    Returns:
        The search count and the items, each with its PMID, its title and its marked passage.
    """
    items: list[dict] = []
    count, total_pages = 0, 1
    page = 1
    while page <= min(max_pages, total_pages):
        payload = common.get_json(client, SEARCH_URL,
                                  {"text": common.pair_query(acc_1, acc_2, relation),
                                   "page": page})
        count = int(payload.get("count", 0) or 0)
        total_pages = int(payload.get("total_pages", 1) or 1)
        for result in payload.get("results", []):
            items.append({"pmid": str(result.get("pmid", "")),
                          "title": result.get("title", "") or "",
                          "text_hl": result.get("text_hl", "") or ""})
        page += 1
        time.sleep(common.DEFAULT_PAUSE_S)
    return {"search_count": count, "items": items}


def save_cache(path: Path, pairs: dict, fetched_at: str, refused: list[str]) -> None:
    """Write the sentence cache.

    The cache is written after every `SAVE_EVERY` pairs, so a run that stops keeps what it read.

    Args:
        path: The cache file.
        pairs: The cached pairs.
        fetched_at: The time of the first fetch.
        refused: The pair keys that the endpoint refused.
    """
    common.write_json(path, {"fetched_at": fetched_at, "endpoint": SEARCH_URL,
                             "refused_pairs": refused, "pairs": pairs})


def fetch_one(spec: tuple[str, str, str, str]) -> tuple[str, dict | None, str]:
    """Retrieve one pair with its own client, for a worker pool.

    Args:
        spec: The cache key, the relation and the two accessions.

    Returns:
        The cache key, the pages, and the error name. The pages are None when the endpoint
        refused the query, and the error name is empty on success.
    """
    key, relation, acc_1, acc_2 = spec
    with common.make_client() as client:
        try:
            return key, fetch_pages(client, relation, acc_1, acc_2, MAX_PAGES_ALL), ""
        except httpx.HTTPError as error:
            return key, None, type(error).__name__


def run_fetch(args: argparse.Namespace) -> int:
    """Fill the sentence cache of every literature edge of the disease node.

    A pair that the endpoint refuses is recorded and skipped. PubTator3 returns a gene accession
    that holds an error string for a small number of records, for example
    `@GENE_ERROR_OCCURRED:_CANNOT_GET_DOCUMENT_SUMMARY`, and the search endpoint answers a query on
    such an accession with HTTP 400. The edge stays in the graph and the run reports it, because a
    resource defect is not a reason to lose the sentences of every other edge.

    Args:
        args: The parsed command line.

    Returns:
        0 on success, 1 when an edge endpoint has no accession at all.
    """
    cache = common.read_json(args.cache) or {}
    pairs = cache.get("pairs", {})
    refused: list[str] = list(cache.get("refused_pairs", []))
    fetched_at = (cache.get("fetched_at")
                  or dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"))
    accessions = accession_by_node(args.entities)
    unresolved: list[str] = []
    todo: list[tuple[str, str, str, str]] = []
    seen: set[str] = set()
    for edge in literature_edges(args.edges, args.seed):
        acc_1 = accession_for(edge["source_id"], accessions)
        acc_2 = accession_for(edge["target_id"], accessions)
        if not acc_1 or not acc_2:
            unresolved.append(f"{edge['source_id']}|{edge['target_id']}")
            continue
        key = f"{edge['relation']}|{acc_1}|{acc_2}|{MAX_PAGES_ALL}"
        if key in pairs or key in refused or key in seen:
            continue
        seen.add(key)
        todo.append((key, edge["relation"], acc_1, acc_2))
    print(f"  {len(todo)} pairs to fetch on {args.workers} workers")
    fetched = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for key, pages, error in pool.map(fetch_one, todo):
            if pages is None:
                refused.append(key)
                print(f"  refused by the endpoint: {key} ({error})")
            else:
                pairs[key] = pages
            fetched += 1
            if fetched % SAVE_EVERY == 0:
                save_cache(args.cache, pairs, fetched_at, refused)
                print(f"    {fetched} of {len(todo)} pairs fetched, {len(pairs)} cached, "
                      f"{len(refused)} refused")
    save_cache(args.cache, pairs, fetched_at, refused)
    sentences = sum(len(pair["items"]) for pair in pairs.values())
    print(f"seed {args.seed}: {len(pairs)} cached pairs, {sentences} sentences, "
          f"{len(refused)} pairs refused by the endpoint")
    if unresolved:
        print(f"no accession for {len(unresolved)} edge endpoints: {', '.join(unresolved[:5])}")
    return 0


def sentence_rows(verdicts: dict, entities: dict, rule: seed_filter.SurfaceRule) -> list[dict]:
    """Build one table row per cached sentence.

    Args:
        verdicts: The records of `seed_filter.pair_verdicts`.
        entities: The resolved entity cache, keyed by accession.
        rule: The rule.

    Returns:
        The rows, in cache order.
    """
    rows: list[dict] = []
    for (relation, acc_1, acc_2), record in sorted(verdicts.items()):
        partner_accession = acc_2 if acc_1 == rule.seed_accession else acc_1
        partner = seed_filter.node_identifier(partner_accession, entities) or partner_accession
        for index, item in enumerate(record["items"]):
            forms = item["forms"] or [""]
            rows.append({"relation": relation, "partner": partner, "item_index": index,
                         "pmid": item["pmid"],
                         "surface_form": " | ".join(form for form in forms if form) or NO_FORM,
                         "rule_accepts": "yes" if item["accepted"] else "no",
                         "reason": item["reason"], "sentence": item["sentence"]})
    return rows


def form_records(verdicts: dict, rule: seed_filter.SurfaceRule) -> list[dict]:
    """Count every observed surface form and give it the verdict of the rule.

    Args:
        verdicts: The records of `seed_filter.pair_verdicts`.
        rule: The rule.

    Returns:
        One record per surface form, in descending order of mention count.
    """
    mentions: Counter[str] = Counter()
    accepted: Counter[str] = Counter()
    for record in verdicts.values():
        for item in record["items"]:
            if not item["forms"]:
                mentions[""] += 1
                continue
            for form in item["forms"]:
                mentions[seed_filter.normalize_form(form)] += 1
            for form in item["accepted_forms"]:
                accepted[seed_filter.normalize_form(form)] += 1
    out: list[dict] = []
    for form, count in mentions.most_common():
        keep = accepted.get(form, 0)
        if not form:
            verdict, generator = "rejected", "no marked mention"
        elif keep == count:
            verdict = "accepted"
            generator = ("mesh_entry_term" if seed_filter.long_form_match(form, rule)
                         else "abbreviation")
        elif keep == 0:
            verdict = "rejected"
            generator = ("abbreviation" if seed_filter.abbreviation_match(form, rule)
                         else "no match")
        else:
            verdict, generator = "conditional", "abbreviation"
        out.append({"form": form, "mentions": count, "accepted": keep, "verdict": verdict,
                    "decided_by": generator})
    return out


def read_reader_verdicts(path: Path) -> dict[tuple[str, str, int], bool]:
    """Read the sentence verdicts that a person wrote.

    Args:
        path: A TSV with the columns `relation`, `partner`, `item_index` and `names_syndrome`.

    Returns:
        The verdict of each sentence, keyed by relation, partner and item index.
    """
    out: dict[tuple[str, str, int], bool] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["relation"], row["partner"], int(row["item_index"]))
            out[key] = row["names_syndrome"].strip().lower() == "yes"
    return out


def agreement(rows: list[dict], reader: dict[tuple[str, str, int], bool]) -> dict:
    """Compare the verdict of the rule with the verdict of the reader.

    The reader is the reference. A false accept is a sentence that the rule accepts and the reader
    rejects. A false reject is a sentence that the rule rejects and the reader accepts.

    Args:
        rows: The sentence rows of `sentence_rows`.
        reader: The verdicts of `read_reader_verdicts`.

    Returns:
        The counts, the three rates and the disagreeing rows.
    """
    matched, agree = 0, 0
    false_accept, false_reject = [], []
    reader_yes = 0
    for row in rows:
        key = (row["relation"], row["partner"], int(row["item_index"]))
        if key not in reader:
            continue
        matched += 1
        truth = reader[key]
        reader_yes += 1 if truth else 0
        call = row["rule_accepts"] == "yes"
        if call == truth:
            agree += 1
        elif call:
            false_accept.append(row)
        else:
            false_reject.append(row)
    reader_no = matched - reader_yes
    return {
        "sentences_compared": matched,
        "reader_accepts": reader_yes,
        "reader_rejects": reader_no,
        "rule_accepts": sum(1 for row in rows
                            if (row["relation"], row["partner"], int(row["item_index"])) in reader
                            and row["rule_accepts"] == "yes"),
        "agreements": agree,
        "agreement_rate": round(agree / matched, 4) if matched else None,
        "false_accepts": len(false_accept),
        "false_accept_rate": round(len(false_accept) / reader_no, 4) if reader_no else None,
        "false_rejects": len(false_reject),
        "false_reject_rate": round(len(false_reject) / reader_yes, 4) if reader_yes else None,
        "false_accept_forms": dict(Counter(row["surface_form"] for row in false_accept)
                                   .most_common()),
        "false_reject_forms": dict(Counter(row["surface_form"] for row in false_reject)
                                   .most_common()),
        "false_reject_reasons": dict(Counter(row["reason"] for row in false_reject).most_common()),
    }


def edge_agreement(verdicts: dict, entities: dict, rule: seed_filter.SurfaceRule,
                   reader: dict[tuple[str, str, int], bool]) -> dict:
    """Compare the edge verdict of the rule with the edge verdict that the reader implies.

    An edge stays when at least one of its sentences is accepted, so the reader verdict of an edge
    is the disjunction of the reader verdicts of its sentences.

    Args:
        verdicts: The records of `seed_filter.pair_verdicts`.
        entities: The resolved entity cache, keyed by accession.
        rule: The rule.
        reader: The verdicts of `read_reader_verdicts`.

    Returns:
        The counts and the two rates.
    """
    compared, agree, false_keep, false_drop = 0, 0, 0, 0
    reader_keeps = 0
    for (relation, acc_1, acc_2), record in verdicts.items():
        partner_accession = acc_2 if acc_1 == rule.seed_accession else acc_1
        partner = seed_filter.node_identifier(partner_accession, entities) or partner_accession
        keys = [(relation, partner, index) for index in range(len(record["items"]))]
        known = [reader[key] for key in keys if key in reader]
        if not known:
            continue
        compared += 1
        truth = any(known)
        reader_keeps += 1 if truth else 0
        call = record["verdict"] == "keep"
        if call == truth:
            agree += 1
        elif call:
            false_keep += 1
        else:
            false_drop += 1
    return {"edges_compared": compared, "reader_keeps": reader_keeps,
            "rule_keeps": sum(1 for record in verdicts.values() if record["verdict"] == "keep"),
            "agreements": agree,
            "agreement_rate": round(agree / compared, 4) if compared else None,
            "false_keeps": false_keep, "false_drops": false_drop}


def write_forms_into_rule(rule_path: Path, forms: list[dict], counts: dict) -> Path:
    """Write the observed surface forms with a count and a verdict back into the rule file.

    Args:
        rule_path: The rule file.
        forms: The records of `form_records`.
        counts: The summary counts of the run.

    Returns:
        The rule file path.
    """
    data = yaml.safe_load(rule_path.read_text()) or {}
    data["audit"] = counts
    data["audit_forms"] = forms
    rule_path.write_text(yaml.safe_dump(data, sort_keys=False, allow_unicode=True, width=98))
    return rule_path


def run_apply(args: argparse.Namespace) -> int:
    """Decide every cached sentence of the disease node and write the tables.

    Args:
        args: The parsed command line.

    Returns:
        0 on success, 1 when the sentence cache is missing.
    """
    cache = common.read_json(args.cache) or {}
    if not cache.get("pairs"):
        print(f"FAIL: no cached sentence at {args.cache}.")
        return 1
    rule = seed_filter.load_rule(args.rule)
    verdicts = seed_filter.pair_verdicts(cache["pairs"], rule)
    entities = (common.read_json(args.entities) or {}).get("entities", {})
    rows = sentence_rows(verdicts, entities, rule)
    common.write_tsv(args.sentences_out, SENTENCE_HEADER, rows)
    verdict_rows = [{"relation": record["relation"],
                     "source_accession": record["source_accession"],
                     "target_accession": record["target_accession"],
                     "sentences": record["sentences"],
                     "accepted_sentences": record["accepted_sentences"],
                     "verdict": record["verdict"],
                     "accepted_forms": " | ".join(record["accepted_forms"]),
                     "rejected_forms": " | ".join(record["rejected_forms"])}
                    for record in sorted(verdicts.values(),
                                         key=lambda r: (r["verdict"], -r["accepted_sentences"],
                                                        r["relation"], r["source_accession"],
                                                        r["target_accession"]))]
    common.write_tsv(args.verdicts_out, seed_filter.VERDICT_HEADER, verdict_rows)
    forms = form_records(verdicts, rule)
    counts = {
        "run_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "sentence_table": seed_filter.repo_path(args.sentences_out),
        "verdict_table": seed_filter.repo_path(args.verdicts_out),
        "edges": len(verdicts),
        "edges_kept": sum(1 for record in verdicts.values() if record["verdict"] == "keep"),
        "edges_dropped": sum(1 for record in verdicts.values() if record["verdict"] == "drop"),
        "sentences": len(rows),
        "sentences_accepted": sum(1 for row in rows if row["rule_accepts"] == "yes"),
        "distinct_surface_forms": len(forms),
        "rejection_reasons": dict(Counter(row["reason"] for row in rows
                                          if row["rule_accepts"] == "no").most_common()),
    }
    if args.update_rule:
        write_forms_into_rule(args.rule, forms, counts)
    report = dict(counts)
    if args.ground_truth and args.ground_truth.is_file():
        report["sentence_agreement"] = agreement(rows, read_reader_verdicts(args.ground_truth))
        report["edge_agreement"] = edge_agreement(verdicts, entities, rule,
                                                  read_reader_verdicts(args.ground_truth))
        report["ground_truth"] = seed_filter.repo_path(args.ground_truth)
    report["rule_file"] = rule.source
    report["rule_sha256_prefix"] = rule.digest[:16]
    report["seed_node"] = rule.seed_node
    if args.report_out:
        common.write_json(args.report_out, report)
    print(f"{rule.seed_node}: {counts['edges']} literature edges, {counts['sentences']} sentences, "
          f"{counts['sentences_accepted']} accepted, {counts['edges_kept']} edges kept")
    if "sentence_agreement" in report:
        block = report["sentence_agreement"]
        print(f"  reader agreement {block['agreement_rate']} over {block['sentences_compared']} "
              f"sentences, false accepts {block['false_accepts']} "
              f"(rate {block['false_accept_rate']}), false rejects {block['false_rejects']} "
              f"(rate {block['false_reject_rate']})")
        print(f"  edge agreement {report['edge_agreement']['agreement_rate']} over "
              f"{report['edge_agreement']['edges_compared']} edges")
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Build the command line parser.

    Returns:
        The parser.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    steps = parser.add_subparsers(dest="step", required=True)
    fetch = steps.add_parser("fetch", help="Cache every sentence of every edge of the node.")
    fetch.add_argument("--seed", required=True, help="The disease node identifier.")
    fetch.add_argument("--edges", type=Path, required=True, help="The graph edges.tsv.")
    fetch.add_argument("--entities", type=Path, required=True, help="The resolved entity cache.")
    fetch.add_argument("--cache", type=Path, required=True, help="The sentence cache to fill.")
    fetch.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                       help="The number of concurrent requests to the search endpoint.")
    apply_step = steps.add_parser("apply", help="Decide every cached sentence.")
    apply_step.add_argument("--rule", type=Path, required=True, help="The rule file.")
    apply_step.add_argument("--cache", type=Path, required=True, help="The sentence cache.")
    apply_step.add_argument("--entities", type=Path, required=True,
                            help="The resolved entity cache.")
    apply_step.add_argument("--sentences-out", type=Path, required=True,
                            help="The sentence table to write.")
    apply_step.add_argument("--verdicts-out", type=Path, required=True,
                            help="The edge verdict table to write.")
    apply_step.add_argument("--report-out", type=Path, default=None,
                            help="The JSON summary to write.")
    apply_step.add_argument("--ground-truth", type=Path, default=None,
                            help="A TSV of reader verdicts to measure the rule against.")
    apply_step.add_argument("--update-rule", action="store_true",
                            help="Write the observed surface forms back into the rule file.")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run one step of the mention audit.

    Args:
        argv: The command line arguments.

    Returns:
        The exit status.
    """
    args = build_parser().parse_args(argv)
    if args.step == "fetch":
        return run_fetch(args)
    return run_apply(args)


if __name__ == "__main__":
    sys.exit(main())
