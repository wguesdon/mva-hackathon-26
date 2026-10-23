"""Build a class A literature pack from Europe PMC for a research board.

The pack holds public metadata and abstracts only. It never reads patient data. Each run
writes a new directory under reports/literature/<pack_id>_<date>/ with:

* papers.tsv: one row per unique paper, with the topics that retrieved it.
* abstracts/<topic>.md: the title, citation and abstract of each paper in the topic.
* index.md: a table of topics with counts and links.
* manifest.json: the queries, the retrieval time, the counts and the SHA-256 of each file.

Usage:
    uv run python scripts/literature_pack.py literature/task1_queries.yaml
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import html
import json
import re
import sys
import time
from pathlib import Path

import httpx
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "harness"))
import privacy_guard  # noqa: E402

SEARCH_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
PAGE_SIZE = 100
PAUSE_S = 0.4
FIELDS = ("paper_id", "pmid", "pmcid", "doi", "year", "journal", "first_author", "title",
          "cited_by", "open_access", "pub_type", "topics")


def clean(text: str | None) -> str:
    """Remove markup and repeated whitespace from a Europe PMC text field.

    Args:
        text: The raw field.

    Returns:
        Plain text on one line.
    """
    if not text:
        return ""
    return " ".join(re.sub(r"<[^>]+>", " ", html.unescape(text)).split())


def search(client: httpx.Client, query: str, limit: int) -> list[dict]:
    """Retrieve up to `limit` core records for one query, in Europe PMC relevance order.

    Args:
        client: The HTTP client.
        query: The Europe PMC query.
        limit: The maximum number of records.

    Returns:
        The raw result records.
    """
    records: list[dict] = []
    cursor = "*"
    while len(records) < limit:
        params = {"query": f"({query}) AND (SRC:MED OR SRC:PMC OR SRC:PPR)", "format": "json",
                  "resultType": "core", "pageSize": min(PAGE_SIZE, limit - len(records)),
                  "cursorMark": cursor}
        response = client.get(SEARCH_URL, params=params)
        response.raise_for_status()
        payload = response.json()
        batch = payload.get("resultList", {}).get("result", [])
        records.extend(batch)
        next_cursor = payload.get("nextCursorMark")
        if not batch or not next_cursor or next_cursor == cursor:
            break
        cursor = next_cursor
        time.sleep(PAUSE_S)
    return records[:limit]


def paper_row(record: dict) -> dict:
    """Convert one Europe PMC core record to a pack row.

    Args:
        record: The raw record.

    Returns:
        The row with an abstract field.
    """
    journal = (record.get("journalInfo") or {}).get("journal", {}).get("title", "")
    authors = (record.get("authorList") or {}).get("author", [])
    first = authors[0].get("fullName", "") if authors else clean(record.get("authorString", "")).split(",")[0]
    types = (record.get("pubTypeList") or {}).get("pubType", [])
    paper_id = f"PMID:{record['pmid']}" if record.get("pmid") else f"{record.get('source')}:{record.get('id')}"
    return {"paper_id": paper_id, "pmid": record.get("pmid", ""), "pmcid": record.get("pmcid", ""),
            "doi": record.get("doi", ""), "year": record.get("pubYear", ""), "journal": clean(journal),
            "first_author": clean(first), "title": clean(record.get("title")),
            "cited_by": str(record.get("citedByCount", 0)), "open_access": record.get("isOpenAccess", ""),
            "pub_type": "; ".join(types[:3]), "abstract": clean(record.get("abstractText"))}


def citation_order(row: dict) -> tuple[int, str]:
    """Return the sort key that puts the most cited papers first.

    Args:
        row: The pack row.

    Returns:
        The negative citation count and the paper identifier.
    """
    return -int(row["cited_by"] or 0), row["paper_id"]


def citation(row: dict) -> str:
    """Format a short citation line.

    Args:
        row: The pack row.

    Returns:
        The citation text with identifiers.
    """
    ids = [f"PMID {row['pmid']}" if row["pmid"] else row["paper_id"]]
    if row["doi"]:
        ids.append(f"doi:{row['doi']}")
    return f"{row['first_author']} ({row['year']}). {row['journal']}. {', '.join(ids)}. Cited by {row['cited_by']}."


def write_pack(out: Path, config: dict, rows: dict[str, dict], topic_hits: dict[str, list[str]]) -> dict:
    """Write the pack files and the manifest.

    Args:
        out: The new pack directory.
        config: The query configuration.
        rows: The unique papers by paper identifier.
        topic_hits: The paper identifiers of each topic, in retrieval order.

    Returns:
        The manifest.
    """
    (out / "abstracts").mkdir(parents=True)
    with (out / "papers.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in sorted(rows.values(), key=citation_order):
            writer.writerow(row)
    index = ["# Literature pack", "", f"Pack `{out.name}` holds public Europe PMC records. Data class A.", "",
             "| Topic | Papers | File |", "| --- | ---: | --- |"]
    for topic in config["topics"]:
        ids = topic_hits[topic["id"]]
        lines = [f"# {topic['id']}", "", f"Query: `{topic['query']}`", "",
                 f"The query returned {len(ids)} papers in Europe PMC relevance order.", ""]
        for paper_id in ids:
            row = rows[paper_id]
            lines += [f"## {row['title']}", "", citation(row), "",
                      row["abstract"] or "The record has no abstract.", ""]
        (out / "abstracts" / f"{topic['id']}.md").write_text("\n".join(lines))
        index.append(f"| {topic['id']} | {len(ids)} | [abstracts/{topic['id']}.md](abstracts/{topic['id']}.md) |")
    index += ["", f"The pack holds {len(rows)} unique papers. Read `papers.tsv` for the full table."]
    (out / "index.md").write_text("\n".join(index) + "\n")
    files = sorted(path for path in out.rglob("*") if path.is_file())
    manifest = {"pack": out.name, "source": SEARCH_URL, "retrieved_utc": dt.datetime.now(dt.UTC).isoformat(),
                "unique_papers": len(rows), "topics": [dict(topic, returned=len(topic_hits[topic["id"]]))
                                                       for topic in config["topics"]],
                "files": {str(path.relative_to(out)): hashlib.sha256(path.read_bytes()).hexdigest() for path in files}}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> int:
    """Build one literature pack.

    Returns:
        Zero on success, two when a file fails the privacy scan.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("queries", type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(args.queries.read_text())
    stamp = dt.datetime.now(dt.UTC).strftime("%Y_%m_%d")
    out = ROOT / "reports/literature" / f"{config['pack_id']}_{stamp}"
    if out.exists():
        print(f"The pack directory exists: {out}. Change pack_id or remove the old pack.", file=sys.stderr)
        return 2
    rows: dict[str, dict] = {}
    topic_hits: dict[str, list[str]] = {}
    with httpx.Client(timeout=httpx.Timeout(60, connect=20), headers={"User-Agent": "mva-hackathon-2026"}) as client:
        for topic in config["topics"]:
            records = search(client, topic["query"], int(topic["max"]))
            topic_hits[topic["id"]] = []
            for record in records:
                row = paper_row(record)
                existing = rows.setdefault(row["paper_id"], dict(row, topics=""))
                existing["topics"] = "; ".join(filter(None, [existing["topics"], topic["id"]]))
                if row["paper_id"] not in topic_hits[topic["id"]]:
                    topic_hits[topic["id"]].append(row["paper_id"])
            print(f"{topic['id']}\t{len(records)}")
            time.sleep(PAUSE_S)
    manifest = write_pack(out, config, rows, topic_hits)
    failed = [name for name in manifest["files"] if privacy_guard.scan((out / name).read_text(), tokens=None)]
    print(f"The pack {out.relative_to(ROOT)} holds {manifest['unique_papers']} unique papers.")
    if failed:
        print(f"These files failed the privacy scan: {failed}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
