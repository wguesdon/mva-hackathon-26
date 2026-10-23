#!/usr/bin/env python3
"""Search PubMed for a query and print the title of each result.

An author who adds a citation to a report needs the PMID of the paper that supports the sentence. A
PMID written from memory points at the wrong paper often enough that this project checks every one.
This script runs one `esearch` query and prints the identifier, the journal, the year and the title
of each result, so the author selects the record and `scripts/build_bibliography.py` then applies the
keyword gate.

Usage:
    uv run python scripts/t2_search_pubmed.py "nicotinamide mononucleotide adenylyltransferase NAD" --retmax 10
"""

from __future__ import annotations

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import httpx  # noqa: F401

sys.path.insert(0, str(Path(__file__).resolve().parent))

import t2_http_retry  # noqa: E402

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def search(query: str, retmax: int) -> list[str]:
    """Run one PubMed search and return the identifiers.

    Args:
        query: The PubMed query string.
        retmax: The maximum number of results.

    Returns:
        The PubMed identifiers, most relevant first.
    """
    params = {"db": "pubmed", "term": query, "retmax": str(retmax), "retmode": "xml",
              "sort": "relevance"}
    response = t2_http_retry.get(ESEARCH, params, timeout=60.0)
    root = ET.fromstring(response.text)
    return [node.text or "" for node in root.findall(".//IdList/Id")]


def titles(pmids: list[str]) -> list[tuple[str, str, str, str]]:
    """Read the title, the journal and the year of each identifier.

    Args:
        pmids: The PubMed identifiers.

    Returns:
        One tuple per record with the PMID, the journal, the year and the title.
    """
    if not pmids:
        return []
    params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}
    response = t2_http_retry.get(EFETCH, params, timeout=60.0)
    root = ET.fromstring(response.text)
    out = []
    for article in root.findall(".//PubmedArticle"):
        pmid = article.findtext(".//PMID") or ""
        node = article.find(".//ArticleTitle")
        title = "".join(node.itertext()) if node is not None else ""
        journal = article.findtext(".//Journal/ISOAbbreviation") or ""
        year = article.findtext(".//JournalIssue/PubDate/Year") or ""
        out.append((pmid, journal, year, title))
    return out


def main() -> int:
    """Print one line per search result.

    Returns:
        Zero when the search returns at least one record.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("query")
    parser.add_argument("--retmax", type=int, default=10)
    args = parser.parse_args()
    records = titles(search(args.query, args.retmax))
    for pmid, journal, year, title in records:
        print(f"{pmid}\t{journal} {year}\t{title}")
    return 0 if records else 1


if __name__ == "__main__":
    sys.exit(main())
