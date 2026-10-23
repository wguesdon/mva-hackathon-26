#!/usr/bin/env python3
"""Print the title and the abstract of each PubMed identifier.

`scripts/t2_search_pubmed.py` returns the title of a search result. An author who has to read a
measured quantity out of a publication needs the abstract as well, and the abstract of a
pharmacokinetic study or of a tumour incidence study carries the number in most cases. This script
reads one `efetch` response and prints the journal, the year, the title and every labelled abstract
section of each record.

Usage:
    uv run python scripts/t2_pubmed_abstracts.py 31685720 36002548
"""

from __future__ import annotations

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import t2_http_retry  # noqa: E402

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def node_text(node: ET.Element) -> str:
    """Return the text of one element and of every element inside it.

    Args:
        node: The element to read.

    Returns:
        The text with no markup.
    """
    return "".join(node.itertext()).strip()


def record_lines(article: ET.Element) -> list[str]:
    """Return the printed lines of one PubMed record.

    Args:
        article: The `PubmedArticle` element.

    Returns:
        One line per field and one line per abstract section.
    """
    pmid = article.findtext(".//PMID") or ""
    journal = article.findtext(".//Journal/ISOAbbreviation") or ""
    year = article.findtext(".//JournalIssue/PubDate/Year") or ""
    title_node = article.find(".//ArticleTitle")
    title = node_text(title_node) if title_node is not None else ""
    lines = [f"=== PMID {pmid} | {journal} {year} ===", title]
    for section in article.findall(".//Abstract/AbstractText"):
        label = section.get("Label") or ""
        text = node_text(section)
        lines.append(f"[{label}] {text}" if label else text)
    if not article.findall(".//Abstract/AbstractText"):
        lines.append("[no abstract in the PubMed record]")
    return lines


def fetch(pmids: list[str]) -> str:
    """Read the PubMed records of the identifiers.

    Args:
        pmids: The PubMed identifiers.

    Returns:
        The raw XML of the `efetch` response.
    """
    params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}
    response = t2_http_retry.get(EFETCH, params, timeout=90.0)
    return response.text


def main() -> int:
    """Print every requested record.

    Returns:
        Zero when every identifier is read.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pmids", nargs="+")
    args = parser.parse_args()
    root = ET.fromstring(fetch(args.pmids))
    for article in root.findall(".//PubmedArticle"):
        print("\n".join(record_lines(article)))
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
