#!/usr/bin/env python3
"""Check that every PMID of an evidence table resolves and carries its recorded keyword.

`scripts/t2_preclinical_protocol.py` applies this gate to the preclinical evidence table. The
pharmacokinetic table and the tumour endpoint table of this project carry the same shape and need
the same gate, so the two checks live here and each table script calls them.

The gate tests the citation and not the claim. A PMID that resolves to a real record whose title
holds the recorded keyword still says nothing about whether that publication supports the sentence
beside it, and every receipt that this module fills states that limit.
"""

from __future__ import annotations

import csv
import hashlib
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import t2_http_retry  # noqa: E402

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NON_WORD = re.compile(r"[^a-z0-9]+")


def normalize(text: str) -> str:
    """Return a title or a keyword in one comparable form.

    Args:
        text: The title or the keyword.

    Returns:
        The lower case text with every run of non alphanumeric characters read as one space.
    """
    return NON_WORD.sub(" ", text.lower()).strip()


def sha256_of_file(path: Path) -> str:
    """Return the SHA-256 of one file.

    Args:
        path: The file path.

    Returns:
        The hexadecimal digest.
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_table(path: Path, rows: list[dict]) -> str:
    """Write one evidence table and return its digest.

    Args:
        path: The output path.
        rows: The rows, each one a mapping with the same keys.

    Returns:
        The SHA-256 of the written file.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return sha256_of_file(path)


def pubmed_titles(pmids: list[str]) -> dict[str, str]:
    """Read the PubMed title of each identifier.

    Args:
        pmids: The PubMed identifiers.

    Returns:
        The title of each identifier that PubMed returns.
    """
    if not pmids:
        return {}
    params = {"db": "pubmed", "id": ",".join(sorted(set(pmids))), "retmode": "xml"}
    root = ET.fromstring(t2_http_retry.get(EFETCH, params, timeout=90.0).text)
    titles = {}
    for article in root.findall(".//PubmedArticle"):
        pmid = article.findtext(".//PMID") or ""
        node = article.find(".//ArticleTitle")
        titles[pmid] = "".join(node.itertext()).strip() if node is not None else ""
    return titles


def citation_findings(rows: list[dict], titles: dict[str, str]) -> list[dict]:
    """Return one finding per row whose citation fails the gate.

    Args:
        rows: The evidence rows, each with a `pmid` field and a `citation_keyword` field.
        titles: The PubMed title of each identifier.

    Returns:
        The findings, empty when every citation passes.
    """
    findings = []
    for row in rows:
        pmid = (row.get("pmid") or "").strip()
        if not pmid:
            continue
        title = titles.get(pmid, "")
        keyword = (row.get("citation_keyword") or "").strip()
        if not title:
            findings.append({"row": row.get("claim_id", ""), "check": "resolve",
                             "detail": f"PMID {pmid} returns no PubMed record"})
        elif not keyword:
            findings.append({"row": row.get("claim_id", ""), "check": "keyword",
                             "detail": f"PMID {pmid} carries no keyword"})
        elif normalize(keyword) not in normalize(title):
            findings.append({"row": row.get("claim_id", ""), "check": "keyword",
                             "detail": f"the keyword {keyword} is absent from the title of "
                                       f"PMID {pmid}"})
    return findings


def repository_commit() -> str:
    """Return the current commit of the repository.

    Returns:
        The short commit hash.
    """
    out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                         capture_output=True, text=True, check=True)
    return out.stdout.strip()
