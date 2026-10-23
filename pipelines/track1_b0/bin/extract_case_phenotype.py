#!/usr/bin/env python3
"""Map a clinical phenotype document to HPO terms using the Python standard library.

The source phenotype document is class C and remains outside model inputs.
The program writes derived terms as class B2 and a separate class A process gate.
Approved B2 processing routes include qualifying subscription sessions under the
organizers' conditions, as described in docs/llm_data_privacy.md.

The mapping uses dictionary lookup over labels and exact or narrow synonyms under
HP:0000118. Matches use normalized text and word boundaries. When matched terms
include an ancestor and its descendant, only the more specific term is kept.
No clinician reviewed the output, and the report states this.

Usage:
    extract_case_phenotype.py --document phenotype.docx --obo hp.obo \
        --out-terms case_hpo_terms.tsv --out-gate phenotype_extraction_gate.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ElementTree
import zipfile
from pathlib import Path

PHENOTYPIC_ABNORMALITY = "HP:0000118"
MIN_LABEL_CHARACTERS = 6
WORD = re.compile(r"[a-z0-9]+")
DOCX_NAMESPACE = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"

# Labels that are common English or that name a document section rather than a finding.
STOP_LABELS = frozenset({
    "abnormality", "abnormal", "all", "phenotypic abnormality", "growth abnormality",
    "clinical modifier", "mode of inheritance", "frequency", "neoplasm", "history", "past history",
    "family history", "physical examination", "laboratory", "abnormal cellular phenotype",
    "constitutional symptom", "abnormality of the cardiovascular system", "abnormal heart morphology",
})


def sha256_of(path: Path) -> str:
    """Hash a file.

    Args:
        path: The file.

    Returns:
        The hexadecimal SHA-256.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def docx_text(path: Path) -> str:
    """Read the visible text of a Word document.

    Args:
        path: The `.docx` file.

    Returns:
        The paragraph text, one paragraph per line.
    """
    lines = []
    with zipfile.ZipFile(path) as archive:
        for name in ("word/document.xml",) + tuple(
                item for item in archive.namelist()
                if item.startswith("word/") and item.endswith(".xml")
                and ("header" in item or "footer" in item)):
            if name not in archive.namelist():
                continue
            root = ElementTree.fromstring(archive.read(name))
            for paragraph in root.iter(f"{DOCX_NAMESPACE}p"):
                text = "".join(node.text or "" for node in paragraph.iter(f"{DOCX_NAMESPACE}t"))
                if text.strip():
                    lines.append(text.strip())
    return "\n".join(lines)


def plain_text(path: Path) -> str:
    """Read a text file.

    Args:
        path: The file.

    Returns:
        The text.
    """
    return path.read_text(errors="replace")


def normalize(text: str) -> str:
    """Lower the case, drop the punctuation and collapse the whitespace.

    Args:
        text: The input text.

    Returns:
        The normalized text, padded with one space at each end.
    """
    return " " + " ".join(WORD.findall(text.lower())) + " "


def parse_obo(path: Path) -> tuple[dict[str, set[str]], dict[str, str], dict[str, list[str]]]:
    """Read the parents, the names and the usable synonyms of every live term.

    Args:
        path: The `hp.obo` file.

    Returns:
        The parents of each term, the name of each term and the labels of each term.
    """
    parents: dict[str, set[str]] = {}
    names: dict[str, str] = {}
    labels: dict[str, list[str]] = {}
    term_id = ""
    in_term = False
    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                in_term, term_id = True, ""
                continue
            if line.startswith("["):
                in_term = False
                continue
            if not in_term or not line:
                continue
            if line.startswith("id: "):
                term_id = line[4:].strip()
                parents.setdefault(term_id, set())
                labels.setdefault(term_id, [])
            elif line.startswith("name: "):
                names[term_id] = line[6:].strip()
                labels.setdefault(term_id, []).append(line[6:].strip())
            elif line.startswith("is_a: "):
                parents.setdefault(term_id, set()).add(line[6:].split("!")[0].strip())
            elif line.startswith("synonym: ") and ('"EXACT"' in line or "EXACT " in line
                                                   or "NARROW " in line or '"NARROW"' in line):
                match = re.search(r'"([^"]+)"\s+(EXACT|NARROW)', line)
                if match:
                    labels.setdefault(term_id, []).append(match.group(1))
            elif line.startswith("is_obsolete: true"):
                parents.pop(term_id, None)
                names.pop(term_id, None)
                labels.pop(term_id, None)
    return parents, names, labels


def ancestor_closure(parents: dict[str, set[str]]) -> dict[str, frozenset[str]]:
    """Build the reflexive ancestor set of every term.

    Args:
        parents: The direct parents of each term.

    Returns:
        Each term mapped to itself and to every ancestor.
    """
    closure: dict[str, frozenset[str]] = {}
    pending = list(parents)
    while pending:
        progressed = False
        remaining = []
        for term in pending:
            live = [parent for parent in parents[term] if parent in parents]
            if all(parent in closure for parent in live):
                collected = {term}
                for parent in live:
                    collected.update(closure[parent])
                closure[term] = frozenset(collected)
                progressed = True
            else:
                remaining.append(term)
        if not progressed:
            for term in remaining:
                closure[term] = frozenset({term})
            break
        pending = remaining
    return closure


def build_lexicon(labels: dict[str, list[str]], closure: dict[str, frozenset[str]]) -> dict[str, str]:
    """Build the label to term map of the phenotypic abnormality branch.

    Args:
        labels: The labels of every term.
        closure: The reflexive ancestor set of every term.

    Returns:
        The normalized label mapped to the most specific term that carries it.
    """
    lexicon: dict[str, str] = {}
    for term, names in labels.items():
        if PHENOTYPIC_ABNORMALITY not in closure.get(term, frozenset()):
            continue
        for name in names:
            key = normalize(name).strip()
            if len(key) < MIN_LABEL_CHARACTERS or key in STOP_LABELS:
                continue
            current = lexicon.get(key)
            if current is None or len(closure[term]) > len(closure[current]):
                lexicon[key] = term
    return lexicon


def match_terms(text: str, lexicon: dict[str, str],
                closure: dict[str, frozenset[str]]) -> list[tuple[str, str]]:
    """Find the most specific HPO terms whose label occurs in the text.

    Args:
        text: The normalized document text.
        lexicon: The normalized label mapped to its term.
        closure: The reflexive ancestor set of every term.

    Returns:
        The term and the label that matched it, sorted by term identifier.
    """
    hits: dict[str, str] = {}
    for label, term in lexicon.items():
        if f" {label} " in text and (term not in hits or len(label) > len(hits[term])):
            hits[term] = label
    ancestors = {ancestor for term in hits for ancestor in closure[term] if ancestor != term}
    return sorted((term, label) for term, label in hits.items() if term not in ancestors)


def main() -> int:
    """Read the document, write the term list and write the gate file.

    Returns:
        Zero when the extraction produced at least one term, and 3 otherwise.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--document", required=True, type=Path)
    parser.add_argument("--obo", required=True, type=Path)
    parser.add_argument("--out-terms", required=True, type=Path)
    parser.add_argument("--out-gate", required=True, type=Path)
    args = parser.parse_args()

    text = docx_text(args.document) if args.document.suffix.lower() == ".docx" \
        else plain_text(args.document)
    parents, names, labels = parse_obo(args.obo)
    closure = ancestor_closure(parents)
    lexicon = build_lexicon(labels, closure)
    matches = match_terms(normalize(text), lexicon, closure)

    lines = ["hpo_id\thpo_label\tmatched_text\n"]
    lines += [f"{term}\t{names.get(term, '')}\t{label}\n" for term, label in matches]
    args.out_terms.parent.mkdir(parents=True, exist_ok=True)
    args.out_terms.write_text("".join(lines))

    gate = {
        "gate_id": "case_phenotype_extraction_v1",
        "data_class": "A",
        "scope": "This file carries a verdict and the checksums of the inputs. It carries no count "
                 "and no term of the child, because both state a fact about the child.",
        "document_sha256": sha256_of(args.document),
        "hp_obo_sha256": sha256_of(args.obo),
        "lexicon_labels": len(lexicon),
        "method": "A dictionary lookup of the HPO label and of the exact and narrow synonyms of "
                  "every term under HP:0000118, on word boundaries, with the ancestors of a matched "
                  "term removed. No clinician reviewed the output.",
        "verdict": "PASS" if matches else "FAIL",
        "condition": "The extraction produced at least one HPO term.",
        "out_terms": str(args.out_terms),
        "out_terms_sha256": sha256_of(args.out_terms),
    }
    args.out_gate.parent.mkdir(parents=True, exist_ok=True)
    args.out_gate.write_text(json.dumps(gate, indent=2) + "\n")
    print(json.dumps({"verdict": gate["verdict"], "lexicon_labels": gate["lexicon_labels"]}))
    return 0 if matches else 3


if __name__ == "__main__":
    sys.exit(main())
