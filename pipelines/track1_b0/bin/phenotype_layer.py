#!/usr/bin/env python3
"""Phenotype similarity layer for the Track 1 B0 ranker. Python standard library only.

The layer scores a gene by the similarity between a case phenotype profile and the phenotype profile
of that gene. The term similarity is the Resnik information content measure, which takes the
information content of the most informative common ancestor of two terms (Resnik 1995, IJCAI
proceedings, pages 448 to 453). The two profiles are combined by the best match average, which is
the mean over both directions of the best match of each term (Pesquita and colleagues, PLoS
Computational Biology 2009, PMID 19649320). The Phenomizer (PMID 19800049) and the Exomiser (PMID
26562621) use the same pair of choices.

The information content corpus is the gene to phenotype annotation of one Human Phenotype Ontology
release, propagated over the `is_a` closure of `hp.obo`. The information content of a term is

    IC(t) = -ln( genes annotated to t or to a descendant of t / genes annotated at all )

The score of a gene is the best match average divided by the largest information content of the
corpus, so the score lies in the closed interval from 0 to 1.

Every input of this module is public. The module reads no patient file. The caller supplies the case
profile, and a case profile of the child is class B2 content under ADR 0005.
"""

from __future__ import annotations

import csv
import hashlib
import math
from collections import defaultdict
from pathlib import Path

PHENOTYPIC_ABNORMALITY = "HP:0000118"
INHERITANCE_ROOT = "HP:0000005"


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


def parse_obo(path: Path) -> tuple[dict[str, set[str]], dict[str, str], dict[str, str]]:
    """Read the term identifiers, the `is_a` parents and the alternate identifiers of `hp.obo`.

    Args:
        path: The `hp.obo` file of one HPO release.

    Returns:
        The parents of each live term, the alternate identifier map and the name of each term.
    """
    parents: dict[str, set[str]] = {}
    alternates: dict[str, str] = {}
    names: dict[str, str] = {}
    term_id = ""
    obsolete = False
    in_term = False
    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                in_term, term_id, obsolete = True, "", False
                continue
            if line.startswith("["):
                in_term = False
                continue
            if not in_term or not line:
                continue
            if line.startswith("id: "):
                term_id = line[4:].strip()
                parents.setdefault(term_id, set())
            elif line.startswith("name: "):
                names[term_id] = line[6:].strip()
            elif line.startswith("is_a: "):
                parents.setdefault(term_id, set()).add(line[6:].split("!")[0].strip())
            elif line.startswith("alt_id: "):
                alternates[line[8:].strip()] = term_id
            elif line.startswith("is_obsolete: true"):
                obsolete = True
                parents.pop(term_id, None)
                names.pop(term_id, None)
            elif line.startswith("replaced_by: ") and obsolete:
                alternates[term_id] = line[13:].strip()
    return parents, alternates, names


def ancestor_closure(parents: dict[str, set[str]]) -> dict[str, frozenset[str]]:
    """Build the reflexive `is_a` ancestor set of every term.

    Args:
        parents: The direct parents of each term.

    Returns:
        Each term mapped to itself and to every ancestor.

    Raises:
        ValueError: When the parent graph holds a cycle.
    """
    closure: dict[str, frozenset[str]] = {}
    order = topological_order(parents)
    for term in order:
        collected = {term}
        for parent in parents.get(term, ()):
            collected.update(closure.get(parent, frozenset({parent})))
        closure[term] = frozenset(collected)
    if len(closure) != len(parents):
        raise ValueError("The parent graph holds a cycle.")
    return closure


def topological_order(parents: dict[str, set[str]]) -> list[str]:
    """Order the terms so that every parent comes before its children.

    Args:
        parents: The direct parents of each term.

    Returns:
        The term identifiers in a parent first order.
    """
    children: dict[str, list[str]] = defaultdict(list)
    pending = {term: 0 for term in parents}
    for term, direct in parents.items():
        live = [parent for parent in direct if parent in parents]
        pending[term] = len(live)
        for parent in live:
            children[parent].append(term)
    queue = [term for term, count in pending.items() if count == 0]
    order = []
    while queue:
        term = queue.pop()
        order.append(term)
        for child in children[term]:
            pending[child] -= 1
            if pending[child] == 0:
                queue.append(child)
    return order


def descendant_map(closure: dict[str, frozenset[str]]) -> dict[str, set[str]]:
    """Invert the ancestor closure into a descendant map.

    Args:
        closure: The reflexive ancestor set of every term.

    Returns:
        Each term mapped to itself and to every descendant.
    """
    descendants: dict[str, set[str]] = defaultdict(set)
    for term, ancestors in closure.items():
        for ancestor in ancestors:
            descendants[ancestor].add(term)
    return descendants


def normalize_term(term: str, closure: dict[str, frozenset[str]], alternates: dict[str, str]) -> str | None:
    """Map a term identifier to the live identifier of the release.

    Args:
        term: The term identifier.
        closure: The ancestor closure, used as the set of live terms.
        alternates: The alternate identifier map.

    Returns:
        The live identifier, or None when the release does not hold the term.
    """
    if term in closure:
        return term
    replacement = alternates.get(term)
    return replacement if replacement in closure else None


def load_gene_phenotypes(path: Path, closure: dict[str, frozenset[str]],
                         alternates: dict[str, str]) -> dict[str, set[str]]:
    """Read the phenotype profile of every gene from `genes_to_phenotype.txt`.

    Args:
        path: The `genes_to_phenotype.txt` file of one HPO release.
        closure: The ancestor closure of the release.
        alternates: The alternate identifier map.

    Returns:
        The gene symbol mapped to its set of live HPO term identifiers.
    """
    profiles: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            term = normalize_term(row["hpo_id"], closure, alternates)
            if term is not None:
                profiles[row["gene_symbol"]].add(term)
    return dict(profiles)


def information_content(profiles: dict[str, set[str]],
                        closure: dict[str, frozenset[str]]) -> dict[str, float]:
    """Compute the information content of every annotated term over the gene corpus.

    Args:
        profiles: The phenotype profile of every gene.
        closure: The ancestor closure of the release.

    Returns:
        Each term that the corpus annotates, mapped to its information content in nats.
    """
    counts: dict[str, int] = defaultdict(int)
    for terms in profiles.values():
        propagated: set[str] = set()
        for term in terms:
            propagated.update(closure[term])
        for term in propagated:
            counts[term] += 1
    total = len(profiles)
    return {term: -math.log(count / total) for term, count in counts.items() if count}


class PhenotypeModel:
    """A Resnik and best match average phenotype model over one HPO release.

    Attributes:
        closure: The reflexive ancestor set of every term.
        descendants: The reflexive descendant set of every term.
        content: The information content of every annotated term.
        profiles: The phenotype profile of every gene.
        max_content: The largest information content of the corpus.
    """

    def __init__(self, closure: dict[str, frozenset[str]], descendants: dict[str, set[str]],
                 content: dict[str, float], profiles: dict[str, set[str]],
                 alternates: dict[str, str]) -> None:
        """Store the ontology, the corpus and the gene profiles.

        Args:
            closure: The reflexive ancestor set of every term.
            descendants: The reflexive descendant set of every term.
            content: The information content of every annotated term.
            profiles: The phenotype profile of every gene.
            alternates: The alternate identifier map of the release.
        """
        self.closure = closure
        self.descendants = descendants
        self.content = content
        self.profiles = profiles
        self.alternates = alternates
        self.max_content = max(content.values()) if content else 0.0
        self._cache: dict[str, dict[str, float]] = {}

    def term_similarity_row(self, query_term: str) -> dict[str, float]:
        """Return the Resnik similarity of one query term against every term of the ontology.

        The similarity of two terms is the largest information content of a common ancestor. The row
        is built by a walk over the ancestors of the query term in descending information content,
        which assigns each descendant of an ancestor the information content of the most informative
        common ancestor that it shares with the query term.

        Args:
            query_term: A live HPO term identifier.

        Returns:
            Each term that shares an annotated ancestor with the query term, mapped to the
            similarity in nats.
        """
        cached = self._cache.get(query_term)
        if cached is not None:
            return cached
        row: dict[str, float] = {}
        ancestors = sorted((term for term in self.closure[query_term] if term in self.content),
                           key=lambda term: -self.content[term])
        for ancestor in ancestors:
            value = self.content[ancestor]
            for term in self.descendants[ancestor]:
                if term not in row:
                    row[term] = value
        self._cache[query_term] = row
        return row

    def normalize(self, terms: list[str]) -> list[str]:
        """Map a term list to the live identifiers of the release and drop the unknown ones.

        Args:
            terms: The term identifiers.

        Returns:
            The live identifiers, without duplicates and in the input order.
        """
        live: list[str] = []
        for term in terms:
            mapped = normalize_term(term, self.closure, self.alternates)
            if mapped is not None and mapped not in live:
                live.append(mapped)
        return live

    def score_gene(self, gene: str, case_terms: list[str]) -> float:
        """Score one gene against a case phenotype profile.

        Args:
            gene: The gene symbol.
            case_terms: The live HPO term identifiers of the case profile.

        Returns:
            The best match average divided by the largest information content, in the closed
            interval from 0 to 1. A gene with no annotation scores 0.0.
        """
        target = self.profiles.get(gene)
        if not target or not case_terms or not self.max_content:
            return 0.0
        rows = [self.term_similarity_row(term) for term in case_terms]
        forward = 0.0
        for row in rows:
            forward += max((row.get(term, 0.0) for term in target), default=0.0)
        backward = 0.0
        for term in target:
            backward += max((row.get(term, 0.0) for row in rows), default=0.0)
        average = 0.5 * (forward / len(rows) + backward / len(target))
        return min(1.0, average / self.max_content)


def build_model(obo: Path, genes_to_phenotype: Path,
                drop_rows: set[tuple[str, str]] | None = None) -> PhenotypeModel:
    """Build a phenotype model from one HPO release.

    Args:
        obo: The `hp.obo` file.
        genes_to_phenotype: The `genes_to_phenotype.txt` file.
        drop_rows: Gene and term pairs to remove from the gene profiles before the model is built.
            The information content corpus is built after the removal, so a held out annotation
            source changes neither the profile nor the corpus.

    Returns:
        The model.
    """
    parents, alternates, _ = parse_obo(obo)
    closure = ancestor_closure(parents)
    descendants = descendant_map(closure)
    profiles = load_gene_phenotypes(genes_to_phenotype, closure, alternates)
    if drop_rows:
        for gene, term in drop_rows:
            if gene in profiles:
                profiles[gene].discard(term)
        profiles = {gene: terms for gene, terms in profiles.items() if terms}
    content = information_content(profiles, closure)
    return PhenotypeModel(closure, descendants, content, profiles, alternates)


def read_case_profile(path: Path) -> list[str]:
    """Read a case phenotype profile from a text file of one HPO term per line.

    A comment line starts with `#`. The caller states the data class of the file. A profile of the
    child is class B2 content under ADR 0005.

    Args:
        path: The profile file.

    Returns:
        The term identifiers in file order.
    """
    terms = []
    for line in path.read_text().splitlines():
        value = line.split("\t")[0].strip()
        if value and not value.startswith("#"):
            terms.append(value)
    return terms
