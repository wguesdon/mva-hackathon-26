"""Derive the accepted surface forms of a disease node from the public MeSH synonym set.

`literature/track2_seed_disease_surface_forms_v1.yaml` holds one list for one node. That
list does not transfer to another disease, so the correction that it carries covers one node of the
graph. This module removes the hand step. It reads the MeSH entry terms of a concept from the NCBI
E-utilities `esummary` endpoint and writes a rule file in the schema that
`scripts/track2_seed_filter.py` reads.

Four generators build the rule from the entry term list.

* G1 keeps the entry term itself, normalized to lower case with single spaces.
* G2 turns a MeSH inverted term such as "Fibrosis, Cystic" into "cystic fibrosis".
* G3 keeps the content tokens of the term, which are the tokens that are neither a generic disease
  word of `GENERIC_WORDS` nor a number. A surface form matches a content token set when one token
  set contains the other and the shared part is large enough to be distinctive.
* G4 builds the abbreviations of the term, which are its initialism and any short all upper case
  token that the term already holds. An abbreviation is ambiguous by construction, so a mention that
  matches only an abbreviation needs a corroboration in its own publication.

The module writes no surface form of its own. Every accepted form comes from the MeSH record, and
the rule file records the MeSH unique identifier, the retrieval time and the entry term list, so a
reader checks the rule against the source.

Usage:
    uv run python scripts/track2_disease_surface_forms.py --mesh-id D003550 \\
        --node DISEASE:MESH:D003550 --accession @DISEASE_Cystic_Fibrosis \\
        --out literature/track2_surface_forms/cystic_fibrosis_v1.yaml
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402

ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
DEFAULT_OUT_DIR = common.ROOT / "literature" / "track2_surface_forms"
RULE_SCHEMA_VERSION = 2

# A word that names a class of disease and not a disease. The list is the same for every disease,
# so no entry here is chosen for one concept. A content token set drops these words, and an
# abbreviation is corroborated when one of them stands beside it.
GENERIC_WORDS = (
    "syndrome", "syndromes", "disease", "diseases", "disorder", "disorders", "condition",
    "conditions", "deficiency", "deficiencies", "defect", "defects", "type", "types", "variant",
    "variants", "form", "forms", "subtype", "subtypes", "autosomal", "recessive", "dominant",
    "linked", "congenital", "hereditary", "inherited", "familial", "idiopathic", "primary",
    "secondary", "classic", "classical", "atypical", "progressive", "nonprogressive", "childhood",
    "juvenile", "infantile", "neonatal", "adult", "early", "late", "onset", "chronic", "acute",
    "severe", "mild", "moderate", "and", "or", "of", "the", "with", "without", "due", "to", "in",
    "a", "an", "non",
)

# The number of shared content tokens that a match needs. One token is enough when the token is
# long, because a long token is distinctive on its own.
MIN_SHARED_TOKENS = 2
LONG_TOKEN_LENGTH = 8

# A term this short is an abbreviation and not a name.
MAX_ABBREVIATION_LENGTH = 6
MIN_ABBREVIATION_LENGTH = 2

_SEPARATORS = re.compile(r"[\-/(),.;:\[\]{}'‘’“”]+")
_UPPER_TOKEN = re.compile(r"^[A-Z][A-Z0-9]{1,5}$")


def mesh_uid(mesh_id: str) -> str:
    """Return the NCBI unique identifier of a MeSH identifier.

    The MeSH database keys a descriptor `D003550` as `68003550` and a supplementary concept record
    `C536987` as `67536987`.

    Args:
        mesh_id: The MeSH identifier, with or without a `MESH:` prefix.

    Returns:
        The unique identifier.

    Raises:
        ValueError: When the identifier is neither a descriptor nor a supplementary record.
    """
    body = mesh_id.split(":")[-1].strip().upper()
    if body.startswith("D"):
        return f"68{body[1:]}"
    if body.startswith("C"):
        return f"67{body[1:]}"
    raise ValueError(f"{mesh_id} is neither a MeSH descriptor nor a supplementary concept record.")


def fetch_entry_terms(mesh_id: str) -> tuple[list[str], str]:
    """Read the entry terms of a MeSH concept.

    Args:
        mesh_id: The MeSH identifier, for example `D003550`.

    Returns:
        The entry terms in the order that the record gives them, and the retrieval time.

    Raises:
        ValueError: When the record holds no entry term.
    """
    uid = mesh_uid(mesh_id)
    with common.make_client() as client:
        payload = common.get_json(client, ESUMMARY_URL,
                                  {"db": "mesh", "id": uid, "retmode": "json"})
    record = (payload.get("result") or {}).get(uid) or {}
    terms = [str(term).strip() for term in record.get("ds_meshterms", []) if str(term).strip()]
    if not terms:
        raise ValueError(f"MeSH record {mesh_id} returned no entry term.")
    return terms, dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def normalize(text: str) -> str:
    """Fold a term or a surface form to one comparable shape.

    Args:
        text: The term or the matched text of a mention.

    Returns:
        The text in lower case, with every separator read as a space and with single spaces.
    """
    return " ".join(_SEPARATORS.sub(" ", text or "").lower().split())


def uninvert(term: str) -> str:
    """Turn a MeSH inverted term into the order that an author writes.

    Args:
        term: The entry term, for example `Fibrosis, Cystic`.

    Returns:
        The uninverted term, or the normalized term when it holds no comma.
    """
    parts = [part.strip() for part in term.split(",") if part.strip()]
    if len(parts) < 2:
        return normalize(term)
    return normalize(" ".join(reversed(parts)))


def content_tokens(text: str) -> list[str]:
    """Return the tokens of a term that name the disease itself.

    Args:
        text: The term or the matched text of a mention.

    Returns:
        The tokens that are neither a generic disease word nor a number, in order and without a
        duplicate.
    """
    kept: list[str] = []
    for token in normalize(text).split():
        if token in GENERIC_WORDS or token.isdigit():
            continue
        if token not in kept:
            kept.append(token)
    return kept


def is_abbreviation(term: str) -> bool:
    """Say whether an entry term is an abbreviation and not a name.

    Args:
        term: The entry term.

    Returns:
        True when the normalized term is short enough to be an abbreviation.
    """
    folded = normalize(term).replace(" ", "")
    return MIN_ABBREVIATION_LENGTH <= len(folded) <= MAX_ABBREVIATION_LENGTH


def initialisms(term: str) -> list[str]:
    """Build the initialisms of one entry term.

    Two families of initialism come from a term. A content initialism takes the first letter of
    every content token, which is how an author builds the abbreviation of a disease name. A full
    initialism takes the first letter of every word token as well, which is how an abbreviation
    such as NBS for Nijmegen breakage syndrome is built. A full initialism of two letters is
    dropped, because a two letter string built from a generic word carries no meaning.

    A MeSH inverted term is read in both orders, so "Muscular Dystrophy, Duchenne" gives DMD as
    well as MDD.

    Args:
        term: The entry term.

    Returns:
        The initialisms, without a duplicate.
    """
    out: list[str] = []
    for surface_form in (normalize(term), uninvert(term)):
        content = "".join(token[0] for token in content_tokens(surface_form))
        full = "".join(token[0] for token in surface_form.split()
                       if token and not token.isdigit())
        for candidate, floor in ((content, MIN_ABBREVIATION_LENGTH), (full, 3)):
            if floor <= len(candidate) <= MAX_ABBREVIATION_LENGTH and candidate not in out:
                out.append(candidate)
    return out


def embedded_abbreviations(term: str) -> list[str]:
    """Read the short upper case tokens that a term already holds.

    Args:
        term: The entry term, for example `MVA Syndrome`.

    Returns:
        The tokens in lower case, without a duplicate.
    """
    out: list[str] = []
    for token in _SEPARATORS.sub(" ", term).split():
        if _UPPER_TOKEN.match(token) and token.lower() not in out:
            out.append(token.lower())
    return out


def token_sets(terms: list[str]) -> list[dict]:
    """Build the content token set of every entry term and of its uninverted form.

    Args:
        terms: The entry terms.

    Returns:
        One record per distinct token set, with the terms that produced it.
    """
    by_key: dict[tuple[str, ...], dict] = {}
    for term in terms:
        if is_abbreviation(term):
            continue
        for surface in (normalize(term), uninvert(term)):
            tokens = tuple(sorted(content_tokens(surface)))
            if not tokens:
                continue
            record = by_key.setdefault(tokens, {"tokens": list(tokens), "terms": []})
            if term not in record["terms"]:
                record["terms"].append(term)
    return list(by_key.values())


def abbreviation_set(terms: list[str]) -> list[dict]:
    """Build the abbreviations of the concept.

    Args:
        terms: The entry terms.

    Returns:
        One record per abbreviation, with the terms that produced it and the generator that found
        it.
    """
    by_form: dict[str, dict] = {}
    for term in terms:
        found = [(form, "embedded") for form in embedded_abbreviations(term)]
        if is_abbreviation(term):
            folded = normalize(term).replace(" ", "")
            found.append((folded, "entry_term"))
        else:
            found.extend((form, "initialism") for form in initialisms(term))
        for form, origin in found:
            record = by_form.setdefault(form, {"form": form, "origin": origin, "terms": []})
            if term not in record["terms"]:
                record["terms"].append(term)
    return [by_form[form] for form in sorted(by_form)]


def build_rule(mesh_id: str, node: str, accession: str, label: str, terms: list[str],
               fetched_at: str) -> dict:
    """Build the rule record of one disease node.

    Args:
        mesh_id: The MeSH identifier of the concept.
        node: The graph node identifier.
        accession: The PubTator3 accession.
        label: The label of the node.
        terms: The MeSH entry terms.
        fetched_at: The retrieval time of the entry terms.

    Returns:
        The rule record, ready to write as YAML.
    """
    return {
        "version": RULE_SCHEMA_VERSION,
        "created": dt.datetime.now(dt.timezone.utc).date().isoformat(),
        "seed_node": node,
        "seed_accession": accession,
        "seed_label": label,
        "source": {
            "database": "MeSH",
            "endpoint": ESUMMARY_URL,
            "mesh_id": mesh_id,
            "mesh_uid": mesh_uid(mesh_id),
            "field": "ds_meshterms",
            "fetched_at": fetched_at,
            "entry_terms": terms,
        },
        "policy": {
            "default": "reject",
            "description": (
                "A mention is accepted when the content tokens of its matched text and the content "
                "tokens of one entry term contain one another and share enough tokens to be "
                "distinctive. A mention whose matched text is only an abbreviation of the concept "
                "needs a corroboration in its own publication."),
            "edge_rule": (
                "An edge of the disease node stays in the graph when at least one of its cached "
                "sentences keeps an accepted mention. An edge with no cached sentence is removed, "
                "because no person can read it."),
            "min_shared_tokens": MIN_SHARED_TOKENS,
            "long_token_length": LONG_TOKEN_LENGTH,
            "generic_words": list(GENERIC_WORDS),
        },
        "token_sets": token_sets(terms),
        "abbreviation_patterns": abbreviation_set(terms),
        "accept_patterns": [],
        "deny_forms": [],
        "other_sense_patterns": [],
        "audit_forms": [],
    }


def main() -> int:
    """Fetch the MeSH entry terms of a disease node and write its rule file.

    Returns:
        0 on success.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mesh-id", required=True, help="The MeSH identifier, for example D003550.")
    parser.add_argument("--node", required=True, help="The graph node identifier.")
    parser.add_argument("--accession", required=True, help="The PubTator3 accession.")
    parser.add_argument("--label", default="", help="The label of the node.")
    parser.add_argument("--out", type=Path, default=None, help="The rule file to write.")
    args = parser.parse_args()

    terms, fetched_at = fetch_entry_terms(args.mesh_id)
    label = args.label or terms[0]
    rule = build_rule(args.mesh_id, args.node, args.accession, label, terms, fetched_at)
    out = args.out or DEFAULT_OUT_DIR / f"{args.mesh_id.split(':')[-1].lower()}_v1.yaml"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(rule, sort_keys=False, allow_unicode=True, width=98))
    print(f"{label}: {len(terms)} entry terms, {len(rule['token_sets'])} token sets, "
          f"{len(rule['abbreviation_patterns'])} abbreviations")
    for record in rule["abbreviation_patterns"]:
        print(f"  abbreviation {record['form']} from {record['origin']}")
    print(f"rule file: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
