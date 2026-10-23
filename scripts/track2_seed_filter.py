"""Decide which literature edge of the Track 2 seed disease node names the syndrome.

PubTator3 maps a stretch of text onto a concept. The seed node of the Track 2 graph is
`@DISEASE_Mosaic_variegated_aneuploidy_syndrome`, MeSH C536987. The recogniser also maps a general
description of a mitotic phenotype onto that concept, for example the words "mitotic abnormalities".
An edge that rests only on such a mention states nothing about the syndrome, and a ranking that
leaves the seed on such an edge ranks a phrase and not a disease.

This module applies the rule of `literature/track2_seed_disease_surface_forms_v1.yaml` to the cached
PubTator3 search results in `reports/track2/raw/edge_provenance_evidence.json`. Every verdict comes
from that list file. No surface form and no reason is written inside a function body here.

The command line writes a verdict table with one row per literature edge of the seed and prints a
summary.

Usage:
    uv run python scripts/track2_seed_filter.py
    uv run python scripts/track2_seed_filter.py --verdicts-out reports/track2/provenance/x.tsv
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402
import track2_disease_surface_forms as surface  # noqa: E402

DEFAULT_RULE = common.ROOT / "literature" / "track2_seed_disease_surface_forms_v1.yaml"
DEFAULT_EVIDENCE = common.RAW_DIR / "edge_provenance_evidence.json"
DEFAULT_VERDICTS = (common.ROOT / "reports" / "track2" / "provenance"
                    / "seed_edge_verdicts_2026_09_15.tsv")
VERDICT_HEADER = ["relation", "source_accession", "target_accession", "sentences",
                  "accepted_sentences", "verdict", "accepted_forms", "rejected_forms"]

# The `text_hl` field of a search result marks an entity as `@<m>BIOTYPE_body</m> @BIOTYPE_id
# @@@surface@@@`. The first group is the accession body and the second is the matched text.
_MARKED_SURFACE = re.compile(r"@<m>([A-Z]+_[^<]*)</m>\s*@[A-Z]+_\S+\s*@@@(.*?)@@@", re.DOTALL)

REASON_NO_ACCEPT_PATTERN = "no accept pattern matches the surface form"
REASON_DENY_FORM = "the surface form is on the deny list"
REASON_OTHER_SENSE = "the publication uses the abbreviation in another sense"
REASON_ACCEPTED = "the surface form names the syndrome"
REASON_NO_MENTION = "the search result marks no mention of the seed"
REASON_NO_EVIDENCE = "the edge has no cached sentence"
REASON_UNCORROBORATED = "the abbreviation carries no long form in the publication"


@dataclass(frozen=True)
class SurfaceRule:
    """The rule that decides whether a seed mention names the syndrome.

    Attributes:
        version: The version of the list file.
        seed_node: The graph node identifier of the seed.
        seed_accession: The PubTator3 accession of the seed.
        accept_patterns: The compiled accept patterns, in file order.
        deny_forms: The normalized surface forms that the deny list holds.
        other_sense_patterns: The compiled patterns of a competing sense.
        source: The path of the list file.
        digest: The SHA-256 of the list file.
        token_sets: The content token set of each MeSH entry term, for a version 2 rule.
        abbreviations: The abbreviations of the concept, for a version 2 rule.
        generic_words: The words that name a class of disease and not a disease.
        min_shared_tokens: The number of shared content tokens that a match needs.
        long_token_length: The length at which one shared token is distinctive on its own.
        preferred_tokens: The content tokens of the preferred term of the concept.
        require_preferred_overlap: True when a matched text must also share one content token
            with the preferred term.
    """

    version: int
    seed_node: str
    seed_accession: str
    accept_patterns: tuple[re.Pattern[str], ...]
    deny_forms: frozenset[str]
    other_sense_patterns: tuple[re.Pattern[str], ...]
    source: str
    digest: str
    token_sets: tuple[frozenset[str], ...] = ()
    abbreviations: frozenset[str] = frozenset()
    generic_words: frozenset[str] = frozenset()
    min_shared_tokens: int = 2
    long_token_length: int = 8
    preferred_tokens: frozenset[str] = frozenset()
    require_preferred_overlap: bool = False


def repo_path(path: Path) -> str:
    """Write a path relative to the repository root when it sits inside it.

    Args:
        path: The path.

    Returns:
        The relative path, or the absolute path when the file sits outside the repository.
    """
    try:
        return str(path.relative_to(common.ROOT))
    except ValueError:
        return str(path)


def normalize_form(text: str) -> str:
    """Fold a surface form to the shape that the deny list uses.

    Args:
        text: The matched text of one mention.

    Returns:
        The text with collapsed white space, in lower case.
    """
    return " ".join((text or "").split()).lower()


def load_rule(path: Path = DEFAULT_RULE) -> SurfaceRule:
    """Read the surface form list file.

    Args:
        path: The path of the list file.

    Returns:
        The rule.

    Raises:
        ValueError: When the file holds no accept pattern.
    """
    data = yaml.safe_load(path.read_text()) or {}
    accept = tuple(re.compile(entry["pattern"], re.IGNORECASE)
                   for entry in data.get("accept_patterns", []))
    sets = tuple(frozenset(entry["tokens"]) for entry in data.get("token_sets", [])
                 if entry.get("tokens"))
    if not accept and not sets:
        raise ValueError(f"{path} holds neither an accept pattern nor a token set.")
    deny = frozenset(normalize_form(entry["form"]) for entry in data.get("deny_forms", []))
    other = tuple(re.compile(entry["pattern"], re.IGNORECASE)
                  for entry in data.get("other_sense_patterns", []))
    policy = data.get("policy", {}) or {}
    return SurfaceRule(version=int(data.get("version", 0)),
                       seed_node=data.get("seed_node", ""),
                       seed_accession=data.get("seed_accession", ""),
                       accept_patterns=accept, deny_forms=deny, other_sense_patterns=other,
                       source=repo_path(path), digest=common.sha256_of(path),
                       token_sets=sets,
                       abbreviations=frozenset(entry["form"] for entry
                                               in data.get("abbreviation_patterns", [])),
                       generic_words=frozenset(policy.get("generic_words",
                                                          surface.GENERIC_WORDS)),
                       min_shared_tokens=int(policy.get("min_shared_tokens",
                                                        surface.MIN_SHARED_TOKENS)),
                       long_token_length=int(policy.get("long_token_length",
                                                        surface.LONG_TOKEN_LENGTH)),
                       preferred_tokens=frozenset(policy.get("preferred_tokens", [])),
                       require_preferred_overlap=bool(
                           policy.get("require_preferred_term_overlap", False)))


def marked_surfaces(text_hl: str) -> dict[str, list[str]]:
    """Read the surface form that the recogniser mapped onto each marked entity.

    Args:
        text_hl: The `text_hl` field of a PubTator3 search result.

    Returns:
        A mapping of accession to the surface forms of its marked mentions, in order.
    """
    found: dict[str, list[str]] = {}
    for body, surface in _MARKED_SURFACE.findall(text_hl or ""):
        found.setdefault(f"@{body}", []).append(" ".join(surface.split()))
    return found


def other_sense(context: str, rule: SurfaceRule) -> bool:
    """Say whether the text of a publication uses the abbreviation in another sense.

    Args:
        context: The title and the sentence of the item, joined.
        rule: The rule.

    Returns:
        True when one other sense pattern matches.
    """
    return any(pattern.search(context or "") for pattern in rule.other_sense_patterns)


def distinctive(shared: set[str], rule: SurfaceRule) -> bool:
    """Say whether a set of shared content tokens is large enough to name the concept.

    Args:
        shared: The tokens that the matched text and an entry term have in common.
        rule: The rule.

    Returns:
        True when the shared part reaches the token count, or when one long token carries it.
    """
    if len(shared) >= rule.min_shared_tokens:
        return True
    return len(shared) == 1 and len(next(iter(shared))) >= rule.long_token_length


def long_form_match(form: str, rule: SurfaceRule) -> bool:
    """Say whether the matched text of a mention is a name of the concept.

    A version 1 rule answers from its accept patterns. A version 2 rule answers from the content
    token set of each MeSH entry term. The matched text matches a term when one token set contains
    the other and the shared part is distinctive.

    Args:
        form: The matched text of the mention.
        rule: The rule.

    Returns:
        True when the matched text names the concept.
    """
    if rule.accept_patterns:
        return any(pattern.search(normalize_form(form)) for pattern in rule.accept_patterns)
    tokens = set(surface.content_tokens(form))
    if not tokens:
        return False
    if rule.require_preferred_overlap and rule.preferred_tokens \
            and not tokens & rule.preferred_tokens:
        return False
    for entry in rule.token_sets:
        if tokens <= entry or entry <= tokens:
            if distinctive(tokens & entry, rule):
                return True
    return False


def abbreviation_match(form: str, rule: SurfaceRule) -> bool:
    """Say whether the matched text of a mention is an abbreviation of the concept.

    Args:
        form: The matched text of the mention.
        rule: The rule.

    Returns:
        True when one token of the matched text is an abbreviation of the concept, with or without
        a trailing plural letter or type number.
    """
    for token in surface.normalize(form).split():
        stripped = token.rstrip("0123456789")
        for candidate in (token, stripped, stripped.rstrip("s")):
            if candidate in rule.abbreviations:
                return True
    return False


def context_has_long_form(context: str, rule: SurfaceRule) -> bool:
    """Say whether the publication of a mention writes a name of the concept.

    Args:
        context: The title and the sentence of the item, joined.
        rule: The rule.

    Returns:
        True when every content token of one entry term appears in the text.
    """
    tokens = set(surface.content_tokens(context))
    if not tokens:
        return False
    for entry in rule.token_sets:
        if entry <= tokens and distinctive(entry, rule):
            return True
    return False


def abbreviation_corroborated(form: str, context: str, rule: SurfaceRule) -> bool:
    """Say whether an abbreviation mention carries enough support to name the concept.

    An abbreviation is ambiguous by construction. Three signals corroborate it. The matched text
    itself carries a generic disease word of the concept, for example "MVA syndrome". The
    abbreviation carries a type number, for example "MVA1". The publication writes a name of the
    concept in its title or in its sentence.

    Args:
        form: The matched text of the mention.
        context: The title and the sentence of the item, joined.
        rule: The rule.

    Returns:
        True when one signal holds.
    """
    tokens = surface.normalize(form).split()
    if any(token in rule.generic_words for token in tokens):
        return True
    if any(token.rstrip("0123456789") in rule.abbreviations and token[-1].isdigit()
           for token in tokens if token):
        return True
    return context_has_long_form(context, rule)


def mention_verdict(form: str, context: str, rule: SurfaceRule) -> tuple[bool, str]:
    """Decide whether one mention of the seed names the disease.

    Args:
        form: The matched text of the mention.
        context: The title and the sentence of the item, joined.
        rule: The rule.

    Returns:
        The verdict and the reason for it.
    """
    normalized = normalize_form(form)
    names = long_form_match(form, rule)
    abbreviated = abbreviation_match(form, rule) if not names else False
    if not names and not abbreviated:
        return False, REASON_NO_ACCEPT_PATTERN
    if normalized in rule.deny_forms:
        return False, REASON_DENY_FORM
    if other_sense(context, rule):
        return False, REASON_OTHER_SENSE
    if abbreviated and not abbreviation_corroborated(form, context, rule):
        return False, REASON_UNCORROBORATED
    return True, REASON_ACCEPTED


def item_verdict(item: dict, rule: SurfaceRule) -> dict:
    """Decide whether one search result names the syndrome.

    Args:
        item: One item of a cached search result. It carries `text_hl`, `title` and `pmid`.
        rule: The rule.

    Returns:
        The PMID, the surface forms, the verdict and the reason.
    """
    forms = marked_surfaces(item.get("text_hl", "")).get(rule.seed_accession, [])
    sentence = common.plain_sentence(item.get("text_hl", ""))
    context = " ".join([item.get("title", "") or "", sentence])
    if not forms:
        return {"pmid": item.get("pmid", ""), "forms": [], "accepted": False,
                "accepted_forms": [], "rejected_forms": [], "reason": REASON_NO_MENTION,
                "sentence": sentence, "title": item.get("title", "") or ""}
    accepted_forms, rejected_forms, reasons = [], [], []
    for form in forms:
        keep, reason = mention_verdict(form, context, rule)
        reasons.append(reason)
        (accepted_forms if keep else rejected_forms).append(form)
    return {"pmid": item.get("pmid", ""), "forms": forms, "accepted": bool(accepted_forms),
            "accepted_forms": accepted_forms, "rejected_forms": rejected_forms,
            "reason": REASON_ACCEPTED if accepted_forms else reasons[0],
            "sentence": sentence, "title": item.get("title", "") or ""}


def split_pair_key(key: str) -> tuple[str, str, str]:
    """Split a key of the provenance evidence cache.

    Args:
        key: A key of the form `relation|accession_1|accession_2|max_pages`.

    Returns:
        The relation, the first accession and the second accession.
    """
    parts = key.split("|")
    return parts[0], parts[1], parts[2]


def pair_verdicts(pairs: dict, rule: SurfaceRule) -> dict[tuple[str, str, str], dict]:
    """Decide every cached entity pair that touches the seed.

    Args:
        pairs: The `pairs` block of `edge_provenance_evidence.json`.
        rule: The rule.

    Returns:
        One record per pair, keyed by the relation and the two accessions in cache order.
    """
    out: dict[tuple[str, str, str], dict] = {}
    for key, pair in pairs.items():
        relation, acc_1, acc_2 = split_pair_key(key)
        if rule.seed_accession not in (acc_1, acc_2):
            continue
        items = [item_verdict(item, rule) for item in pair.get("items", [])]
        accepted = [item for item in items if item["accepted"]]
        record = {
            "relation": relation, "source_accession": acc_1, "target_accession": acc_2,
            "sentences": len(items), "accepted_sentences": len(accepted),
            "verdict": "keep" if accepted else "drop",
            "accepted_forms": sorted({form for item in accepted for form in item["accepted_forms"]}),
            "rejected_forms": sorted({form for item in items for form in item["rejected_forms"]}),
            "items": items,
        }
        out[(relation, acc_1, acc_2)] = record
    return out


def node_key(relation: str, acc_1: str, acc_2: str, entities: dict) -> tuple[str, str, str] | None:
    """Map a cached pair key to the node key of the graph edge.

    Args:
        relation: The relation name.
        acc_1: The accession of the first entity.
        acc_2: The accession of the second entity.
        entities: The resolved entity cache, keyed by accession.

    Returns:
        The source identifier, the target identifier and the relation, or None when either
        accession has no node identifier.
    """
    first = node_identifier(acc_1, entities)
    second = node_identifier(acc_2, entities)
    if not first or not second:
        return None
    return first, second, relation


def node_identifier(accession: str, entities: dict) -> str:
    """Return the graph node identifier of one accession.

    Args:
        accession: The PubTator3 accession.
        entities: The resolved entity cache, keyed by accession.

    Returns:
        The node identifier, or an empty string when the entity does not resolve.
    """
    if accession.startswith("@GENE_"):
        return f"GENE:{accession.split('_', 1)[1]}"
    entry = entities.get(accession, {})
    return common.node_id_for_accession(accession, entry.get("identifier", ""),
                                        entry.get("database", ""))


def dropped_edge_keys(verdicts: dict[tuple[str, str, str], dict],
                      entities: dict) -> set[tuple[str, str, str]]:
    """List the graph edge keys that the rule removes.

    An edge key is unordered, because the graph stores one direction and the cache may store the
    other. The set holds both orders so that a caller can test an edge row directly.

    Args:
        verdicts: The records of `pair_verdicts`.
        entities: The resolved entity cache, keyed by accession.

    Returns:
        The source, target and relation of every removed edge, in both orders.
    """
    keys: set[tuple[str, str, str]] = set()
    for (relation, acc_1, acc_2), record in verdicts.items():
        if record["verdict"] == "keep":
            continue
        mapped = node_key(relation, acc_1, acc_2, entities)
        if not mapped:
            continue
        keys.add(mapped)
        keys.add((mapped[1], mapped[0], mapped[2]))
    return keys


def unlisted_forms(verdicts: dict[tuple[str, str, str], dict], listed: set[str]) -> list[str]:
    """List the surface forms that the cache holds and the list file does not record.

    Args:
        verdicts: The records of `pair_verdicts`.
        listed: The normalized forms of the `audit_forms` block.

    Returns:
        The unlisted forms, sorted.
    """
    seen = {normalize_form(form) for record in verdicts.values()
            for item in record["items"] for form in item["forms"]}
    return sorted(seen - listed)


def summary(verdicts: dict[tuple[str, str, str], dict], rule: SurfaceRule,
            listed: set[str]) -> dict:
    """Summarize the decision over every cached seed pair.

    Args:
        verdicts: The records of `pair_verdicts`.
        rule: The rule.
        listed: The normalized forms of the `audit_forms` block.

    Returns:
        The counts of pairs, sentences and forms, for a receipt.
    """
    items = [item for record in verdicts.values() for item in record["items"]]
    return {
        "rule_file": rule.source,
        "rule_sha256_prefix": rule.digest[:16],
        "rule_version": rule.version,
        "seed_node": rule.seed_node,
        "seed_accession": rule.seed_accession,
        "pairs": len(verdicts),
        "pairs_kept": sum(1 for r in verdicts.values() if r["verdict"] == "keep"),
        "pairs_dropped": sum(1 for r in verdicts.values() if r["verdict"] == "drop"),
        "sentences": len(items),
        "sentences_accepted": sum(1 for item in items if item["accepted"]),
        "rejection_reasons": dict(Counter(item["reason"] for item in items
                                          if not item["accepted"]).most_common()),
        "accepted_forms": dict(Counter(normalize_form(form) for item in items
                                       for form in item["accepted_forms"]).most_common()),
        "rejected_forms": dict(Counter(normalize_form(form) for item in items
                                       for form in item["rejected_forms"]).most_common()),
        "forms_not_in_the_list_file": unlisted_forms(verdicts, listed),
    }


def listed_forms(path: Path = DEFAULT_RULE) -> set[str]:
    """Read the normalized surface forms of the `audit_forms` block.

    Args:
        path: The path of the list file.

    Returns:
        The normalized forms.
    """
    data = yaml.safe_load(path.read_text()) or {}
    return {normalize_form(entry["form"]) for entry in data.get("audit_forms", [])}


def main() -> int:
    """Apply the rule to the cached seed evidence and write the verdict table.

    Returns:
        0 on success, 1 when the evidence cache is missing.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--rule", type=Path, default=DEFAULT_RULE,
                        help="The surface form list file.")
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE,
                        help="The cached PubTator3 search results of the seed edges.")
    parser.add_argument("--entities", type=Path, default=common.RAW_DIR / "pubtator_entities.json",
                        help="The resolved entity cache.")
    parser.add_argument("--verdicts-out", type=Path, default=DEFAULT_VERDICTS,
                        help="The verdict table to write.")
    args = parser.parse_args()

    cache = common.read_json(args.evidence) or {}
    if not cache.get("pairs"):
        print(f"FAIL: no cached evidence at {args.evidence}.")
        return 1
    rule = load_rule(args.rule)
    verdicts = pair_verdicts(cache["pairs"], rule)
    entities = (common.read_json(args.entities) or {}).get("entities", {})
    rows = [{"relation": record["relation"], "source_accession": record["source_accession"],
             "target_accession": record["target_accession"], "sentences": record["sentences"],
             "accepted_sentences": record["accepted_sentences"], "verdict": record["verdict"],
             "accepted_forms": " | ".join(record["accepted_forms"]),
             "rejected_forms": " | ".join(record["rejected_forms"])}
            for record in sorted(verdicts.values(),
                                 key=lambda r: (r["verdict"], -r["accepted_sentences"],
                                                r["relation"], r["source_accession"],
                                                r["target_accession"]))]
    path = common.write_tsv(args.verdicts_out, VERDICT_HEADER, rows)
    report = summary(verdicts, rule, listed_forms(args.rule))
    print(f"rule {report['rule_file']} version {report['rule_version']}")
    print(f"pairs {report['pairs']}, kept {report['pairs_kept']}, "
          f"dropped {report['pairs_dropped']}")
    print(f"sentences {report['sentences']}, accepted {report['sentences_accepted']}")
    print(f"rejection reasons: {report['rejection_reasons']}")
    print(f"removed edge keys: {len(dropped_edge_keys(verdicts, entities)) // 2}")
    if report["forms_not_in_the_list_file"]:
        print(f"WARNING: forms not in the list file: {report['forms_not_in_the_list_file']}")
    print(f"verdict table: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
