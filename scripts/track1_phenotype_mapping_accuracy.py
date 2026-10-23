#!/usr/bin/env python3
"""Measure the accuracy of the clinical document to HPO term mapping. Data class A.

Defect R2-8 of the second Track 1 jury round states that the accuracy of the clinical document to
HPO mapping is not measured. The mapping is `pipelines/track1_b0/bin/extract_case_phenotype.py`. It
builds a lexicon from the label and the exact and narrow synonyms of every term under `HP:0000118`,
it matches a label on word boundaries, and it removes the ancestors of a matched term. It handles no
negation and no family history.

This program imports that module and calls its own functions. It changes no line of the mapping. It
scores the mapping on four arms.

    gsc_plus            The GSC+ gold standard corpus for HPO concept recognition, 228 PubMed
                        abstracts with curated HPO spans. Precision, recall and F1.
    orphanet            The Orphanet clinical definition text of a rare disease against the curated
                        Orphanet HPO annotation set of the same disease. Precision and text
                        coverage.
    synthetic           An invented corpus of 20 short clinical notes. Each sentence carries a
                        design tag, so the count of each error mode is exact.
    label_control       A positive control. One document that holds the label of one term under
                        `HP:0000118` must return that term.

Every input of this run is public. No patient file and no term of the child enters any step. The
synthetic notes describe no real person.

The gate conditions are in `docs/benchmarks/t1_phenotype_mapping_accuracy_v1.md`, recorded before
the result.

Usage:
    uv run python scripts/track1_phenotype_mapping_accuracy.py \\
        --hpo-dir tracking/cache/hpo_v2026_09_01 \\
        --cache-dir tracking/cache/t1_hpo_mapping \\
        --receipt tracking/receipts/t1_phenotype_mapping_accuracy_v1.json
"""

from __future__ import annotations

import argparse
import dataclasses
import datetime as dt
import hashlib
import json
import math
import random
import re
import subprocess
import sys
import urllib.request
import xml.etree.ElementTree as ElementTree
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "pipelines" / "track1_b0" / "bin"
sys.path.insert(0, str(BIN))

import extract_case_phenotype as extract  # noqa: E402

GSC_PLUS_URL = "https://raw.githubusercontent.com/lasigeBioTM/IHP/master/GSC%2B.rar"
GSC_PLUS_COMMIT = "05872b11d5f52538d6d5de117dbcf0e9065ba9e9"
ORPHANET_PRODUCT4_URL = "https://www.orphadata.com/data/xml/en_product4.xml"
ORPHANET_PRODUCT1_URL = "https://www.orphadata.com/data/xml/en_product1.xml"

LABEL_CONTROL_SAMPLE = 2000
LABEL_CONTROL_SEED = 20260917

ANNOTATION_ROW = re.compile(r"^\[(\d+)::(\d+)\]\t(HP[_:]\d+)\s*\|\s*(.*)$")
SENTENCE_SPLIT = re.compile(r"(?<=[.;:!?])\s+|\n+")
TAG = re.compile(r"<[^>]+>")

# A negation cue applies only inside a window of tokens before the matched label. "absence of" is
# not a cue, because the HPO itself uses that wording in a label of a present finding.
NEGATION_CUES = (
    "no", "not", "never", "without", "denies", "denied", "negative", "neither", "nor",
    "lacks", "unremarkable", "excludes", "excluded", "ruled", "rules", "absent",
)
NEGATION_WINDOW = 6

# A family history cue and a differential diagnosis cue apply to the whole sentence, because the
# scope of both constructs is the sentence.
FAMILY_HISTORY_CUES = (
    "family history", "mother", "father", "maternal", "paternal", "sibling", "siblings",
    "brother", "sister", "grandmother", "grandfather", "cousin", "aunt", "uncle", "pedigree",
    "kindred", "relative", "relatives", "consanguineous",
)
DIFFERENTIAL_CUES = (
    "differential diagnosis", "differential", "differentiated from", "distinguished from",
    "distinguish", "mimic", "mimics", "resembles", "suggestive of", "considered", "consider",
    "in contrast to", "as opposed to", "unlike", "versus",
)


@dataclasses.dataclass(frozen=True)
class SyntheticSentence:
    """One sentence of the synthetic corpus.

    Attributes:
        text: The sentence.
        construct: The design tag. One of `assert`, `negated`, `family_history`, `differential`
            and `paraphrase`.
        term: The HPO identifier that the sentence names.
    """

    text: str
    construct: str
    term: str


@dataclasses.dataclass(frozen=True)
class SyntheticDocument:
    """One note of the synthetic corpus.

    Attributes:
        doc_id: The note identifier.
        sentences: The sentences of the note.
    """

    doc_id: str
    sentences: tuple[SyntheticSentence, ...]


# The 20 notes below are invented for this measurement. They describe no real person. Each HPO
# identifier occurs in exactly one sentence of the whole corpus, so a false positive attributes to
# one construct with no ambiguity. `build_synthetic_corpus` checks that property.
SYNTHETIC_CORPUS = (
    SyntheticDocument("s01", (
        SyntheticSentence("The child shows microcephaly on serial measurement.",
                          "assert", "HP:0000252"),
        SyntheticSentence("There is no splenomegaly.", "negated", "HP:0001744"),
    )),
    SyntheticDocument("s02", (
        SyntheticSentence("Short stature was recorded at each visit.", "assert", "HP:0004322"),
        SyntheticSentence("The mother of the child has hearing impairment.",
                          "family_history", "HP:0000365"),
    )),
    SyntheticDocument("s03", (
        SyntheticSentence("Seizure activity began in the second year.", "assert", "HP:0001250"),
        SyntheticSentence("The differential diagnosis included atrial septal defect.",
                          "differential", "HP:0001631"),
    )),
    SyntheticDocument("s04", (
        SyntheticSentence("Intellectual disability is documented in the school report.",
                          "assert", "HP:0001249"),
        SyntheticSentence("The examination found no ptosis.", "negated", "HP:0000508"),
    )),
    SyntheticDocument("s05", (
        SyntheticSentence("The child cannot walk more than ten steps without holding a hand.",
                          "paraphrase", "HP:0001288"),
        SyntheticSentence("Brachydactyly of both hands is present.", "assert", "HP:0001156"),
    )),
    SyntheticDocument("s06", (
        SyntheticSentence("A brother has adrenal insufficiency.",
                          "family_history", "HP:0000846"),
        SyntheticSentence("Unsteady gait is present from the age of three years.",
                          "assert", "HP:0002317"),
    )),
    SyntheticDocument("s07", (
        SyntheticSentence("Hepatomegaly is not present on palpation.", "negated", "HP:0002240"),
        SyntheticSentence("Anemia was found on the blood count.", "assert", "HP:0001903"),
    )),
    SyntheticDocument("s08", (
        SyntheticSentence("The differential diagnosis listed recurrent respiratory infections.",
                          "differential", "HP:0002205"),
        SyntheticSentence("Joint hypermobility is present at the knees and the elbows.",
                          "assert", "HP:0001382"),
    )),
    SyntheticDocument("s09", (
        SyntheticSentence("Dysphagia started at the age of four years.", "assert", "HP:0002015"),
        SyntheticSentence("The father has thrombocytopenia.", "family_history", "HP:0001873"),
    )),
    SyntheticDocument("s10", (
        SyntheticSentence("No cataract is seen on slit lamp examination.",
                          "negated", "HP:0000518"),
        SyntheticSentence("Nystagmus is present on lateral gaze.", "assert", "HP:0000639"),
    )),
    SyntheticDocument("s11", (
        SyntheticSentence("The child speaks fewer than ten words at the age of three years.",
                          "paraphrase", "HP:0000750"),
        SyntheticSentence("Hypotonia is present in the trunk.", "assert", "HP:0001252"),
    )),
    SyntheticDocument("s12", (
        SyntheticSentence("The team considered hypertrophic cardiomyopathy.",
                          "differential", "HP:0001639"),
        SyntheticSentence("Scoliosis of 25 degrees is recorded.", "assert", "HP:0002650"),
    )),
    SyntheticDocument("s13", (
        SyntheticSentence("A maternal cousin has optic atrophy.",
                          "family_history", "HP:0000648"),
        SyntheticSentence("Failure to thrive was recorded in the first year.",
                          "assert", "HP:0001508"),
    )),
    SyntheticDocument("s14", (
        SyntheticSentence("There is no polydactyly of the hands or the feet.",
                          "negated", "HP:0010442"),
        SyntheticSentence("Cleft palate was repaired at the age of one year.",
                          "assert", "HP:0000175"),
    )),
    SyntheticDocument("s15", (
        SyntheticSentence("Macrocephaly is present and the centile is above 97.",
                          "assert", "HP:0000256"),
        SyntheticSentence("The differential diagnosis included hydrocephalus.",
                          "differential", "HP:0000238"),
    )),
    SyntheticDocument("s16", (
        SyntheticSentence("The child cannot lift a cup with one hand and the grip gives way.",
                          "paraphrase", "HP:0001324"),
        SyntheticSentence("Spasticity is present in both legs.", "assert", "HP:0001257"),
    )),
    SyntheticDocument("s17", (
        SyntheticSentence("Strabismus was not observed at any visit.", "negated", "HP:0000486"),
        SyntheticSentence("Hyperreflexia is present at the knee.", "assert", "HP:0001347"),
    )),
    SyntheticDocument("s18", (
        SyntheticSentence("The sister of the child has a renal cyst.",
                          "family_history", "HP:0000107"),
        SyntheticSentence("Bruising susceptibility is reported after a minor knock.",
                          "assert", "HP:0000978"),
    )),
    SyntheticDocument("s19", (
        SyntheticSentence("Ataxia is present on the finger to nose test.",
                          "assert", "HP:0001251"),
        SyntheticSentence("The differential diagnosis listed retinal dystrophy.",
                          "differential", "HP:0000556"),
    )),
    SyntheticDocument("s20", (
        SyntheticSentence("The child was admitted four times in one year for a chest illness "
                          "that needed antibiotics.", "paraphrase", "HP:0002719"),
        SyntheticSentence("Pes planus is present on both feet.", "assert", "HP:0001763"),
    )),
)

PRESENT_CONSTRUCTS = frozenset({"assert", "paraphrase"})
ABSENT_CONSTRUCTS = frozenset({"negated", "family_history", "differential"})


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


def fetch_to_cache(url: str, path: Path) -> dict[str, str]:
    """Download a file once and record its provenance.

    Args:
        url: The download address.
        path: The cache file.

    Returns:
        The url, the cache path, the SHA-256 and the download date of the cache file.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fetched = False
    if not path.exists():
        with urllib.request.urlopen(url, timeout=600) as response, path.open("wb") as handle:
            while block := response.read(1 << 20):
                handle.write(block)
        fetched = True
    stamp = dt.datetime.fromtimestamp(path.stat().st_mtime, dt.timezone.utc)
    return {
        "url": url,
        "path": str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
        "sha256": sha256_of(path),
        "cached_utc": stamp.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "downloaded_in_this_run": fetched,
    }


def wilson_interval(successes: int, total: int) -> tuple[float, float]:
    """Compute the Wilson score interval of a proportion at 95 percent.

    Args:
        successes: The count of successes.
        total: The count of trials.

    Returns:
        The lower bound and the upper bound. Both are zero when the total is zero.
    """
    if total == 0:
        return 0.0, 0.0
    z = 1.959963984540054
    proportion = successes / total
    centre = proportion + z * z / (2 * total)
    spread = z * math.sqrt(proportion * (1 - proportion) / total + z * z / (4 * total * total))
    scale = 1 + z * z / total
    return max(0.0, (centre - spread) / scale), min(1.0, (centre + spread) / scale)


def proportion_record(successes: int, total: int) -> dict[str, float | int]:
    """Describe a proportion with its denominator and its interval.

    Args:
        successes: The numerator.
        total: The denominator.

    Returns:
        The numerator, the denominator, the value and the Wilson bounds.
    """
    low, high = wilson_interval(successes, total)
    return {
        "numerator": successes,
        "denominator": total,
        "value": round(successes / total, 4) if total else 0.0,
        "ci95_low": round(low, 4),
        "ci95_high": round(high, 4),
    }


def parse_obo_identifier_history(path: Path) -> tuple[dict[str, str], dict[str, str], set[str]]:
    """Read the secondary identifiers, the replacements and the obsolete terms of the release.

    Args:
        path: The `hp.obo` file.

    Returns:
        The secondary identifier mapped to its primary identifier, the obsolete identifier mapped to
        its replacement, and the set of obsolete identifiers.
    """
    alternates: dict[str, str] = {}
    replacements: dict[str, str] = {}
    obsolete: set[str] = set()
    term_id = ""
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if line == "[Term]":
                term_id = ""
            elif line.startswith("id: "):
                term_id = line[4:].strip()
            elif line.startswith("alt_id: "):
                alternates[line[8:].strip()] = term_id
            elif line.startswith("replaced_by: "):
                replacements[term_id] = line[13:].strip()
            elif line == "is_obsolete: true":
                obsolete.add(term_id)
    return alternates, replacements, obsolete


def resolve_identifier(term: str, live: frozenset[str], alternates: dict[str, str],
                       replacements: dict[str, str]) -> str:
    """Map an identifier of an older release to a live identifier of the current release.

    Args:
        term: The identifier as the ground truth writes it.
        live: The live identifiers of the current release.
        alternates: The secondary identifier mapped to its primary identifier.
        replacements: The obsolete identifier mapped to its replacement.

    Returns:
        The live identifier, or the empty string when no live identifier exists.
    """
    if term in live:
        return term
    if term in alternates and alternates[term] in live:
        return alternates[term]
    if term in replacements:
        target = replacements[term]
        if target in live:
            return target
        if target in alternates and alternates[target] in live:
            return alternates[target]
    return ""


def raw_hits(text: str, lexicon: dict[str, str]) -> dict[str, str]:
    """Find every term whose label occurs in the text, before the ancestor removal.

    This mirrors the first loop of `extract.match_terms`. `verify_mirror` checks the mirror against
    the shipped function on every document of every arm.

    Args:
        text: The normalized document text.
        lexicon: The normalized label mapped to its term.

    Returns:
        Each matched term mapped to the longest label that matched it.
    """
    hits: dict[str, str] = {}
    for label, term in lexicon.items():
        if f" {label} " in text and (term not in hits or len(label) > len(hits[term])):
            hits[term] = label
    return hits


def reduce_to_specific(hits: dict[str, str],
                       closure: dict[str, frozenset[str]]) -> list[tuple[str, str]]:
    """Drop every matched term that is an ancestor of another matched term.

    Args:
        hits: Each matched term mapped to the label that matched it.
        closure: The reflexive ancestor set of every term.

    Returns:
        The term and the label that matched it, sorted by term identifier.
    """
    ancestors = {ancestor for term in hits for ancestor in closure[term] if ancestor != term}
    return sorted((term, label) for term, label in hits.items() if term not in ancestors)


def verify_mirror(text: str, lexicon: dict[str, str],
                  closure: dict[str, frozenset[str]]) -> bool:
    """Check that the local mirror returns what the shipped mapping returns.

    Args:
        text: The normalized document text.
        lexicon: The normalized label mapped to its term.
        closure: The reflexive ancestor set of every term.

    Returns:
        True when the two results are equal.
    """
    return reduce_to_specific(raw_hits(text, lexicon), closure) == extract.match_terms(
        text, lexicon, closure)


def most_specific(terms: set[str], closure: dict[str, frozenset[str]]) -> frozenset[str]:
    """Keep the terms of a set that are not an ancestor of another term of the same set.

    Args:
        terms: The identifiers.
        closure: The reflexive ancestor set of every term.

    Returns:
        The most specific identifiers of the set.
    """
    ancestors = {ancestor for term in terms
                 for ancestor in closure.get(term, frozenset()) if ancestor != term}
    return frozenset(term for term in terms if term not in ancestors)


def is_strict_ancestor(candidate: str, target: str,
                       closure: dict[str, frozenset[str]]) -> bool:
    """Test whether one term is a proper ancestor of another term.

    Args:
        candidate: The possible ancestor.
        target: The possible descendant.
        closure: The reflexive ancestor set of every term.

    Returns:
        True when the candidate is an ancestor of the target and not the target.
    """
    return candidate != target and candidate in closure.get(target, frozenset())


def hierarchy_relation(predicted: str, gold: frozenset[str],
                       closure: dict[str, frozenset[str]]) -> str:
    """Name the relation between a predicted term and a gold term set.

    Args:
        predicted: The predicted identifier.
        gold: The gold identifiers.
        closure: The reflexive ancestor set of every term.

    Returns:
        `exact`, `over_general`, `under_general` or `unrelated`.
    """
    if predicted in gold:
        return "exact"
    if any(is_strict_ancestor(predicted, term, closure) for term in gold):
        return "over_general"
    if any(is_strict_ancestor(term, predicted, closure) for term in gold):
        return "under_general"
    return "unrelated"


def split_sentences(text: str) -> list[str]:
    """Split a document into sentences on a full stop, a semicolon, a colon or a line break.

    Args:
        text: The raw document text.

    Returns:
        The sentences with no empty entry.
    """
    return [part.strip() for part in SENTENCE_SPLIT.split(text) if part.strip()]


def label_is_negated(sentence: str, label: str) -> bool:
    """Test whether a negation cue stands in the window before a label in a sentence.

    Args:
        sentence: The raw sentence.
        label: The normalized label that matched.

    Returns:
        True when a cue of `NEGATION_CUES` stands within `NEGATION_WINDOW` tokens before the label.
    """
    tokens = extract.normalize(sentence).split()
    target = label.split()
    for index in range(len(tokens) - len(target) + 1):
        if tokens[index:index + len(target)] != target:
            continue
        window = tokens[max(0, index - NEGATION_WINDOW):index]
        if any(token in NEGATION_CUES for token in window):
            return True
    return False


def classify_context(text: str, label: str) -> str:
    """Name the construct that carries a label in a document.

    A label that stands in one plain sentence is reported as `assertion`, even when it also stands
    in a negated sentence. The rule does not blame a construct when a plain statement of the same
    finding exists.

    Args:
        text: The raw document text.
        label: The normalized label that matched.

    Returns:
        `negation`, `family_history`, `differential_diagnosis`, `assertion` or `no_sentence_found`.
    """
    padded = f" {label} "
    carriers = [sentence for sentence in split_sentences(text)
                if padded in extract.normalize(sentence)]
    if not carriers:
        return "no_sentence_found"
    verdicts = []
    for sentence in carriers:
        lowered = " " + " ".join(extract.WORD.findall(sentence.lower())) + " "
        if label_is_negated(sentence, label):
            verdicts.append("negation")
        elif any(f" {cue} " in lowered for cue in FAMILY_HISTORY_CUES):
            verdicts.append("family_history")
        elif any(f" {cue} " in lowered for cue in DIFFERENTIAL_CUES):
            verdicts.append("differential_diagnosis")
        else:
            verdicts.append("assertion")
    if "assertion" in verdicts:
        return "assertion"
    return verdicts[0]


def classify_false_negative(term: str, mentions: list[str], text: str,
                            hits: dict[str, str], lexicon: dict[str, str],
                            labels: dict[str, list[str]]) -> str:
    """Name the reason why a gold term is absent from the output.

    Args:
        term: The gold identifier.
        mentions: The gold mention strings of the term in the document.
        text: The normalized document text.
        hits: The matched terms of the document before the ancestor removal.
        lexicon: The normalized label mapped to its term.
        labels: The labels of every term.

    Returns:
        One of `removed_as_an_ancestor`, `wording_absent_from_the_lexicon`,
        `label_maps_to_another_term`, `label_below_length_or_on_the_stop_list` and `other`.
    """
    if term in hits:
        return "removed_as_an_ancestor"
    own = [extract.normalize(name).strip() for name in labels.get(term, [])]
    kept = [key for key in own
            if len(key) >= extract.MIN_LABEL_CHARACTERS and key not in extract.STOP_LABELS]
    if own and not kept:
        return "label_below_length_or_on_the_stop_list"
    if any(f" {key} " in text for key in kept):
        return "other"
    normalized_mentions = [extract.normalize(mention).strip() for mention in mentions]
    in_lexicon = [key for key in normalized_mentions if key in lexicon]
    if not in_lexicon:
        return "wording_absent_from_the_lexicon"
    if all(lexicon[key] != term for key in in_lexicon):
        return "label_maps_to_another_term"
    return "other"


def read_gsc_plus(root: Path) -> list[dict[str, object]]:
    """Read the GSC+ corpus text and annotations.

    Args:
        root: The directory that holds `Text` and `Annotations`.

    Returns:
        One record per document with the identifier, the text, the mentions of each term and the
        count of annotation rows.
    """
    documents = []
    for text_path in sorted((root / "Text").iterdir()):
        annotation_path = root / "Annotations" / text_path.name
        mentions: dict[str, list[str]] = {}
        rows = 0
        if annotation_path.exists():
            for raw in annotation_path.read_text(errors="replace").splitlines():
                line = raw.rstrip("\r").strip()
                if not line:
                    continue
                match = ANNOTATION_ROW.match(line)
                if not match:
                    continue
                rows += 1
                term = match.group(3).replace("_", ":")
                mentions.setdefault(term, []).append(match.group(4).strip())
        documents.append({
            "doc_id": text_path.name,
            "text": text_path.read_text(errors="replace"),
            "mentions": mentions,
            "annotation_rows": rows,
        })
    return documents


def parse_orphanet_annotations(path: Path) -> dict[str, set[str]]:
    """Read the curated HPO annotation set of every Orphanet disorder.

    Args:
        path: The `en_product4.xml` file.

    Returns:
        The Orpha code mapped to its annotated HPO identifiers.
    """
    annotations: dict[str, set[str]] = {}
    for _, element in ElementTree.iterparse(str(path), events=("end",)):
        if element.tag != "Disorder":
            continue
        code = element.findtext("OrphaCode", default="").strip()
        terms = {node.text.strip() for node in element.iter("HPOId") if node.text}
        if code and terms:
            annotations[code] = terms
        element.clear()
    return annotations


def parse_orphanet_definitions(path: Path) -> dict[str, dict[str, str]]:
    """Read the clinical definition text of every Orphanet disorder.

    Args:
        path: The `en_product1.xml` file.

    Returns:
        The Orpha code mapped to the disorder name and the definition text.
    """
    definitions: dict[str, dict[str, str]] = {}
    for _, element in ElementTree.iterparse(str(path), events=("end",)):
        if element.tag != "Disorder":
            continue
        code = element.findtext("OrphaCode", default="").strip()
        name = element.findtext("Name", default="").strip()
        text = ""
        for section in element.iter("TextSection"):
            kind = section.find("TextSectionType")
            label = kind.findtext("Name", default="") if kind is not None else ""
            if label.strip() == "Definition":
                text = TAG.sub(" ", section.findtext("Contents", default="")).strip()
                break
        if code and text:
            definitions[code] = {"name": name, "definition": text}
        element.clear()
    return definitions


def build_synthetic_corpus() -> tuple[SyntheticDocument, ...]:
    """Return the synthetic corpus and check that each term occurs in one sentence only.

    Returns:
        The synthetic documents.

    Raises:
        ValueError: When one HPO identifier occurs in more than one sentence, because a false
            positive would then attribute to more than one construct.
    """
    seen = Counter(sentence.term for document in SYNTHETIC_CORPUS
                   for sentence in document.sentences)
    repeated = sorted(term for term, count in seen.items() if count > 1)
    if repeated:
        raise ValueError(f"These identifiers occur in more than one sentence: {repeated}")
    return SYNTHETIC_CORPUS


def synthetic_text(document: SyntheticDocument) -> str:
    """Join the sentences of a synthetic note.

    Args:
        document: The note.

    Returns:
        The note text, one sentence per line.
    """
    return "\n".join(sentence.text for sentence in document.sentences)


def score_synthetic(documents: tuple[SyntheticDocument, ...], lexicon: dict[str, str],
                    closure: dict[str, frozenset[str]]) -> dict[str, object]:
    """Score the mapping on the synthetic corpus and count each designed error mode.

    Args:
        documents: The synthetic notes.
        lexicon: The normalized label mapped to its term.
        closure: The reflexive ancestor set of every term.

    Returns:
        The precision, the recall, the F1, the count of each construct and its denominator, and the
        terms that the lexicon matched by accident.
    """
    true_positive = 0
    false_negative_by_construct: Counter[str] = Counter()
    false_positive_by_construct: Counter[str] = Counter()
    denominators: Counter[str] = Counter()
    unintended: list[dict[str, str]] = []
    mirror_agrees = 0
    for document in documents:
        text = synthetic_text(document)
        normalized = extract.normalize(text)
        hits = raw_hits(normalized, lexicon)
        predicted = {term for term, _ in reduce_to_specific(hits, closure)}
        mirror_agrees += int(verify_mirror(normalized, lexicon, closure))
        by_term = {sentence.term: sentence.construct for sentence in document.sentences}
        for sentence in document.sentences:
            denominators[sentence.construct] += 1
            if sentence.construct in PRESENT_CONSTRUCTS:
                if sentence.term in predicted:
                    true_positive += 1
                else:
                    false_negative_by_construct[sentence.construct] += 1
            elif sentence.term in predicted:
                false_positive_by_construct[sentence.construct] += 1
        for term in sorted(predicted):
            if term not in by_term:
                unintended.append({"doc_id": document.doc_id, "hpo_id": term,
                                   "matched_text": hits[term]})
    designed_false_positive = sum(false_positive_by_construct.values())
    gold_total = sum(denominators[name] for name in sorted(PRESENT_CONSTRUCTS))
    predicted_total = true_positive + designed_false_positive + len(unintended)
    precision = proportion_record(true_positive, predicted_total)
    recall = proportion_record(true_positive, gold_total)
    return {
        "documents": len(documents),
        "sentences": sum(denominators.values()),
        "mirror_agreements": mirror_agrees,
        "gold_terms": gold_total,
        "predicted_terms": predicted_total,
        "true_positive": true_positive,
        "precision": precision,
        "recall": recall,
        "f1": round(f1_from_counts(true_positive, predicted_total, gold_total), 4),
        "false_positive_by_construct": dict(sorted(false_positive_by_construct.items())),
        "false_negative_by_construct": dict(sorted(false_negative_by_construct.items())),
        "construct_denominators": dict(sorted(denominators.items())),
        "negated_false_positive_rate": proportion_record(
            false_positive_by_construct["negated"], denominators["negated"]),
        "family_history_false_positive_rate": proportion_record(
            false_positive_by_construct["family_history"], denominators["family_history"]),
        "differential_false_positive_rate": proportion_record(
            false_positive_by_construct["differential"], denominators["differential"]),
        "paraphrase_false_negative_rate": proportion_record(
            false_negative_by_construct["paraphrase"], denominators["paraphrase"]),
        "unintended_lexicon_matches": unintended,
    }


def f1_from_counts(true_positive: int, predicted: int, gold: int) -> float:
    """Compute the F1 from the integer counts and not from the rounded rates.

    The third jury round of the Track 1 report recorded that an F1 taken from a precision and a
    recall that are already rounded to four decimals can differ from the count based value in the
    last digit. This function uses the counts.

    Args:
        true_positive: The count of correct predictions.
        predicted: The count of predictions.
        gold: The count of gold items.

    Returns:
        The F1, and zero when both denominators are zero.
    """
    if predicted + gold == 0:
        return 0.0
    return 2 * true_positive / (predicted + gold)


def score_gsc_plus(documents: list[dict[str, object]], lexicon: dict[str, str],
                   closure: dict[str, frozenset[str]], labels: dict[str, list[str]],
                   names: dict[str, str], branch: frozenset[str], live: frozenset[str],
                   alternates: dict[str, str],
                   replacements: dict[str, str]) -> dict[str, object]:
    """Score the mapping on the GSC+ gold standard corpus at the document level.

    Args:
        documents: The corpus records from `read_gsc_plus`.
        lexicon: The normalized label mapped to its term.
        closure: The reflexive ancestor set of every term.
        labels: The labels of every term.
        names: The name of every term.
        branch: The identifiers under `HP:0000118`.
        live: The live identifiers of the release.
        alternates: The secondary identifier mapped to its primary identifier.
        replacements: The obsolete identifier mapped to its replacement.

    Returns:
        The precision, the recall and the F1 against the ancestor reduced gold set and against the
        full gold set, the error mode counts and the identifier mapping counts.
    """
    strict = Counter()
    full = Counter()
    relation_counts: Counter[str] = Counter()
    false_positive_context: Counter[str] = Counter()
    false_negative_reason: Counter[str] = Counter()
    mapping = Counter()
    mirror_agrees = 0
    examples: list[dict[str, str]] = []
    for record in documents:
        text = str(record["text"])
        normalized = extract.normalize(text)
        hits = raw_hits(normalized, lexicon)
        matches = reduce_to_specific(hits, closure)
        mirror_agrees += int(verify_mirror(normalized, lexicon, closure))
        predicted = {term for term, _ in matches}
        matched_label = dict(matches)

        gold_mentions: dict[str, list[str]] = {}
        for annotated, mentions in dict(record["mentions"]).items():  # type: ignore[arg-type]
            mapping["annotated_identifiers"] += 1
            resolved = resolve_identifier(annotated, live, alternates, replacements)
            if not resolved:
                mapping["identifiers_with_no_live_target"] += 1
                continue
            if resolved != annotated:
                mapping["identifiers_mapped_to_a_new_identifier"] += 1
            if resolved not in branch:
                mapping["identifiers_outside_the_phenotypic_abnormality_branch"] += 1
                continue
            mapping["identifiers_kept"] += 1
            gold_mentions.setdefault(resolved, []).extend(mentions)

        gold_all = frozenset(gold_mentions)
        gold_specific = most_specific(set(gold_all), closure)

        strict["true_positive"] += len(predicted & gold_specific)
        strict["false_positive"] += len(predicted - gold_specific)
        strict["false_negative"] += len(gold_specific - predicted)
        full["true_positive"] += len(predicted & gold_all)
        full["false_positive"] += len(predicted - gold_all)
        full["false_negative"] += len(gold_all - predicted)

        for term in sorted(predicted - gold_specific):
            relation = hierarchy_relation(term, gold_all, closure)
            relation_counts[relation] += 1
            if relation == "unrelated":
                context = classify_context(text, matched_label[term])
                false_positive_context[context] += 1
                if context != "assertion" and len(examples) < 25:
                    examples.append({"doc_id": str(record["doc_id"]), "hpo_id": term,
                                     "hpo_label": names.get(term, ""),
                                     "matched_text": matched_label[term], "context": context})
        for term in sorted(gold_specific - predicted):
            false_negative_reason[classify_false_negative(
                term, gold_mentions[term], normalized, hits, lexicon, labels)] += 1

    precision = proportion_record(strict["true_positive"],
                                  strict["true_positive"] + strict["false_positive"])
    recall = proportion_record(strict["true_positive"],
                               strict["true_positive"] + strict["false_negative"])
    full_precision = proportion_record(full["true_positive"],
                                       full["true_positive"] + full["false_positive"])
    full_recall = proportion_record(full["true_positive"],
                                    full["true_positive"] + full["false_negative"])
    tolerant_true = strict["true_positive"] + relation_counts["over_general"] \
        + relation_counts["under_general"]
    tolerant_precision = proportion_record(
        tolerant_true, strict["true_positive"] + strict["false_positive"])
    return {
        "documents": len(documents),
        "annotation_rows": sum(int(record["annotation_rows"]) for record in documents),
        "mirror_agreements": mirror_agrees,
        "identifier_mapping": dict(sorted(mapping.items())),
        "primary_scoring": {
            "gold": "the annotated terms under HP:0000118, with the ancestors of another gold "
                    "term of the same document removed",
            "true_positive": strict["true_positive"],
            "false_positive": strict["false_positive"],
            "false_negative": strict["false_negative"],
            "precision": precision,
            "recall": recall,
            "f1": round(f1_from_counts(strict["true_positive"], precision["denominator"],
                                       recall["denominator"]), 4),
        },
        "full_gold_scoring": {
            "gold": "every annotated term under HP:0000118, with no ancestor removal",
            "true_positive": full["true_positive"],
            "false_positive": full["false_positive"],
            "false_negative": full["false_negative"],
            "precision": full_precision,
            "recall": full_recall,
            "f1": round(f1_from_counts(full["true_positive"], full_precision["denominator"],
                                       full_recall["denominator"]), 4),
        },
        "hierarchy_tolerant_precision": tolerant_precision,
        "false_positive_relation": dict(sorted(relation_counts.items())),
        "false_positive_context_of_an_unrelated_term": dict(sorted(
            false_positive_context.items())),
        "false_negative_reason": dict(sorted(false_negative_reason.items())),
        "false_positive_examples": examples,
    }


def score_orphanet(annotations: dict[str, set[str]], definitions: dict[str, dict[str, str]],
                   lexicon: dict[str, str], closure: dict[str, frozenset[str]],
                   branch: frozenset[str], live: frozenset[str], alternates: dict[str, str],
                   replacements: dict[str, str]) -> dict[str, object]:
    """Score the mapping on the Orphanet clinical definition text of a rare disease.

    The Orphanet annotation set describes the disease and not the text, so the fraction of the
    annotation set that the text names is a property of the text. That fraction is reported as text
    coverage and not as a recall of the mapping.

    Args:
        annotations: The Orpha code mapped to its curated HPO identifiers.
        definitions: The Orpha code mapped to its name and its definition text.
        lexicon: The normalized label mapped to its term.
        closure: The reflexive ancestor set of every term.
        branch: The identifiers under `HP:0000118`.
        live: The live identifiers of the release.
        alternates: The secondary identifier mapped to its primary identifier.
        replacements: The obsolete identifier mapped to its replacement.

    Returns:
        The strict precision, the hierarchy tolerant precision, the text coverage and the context of
        every unrelated false positive.
    """
    codes = sorted(set(annotations) & set(definitions), key=int)
    relation_counts: Counter[str] = Counter()
    context_counts: Counter[str] = Counter()
    covered = 0
    gold_total = 0
    predicted_total = 0
    exact = 0
    mirror_agrees = 0
    examples: list[dict[str, str]] = []
    for code in codes:
        text = definitions[code]["definition"]
        normalized = extract.normalize(text)
        hits = raw_hits(normalized, lexicon)
        matches = reduce_to_specific(hits, closure)
        mirror_agrees += int(verify_mirror(normalized, lexicon, closure))
        predicted = {term for term, _ in matches}
        matched_label = dict(matches)
        gold_all = frozenset(
            resolved for resolved in
            (resolve_identifier(term, live, alternates, replacements)
             for term in annotations[code])
            if resolved and resolved in branch)
        gold_specific = most_specific(set(gold_all), closure)
        gold_total += len(gold_specific)
        covered += len(predicted & gold_specific)
        predicted_total += len(predicted)
        exact += len(predicted & gold_all)
        for term in sorted(predicted - gold_all):
            relation = hierarchy_relation(term, gold_all, closure)
            relation_counts[relation] += 1
            if relation == "unrelated":
                context = classify_context(text, matched_label[term])
                context_counts[context] += 1
                if context != "assertion" and len(examples) < 25:
                    examples.append({"orpha_code": code, "hpo_id": term,
                                     "matched_text": matched_label[term], "context": context})
    tolerant = exact + relation_counts["over_general"] + relation_counts["under_general"]
    return {
        "disorders_with_an_annotation_set": len(annotations),
        "disorders_with_a_definition": len(definitions),
        "disorders_scored": len(codes),
        "mirror_agreements": mirror_agrees,
        "predicted_terms": predicted_total,
        "gold_terms_of_the_scored_disorders": gold_total,
        "strict_precision": proportion_record(exact, predicted_total),
        "hierarchy_tolerant_precision": proportion_record(tolerant, predicted_total),
        "text_coverage_of_the_curated_set": proportion_record(covered, gold_total),
        "false_positive_relation": dict(sorted(relation_counts.items())),
        "false_positive_context_of_an_unrelated_term": dict(sorted(context_counts.items())),
        "false_positive_examples": examples,
    }


def score_label_control(lexicon: dict[str, str], closure: dict[str, frozenset[str]],
                        branch: frozenset[str], names: dict[str, str]) -> dict[str, object]:
    """Run the positive control. A document that holds one label must return that term.

    Args:
        lexicon: The normalized label mapped to its term.
        closure: The reflexive ancestor set of every term.
        branch: The identifiers under `HP:0000118`.
        names: The name of every term.

    Returns:
        The sample size, the recovery rate of the term and the mean count of extra terms.
    """
    candidates = sorted(term for term in branch if names.get(term))
    sample = random.Random(LABEL_CONTROL_SEED).sample(
        candidates, min(LABEL_CONTROL_SAMPLE, len(candidates)))
    recovered = 0
    extra = 0
    missed: list[dict[str, str]] = []
    for term in sample:
        normalized = extract.normalize(f"The patient shows {names[term]}.")
        predicted = {found for found, _ in extract.match_terms(normalized, lexicon, closure)}
        if term in predicted:
            recovered += 1
        elif len(missed) < 25:
            missed.append({"hpo_id": term, "hpo_label": names[term],
                           "returned": ",".join(sorted(predicted))})
        extra += len(predicted - {term})
    return {
        "sample_size": len(sample),
        "seed": LABEL_CONTROL_SEED,
        "term_recovered": proportion_record(recovered, len(sample)),
        "mean_extra_terms_per_document": round(extra / len(sample), 4) if sample else 0.0,
        "missed_examples": missed,
    }


def tool_versions() -> dict[str, str]:
    """Collect the version of every tool of this run.

    Returns:
        The interpreter version, the `uv` version and the checksum of the two code files.
    """
    try:
        uv_version = subprocess.run(["uv", "--version"], capture_output=True, text=True,
                                    check=True, timeout=60).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        uv_version = "not available"
    return {
        "python": sys.version.split()[0],
        "uv": uv_version,
        "measurement_script_sha256": sha256_of(Path(__file__).resolve()),
        "mapping_under_test_sha256": sha256_of(BIN / "extract_case_phenotype.py"),
        "mapping_under_test": "pipelines/track1_b0/bin/extract_case_phenotype.py",
    }


def evaluate_gate(gsc: dict[str, object], synthetic: dict[str, object],
                  control: dict[str, object], orphanet: dict[str, object],
                  sources: dict[str, object]) -> dict[str, object]:
    """Apply the gate conditions of the record to the result.

    Args:
        gsc: The GSC+ arm result.
        synthetic: The synthetic arm result.
        control: The positive control result.
        orphanet: The Orphanet arm result.
        sources: The provenance of every downloaded input.

    Returns:
        Each condition with its statement, its value, its verdict, and the outcome of the gate.
    """
    primary_f1 = dict(gsc["primary_scoring"])["f1"]  # type: ignore[index]
    full_f1 = dict(gsc["full_gold_scoring"])["f1"]  # type: ignore[index]
    negated = dict(synthetic["negated_false_positive_rate"])["value"]
    family = dict(synthetic["family_history_false_positive_rate"])["value"]
    differential = dict(synthetic["differential_false_positive_rate"])["value"]
    recovered = dict(control["term_recovered"])["value"]
    mirror_ok = int(gsc["mirror_agreements"]) == int(gsc["documents"]) \
        and int(synthetic["mirror_agreements"]) == int(synthetic["documents"])
    mirror_value = (f"{gsc['mirror_agreements']} of {gsc['documents']} GSC+ documents and "
                    f"{synthetic['mirror_agreements']} of {synthetic['documents']} "
                    f"synthetic notes")
    if orphanet.get("ran"):
        mirror_ok = mirror_ok and int(orphanet["mirror_agreements"]) == int(
            orphanet["disorders_scored"])
        mirror_value += (f" and {orphanet['mirror_agreements']} of "
                         f"{orphanet['disorders_scored']} Orphanet definition texts")
    conditions = [
        {"id": "B1",
         "statement": "The ground truth is a published gold standard corpus for HPO concept "
                      "recognition with a verified open licence, and the measurement does not "
                      "rest on a synthetic corpus alone.",
         "value": f"{sources['gsc_plus']['url']} at {sources['gsc_plus']['cached_utc']}",
         "prospective": True,
         "holds": bool(gsc["documents"])},
        {"id": "B2",
         "statement": "The local mirror of the matching loop returns what the shipped "
                      "extract_case_phenotype.match_terms returns on every document of every arm.",
         "value": mirror_value,
         "prospective": True,
         "holds": mirror_ok},
        {"id": "B3",
         "statement": "The document level F1 on the published corpus reaches 0.56 under both "
                      "scorings, the ancestor reduced gold set and the full gold set. The "
                      "threshold 0.56 is the F-measure that Bio-LarK CR reached on the GSC, as "
                      "Lobo et al. 2017 report it, PMID 29250549. That value is a span level "
                      "F-measure and this value is a document level term set F1, so the two "
                      "units differ.",
         "value": {"ancestor_reduced_gold_f1": primary_f1, "full_gold_f1": full_f1,
                   "threshold": 0.56},
         "prospective": False,
         "prospective_note": "The GSC+ arm ran before this record entered the repository. The "
                             "threshold comes from an external publication and not from the "
                             "measured value.",
         "holds": primary_f1 >= 0.56 and full_f1 >= 0.56},
        {"id": "B4",
         "statement": "The false positive rate of a negated statement on the synthetic corpus "
                      "is 0.",
         "value": negated,
         "prospective": True,
         "prospective_note": "The synthetic corpus had not run. The source of the mapping holds "
                             "no negation rule, so the outcome was predictable from the code.",
         "holds": negated == 0.0},
        {"id": "B5",
         "statement": "The false positive rate of a family history statement on the synthetic "
                      "corpus is 0.",
         "value": family,
         "prospective": True,
         "prospective_note": "The synthetic corpus had not run. The source of the mapping holds "
                             "no family history rule, so the outcome was predictable from the "
                             "code.",
         "holds": family == 0.0},
        {"id": "B6",
         "statement": "The false positive rate of a differential diagnosis statement on the "
                      "synthetic corpus is 0.",
         "value": differential,
         "prospective": True,
         "holds": differential == 0.0},
        {"id": "B7",
         "statement": "The positive control recovers the term of its own label in at least 0.95 "
                      "of the sampled terms.",
         "value": recovered,
         "prospective": True,
         "holds": recovered >= 0.95},
    ]
    met = all(item["holds"] for item in conditions)
    return {
        "conditions": conditions,
        "outcome": "GATE MET" if met else "GATE NOT MET",
        "external_anchor": {
            "bio_lark_cr_f_measure_on_the_gsc": 0.56,
            "ihp_f_measure_on_the_gsc": 0.65,
            "ihp_f_measure_on_the_gsc_plus": 0.863,
            "source": "Lobo et al. 2017, PMID 29250549, DOI 10.1155/2017/8565739",
            "unit": "A span level F-measure on the corpus test split. The unit differs from the "
                    "document level term set F1 of this measurement.",
        },
        "phenotype_weight_decision": {
            "current_weight": 0.0,
            "source_of_the_current_weight": "docs/benchmarks/t1_phenotype_mismatch_v1.md, "
                                            "condition A3",
            "rule_recorded_before_the_result":
                "The weight stays at 0.0 when the gate is not met, and this measurement then adds "
                "a second and independent reason for weight 0.0. The weight returns to the four "
                "conditions A1 to A4 of the mismatch record when the gate is met. A1 to A4 are "
                "GATE NOT MET, because A3 failed at 93.0 rank points against 94.0. The weight "
                "therefore stays at 0.0 in both states.",
            "effect": "unchanged",
        },
    }


def main() -> int:
    """Run the four arms, write the receipt and print the summary.

    Returns:
        Zero when every arm ran.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hpo-dir", type=Path, default=ROOT / "tracking/cache/hpo_v2026_09_01")
    parser.add_argument("--hpo-release", default="v2026-09-01")
    parser.add_argument("--cache-dir", type=Path, default=ROOT / "tracking/cache/t1_hpo_mapping")
    parser.add_argument("--receipt", type=Path,
                        default=ROOT / "tracking/receipts/t1_phenotype_mapping_accuracy_v1.json")
    parser.add_argument("--skip-orphanet", action="store_true",
                        help="Run the GSC+ arm, the synthetic arm and the control only.")
    args = parser.parse_args()

    obo = args.hpo_dir / "hp.obo"
    parents, names, labels = extract.parse_obo(obo)
    closure = extract.ancestor_closure(parents)
    lexicon = extract.build_lexicon(labels, closure)
    live = frozenset(parents)
    branch = frozenset(term for term in live
                       if extract.PHENOTYPIC_ABNORMALITY in closure.get(term, frozenset()))
    alternates, replacements, obsolete = parse_obo_identifier_history(obo)

    sources: dict[str, object] = {}
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    archive = args.cache_dir / "gsc_plus.rar"
    sources["gsc_plus"] = fetch_to_cache(GSC_PLUS_URL, archive)
    corpus_root = args.cache_dir / "gsc_plus"
    if not (corpus_root / "Text").is_dir():
        corpus_root.mkdir(parents=True, exist_ok=True)
        subprocess.run(["bsdtar", "-xf", str(archive), "-C", str(corpus_root)], check=True)
    gsc_documents = read_gsc_plus(corpus_root)
    gsc = score_gsc_plus(gsc_documents, lexicon, closure, labels, names, branch, live,
                         alternates, replacements)

    synthetic = score_synthetic(build_synthetic_corpus(), lexicon, closure)
    control = score_label_control(lexicon, closure, branch, names)

    orphanet: dict[str, object] = {"ran": False,
                                   "reason": "The run passed --skip-orphanet."}
    if not args.skip_orphanet:
        product4 = args.cache_dir / "en_product4.xml"
        product1 = args.cache_dir / "en_product1.xml"
        sources["orphanet_en_product4"] = fetch_to_cache(ORPHANET_PRODUCT4_URL, product4)
        sources["orphanet_en_product1"] = fetch_to_cache(ORPHANET_PRODUCT1_URL, product1)
        orphanet = score_orphanet(parse_orphanet_annotations(product4),
                                  parse_orphanet_definitions(product1), lexicon, closure,
                                  branch, live, alternates, replacements)
        orphanet["ran"] = True

    gate = evaluate_gate(gsc, synthetic, control, orphanet, sources)
    receipt = {
        "receipt_id": "t1_phenotype_mapping_accuracy_v1",
        "created_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "track": "t1_phenotype",
        "defect": "R2-8 of reports/jury/defect_ledger_track1_round2_2026-09-17.md",
        "data_class": "A",
        "gate_record": "docs/benchmarks/t1_phenotype_mapping_accuracy_v1.md",
        "receipt_scope": "Public ground truth only. No patient file and no term of the child "
                         "enters any step. The synthetic notes describe no real person.",
        "tool_versions": tool_versions(),
        "hpo": {
            "release": args.hpo_release,
            "hp_obo_sha256": sha256_of(obo),
            "terms": len(live),
            "terms_under_HP_0000118": len(branch),
            "obsolete_terms": len(obsolete),
            "lexicon_labels": len(lexicon),
        },
        "parameters": {
            "phenotypic_abnormality_root": extract.PHENOTYPIC_ABNORMALITY,
            "minimum_label_characters": extract.MIN_LABEL_CHARACTERS,
            "stop_labels": sorted(extract.STOP_LABELS),
            "unit_of_scoring": "one HPO identifier per document, because the mapping returns a "
                               "term set per document and not a span",
            "negation_cues": list(NEGATION_CUES),
            "negation_window_tokens": NEGATION_WINDOW,
            "family_history_cues": list(FAMILY_HISTORY_CUES),
            "differential_cues": list(DIFFERENTIAL_CUES),
            "context_rule": "A label that stands in one plain sentence of a document is reported "
                            "as an assertion, even when it also stands in a negated sentence.",
            "interval": "Wilson score interval at 95 percent on every proportion.",
            "label_control_sample": LABEL_CONTROL_SAMPLE,
            "label_control_seed": LABEL_CONTROL_SEED,
        },
        "ground_truth": {
            "gsc_plus": {
                "name": "GSC+, the extended gold standard corpus for HPO concept recognition",
                "documents": len(gsc_documents),
                "origin": "Groza et al. 2015, Database, PMID 25725061, DOI "
                          "10.1093/database/bav005, licence CC BY. Extended by Lobo et al. 2017, "
                          "Computational and Mathematical Methods in Medicine, PMID 29250549, DOI "
                          "10.1155/2017/8565739, licence CC BY.",
                "distribution": "The GitHub repository lasigeBioTM/IHP, MIT licence. The original "
                                "address bio-lark.org/hpo_res.html no longer serves the corpus.",
                "distribution_commit": GSC_PLUS_COMMIT,
                "distribution_commit_date": "2018-11-02T11:08:51Z",
                "note": "That repository holds the same archive twice, as GSC+.rar and as "
                        "GSC_v2.rar. Both names carry the git blob "
                        "cb94b6ff5a3e2be46a465ac3f387a2f948854c00.",
            },
            "orphanet": {
                "name": "Orphanet en_product4.xml and en_product1.xml",
                "licence": "CC-BY-4.0, stated in the Availability element of both files",
                "database_version": "1.3.42 / 4.1.8, export date 2026-06-23",
            },
            "synthetic": {
                "name": "The synthetic clinical note corpus of this measurement",
                "documents": len(SYNTHETIC_CORPUS),
                "licence": "Invented for this measurement. It describes no real person.",
            },
        },
        "sources": sources,
        "cache_policy": "Every download stays under tracking/cache/t1_hpo_mapping/. The "
                        "repository .gitignore holds the line tracking/cache/, so no cache file "
                        "enters a commit. The receipt carries the url, the date and the SHA-256 "
                        "of each download.",
        "arms": {
            "gsc_plus": gsc,
            "orphanet": orphanet,
            "synthetic": synthetic,
            "label_control": control,
        },
        "decision_gate": gate,
        "verdict": "MEASURED",
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n")

    summary = {
        "gsc_plus_primary": gsc["primary_scoring"],
        "gsc_plus_full_gold": gsc["full_gold_scoring"],
        "gsc_plus_false_positive_relation": gsc["false_positive_relation"],
        "gsc_plus_false_positive_context": gsc["false_positive_context_of_an_unrelated_term"],
        "gsc_plus_false_negative_reason": gsc["false_negative_reason"],
        "orphanet": {key: value for key, value in orphanet.items()
                     if key != "false_positive_examples"},
        "synthetic": {key: value for key, value in synthetic.items()
                      if key != "unintended_lexicon_matches"},
        "synthetic_unintended_matches": synthetic.get("unintended_lexicon_matches"),
        "label_control": {key: value for key, value in control.items()
                          if key != "missed_examples"},
        "gate": gate,
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
