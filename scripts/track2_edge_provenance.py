"""Extract the evidence sentences behind the graph edges that support a Track 2 candidate.

Publication counts do not establish whether a text-mined relationship has the intended meaning.
This script retrieves source sentences for a configured set of public entity pairs and checks
caller-supplied judgments against those sentences. It includes no completed judgments or rankings.

``extract`` fetches the PubTator3 sentences of every text-mined edge in `EDGE_SPEC` and writes one
row per sentence to a TSV. Each row keeps the surface form that the recogniser mapped onto each of
the two entities, because a wrong surface form is the signature of a mapping error. An Open
Targets edge carries no sentence, so its row holds the curated mechanism of action instead.

* An edge in `all` mode is read in full. The script fetches every result page of the pair.
* An edge in `sample` mode has hundreds or thousands of publications. The script fetches the first
  `SAMPLE_PAGES` pages in endpoint order and draws `SAMPLE_SIZE` of them with the fixed seed
  `SAMPLE_SEED`. The `in_sample` column marks the drawn rows.

``audit-seed`` is a second, mechanical check. It reads every sentence of every literature edge of
the seed node and reports whether the text that the recogniser mapped onto the seed names the
syndrome. It also counts, for the head of a ranking, how often the best path leaves the seed on an
edge that no sentence supports.

``receipt`` reads the per-edge judgments that the reading produced, checks that every
edge of the spec has one judgment with an allowed class, that each judgment quotes a sentence that
was read, and writes a receipt with the counts and the digest of each input.

The script sends public entity accessions to the public PubTator3 search endpoint. No patient data
enters it. The API responses are cached under `reports/track2/raw/`, which git ignores.

Usage:
    uv run python scripts/track2_edge_provenance.py extract
    uv run python scripts/track2_edge_provenance.py audit-seed
    uv run python scripts/track2_edge_provenance.py receipt
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import random
import re
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402

SEARCH_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/"
CACHE_PATH = common.RAW_DIR / "edge_provenance_evidence.json"
SENTENCES_PATH = common.ROOT / "reports" / "track2" / "provenance" / "edge_sentences_2026_09_15.tsv"
JUDGMENTS_PATH = common.ROOT / "tracking" / "track2_edge_provenance_2026_09_15.tsv"
RECEIPT_PATH = common.ROOT / "tracking" / "receipts" / "t2_edge_provenance_v1.json"
DRUG_EVIDENCE_PATH = common.GRAPH_DIR / "drug_evidence.tsv"
EDGES_PATH = common.GRAPH_DIR / "edges.tsv"
NODES_PATH = common.GRAPH_DIR / "nodes.tsv"
RANKING_PATH = (common.ROOT / "reports" / "track2" / "provenance"
                / "ranking_run_f_2026_09_15.tsv")
UNCONSTRAINED_RANKING_PATH = (common.ROOT / "reports" / "track2" / "provenance"
                              / "ranking_run_a_2026_09_15.tsv")
NAD_LABELS = ("Nicotinamide Mononucleotide", "NAD", "Niacinamide")

SAMPLE_SEED = 20260915
SAMPLE_SIZE = 10
SAMPLE_PAGES = 3
MAX_PAGES_ALL = 20

CLASSES = ("SUPPORTED", "WRONG DIRECTION", "NOT AN INDICATION", "OFF TOPIC", "UNCLEAR")

SENTENCE_HEADER = ["edge_id", "group", "supports", "evidence_kind", "relation", "node_1", "node_2",
                   "graph_evidence_count", "graph_specificity_npmi", "search_count", "pages_read",
                   "read_mode", "item_index", "in_sample", "pmid", "year", "journal", "title",
                   "surface_1", "surface_2", "sentence"]

MVA = "@DISEASE_Mosaic_variegated_aneuploidy_syndrome"
CIN = "@DISEASE_Chromosomal_Instability"

# The groups organize public entity pairs for evidence retrieval. Group labels and edge IDs
# identify configuration entries; they do not record candidate ranks or analysis verdicts.
EDGE_SPEC: list[dict[str, str]] = [
    {"edge_id": "A1", "group": "A", "supports": "Configured comparison edge",
     "relation": "associate", "acc_1": MVA, "acc_2": "@GENE_BUB1B",
     "node_1": "DISEASE:MESH:C536987", "node_2": "GENE:BUB1B", "mode": "all"},
    {"edge_id": "A2", "group": "A", "supports": "Configured comparison edge",
     "relation": "associate", "acc_1": "@GENE_BUB1B", "acc_2": "@GENE_SIRT2",
     "node_1": "GENE:BUB1B", "node_2": "GENE:SIRT2", "mode": "all"},
    {"edge_id": "A3", "group": "A", "supports": "Configured comparison edge",
     "relation": "positive_correlate", "acc_1": "@GENE_BUB1B", "acc_2": "@GENE_SIRT2",
     "node_1": "GENE:BUB1B", "node_2": "GENE:SIRT2", "mode": "all"},
    {"edge_id": "A4", "group": "A", "supports": "Configured comparison edge",
     "relation": "positive_correlate", "acc_1": "@CHEMICAL_Nicotinamide_Mononucleotide",
     "acc_2": "@GENE_BUB1B", "node_1": "CHEMICAL:MESH:D009537", "node_2": "GENE:BUB1B",
     "mode": "all"},
    {"edge_id": "A5", "group": "A", "supports": "Configured comparison edge",
     "relation": "associate", "acc_1": "@CHEMICAL_Nicotinamide_Mononucleotide",
     "acc_2": "@GENE_SIRT2", "node_1": "CHEMICAL:MESH:D009537", "node_2": "GENE:SIRT2",
     "mode": "all"},
    {"edge_id": "A6", "group": "A", "supports": "Configured comparison edge",
     "relation": "positive_correlate", "acc_1": "@CHEMICAL_Nicotinamide_Mononucleotide",
     "acc_2": "@GENE_SIRT2", "node_1": "CHEMICAL:MESH:D009537", "node_2": "GENE:SIRT2",
     "mode": "all"},
    {"edge_id": "A7", "group": "A", "supports": "Configured comparison edge",
     "relation": "associate", "acc_1": "@CHEMICAL_NAD", "acc_2": "@GENE_BUB1B",
     "node_1": "CHEMICAL:MESH:D009243", "node_2": "GENE:BUB1B", "mode": "all"},
    {"edge_id": "A8", "group": "A", "supports": "Configured comparison edge",
     "relation": "negative_correlate", "acc_1": "@CHEMICAL_NAD", "acc_2": "@GENE_BUB1B",
     "node_1": "CHEMICAL:MESH:D009243", "node_2": "GENE:BUB1B", "mode": "all"},
    {"edge_id": "A9", "group": "A", "supports": "Configured comparison edge",
     "relation": "associate", "acc_1": "@CHEMICAL_NAD", "acc_2": "@GENE_SIRT2",
     "node_1": "CHEMICAL:MESH:D009243", "node_2": "GENE:SIRT2", "mode": "all"},
    {"edge_id": "A10", "group": "A", "supports": "Configured comparison edge",
     "relation": "positive_correlate", "acc_1": "@CHEMICAL_NAD", "acc_2": "@GENE_SIRT2",
     "node_1": "CHEMICAL:MESH:D009243", "node_2": "GENE:SIRT2", "mode": "all"},
    {"edge_id": "A11", "group": "A", "supports": "Configured comparison edge",
     "relation": "negative_correlate", "acc_1": "@CHEMICAL_NAD", "acc_2": "@GENE_SIRT2",
     "node_1": "CHEMICAL:MESH:D009243", "node_2": "GENE:SIRT2", "mode": "all"},
    {"edge_id": "B1", "group": "B", "supports": "Configured comparison edge",
     "relation": "negative_correlate", "acc_1": "@CHEMICAL_Niacinamide", "acc_2": "@GENE_SIRT2",
     "node_1": "CHEMICAL:MESH:D009536", "node_2": "GENE:SIRT2", "mode": "all"},
    {"edge_id": "B2", "group": "B", "supports": "Configured comparison edge",
     "relation": "associate", "acc_1": "@CHEMICAL_Niacinamide", "acc_2": "@GENE_SIRT2",
     "node_1": "CHEMICAL:MESH:D009536", "node_2": "GENE:SIRT2", "mode": "all"},
    {"edge_id": "C1", "group": "C", "supports": "Configured comparison edge",
     "relation": "associate", "acc_1": MVA, "acc_2": "@GENE_APC",
     "node_1": "DISEASE:MESH:C536987", "node_2": "GENE:APC", "mode": "all"},
    {"edge_id": "C2", "group": "C", "supports": "Configured comparison edge",
     "relation": "associate", "acc_1": "@GENE_APC", "acc_2": "@GENE_CTNNB1",
     "node_1": "GENE:APC", "node_2": "GENE:CTNNB1", "mode": "sample"},
    {"edge_id": "C3", "group": "C", "supports": "Configured comparison edge", "kind": "curated",
     "relation": "inhibitor", "node_1": "CHEMBL:CHEMBL4297813", "node_2": "GENE:CTNNB1"},
    {"edge_id": "C4", "group": "C",
     "supports": "Configured comparison edge",
     "relation": "associate", "acc_1": MVA, "acc_2": "@GENE_EGFR",
     "node_1": "DISEASE:MESH:C536987", "node_2": "GENE:EGFR", "mode": "all"},
    {"edge_id": "C5", "group": "C", "supports": "Configured comparison edge",
     "relation": "positive_correlate", "acc_1": "@GENE_EGFR", "acc_2": "@GENE_AKT1",
     "node_1": "GENE:EGFR", "node_2": "GENE:AKT1", "mode": "sample"},
    {"edge_id": "C6", "group": "C", "supports": "Configured comparison edge", "kind": "curated",
     "relation": "inhibitor", "node_1": "CHEMBL:CHEMBL1079175", "node_2": "GENE:AKT1"},
    {"edge_id": "C7", "group": "C", "supports": "Configured comparison edge", "kind": "curated",
     "relation": "inhibitor", "node_1": "CHEMBL:CHEMBL2016893", "node_2": "GENE:AKT1"},
    {"edge_id": "C8", "group": "C", "supports": "Configured comparison edge", "kind": "curated",
     "relation": "inhibitor", "node_1": "CHEMBL:CHEMBL2177390", "node_2": "GENE:AKT1"},
    {"edge_id": "C9", "group": "C", "supports": "Configured comparison edge", "kind": "curated",
     "relation": "inhibitor", "node_1": "CHEMBL:CHEMBL2219422", "node_2": "GENE:AKT1"},
    {"edge_id": "C10", "group": "C", "supports": "Configured comparison edge", "kind": "curated",
     "relation": "inhibitor", "node_1": "CHEMBL:CHEMBL2325741", "node_2": "GENE:AKT1"},
    {"edge_id": "C11", "group": "C", "supports": "Configured comparison edge", "kind": "curated",
     "relation": "inhibitor", "node_1": "CHEMBL:CHEMBL245030", "node_2": "GENE:AKT1"},
    {"edge_id": "C12", "group": "C", "supports": "Configured comparison edge", "kind": "curated",
     "relation": "inhibitor", "node_1": "CHEMBL:CHEMBL3137336", "node_2": "GENE:AKT1"},
    {"edge_id": "C13", "group": "C", "supports": "Configured comparison edge", "kind": "curated",
     "relation": "inhibitor", "node_1": "CHEMBL:CHEMBL3545134", "node_2": "GENE:AKT1"},
    {"edge_id": "C14", "group": "C", "supports": "Configured comparison edge", "kind": "curated",
     "relation": "inhibitor", "node_1": "CHEMBL:CHEMBL4297188", "node_2": "GENE:AKT1"},
    {"edge_id": "D1", "group": "D", "supports": "Configured comparison edge",
     "relation": "associate", "acc_1": MVA, "acc_2": "@GENE_PLK1",
     "node_1": "DISEASE:MESH:C536987", "node_2": "GENE:PLK1", "mode": "all"},
    {"edge_id": "D2", "group": "D", "supports": "Configured comparison edge", "kind": "curated",
     "relation": "inhibitor", "node_1": "CHEMBL:CHEMBL1233528", "node_2": "GENE:PLK1"},
    {"edge_id": "D3", "group": "D", "supports": "Configured comparison edge",
     "relation": "associate", "acc_1": MVA, "acc_2": "@GENE_AURKA",
     "node_1": "DISEASE:MESH:C536987", "node_2": "GENE:AURKA", "mode": "all"},
    {"edge_id": "D4", "group": "D", "supports": "Configured comparison edge", "kind": "curated",
     "relation": "inhibitor", "node_1": "CHEMBL:CHEMBL483158", "node_2": "GENE:AURKA"},
    {"edge_id": "E1", "group": "E", "supports": "Configured comparison edge",
     "relation": "stimulate", "acc_1": MVA, "acc_2": "@GENE_CASP3",
     "node_1": "DISEASE:MESH:C536987", "node_2": "GENE:CASP3", "mode": "all"},
    {"edge_id": "E2", "group": "E", "supports": "Configured comparison edge",
     "relation": "positive_correlate", "acc_1": "@CHEMICAL_Lipopolysaccharides",
     "acc_2": "@GENE_CASP3", "node_1": "CHEMICAL:MESH:D008070", "node_2": "GENE:CASP3",
     "mode": "sample"},
    {"edge_id": "E3", "group": "E", "supports": "Configured comparison edge",
     "relation": "negative_correlate", "acc_1": "@CHEMICAL_Gefitinib", "acc_2": "@GENE_EGFR",
     "node_1": "CHEMICAL:MESH:D000077156", "node_2": "GENE:EGFR", "mode": "sample"},
    {"edge_id": "E4", "group": "E", "supports": "Configured comparison edge",
     "relation": "associate", "acc_1": "@CHEMICAL_Docetaxel", "acc_2": "@GENE_BUB1B",
     "node_1": "CHEMICAL:MESH:D000077143", "node_2": "GENE:BUB1B", "mode": "all"},
    {"edge_id": "E5", "group": "E", "supports": "Configured comparison edge",
     "relation": "negative_correlate", "acc_1": "@CHEMICAL_Acetylcysteine", "acc_2": "@GENE_CASP3",
     "node_1": "CHEMICAL:MESH:D000111", "node_2": "GENE:CASP3", "mode": "sample"},
    {"edge_id": "E6", "group": "E", "supports": "Configured comparison edge",
     "relation": "stimulate", "acc_1": MVA, "acc_2": "@GENE_FGF2",
     "node_1": "DISEASE:MESH:C536987", "node_2": "GENE:FGF2", "mode": "all"},
    {"edge_id": "E7", "group": "E", "supports": "Configured comparison edge",
     "relation": "associate", "acc_1": "@CHEMICAL_Heparin", "acc_2": "@GENE_FGF2",
     "node_1": "CHEMICAL:MESH:D006493", "node_2": "GENE:FGF2", "mode": "sample"},
    {"edge_id": "F1", "group": "F", "supports": "Configured comparison edge",
     "relation": "treat", "acc_1": "@CHEMICAL_oxybutynin", "acc_2": CIN,
     "node_1": "CHEMICAL:MESH:C005419", "node_2": "DISEASE:MESH:D043171", "mode": "all"},
    {"edge_id": "F2", "group": "F", "supports": "Configured comparison edge",
     "relation": "treat", "acc_1": "@CHEMICAL_Steroids", "acc_2": CIN,
     "node_1": "CHEMICAL:MESH:D013256", "node_2": "DISEASE:MESH:D043171", "mode": "all"},
    {"edge_id": "F3", "group": "F", "supports": "Configured comparison edge",
     "relation": "treat", "acc_1": "@CHEMICAL_Methylprednisolone", "acc_2": CIN,
     "node_1": "CHEMICAL:MESH:D008775", "node_2": "DISEASE:MESH:D043171", "mode": "all"},
    {"edge_id": "F4", "group": "F", "supports": "Configured comparison edge",
     "relation": "treat", "acc_1": "@CHEMICAL_Propofol", "acc_2": CIN,
     "node_1": "CHEMICAL:MESH:D015742", "node_2": "DISEASE:MESH:D043171", "mode": "all"},
    {"edge_id": "F5", "group": "F", "supports": "Configured comparison edge",
     "relation": "treat", "acc_1": "@CHEMICAL_Rituximab", "acc_2": CIN,
     "node_1": "CHEMICAL:MESH:D000069283", "node_2": "DISEASE:MESH:D043171", "mode": "all"},
]

SEED = "DISEASE:MESH:C536987"
ENTITIES_PATH = common.RAW_DIR / "pubtator_entities.json"
AUDIT_PATH = common.ROOT / "reports" / "track2" / "provenance" / "seed_edge_surfaces_2026_09_15.tsv"
AUDIT_HEADER = ["relation", "partner", "graph_evidence_count", "search_count", "item_index", "pmid",
                "seed_surface", "rule_names_syndrome", "reader_override", "names_syndrome",
                "sentence"]
# These configured exceptions disambiguate "MVA" as motor vehicle accident or modified vaccinia
# Ankara. Each key contains a public partner accession and publication identifier.
OTHER_MEANING = {
    ("CHEMICAL:MESH:D000438", "36268435"): "motor vehicle accident",
    ("CHEMICAL:MESH:D000438", "23434842"): "motor vehicle accident",
    ("CHEMICAL:MESH:D000438", "18333825"): "motor vehicle accident",
    ("CHEMICAL:MESH:D003042", "23434842"): "motor vehicle accident",
    ("CHEMICAL:MESH:D008694", "23434842"): "motor vehicle accident",
    ("CHEMICAL:MESH:D000597", "22690067"): "motor vehicle accident",
    ("GENE:CCL2", "26992431"): "modified vaccinia Ankara",
    ("GENE:IL6", "17537851"): "modified vaccinia Ankara",
}
SYNDROME_NAME = re.compile(r"variegated|\bMVAS?\d?\b|\bPCS\b|premature chromatid separation|"
                           r"premature centromere division", re.IGNORECASE)

_MARKED_SURFACE = re.compile(r"@<m>([A-Z]+_[^<]*)</m>\s*@[A-Z]+_\S+\s*@@@(.*?)@@@", re.DOTALL)


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


def fetch_pages(client, spec: dict[str, str], max_pages: int) -> dict:
    """Fetch the search results of one entity pair, page by page.

    Args:
        client: The HTTP client.
        spec: One entry of `EDGE_SPEC`.
        max_pages: The largest number of pages to read.

    Returns:
        The query, the reported count, the number of pages read and one record per result.
    """
    text = common.pair_query(spec["acc_1"], spec["acc_2"], spec["relation"])
    items: list[dict] = []
    total, pages, page = 0, 1, 1
    while page <= min(max_pages, pages):
        payload = common.get_json(client, SEARCH_URL, {"text": text, "page": page})
        if not isinstance(payload, dict):
            break
        total = int(payload.get("count", 0) or 0)
        pages = int(payload.get("total_pages", 1) or 1)
        for result in payload.get("results", []) or []:
            items.append({"pmid": str(result.get("pmid", "") or ""),
                          "year": str(result.get("date", "") or "")[:4],
                          "journal": result.get("journal", "") or "",
                          "title": result.get("title", "") or "",
                          "text_hl": result.get("text_hl", "") or ""})
        page += 1
        time.sleep(common.DEFAULT_PAUSE_S)
    return {"query": text, "search_count": total, "total_pages": pages, "pages_read": page - 1,
            "items": items}


def graph_rows(edges_path: Path) -> dict[tuple[str, str, str], dict]:
    """Index the graph edges by their unordered node pair and relation.

    Args:
        edges_path: The `edges.tsv` of the graph.

    Returns:
        The summed evidence count, the specificity weight and the sources of each edge key.
    """
    per_source: dict[tuple[str, str, str], Counter] = {}
    npmi: dict[tuple[str, str, str], str] = {}
    with edges_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            first, second = sorted((row["source_id"], row["target_id"]))
            key = (first, second, row["relation"])
            per_source.setdefault(key, Counter())[row["source"]] += int(row["evidence_count"] or 0)
            npmi.setdefault(key, row["specificity_npmi"])
    # The consequence layer repeats a literature edge that the base layer already holds, so the
    # count of an edge is the largest count of one source, summed over both orientations.
    return {key: {"count": max(counts.values()), "npmi": npmi[key],
                  "sources": ",".join(sorted(counts))}
            for key, counts in per_source.items()}


def curated_mechanisms(path: Path) -> dict[tuple[str, str], dict]:
    """Read the Open Targets mechanism of each compound and target gene.

    Args:
        path: The `drug_evidence.tsv` of the graph.

    Returns:
        The drug record, keyed by the compound node and the target node.
    """
    with path.open() as handle:
        return {(row["drug_id"], row["target_gene_id"]): row
                for row in csv.DictReader(handle, delimiter="\t")}


def sample_indices(count: int, mode: str) -> set[int]:
    """Choose the result indices that the reader reads.

    Args:
        count: The number of fetched results.
        mode: `all` or `sample`.

    Returns:
        Every index in `all` mode, and `SAMPLE_SIZE` indices drawn with `SAMPLE_SEED` otherwise.
    """
    if mode == "all" or count <= SAMPLE_SIZE:
        return set(range(count))
    return set(random.Random(SAMPLE_SEED).sample(range(count), SAMPLE_SIZE))


def text_rows(spec: dict[str, str], fetched: dict, graph: dict) -> list[dict]:
    """Build the sentence rows of one text-mined edge.

    Args:
        spec: One entry of `EDGE_SPEC`.
        fetched: The cached search result of the pair.
        graph: The graph edge index from `graph_rows`.

    Returns:
        One row per fetched result.
    """
    first, second = sorted((spec["node_1"], spec["node_2"]))
    entry = graph.get((first, second, spec["relation"]), {})
    chosen = sample_indices(len(fetched["items"]), spec["mode"])
    rows = []
    for index, item in enumerate(fetched["items"]):
        surfaces = marked_surfaces(item["text_hl"])
        rows.append({
            "edge_id": spec["edge_id"], "group": spec["group"], "supports": spec["supports"],
            "evidence_kind": "pubtator3_sentence", "relation": spec["relation"],
            "node_1": spec["node_1"], "node_2": spec["node_2"],
            "graph_evidence_count": entry.get("count", ""),
            "graph_specificity_npmi": entry.get("npmi", ""),
            "search_count": fetched["search_count"], "pages_read": fetched["pages_read"],
            "read_mode": spec["mode"], "item_index": index,
            "in_sample": "yes" if index in chosen else "no",
            "pmid": item["pmid"], "year": item["year"], "journal": item["journal"],
            "title": item["title"],
            "surface_1": " | ".join(surfaces.get(spec["acc_1"], [])),
            "surface_2": " | ".join(surfaces.get(spec["acc_2"], [])),
            "sentence": common.plain_sentence(item["text_hl"]),
        })
    return rows


def curated_row(spec: dict[str, str], drugs: dict, graph: dict) -> dict:
    """Build the single row of one Open Targets edge, which carries no sentence.

    Args:
        spec: One entry of `EDGE_SPEC`.
        drugs: The mechanism index from `curated_mechanisms`.
        graph: The graph edge index from `graph_rows`.

    Returns:
        One row whose sentence column holds the curated mechanism of action.
    """
    record = drugs.get((spec["node_1"], spec["node_2"]), {})
    first, second = sorted((spec["node_1"], spec["node_2"]))
    entry = graph.get((first, second, spec["relation"]), {})
    return {
        "edge_id": spec["edge_id"], "group": spec["group"], "supports": spec["supports"],
        "evidence_kind": "open_targets_mechanism", "relation": spec["relation"],
        "node_1": spec["node_1"], "node_2": spec["node_2"],
        "graph_evidence_count": entry.get("count", ""), "graph_specificity_npmi": "",
        "search_count": "", "pages_read": "", "read_mode": "all", "item_index": 0,
        "in_sample": "yes", "pmid": "", "year": "", "journal": "",
        "title": record.get("drug_name", ""), "surface_1": record.get("drug_name", ""),
        "surface_2": spec["node_2"].split(":", 1)[-1],
        "sentence": (f"{record.get('action_type', '')}; {record.get('mechanism_of_action', '')}; "
                     f"stage on target {record.get('max_clinical_stage_on_target', '')}"),
    }


def run_extract(args: argparse.Namespace) -> int:
    """Fetch or reuse the sentences of every edge of the spec and write the sentence table.

    Args:
        args: The parsed command line.

    Returns:
        0 on success, 1 when an edge of the spec is absent from the graph.
    """
    cache = {} if args.refresh else (common.read_json(args.cache) or {})
    pairs = cache.get("pairs", {})
    graph = graph_rows(args.edges)
    drugs = curated_mechanisms(args.drug_evidence)
    client = common.make_client()
    rows: list[dict] = []
    missing = []
    for spec in EDGE_SPEC:
        first, second = sorted((spec["node_1"], spec["node_2"]))
        if (first, second, spec["relation"]) not in graph:
            missing.append(spec["edge_id"])
        if spec.get("kind") == "curated":
            rows.append(curated_row(spec, drugs, graph))
            continue
        max_pages = MAX_PAGES_ALL if spec["mode"] == "all" else SAMPLE_PAGES
        key = f"{spec['relation']}|{spec['acc_1']}|{spec['acc_2']}|{max_pages}"
        if key not in pairs:
            pairs[key] = fetch_pages(client, spec, max_pages)
            print(f"fetched {spec['edge_id']}: {pairs[key]['search_count']} results, "
                  f"{pairs[key]['pages_read']} pages")
        fetched = pairs[key]
        if spec["mode"] == "all" and fetched["pages_read"] < fetched["total_pages"]:
            print(f"WARNING: {spec['edge_id']} read {fetched['pages_read']} of "
                  f"{fetched['total_pages']} pages")
        rows.extend(text_rows(spec, fetched, graph))
    fetched_at = cache.get("fetched_at") or dt.datetime.now(dt.timezone.utc).isoformat(
        timespec="seconds")
    common.write_json(args.cache, {"fetched_at": fetched_at, "endpoint": SEARCH_URL,
                                   "pairs": pairs})
    common.write_tsv(args.out, SENTENCE_HEADER, rows)
    read = sum(1 for row in rows if row["in_sample"] == "yes")
    print(f"edges: {len(EDGE_SPEC)}, rows: {len(rows)}, rows to read: {read}")
    print(f"output: {args.out}")
    if missing:
        print(f"FAIL: edges absent from the graph: {', '.join(missing)}")
        return 1
    return 0


def accession_by_node(entities_path: Path) -> dict[str, str]:
    """Map each graph node identifier back to its PubTator3 accession.

    Args:
        entities_path: The resolved entity cache, `pubtator_entities.json`.

    Returns:
        The accession of each node identifier that the entity cache resolves.
    """
    entities = (common.read_json(entities_path) or {}).get("entities", {})
    mapping: dict[str, str] = {}
    for accession, entry in entities.items():
        node_id = common.node_id_for_accession(accession, entry.get("identifier", ""),
                                               entry.get("database", ""))
        if node_id:
            mapping.setdefault(node_id, accession)
    return mapping


def accession_for(node_id: str, mapping: dict[str, str]) -> str:
    """Return the PubTator3 accession of one graph node.

    Args:
        node_id: The graph node identifier, for example `GENE:BUB1B`.
        mapping: The index from `accession_by_node`.

    Returns:
        The accession. A gene node keeps its symbol, so its accession needs no lookup.
    """
    if node_id.startswith("GENE:"):
        return f"@GENE_{node_id.split(':', 1)[1]}"
    return mapping.get(node_id, "")


def names_the_syndrome(surface: str) -> bool:
    """Say whether a surface form names mosaic variegated aneuploidy syndrome.

    Args:
        surface: The surface forms that the recogniser mapped onto the seed, joined by ` | `.

    Returns:
        True when one surface form is a name or an abbreviation of the syndrome.
    """
    return bool(SYNDROME_NAME.search(surface or ""))


def seed_edges(edges_path: Path, seed: str) -> list[dict]:
    """List the distinct literature edges that touch the seed node.

    Args:
        edges_path: The graph `edges.tsv`.
        seed: The seed node identifier.

    Returns:
        One record per distinct source, target and relation, with its evidence count.
    """
    found: dict[tuple[str, str, str], dict] = {}
    with edges_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if seed not in (row["source_id"], row["target_id"]):
                continue
            if not row["source"].startswith("pubtator3"):
                continue
            key = (row["source_id"], row["target_id"], row["relation"])
            found.setdefault(key, {"source_id": key[0], "target_id": key[1], "relation": key[2],
                                   "evidence_count": row["evidence_count"]})
    return sorted(found.values(), key=lambda r: (-int(r["evidence_count"] or 0), r["target_id"],
                                                 r["relation"]))


def run_audit_seed(args: argparse.Namespace) -> int:
    """Read the surface form of the seed in every sentence of every literature edge of the seed.

    The check is mechanical. A sentence passes when the text that the recogniser mapped onto the
    seed names the syndrome. The distinct surface forms are printed so that a person can check
    the rule itself.

    Args:
        args: The parsed command line.

    Returns:
        0 on success, 1 when a seed edge has no accession to query.
    """
    cache = common.read_json(args.cache) or {}
    pairs = cache.get("pairs", {})
    accessions = accession_by_node(args.entities)
    seed_accession = accessions.get(args.seed, "")
    client = common.make_client()
    rows: list[dict] = []
    unresolved = []
    for edge in seed_edges(args.edges, args.seed):
        acc_1 = accession_for(edge["source_id"], accessions)
        acc_2 = accession_for(edge["target_id"], accessions)
        if not acc_1 or not acc_2:
            unresolved.append(f"{edge['source_id']}|{edge['target_id']}")
            continue
        spec = {"relation": edge["relation"], "acc_1": acc_1, "acc_2": acc_2}
        key = f"{edge['relation']}|{acc_1}|{acc_2}|{MAX_PAGES_ALL}"
        if key not in pairs:
            pairs[key] = fetch_pages(client, spec, MAX_PAGES_ALL)
        partner = edge["target_id"] if edge["source_id"] == args.seed else edge["source_id"]
        for index, item in enumerate(pairs[key]["items"]):
            surfaces = marked_surfaces(item["text_hl"])
            seed_surface = " | ".join(surfaces.get(seed_accession, []))
            rows.append({"relation": edge["relation"], "partner": partner,
                         "graph_evidence_count": edge["evidence_count"],
                         "search_count": pairs[key]["search_count"], "item_index": index,
                         "pmid": item["pmid"], "seed_surface": seed_surface,
                         "rule_names_syndrome": "yes" if names_the_syndrome(seed_surface) else "no",
                         "reader_override": OTHER_MEANING.get((partner, item["pmid"]), ""),
                         "names_syndrome": "yes" if names_the_syndrome(seed_surface) and (
                             partner, item["pmid"]) not in OTHER_MEANING else "no",
                         "sentence": common.plain_sentence(item["text_hl"])})
    common.write_json(args.cache, {"fetched_at": cache.get("fetched_at") or dt.datetime.now(
        dt.timezone.utc).isoformat(timespec="seconds"), "endpoint": SEARCH_URL, "pairs": pairs})
    common.write_tsv(args.audit_out, AUDIT_HEADER, rows)
    surfaces = Counter(row["seed_surface"].lower() for row in rows)
    for surface, count in surfaces.most_common():
        print(f"{count:4d}  {'yes' if names_the_syndrome(surface) else 'no ':3s}  {surface!r}")
    named = sum(1 for row in rows if row["names_syndrome"] == "yes")
    print(f"seed edges: {len(seed_edges(args.edges, args.seed))}, sentences: {len(rows)}, "
          f"naming the syndrome: {named}")
    print(f"output: {args.audit_out}")
    if unresolved:
        print(f"FAIL: no accession for {', '.join(unresolved)}")
        return 1
    return 0


def audit_summary(path: Path) -> dict:
    """Summarize the seed surface audit table for the receipt.

    Args:
        path: The audit TSV written by `run_audit_seed`.

    Returns:
        The counts of sentences, of edges, and of edges with no sentence that names the syndrome.
    """
    if not path.is_file():
        return {}
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_edge: dict[tuple[str, str], list[str]] = {}
    for row in rows:
        by_edge.setdefault((row["relation"], row["partner"]), []).append(row["names_syndrome"])
    return {
        "path": str(path.relative_to(common.ROOT)),
        "sha256_prefix": common.sha256_of(path)[:16],
        "rule": SYNDROME_NAME.pattern,
        "sentences": len(rows),
        "sentences_naming_the_syndrome": sum(1 for row in rows if row["names_syndrome"] == "yes"),
        "edges_with_a_sentence": len(by_edge),
        "edges_with_no_sentence_naming_the_syndrome": sum(
            1 for flags in by_edge.values() if "yes" not in flags),
    }


def first_hop_summary(ranking_path: Path, audit_path: Path, nodes_path: Path) -> dict:
    """Count how many top ranked paths leave the seed on an edge that never names the syndrome.

    The seed surface audit says which neighbour of the seed carries at least one sentence that
    names the syndrome. The best path of a candidate leaves the seed through one neighbour. This
    function reports, for the head of the ranking, how often that neighbour is one of them.

    Args:
        ranking_path: A ranking TSV written by `scripts/track2_rank.py`.
        audit_path: The audit TSV written by `run_audit_seed`.
        nodes_path: The graph `nodes.tsv`, which gives the label of each node.

    Returns:
        The counts at the top 10, the top 50 and the whole ranking, and the digest of the ranking.
    """
    if not ranking_path.is_file() or not audit_path.is_file():
        return {}
    with nodes_path.open() as handle:
        labels = {row["node_id"]: row["label"] for row in csv.DictReader(handle, delimiter="\t")}
    with audit_path.open() as handle:
        named = {labels.get(row["partner"], row["partner"])
                 for row in csv.DictReader(handle, delimiter="\t")
                 if row["names_syndrome"] == "yes"}
    with ranking_path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    counts = {}
    for top in (10, 50, len(rows)):
        head = rows[:top]
        first_hops = [row["path"].split(" -> ")[1] for row in head if " -> " in row["path"]]
        counts[f"top_{top}"] = {
            "candidates": len(head),
            "first_hop_names_the_syndrome": sum(1 for hop in first_hops if hop in named),
            "first_hop_never_names_the_syndrome": sum(1 for hop in first_hops if hop not in named),
            "first_hops": dict(Counter(first_hops).most_common(6)),
        }
    return {"path": str(ranking_path.relative_to(common.ROOT)),
            "sha256_prefix": common.sha256_of(ranking_path)[:16],
            "neighbours_naming_the_syndrome": sorted(named), "counts": counts}


def compound_ranks(ranking_path: Path, labels: tuple[str, ...]) -> dict:
    """Report the rank and the best path of named compounds in one ranking.

    Args:
        ranking_path: A ranking TSV written by `scripts/track2_rank.py`.
        labels: The compound labels to look up.

    Returns:
        The rank, the path and the relations of each label, with `absent` when the ranking lacks
        it, and the digest of the ranking.
    """
    if not ranking_path.is_file():
        return {}
    with ranking_path.open() as handle:
        rows = {row["label"]: row for row in csv.DictReader(handle, delimiter="\t")}
    found = {}
    for label in labels:
        row = rows.get(label)
        found[label] = ({"rank": int(row["rank"]), "of": len(rows), "path": row["path"],
                         "relations": row["relations"]} if row else "absent")
    return {"path": str(ranking_path.relative_to(common.ROOT)),
            "sha256_prefix": common.sha256_of(ranking_path)[:16], "compounds": found}


def read_judgments(path: Path) -> list[dict]:
    """Read the per-edge judgment table.

    Args:
        path: The judgment TSV.

    Returns:
        One record per edge.
    """
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_receipt(args: argparse.Namespace) -> int:
    """Check the judgment table against the spec and the sentence table, and write the receipt.

    Args:
        args: The parsed command line.

    Returns:
        0 when every edge has one valid judgment, 1 otherwise.
    """
    judgments = read_judgments(args.judgments)
    with args.out.open() as handle:
        sentences = list(csv.DictReader(handle, delimiter="\t"))
    spec_ids = [spec["edge_id"] for spec in EDGE_SPEC]
    judged_ids = [row["edge_id"] for row in judgments]
    problems = []
    if sorted(spec_ids) != sorted(judged_ids):
        problems.append(f"edge ids differ: spec {len(spec_ids)}, judged {len(judged_ids)}")
    for row in judgments:
        if row["classification"] not in CLASSES:
            problems.append(f"{row['edge_id']}: class {row['classification']!r}")
    read_by_edge = Counter(row["edge_id"] for row in sentences if row["in_sample"] == "yes")
    quoted = [row for row in sentences if row["in_sample"] == "yes"]
    for row in judgments:
        if not any(row["quote"] in item["sentence"] and item["edge_id"] == row["edge_id"]
                   and (not row["quote_pmid"] or item["pmid"] == row["quote_pmid"])
                   for item in quoted):
            problems.append(f"{row['edge_id']}: the quote is in no sentence that was read")
    for row in judgments:
        if int(row["sentences_read"]) != read_by_edge.get(row["edge_id"], 0):
            problems.append(f"{row['edge_id']}: read {row['sentences_read']}, table "
                            f"{read_by_edge.get(row['edge_id'], 0)}")
    by_class = Counter(row["classification"] for row in judgments)
    text_mined = [row for row in judgments if row["evidence_kind"] == "pubtator3_sentence"]
    receipt = {
        "receipt_id": "t2_edge_provenance_v1",
        "created_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "track": "t2_rank",
        "control": "agent_reading_of_edge_evidence",
        "queue_item": "R3",
        "data_class_of_run": "A",
        "data_class_of_receipt": "A",
        "receipt_scope": "Public entity accessions, PMIDs and published sentences only. No "
                         "patient data enters this run.",
        "script": "scripts/track2_edge_provenance.py",
        "reader": "Agent reading on the owner's behalf; the owner has not read the sentences",
        "inputs": {
            "edges_tsv": {"path": str(args.edges.relative_to(common.ROOT)),
                          "sha256_prefix": common.sha256_of(args.edges)[:16]},
            "drug_evidence_tsv": {"path": str(args.drug_evidence.relative_to(common.ROOT)),
                                  "sha256_prefix": common.sha256_of(args.drug_evidence)[:16]},
            "sentence_cache": {"path": str(args.cache.relative_to(common.ROOT)),
                               "sha256_prefix": common.sha256_of(args.cache)[:16]},
        },
        "outputs": {
            "sentences_tsv": {"path": str(args.out.relative_to(common.ROOT)),
                              "sha256_prefix": common.sha256_of(args.out)[:16],
                              "rows": len(sentences)},
            "judgments_tsv": {"path": str(args.judgments.relative_to(common.ROOT)),
                              "sha256_prefix": common.sha256_of(args.judgments)[:16]},
        },
        "sampling": {"seed": SAMPLE_SEED, "sample_size": SAMPLE_SIZE,
                     "pages_fetched_for_a_sample": SAMPLE_PAGES,
                     "rule": "An edge in all mode is read in full. An edge in sample mode is read "
                             "on a seeded draw from the first pages in endpoint order."},
        "result": {
            "edges": len(judgments),
            "text_mined_edges": len(text_mined),
            "curated_edges": len(judgments) - len(text_mined),
            "sentences_read": sum(int(row["sentences_read"]) for row in judgments),
            "by_class": {name: by_class.get(name, 0) for name in CLASSES},
            "by_class_text_mined": {name: sum(1 for row in text_mined
                                              if row["classification"] == name)
                                    for name in CLASSES},
            "problems": problems,
            "verdict": "PASS" if not problems else "FAIL",
        },
        "seed_surface_audit": audit_summary(args.audit_out),
        "first_hop_of_the_ranking": first_hop_summary(args.ranking, args.audit_out, args.nodes),
        "nad_precursors_in_the_rankings": {
            "run_f_constrained_with_defect": compound_ranks(args.ranking, NAD_LABELS),
            "run_a_unconstrained": compound_ranks(args.unconstrained_ranking, NAD_LABELS),
        },
    }
    common.write_json(args.receipt, receipt)
    print(json.dumps(receipt["result"], indent=1))
    print(f"receipt: {args.receipt}")
    return 0 if not problems else 1


def main() -> int:
    """Parse the command line and run one step.

    Returns:
        The exit status of the step.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--edges", type=Path, default=EDGES_PATH, help="The graph edges.tsv.")
    parser.add_argument("--drug-evidence", type=Path, default=DRUG_EVIDENCE_PATH,
                        help="The graph drug_evidence.tsv.")
    parser.add_argument("--cache", type=Path, default=CACHE_PATH,
                        help="The cache of the search responses.")
    parser.add_argument("--out", type=Path, default=SENTENCES_PATH, help="The sentence table.")
    parser.add_argument("--audit-out", type=Path, default=AUDIT_PATH,
                        help="The seed surface audit table.")
    steps = parser.add_subparsers(dest="step", required=True)
    extract = steps.add_parser("extract", help="Fetch the sentences and write the sentence table.")
    extract.add_argument("--refresh", action="store_true",
                         help="Fetch every pair again instead of reading the cache.")
    receipt = steps.add_parser("receipt", help="Check the judgments and write the receipt.")
    receipt.add_argument("--judgments", type=Path, default=JUDGMENTS_PATH,
                         help="The per-edge judgment table.")
    receipt.add_argument("--receipt", type=Path, default=RECEIPT_PATH, help="The receipt path.")
    receipt.add_argument("--ranking", type=Path, default=RANKING_PATH,
                         help="The ranking whose first hops are checked against the audit.")
    receipt.add_argument("--nodes", type=Path, default=NODES_PATH, help="The graph nodes.tsv.")
    receipt.add_argument("--unconstrained-ranking", type=Path, default=UNCONSTRAINED_RANKING_PATH,
                         help="The ranking with every constraint off.")
    audit = steps.add_parser(
        "audit-seed", help="Check the seed surface form in every sentence of every seed edge.")
    audit.add_argument("--seed", default=SEED, help="The seed node identifier.")
    audit.add_argument("--entities", type=Path, default=ENTITIES_PATH,
                       help="The resolved entity cache.")
    args = parser.parse_args()
    if args.step == "extract":
        return run_extract(args)
    if args.step == "audit-seed":
        return run_audit_seed(args)
    return run_receipt(args)


if __name__ == "__main__":
    sys.exit(main())
