"""Write a synthetic Track 2 graph that follows the input contract of the real graph.

The real graph comes from PubTator3, Open Targets, Reactome and ChEMBL. This module writes a small
graph of the same shape, so the ranker and the two controls run and are tested before a change
reaches the real graph. No identifier and no PMID of this graph is real. Every gene, disease,
chemical and PMID carries a `SYN` prefix.

The fixture plants the structures that the ranker must separate.

1. One reagent. `CHEMICAL:SYNR01` carries the most PMIDs of the graph on a co-mention edge to the
   defect gene, and it has no drug record. It wins a ranking that reads PMIDs alone, and both the
   terminal constraint and the eligibility constraint remove it.
2. One restorer. `CHEMBL:SYNC01` is an agonist of the defect gene. For a loss of function defect it
   pushes against the defect.
3. One blocker. `CHEMBL:SYNC02` is an inhibitor of the defect gene. It reinforces a loss of
   function defect and must sink to the bottom of the table with its reason.
4. One pathway blocker. `CHEMBL:SYNC03` inhibits a gene that shares a pathway with the defect gene.
5. One unrelated drug. `CHEMBL:SYNC04` inhibits a gene that shares nothing with the defect gene, so
   its direction does not compose.
6. One second disease. `DISEASE:SYNDIS20` co-mentions the defect gene and carries its own drugs. A
   path that hops through it is a co-mention of two diseases and the ranker forbids it.
7. Three degree hubs and nine held-out cases, for the two controls.

Usage:
    uv run python scripts/track2_make_graph_fixture.py --out-dir reports/track2/fixture_graph
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_DIR = ROOT / "reports" / "track2" / "fixture_graph"

NODE_COLUMNS = ("node_id", "kind", "label", "source")
EDGE_COLUMNS = ("source_id", "target_id", "relation", "evidence_count", "pmids", "source")
DRUG_COLUMNS = ("drug_id", "drug_name", "drug_type", "target_gene_id", "action_type",
                "max_clinical_stage_on_target", "max_clinical_stage_overall", "mechanism_of_action")

GENES = 40
CHEMICALS = 30
PATHWAYS = 6
HUBS = ("SYNH1", "SYNH2", "SYNH3")
HUB_DEGREES = (18, 20, 17)
HELD_OUT_CASES = 9
FIRST_HELD_OUT_DRUG = 11

SEED_DISEASE = "DISEASE:SYNDIS00"
OTHER_DISEASE = "DISEASE:SYNDIS20"
DEFECT_GENE = "GENE:SYNG01"
PATHWAY_GENE = "GENE:SYNG02"
FAR_GENE = "GENE:SYNG03"
SHARED_PATHWAY = "PATHWAY:SYNP1"
REAGENT = "CHEMICAL:SYNR01"
RESTORER = "CHEMBL:SYNC01"
BLOCKER = "CHEMBL:SYNC02"
PATHWAY_BLOCKER = "CHEMBL:SYNC03"
UNRELATED = "CHEMBL:SYNC04"

CO_MENTION_RELATIONS = ("associate", "positive_correlate", "negative_correlate", "interact")


def gene_id(number: int) -> str:
    """Return the node identifier of one synthetic gene.

    Args:
        number: The gene number, starting at 1.

    Returns:
        The identifier, for example `GENE:SYNG01`.
    """
    return f"GENE:SYNG{number:02d}"


def chemical_id(number: int) -> str:
    """Return the node identifier of one synthetic drug.

    Args:
        number: The chemical number, starting at 1.

    Returns:
        The identifier, for example `CHEMBL:SYNC01`.
    """
    return f"CHEMBL:SYNC{number:02d}"


def disease_id(number: int) -> str:
    """Return the node identifier of one synthetic disease.

    Args:
        number: The disease number, starting at 0.

    Returns:
        The identifier, for example `DISEASE:SYNDIS01`.
    """
    return f"DISEASE:SYNDIS{number:02d}"


def pmid_list(generator: random.Random, count: int) -> str:
    """Return a comma separated list of synthetic PMIDs.

    Args:
        generator: The random generator.
        count: The number of PMIDs.

    Returns:
        The field value of the edge table. An empty string means no publication.
    """
    return ",".join(sorted(f"SYN{generator.randrange(100000, 999999)}" for _ in range(count)))


def build_nodes() -> list[dict[str, str]]:
    """Build the node table of the fixture.

    Returns:
        One row per node.
    """
    rows = [{"node_id": gene_id(number), "kind": "gene", "label": f"syn_gene_{number:02d}",
             "source": "synthetic_fixture"} for number in range(1, GENES + 1)]
    rows += [{"node_id": chemical_id(number), "kind": "chemical",
              "label": f"syn_drug_{number:02d}", "source": "synthetic_fixture"}
             for number in range(1, CHEMICALS + 1)]
    rows += [{"node_id": REAGENT, "kind": "chemical", "label": "syn_reagent",
              "source": "synthetic_fixture"}]
    rows += [{"node_id": f"CHEMBL:{name}", "kind": "chemical",
              "label": f"syn_drug_hub_{position + 1}", "source": "synthetic_fixture"}
             for position, name in enumerate(HUBS)]
    rows += [{"node_id": disease_id(number), "kind": "disease",
              "label": f"syn_disease_{number:02d}", "source": "synthetic_fixture"}
             for number in list(range(0, HELD_OUT_CASES + 3)) + [20]]
    rows += [{"node_id": f"PATHWAY:SYNP{number}", "kind": "pathway",
              "label": f"syn_pathway_{number}", "source": "synthetic_fixture"}
             for number in range(1, PATHWAYS + 1)]
    return rows


def edge_row(source_id: str, target_id: str, relation: str, evidence: int, pmids: str,
             source: str = "synthetic_fixture") -> dict:
    """Build one edge row.

    Args:
        source_id: The first endpoint.
        target_id: The second endpoint.
        relation: The relation name.
        evidence: The evidence count.
        pmids: The comma separated PMIDs.
        source: The database name.

    Returns:
        The row of the edge table.
    """
    return {"source_id": source_id, "target_id": target_id, "relation": relation,
            "evidence_count": str(evidence), "pmids": pmids, "source": source}


def drug_row(drug_id: str, gene: str, action: str, stage: str, mechanism: str) -> dict:
    """Build one drug evidence row.

    Args:
        drug_id: The chemical node identifier.
        gene: The target gene node identifier.
        action: The Open Targets action type.
        stage: The highest clinical stage.
        mechanism: The mechanism text.

    Returns:
        The row of the drug evidence table.
    """
    return {"drug_id": drug_id, "drug_name": drug_id.split(":")[-1], "drug_type": "Small molecule",
            "target_gene_id": gene, "action_type": action, "max_clinical_stage_on_target": stage,
            "max_clinical_stage_overall": stage, "mechanism_of_action": mechanism}


def seed_case(generator: random.Random) -> tuple[list[dict], list[dict]]:
    """Build the seed disease, the reagent and the four drugs of the direction test.

    Args:
        generator: The random generator.

    Returns:
        The edges and the drug evidence rows of the seed case.
    """
    edges = [
        edge_row(SEED_DISEASE, DEFECT_GENE, "associate", 27, pmid_list(generator, 5), "pubtator3"),
        edge_row(SEED_DISEASE, PATHWAY_GENE, "associate", 12, pmid_list(generator, 3), "pubtator3"),
        edge_row(SEED_DISEASE, FAR_GENE, "associate", 6, pmid_list(generator, 2), "pubtator3"),
        edge_row(REAGENT, DEFECT_GENE, "associate", 90, pmid_list(generator, 40), "pubtator3"),
        edge_row(RESTORER, DEFECT_GENE, "agonist", 30, pmid_list(generator, 6), "open_targets"),
        edge_row(BLOCKER, DEFECT_GENE, "inhibitor", 30, pmid_list(generator, 6), "open_targets"),
        edge_row(PATHWAY_BLOCKER, PATHWAY_GENE, "inhibitor", 25, pmid_list(generator, 5),
                 "open_targets"),
        edge_row(UNRELATED, FAR_GENE, "inhibitor", 25, pmid_list(generator, 5), "open_targets"),
        edge_row(DEFECT_GENE, SHARED_PATHWAY, "participates_in", 1, "", "open_targets_reactome"),
        edge_row(PATHWAY_GENE, SHARED_PATHWAY, "participates_in", 1, "", "open_targets_reactome"),
        edge_row(OTHER_DISEASE, DEFECT_GENE, "associate", 40, pmid_list(generator, 8), "pubtator3"),
        edge_row(OTHER_DISEASE, chemical_id(5), "indication", 200, "", "chembl"),
        edge_row(OTHER_DISEASE, chemical_id(6), "indication", 200, "", "chembl"),
    ]
    drugs = [
        drug_row(RESTORER, DEFECT_GENE, "agonist", "PHASE_2", "syn_gene_01 agonist"),
        drug_row(BLOCKER, DEFECT_GENE, "inhibitor", "PHASE_3", "syn_gene_01 inhibitor"),
        drug_row(PATHWAY_BLOCKER, PATHWAY_GENE, "inhibitor", "PHASE_2", "syn_gene_02 inhibitor"),
        drug_row(UNRELATED, FAR_GENE, "inhibitor", "APPROVAL", "syn_gene_03 inhibitor"),
        drug_row(chemical_id(5), gene_id(30), "inhibitor", "APPROVAL", "syn_gene_30 inhibitor"),
        drug_row(chemical_id(6), gene_id(31), "inhibitor", "APPROVAL", "syn_gene_31 inhibitor"),
    ]
    return edges, drugs


def held_out_cases(generator: random.Random) -> tuple[list[dict], list[dict]]:
    """Build the diseases that have an approved drug.

    Each disease keeps one gene that the drug acts on, so the drug stays reachable at two edges
    after the direct edge is removed.

    Args:
        generator: The random generator.

    Returns:
        The edges and the drug evidence rows of the held-out cases.
    """
    edges = []
    drugs = []
    for case in range(1, HELD_OUT_CASES + 1):
        disease = disease_id(case)
        drug = chemical_id(case + FIRST_HELD_OUT_DRUG)
        gene = gene_id(case + 4)
        action = "agonist" if case % 3 == 0 else "inhibitor"
        edges.append(edge_row(disease, drug, "indication", generator.randrange(3, 15), "", "chembl"))
        edges.append(edge_row(disease, gene, "associate", generator.randrange(5, 25),
                              pmid_list(generator, generator.randrange(2, 6)), "pubtator3"))
        edges.append(edge_row(drug, gene, action, generator.randrange(10, 45),
                              pmid_list(generator, generator.randrange(2, 7)), "open_targets"))
        drugs.append(drug_row(drug, gene, action, "APPROVAL", f"{gene} {action}"))
        for other in generator.sample(range(1, GENES + 1), 2):
            edges.append(edge_row(disease, gene_id(other), "associate", generator.randrange(1, 6),
                                  pmid_list(generator, generator.randrange(0, 2)), "pubtator3"))
    for case in range(HELD_OUT_CASES + 1, HELD_OUT_CASES + 3):
        for gene in generator.sample(range(1, GENES + 1), 3):
            edges.append(edge_row(disease_id(case), gene_id(gene), "associate",
                                  generator.randrange(1, 10),
                                  pmid_list(generator, generator.randrange(0, 3)), "pubtator3"))
    return edges, drugs


def gene_background(generator: random.Random) -> list[dict]:
    """Build the gene to gene and gene to pathway background.

    Args:
        generator: The random generator.

    Returns:
        The background edges.
    """
    rows = []
    for _ in range(70):
        first, second = generator.sample(range(1, GENES + 1), 2)
        rows.append(edge_row(gene_id(first), gene_id(second),
                             generator.choice(CO_MENTION_RELATIONS), generator.randrange(1, 9),
                             pmid_list(generator, generator.randrange(0, 4)), "pubtator3"))
    for number in range(4, GENES + 1):
        for pathway in generator.sample(range(2, PATHWAYS + 1), generator.randrange(1, 3)):
            rows.append(edge_row(gene_id(number), f"PATHWAY:SYNP{pathway}", "participates_in", 1, "",
                                 "open_targets_reactome"))
    return rows


def chemical_background(generator: random.Random) -> tuple[list[dict], list[dict]]:
    """Build the chemical to gene background and the three degree hubs.

    Args:
        generator: The random generator.

    Returns:
        The chemical edges and their drug evidence rows.
    """
    edges = []
    drugs = []
    for number in range(7, CHEMICALS + 1):
        if number in range(FIRST_HELD_OUT_DRUG + 1, FIRST_HELD_OUT_DRUG + HELD_OUT_CASES + 1):
            continue
        for gene in generator.sample(range(4, GENES + 1), generator.randrange(1, 4)):
            action = generator.choice(("inhibitor", "agonist", "modulator"))
            edges.append(edge_row(chemical_id(number), gene_id(gene), action,
                                  generator.randrange(1, 40),
                                  pmid_list(generator, generator.randrange(0, 5)), "open_targets"))
            drugs.append(drug_row(chemical_id(number), gene_id(gene), action,
                                  generator.choice(("PHASE_1", "PHASE_2", "APPROVAL")),
                                  f"{gene_id(gene)} {action}"))
        edges.append(edge_row(chemical_id(number), gene_id(generator.randrange(4, GENES + 1)),
                              generator.choice(CO_MENTION_RELATIONS), generator.randrange(1, 30),
                              pmid_list(generator, generator.randrange(0, 6)), "pubtator3"))
    for name, degree in zip(HUBS, HUB_DEGREES):
        for gene in generator.sample(range(1, GENES + 1), degree):
            edges.append(edge_row(f"CHEMBL:{name}", gene_id(gene), "associate", 1, "", "pubtator3"))
    return edges, drugs


def build_fixture(random_seed: int = 20260915) -> tuple[list[dict], list[dict], list[dict]]:
    """Build the whole fixture.

    Args:
        random_seed: The seed of the random generator, so the fixture repeats.

    Returns:
        The node rows, the edge rows and the drug evidence rows.
    """
    generator = random.Random(random_seed)
    seed_edges, seed_drugs = seed_case(generator)
    case_edges, case_drugs = held_out_cases(generator)
    chemical_edges, chemical_drugs = chemical_background(generator)
    edges = seed_edges + case_edges + gene_background(generator) + chemical_edges
    return build_nodes(), edges, seed_drugs + case_drugs + chemical_drugs


def write_table(path: Path, columns: tuple[str, ...], rows: list[dict]) -> None:
    """Write one tab separated table.

    Args:
        path: The output path.
        columns: The column names.
        rows: The rows.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_fixture(out_dir: Path, random_seed: int = 20260915) -> tuple[Path, Path, Path]:
    """Write the fixture tables.

    Args:
        out_dir: The output directory.
        random_seed: The seed of the random generator.

    Returns:
        The nodes path, the edges path and the drug evidence path.
    """
    nodes, edges, drugs = build_fixture(random_seed)
    nodes_path = out_dir / "nodes.tsv"
    edges_path = out_dir / "edges.tsv"
    drugs_path = out_dir / "drug_evidence.tsv"
    write_table(nodes_path, NODE_COLUMNS, nodes)
    write_table(edges_path, EDGE_COLUMNS, edges)
    write_table(drugs_path, DRUG_COLUMNS, drugs)
    return nodes_path, edges_path, drugs_path


def main(argv: list[str] | None = None) -> int:
    """Write the fixture from the command line.

    Args:
        argv: The arguments. `sys.argv[1:]` is used when it is absent.

    Returns:
        0 on success.
    """
    parser = argparse.ArgumentParser(description="Write a synthetic Track 2 graph fixture.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                        help="The output directory.")
    parser.add_argument("--random-seed", type=int, default=20260915,
                        help="The seed of the random generator.")
    args = parser.parse_args(argv)
    nodes_path, edges_path, drugs_path = write_fixture(args.out_dir, args.random_seed)
    nodes, edges, drugs = build_fixture(args.random_seed)
    print(f"nodes: {len(nodes)} -> {nodes_path}")
    print(f"edges: {len(edges)} -> {edges_path}")
    print(f"drug records: {len(drugs)} -> {drugs_path}")
    print(f"seed disease: {SEED_DISEASE}, defect gene: {DEFECT_GENE}, reagent: {REAGENT}, "
          f"restorer: {RESTORER}, blocker: {BLOCKER}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
