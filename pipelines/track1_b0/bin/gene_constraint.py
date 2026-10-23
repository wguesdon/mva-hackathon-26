#!/usr/bin/env python3
"""Gene level constraint from the gnomAD v4.1 constraint metrics table. Python standard library only.

The table is `gnomad.v4.1.constraint_metrics.tsv` of the gnomAD v4.1 release. The gnomAD policy page
states that the primary data of the gnomAD exomes and genomes are available free of restrictions
under the Creative Commons Zero Public Domain Dedication, and that no embargo applies to a result
derived from the data. The same page states that an annotation of a third party can carry another
licence, and it names the SpliceAI annotation, which Illumina supplies under CC BY-NC 4.0. This
module reads the constraint table alone, which is primary gnomAD content under CC0.

The file holds one row per transcript. This module keeps one row per gene symbol and prefers the
MANE Select row, then the canonical row, then the first row.

Three values leave this module.

    loeuf   `lof.oe_ci.upper`, the upper bound of the 90 percent confidence interval of the
            observed over expected ratio of predicted loss of function variants. A small value
            means that the gene tolerates few loss of function variants.
    pli     `lof.pLI`, the probability that the gene is intolerant of a single loss of function
            allele.
    mis_z   `mis.z_score`, the z score of the observed against the expected missense count.

The table holds no variant of any benchmark suite and no clinical significance term. It is a gene
level summary of population sequencing, so it carries no label of a test allele.
"""

from __future__ import annotations

import csv
from pathlib import Path

LOEUF_COLUMN = "lof.oe_ci.upper"
PLI_COLUMN = "lof.pLI"
MIS_Z_COLUMN = "mis.z_score"


def parse_number(row: dict, name: str) -> float | None:
    """Read one numeric column of a constraint row.

    Args:
        row: The table row.
        name: The column name.

    Returns:
        The value, or None when the column is empty or is not a number.
    """
    try:
        return float(row.get(name, ""))
    except (TypeError, ValueError):
        return None


def row_rank(row: dict) -> int:
    """Rank a constraint row by the transcript that it describes.

    Args:
        row: The table row.

    Returns:
        Zero for the MANE Select row, one for the canonical row, two otherwise.
    """
    if row.get("mane_select") == "true":
        return 0
    if row.get("canonical") == "true":
        return 1
    return 2


def load(path: Path) -> dict[str, dict]:
    """Read the constraint values of every gene symbol.

    Args:
        path: The gnomAD v4.1 constraint metrics table.

    Returns:
        The gene symbol mapped to its LOEUF, pLI and missense z score.
    """
    genes: dict[str, dict] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = row.get("gene") or ""
            if not gene:
                continue
            rank = row_rank(row)
            held = genes.get(gene)
            if held is not None and held["row_rank"] <= rank:
                continue
            genes[gene] = {"row_rank": rank, "loeuf": parse_number(row, LOEUF_COLUMN),
                           "pli": parse_number(row, PLI_COLUMN),
                           "mis_z": parse_number(row, MIS_Z_COLUMN)}
    return genes


def loeuf_points(gene: str, table: dict[str, dict], ceiling: float = 1.0) -> float:
    """Score one gene by its loss of function constraint.

    The score is `ceiling - LOEUF` clamped to the range 0 to 1, and it is 0.0 for a gene that the
    table does not hold. A gene with a LOEUF at or above the ceiling therefore earns nothing, and a
    fully constrained gene earns close to the ceiling.

    Args:
        gene: The gene symbol.
        table: The output of `load`.
        ceiling: The LOEUF value at which the score reaches zero.

    Returns:
        The constraint score between 0.0 and 1.0.
    """
    value = table.get(gene, {}).get("loeuf")
    if value is None:
        return 0.0
    return max(0.0, min(1.0, ceiling - value))


def pli_points(gene: str, table: dict[str, dict]) -> float:
    """Score one gene by pLI.

    Args:
        gene: The gene symbol.
        table: The output of `load`.

    Returns:
        The pLI value, or 0.0 for a gene that the table does not hold.
    """
    value = table.get(gene, {}).get("pli")
    return 0.0 if value is None else max(0.0, min(1.0, value))


def mis_z_points(gene: str, table: dict[str, dict], ceiling: float = 4.0) -> float:
    """Score one gene by its missense constraint z score.

    Args:
        gene: The gene symbol.
        table: The output of `load`.
        ceiling: The z score at which the score reaches one.

    Returns:
        The score between 0.0 and 1.0.
    """
    value = table.get(gene, {}).get("mis_z")
    if value is None:
        return 0.0
    return max(0.0, min(1.0, value / ceiling))
