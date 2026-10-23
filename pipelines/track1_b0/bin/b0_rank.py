#!/usr/bin/env python3
"""B0 baseline ranker for Track 1 (ADR 0003). Python standard library only.

B0 ranks causal hypotheses from a VEP-annotated VCF. It uses consequence impact, population
frequency, ClinVar significance, zygosity and an optional list of resolved literature alleles. It
uses no model. A phenotype layer is available and it is off by default.

Variant score:
    impact      HIGH 3, MODERATE 2, LOW 1, MODIFIER 0; a literature allele counts as HIGH
    clinvar     +2 when CLIN_SIG holds pathogenic or likely_pathogenic and no benign or conflicting term
    rarity      +1 when absent from gnomAD, +0.5 when the highest gnomAD AF is 0.001 or less
    A variant is a candidate when its highest gnomAD AF is 0.01 or less and its score is 2 or more.

Hypotheses:
    pair        two heterozygous candidates in one gene at distinct positions; score = sum of the two
    homozygous  one homozygous candidate; score = twice its variant score
    phenotype   + PHENOTYPE_WEIGHT * phenotype_score(gene), which is 0 when the layer is off
The top hypotheses by score become rows. The row budget is 10 and the EPCR falls from 0.95 by 0.09
per rank. `--max-rows`, `--epcr-top` and `--epcr-step` change the three values, and their defaults
reproduce every receipt that this ranker wrote before the options existed.

Scoring profile:
    v1          the rule above, which is the rule of every receipt before the option existed
    v2          the inheritance aware frequency ladder of `rarity_points`, which gives the full
                rarity point to every allele at or below RECESSIVE_AF and half a point up to
                MAX_AF, because a recessive causal allele may sit in gnomAD at a low frequency
`--gene-constraint` adds a gene term of CONSTRAINT_WEIGHT times the gnomAD v4.1 pLI of the gene.
The v2 profile and that term are the two changes that the development split selected. Four other
candidate terms failed on the same split and each one stays reachable as an option, which are the
transcript selection rule of `select_v2`, the LOEUF and missense z forms of the constraint term,
and the phenotype layer.

`--profile` is v1 by default and `--gene-constraint` is unset, so a run with neither option
reproduces every receipt that this ranker wrote before the profile existed. The design, the
evidence per feature and the decision gate of v2 are in
`docs/benchmarks/t1_ranker_v2_heldout_v1.md`.

The phenotype layer scores a gene by the Resnik information content similarity of the case profile
against the phenotype profile of the gene, combined by the best match average. The module
`phenotype_layer.py` holds the measure and its sources. `--phenotype-weight` is 0.0 by default, so a
run without the three phenotype options reproduces every receipt that this ranker wrote before the
layer existed. A case profile of the child is class B2 content under ADR 0005.

Usage:
    b0_rank.py --vcf annotated.vcf.gz [--extra-vcf spikes_annotated.vcf.gz --ids id1,id2]
               [--literature known_alleles.tsv] --out submission.csv [--summary summary.json]
               [--hpo-obo hp.obo --genes-to-phenotype genes_to_phenotype.txt
                --case-phenotype terms.txt --phenotype-weight 2.0]
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

IMPACT = {"HIGH": 3.0, "MODERATE": 2.0, "LOW": 1.0, "MODIFIER": 0.0}
# The scoring profile. "v1" is the rule of every receipt that this ranker wrote before the profile
# option existed. "v2" adds the transcript selection rule, the inheritance aware frequency ladder
# and the gnomAD v4.1 gene constraint term of docs/benchmarks/t1_ranker_v2_heldout_v1.md.
PROFILE = "v1"
# The two switches that the profile sets. They are separate so that a control can measure one
# without the other. TRANSCRIPT_RULE picks `select_v1` or `select_v2`, and AF_LADDER picks the v1
# or the v2 branch of `rarity_points`.
TRANSCRIPT_RULE = "v1"
AF_LADDER = "v1"
# The points that the v2 frequency ladder gives. A recessive causal allele may sit in gnomAD at a
# low frequency, so v2 gives the full rarity point to every allele at or below RECESSIVE_AF and half
# a point between RECESSIVE_AF and MAX_AF. v1 gives the full point only to an absent allele.
RECESSIVE_AF = 0.001
# The gene constraint metric of the v2 gene term. The development split selected pLI. LOEUF and the
# missense z score are reachable as options and both failed on that split.
CONSTRAINT_METRIC = "pli"
# The weight of the gene constraint term in the hypothesis score of v2. The development split holds
# the result flat from 0.25 to 1.0 and it falls at 1.5, so the shipped value is the middle of the
# flat interval on a logarithmic scale.
CONSTRAINT_WEIGHT = 0.5
# The LOEUF value at which the LOEUF form of the constraint term reaches zero. It has no effect on
# the shipped pLI form.
CONSTRAINT_CEILING = 1.0
# The points that a pathogenic or likely pathogenic ClinVar term adds. The shipped value is 2.0.
# `scripts/track1_clinvar_ablation.py` sets it to 0.0 to measure how much of the spike-in suite score
# rests on the same ClinVar status that selected the spiked alleles.
CLINVAR_WEIGHT = 2.0
# The number of heterozygous candidates per gene that enter the compound heterozygous pair search.
# The shipped value is 4. `scripts/track1_pair_pruning_control.py` raises it to measure whether the
# cap removes a true pair.
MAX_HET_PER_GENE = 4
MAX_AF = 0.01
MIN_SCORE = 2.0
MAX_ROWS = 10
EPCR_TOP = 0.95
EPCR_STEP = 0.09
EPCR_DIGITS = 2
HEADER = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2", "alt_2",
          "epcr", "finding_type", "notes"]


def open_text(path: Path):
    """Open a plain or gzipped text file.

    Args:
        path: The file path.

    Returns:
        A text handle.
    """
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def csq_format(header_line: str) -> list[str]:
    """Read the CSQ field names from the VEP INFO header line.

    Args:
        header_line: The ##INFO=<ID=CSQ,...> line.

    Returns:
        The field names in order.
    """
    return header_line.split("Format: ")[1].rstrip('">\n').split("|")


def to_float(value: str) -> float | None:
    """Parse a number, or return None for an empty value.

    Args:
        value: The text.

    Returns:
        The number or None.
    """
    try:
        return float(value) if value not in ("", ".") else None
    except ValueError:
        return None


def clinvar_points(clin_sig: str) -> float:
    """Score a VEP CLIN_SIG value.

    Args:
        clin_sig: The CLIN_SIG text, with & separators.

    Returns:
        CLINVAR_WEIGHT for pathogenic evidence without benign or conflicting terms, else 0.0. The
        shipped weight is 2.0.
    """
    terms = clin_sig.lower().replace("&", ",").split(",")
    pathogenic = any(term.strip() in ("pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic") for term in terms)
    negative = any("benign" in term or "conflicting" in term for term in terms)
    return CLINVAR_WEIGHT if pathogenic and not negative else 0.0


def set_profile(profile: str) -> None:
    """Set the scoring switches from a profile name.

    Profile v2 turns on the inheritance aware frequency ladder and leaves the transcript selection
    rule at v1. The transcript rule ran as a separate arm on the development split and it changed
    no rank 1 outcome while it cost 0.91 mean rank points, so v2 does not carry it.
    `--transcript-rule` reaches it for a control.

    Args:
        profile: The profile name, "v1" or "v2".
    """
    globals()["PROFILE"] = profile
    globals()["TRANSCRIPT_RULE"] = "v1"
    globals()["AF_LADDER"] = profile


def entry_impact(values: dict) -> float:
    """Read the impact points of one CSQ entry.

    Args:
        values: The CSQ entry.

    Returns:
        The points of its IMPACT term.
    """
    return IMPACT.get(values.get("IMPACT", "MODIFIER"), 0.0)


def select_v1(entries: list[dict]) -> dict:
    """Choose the CSQ entry of one gene under the v1 rule.

    The v1 rule walks the entries in file order and keeps an entry when it carries the PICK flag or
    when its IMPACT is higher than the IMPACT of the entry that it holds. A later entry therefore
    overrides an earlier PICK entry when its IMPACT is higher.

    Args:
        entries: The CSQ entries of one gene, in file order.

    Returns:
        The chosen entry with its impact points.
    """
    current: dict | None = None
    for values in entries:
        impact = entry_impact(values)
        if current is None or values.get("PICK") == "1" or impact > current["impact"]:
            current = {"impact": impact, "values": values}
    return current or {"impact": 0.0, "values": {}}


def select_v2(entries: list[dict]) -> dict:
    """Choose the CSQ entry of one gene under the v2 transcript selection rule.

    The rule is stated and it does not depend on the order of the entries. The first entry that
    names a MANE Select transcript wins. When no entry names one, the first canonical entry wins.
    When no entry is canonical, the first entry that carries the PICK flag wins. When no entry
    carries the flag, the entry with the highest IMPACT wins, and a tie keeps the first entry.

    Args:
        entries: The CSQ entries of one gene, in file order.

    Returns:
        The chosen entry with its impact points.
    """
    chosen = next((values for values in entries if values.get("MANE_SELECT")), None)
    if chosen is None:
        chosen = next((values for values in entries if values.get("CANONICAL") == "YES"), None)
    if chosen is None:
        chosen = next((values for values in entries if values.get("PICK") == "1"), None)
    if chosen is None and entries:
        chosen = max(entries, key=entry_impact)
    if chosen is None:
        return {"impact": 0.0, "values": {}}
    return {"impact": entry_impact(chosen), "values": chosen}


def rarity_points(max_af: float | None, profile: str = "v1") -> float:
    """Score the population frequency of one allele.

    Args:
        max_af: The highest gnomAD allele frequency, or None when gnomAD holds no frequency.
        profile: The scoring profile, "v1" or "v2".

    Returns:
        The rarity points. v1 gives 1.0 to an absent allele, 0.5 at or below 0.001 and 0.0 above.
        v2 gives 1.0 at or below RECESSIVE_AF, which covers an absent allele, 0.5 between
        RECESSIVE_AF and MAX_AF, and 0.0 above MAX_AF.
    """
    if profile == "v2":
        if max_af is None or max_af <= RECESSIVE_AF:
            return 1.0
        return 0.5 if max_af <= MAX_AF else 0.0
    if max_af is None or max_af == 0:
        return 1.0
    return 0.5 if max_af <= 0.001 else 0.0


def constraint_points(gene: str, table: dict[str, dict], metric: str, ceiling: float) -> float:
    """Score one gene by one gnomAD constraint metric.

    Args:
        gene: The gene symbol.
        table: The output of `gene_constraint.load`.
        metric: The metric name, "pli", "loeuf" or "mis_z".
        ceiling: The LOEUF value at which the LOEUF form reaches zero.

    Returns:
        The score between 0.0 and 1.0.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import gene_constraint

    if metric == "loeuf":
        return gene_constraint.loeuf_points(gene, table, ceiling)
    if metric == "mis_z":
        return gene_constraint.mis_z_points(gene, table)
    return gene_constraint.pli_points(gene, table)


def constraint_bonus(records: list[dict], table: dict[str, dict], weight: float = CONSTRAINT_WEIGHT,
                     ceiling: float = CONSTRAINT_CEILING,
                     metric: str = CONSTRAINT_METRIC) -> tuple[dict[str, float], dict]:
    """Score every gene of the candidate records by its gnomAD constraint.

    Args:
        records: Candidate variant records.
        table: The output of `gene_constraint.load`.
        weight: The weight of the constraint score in the hypothesis score.
        ceiling: The LOEUF value at which the LOEUF form of the score reaches zero.
        metric: The metric name, "pli", "loeuf" or "mis_z".

    Returns:
        The bonus of each gene, already multiplied by the weight, and a record of the layer.
    """
    column = {"pli": "lof.pLI", "loeuf": "lof.oe_ci.upper", "mis_z": "mis.z_score"}[metric]
    genes = sorted({record["gene"] for record in records})
    points = {gene: constraint_points(gene, table, metric, ceiling) for gene in genes}
    report = {"source": f"gnomAD v4.1 constraint metrics, column {column}",
              "licence": "CC0, the gnomAD policy page",
              "metric": metric, "weight": weight, "ceiling": ceiling, "genes_scored": len(genes),
              "genes_with_the_metric": sum(1 for gene in genes if gene in table
                                           and table[gene].get(metric) is not None)}
    return {gene: weight * value for gene, value in points.items()}, report


def merge_bonus(first: dict[str, float] | None, second: dict[str, float] | None) -> dict[str, float] | None:
    """Add two gene bonus maps.

    Args:
        first: The first map, or None.
        second: The second map, or None.

    Returns:
        The sum of the two maps, or None when both are None.
    """
    if first is None and second is None:
        return None
    merged = dict(first or {})
    for gene, value in (second or {}).items():
        merged[gene] = merged.get(gene, 0.0) + value
    return merged


def load_literature(path: Path | None) -> dict[tuple[str, int, str, str], str]:
    """Read resolved literature alleles as scorer tuples with their gene.

    Args:
        path: The S3 allele table, or None.

    Returns:
        (chrom, pos, ref, alt) with a chr prefix, mapped to the gene symbol.
    """
    if path is None:
        return {}
    keys = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["status"] == "resolved" and row["pos"]:
                key = ("chr" + row["chrom"].removeprefix("chr"), int(row["pos"]), row["ref"].upper(), row["alt"].upper())
                keys[key] = row["gene"]
    return keys


def parse_vcf(path: Path, literature: dict, keep_ids: set[str] | None = None) -> list[dict]:
    """Read candidate variants from a VEP-annotated single-sample VCF.

    Args:
        path: The annotated VCF.
        literature: Literature allele keys.
        keep_ids: When set, keep only records whose ID is in this set.

    Returns:
        One record per variant allele and gene with its score.
    """
    fields: list[str] = []
    records = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                fields = csq_format(line)
                continue
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt, _, filt, info = columns[:8]
            if keep_ids is not None and vid not in keep_ids:
                continue
            if filt not in ("PASS", ".") or "," in alt:
                continue
            sample = dict(zip(columns[8].split(":"), columns[9].split(":"))) if len(columns) > 9 else {}
            gt = sample.get("GT", "").replace("|", "/")
            if gt not in ("0/1", "1/0", "1/1"):
                continue
            gq = to_float(sample.get("GQ", ""))
            if gq is not None and gq < 20:
                continue
            key = (chrom if chrom.startswith("chr") else "chr" + chrom, int(pos), ref.upper(), alt.upper())
            csq = next((item[4:] for item in info.split(";") if item.startswith("CSQ=")), "")
            entries_by_gene: dict[str, list[dict]] = defaultdict(list)
            for entry in csq.split(","):
                values = dict(zip(fields, entry.split("|")))
                gene = values.get("SYMBOL") or values.get("Gene") or ""
                if gene:
                    entries_by_gene[gene].append(values)
            select = select_v2 if TRANSCRIPT_RULE == "v2" else select_v1
            by_gene: dict[str, dict] = {gene: select(entries)
                                        for gene, entries in entries_by_gene.items()}
            # A literature allele joins its published gene, also when VEP names no transcript of that gene.
            # The BUB1B upstream allele lies 44 kb before the gene, outside the VEP upstream window.
            literature_gene = literature.get(key)
            if literature_gene and literature_gene not in by_gene:
                template = next(iter(by_gene.values()))["values"] if by_gene else {}
                by_gene[literature_gene] = {"impact": 0.0, "values": template}
            for gene, chosen in by_gene.items():
                values = chosen["values"]
                frequencies = [to_float(values.get(name, "")) for name in ("gnomADg_AF", "gnomADe_AF")]
                present = [value for value in frequencies if value is not None]
                max_af = max(present) if present else None
                if max_af is not None and max_af > MAX_AF:
                    continue
                impact = 3.0 if key in literature else chosen["impact"]
                rarity = rarity_points(max_af, AF_LADDER)
                score = impact + clinvar_points(values.get("CLIN_SIG", "")) + rarity
                if score < MIN_SCORE:
                    continue
                records.append({"key": key, "gene": gene, "gt": gt, "score": score, "impact": impact,
                                "max_af": max_af, "literature": key in literature})
    return records


def hypotheses(records: list[dict], gene_bonus: dict[str, float] | None = None) -> list[dict]:
    """Build pair and homozygous hypotheses per gene.

    Args:
        records: Candidate variant records.
        gene_bonus: The phenotype bonus of each gene, already multiplied by the weight. The bonus is
            0.0 for a gene that the map does not hold, and the map is None when the layer is off.

    Returns:
        Hypotheses sorted by score, then by genomic position.
    """
    bonus = gene_bonus or {}
    by_gene: dict[str, dict] = defaultdict(dict)
    for record in records:
        best = by_gene[record["gene"]].get(record["key"])
        if best is None or record["score"] > best["score"]:
            by_gene[record["gene"]][record["key"]] = record
    found = []
    for gene, variants in by_gene.items():
        gene_points = bonus.get(gene, 0.0)
        ordered = sorted(variants.values(), key=variant_order)
        for record in ordered:
            if record["gt"] == "1/1":
                found.append({"gene": gene, "variants": [record["key"]],
                              "score": 2 * record["score"] + gene_points,
                              "variant_score": 2 * record["score"], "phenotype_points": gene_points})
        hets = [record for record in ordered if record["gt"] != "1/1"][:MAX_HET_PER_GENE]
        for first_index, first in enumerate(hets):
            for second in hets[first_index + 1:]:
                if first["key"][:2] != second["key"][:2]:
                    found.append({"gene": gene, "variants": sorted([first["key"], second["key"]], key=position_order),
                                  "score": first["score"] + second["score"] + gene_points,
                                  "variant_score": first["score"] + second["score"],
                                  "phenotype_points": gene_points})
    return sorted(found, key=hypothesis_order)


def variant_order(record: dict) -> tuple:
    """Sort variants by descending score, then by position.

    Args:
        record: A variant record.

    Returns:
        The sort key.
    """
    return -record["score"], position_order(record["key"])


def position_order(key: tuple) -> tuple:
    """Sort scorer tuples by chromosome number and position.

    Args:
        key: (chrom, pos, ref, alt).

    Returns:
        The sort key.
    """
    name = key[0].removeprefix("chr")
    rank = int(name) if name.isdigit() else {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(name, 26)
    return rank, key[1], key[2], key[3]


def hypothesis_order(item: dict) -> tuple:
    """Sort hypotheses by descending score, then by the position of the first variant.

    Args:
        item: A hypothesis.

    Returns:
        The sort key.
    """
    return -item["score"], position_order(item["variants"][0])


def rows_for(found: list[dict], proband: str, max_rows: int = MAX_ROWS,
             epcr_top: float = EPCR_TOP, epcr_step: float = EPCR_STEP,
             epcr_digits: int = EPCR_DIGITS) -> list[dict]:
    """Convert the top hypotheses to submission rows without reusing a pair.

    Args:
        found: Sorted hypotheses.
        proband: The proband identifier.
        max_rows: The row budget of this writer.
        epcr_top: The EPCR of the first row.
        epcr_step: The fall of the EPCR per rank.
        epcr_digits: The number of decimal digits of the EPCR.

    Returns:
        At most `max_rows` rows with distinct EPCR values.
    """
    rows, seen = [], set()
    for item in found:
        signature = tuple(item["variants"])
        if signature in seen:
            continue
        seen.add(signature)
        rank = len(rows) + 1
        first = item["variants"][0]
        second = item["variants"][1] if len(item["variants"]) > 1 else ("", "", "", "")
        value = round(epcr_top - epcr_step * (rank - 1), epcr_digits)
        rows.append({"proband_id": proband, "chrom_1": first[0], "pos_1": first[1], "ref_1": first[2], "alt_1": first[3],
                     "chrom_2": second[0], "pos_2": second[1], "ref_2": second[2], "alt_2": second[3],
                     "epcr": f"{value:.{epcr_digits}f}", "finding_type": "primary",
                     "notes": f"B0 {item['gene']} score {item['score']:.1f}"})
        if len(rows) == max_rows:
            break
    return rows


def phenotype_bonus(records: list[dict], obo: Path, genes_to_phenotype: Path, case_profile: Path,
                    weight: float) -> tuple[dict[str, float], dict]:
    """Score every gene of the candidate records against a case phenotype profile.

    Args:
        records: Candidate variant records.
        obo: The `hp.obo` file of one HPO release.
        genes_to_phenotype: The `genes_to_phenotype.txt` file of the same release.
        case_profile: The file that holds one HPO term of the case per line.
        weight: The weight of the phenotype score in the hypothesis score.

    Returns:
        The bonus of each gene, already multiplied by the weight, and a record of the layer that
        holds no term of the case.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import phenotype_layer

    model = phenotype_layer.build_model(obo, genes_to_phenotype)
    terms = model.normalize(phenotype_layer.read_case_profile(case_profile))
    genes = sorted({record["gene"] for record in records})
    scores = {gene: model.score_gene(gene, terms) for gene in genes}
    report = {"hpo_obo_sha256": phenotype_layer.sha256_of(obo),
              "genes_to_phenotype_sha256": phenotype_layer.sha256_of(genes_to_phenotype),
              "case_profile_sha256": phenotype_layer.sha256_of(case_profile),
              "case_terms_used": len(terms), "genes_scored": len(genes),
              "genes_with_an_annotation": sum(1 for gene in genes if gene in model.profiles),
              "weight": weight, "measure": "Resnik information content, best match average",
              "max_information_content_nats": round(model.max_content, 4)}
    return {gene: weight * score for gene, score in scores.items()}, report


def normalization_report(rows: list[dict], args) -> dict:
    """Normalize the rows of a submission against GRCh38 and record every change.

    The check does not decide anything. The writer emits the normalized form and the record names
    every variant that differed, so the owner reads the record and decides what goes out.

    Args:
        rows: The submission rows.
        args: The parsed command line, which holds the reference and the backend options.

    Returns:
        The report of submission_norm.check_submission.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import submission_norm

    reference = submission_norm.Reference(args.reference)
    bcftools = args.bcftools if args.norm_backend == "bcftools" else None
    report = submission_norm.check_submission(rows, reference, bcftools, args.norm_container_image)
    reference.close()
    print(f"Normalization check: {report['variants']} variants, {report['changed']} changed, "
          f"{report['with_a_second_representation']} carry a second representation, "
          f"backend {report['backend']}, {report['bcftools']}", file=sys.stderr)
    return report


def main() -> int:
    """Rank one sample and write a submission CSV.

    Returns:
        Zero on success.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--extra-vcf", type=Path)
    parser.add_argument("--ids", default="")
    parser.add_argument("--literature", type=Path)
    parser.add_argument("--proband", default="PROBAND01")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--reference", type=Path,
                        help="An indexed GRCh38 FASTA. It turns on the normalization check.")
    parser.add_argument("--norm-report", type=Path,
                        help="The record of every row that changed under normalization.")
    parser.add_argument("--norm-backend", choices=("bcftools", "internal"), default="bcftools")
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--norm-container-image")
    parser.add_argument("--hpo-obo", type=Path, help="The hp.obo file of one HPO release.")
    parser.add_argument("--genes-to-phenotype", type=Path,
                        help="The genes_to_phenotype.txt file of the same HPO release.")
    parser.add_argument("--case-phenotype", type=Path,
                        help="One HPO term of the case per line. Class B2 for the child.")
    parser.add_argument("--profile", choices=("v1", "v2"), default=PROFILE,
                        help="The scoring profile. v1 reproduces every receipt before the option.")
    parser.add_argument("--gene-constraint", type=Path,
                        help="The gnomAD v4.1 constraint metrics table. It turns on the v2 "
                             "constraint term.")
    parser.add_argument("--constraint-weight", type=float, default=CONSTRAINT_WEIGHT,
                        help="The weight of the gene constraint score. 0.0 removes the feature.")
    parser.add_argument("--constraint-ceiling", type=float, default=CONSTRAINT_CEILING,
                        help="The LOEUF value at which the LOEUF constraint score reaches zero.")
    parser.add_argument("--constraint-metric", choices=("pli", "loeuf", "mis_z"),
                        default=CONSTRAINT_METRIC,
                        help="The gnomAD constraint metric of the gene term.")
    parser.add_argument("--transcript-rule", choices=("v1", "v2"),
                        help="Override the transcript selection rule that the profile sets.")
    parser.add_argument("--constraint-report", type=Path,
                        help="Where to write the record of the gene constraint layer.")
    parser.add_argument("--clinvar-weight", type=float, default=CLINVAR_WEIGHT,
                        help="Points for a pathogenic ClinVar term. 0.0 removes the feature.")
    parser.add_argument("--max-het-per-gene", type=int, default=MAX_HET_PER_GENE,
                        help="Heterozygous candidates per gene in the pair search.")
    parser.add_argument("--phenotype-weight", type=float, default=0.0,
                        help="The weight of the phenotype score. 0.0 turns the layer off.")
    parser.add_argument("--phenotype-report", type=Path,
                        help="Where to write the record of the phenotype layer.")
    parser.add_argument("--max-rows", type=int, default=MAX_ROWS,
                        help="The row budget of this writer.")
    parser.add_argument("--epcr-top", type=float, default=EPCR_TOP,
                        help="The EPCR of the first row.")
    parser.add_argument("--epcr-step", type=float, default=EPCR_STEP,
                        help="The fall of the EPCR per rank.")
    parser.add_argument("--epcr-digits", type=int, default=EPCR_DIGITS,
                        help="The number of decimal digits of the EPCR.")
    args = parser.parse_args()
    literature = load_literature(args.literature)
    globals()["CLINVAR_WEIGHT"] = args.clinvar_weight
    globals()["MAX_HET_PER_GENE"] = args.max_het_per_gene
    set_profile(args.profile)
    if args.transcript_rule:
        globals()["TRANSCRIPT_RULE"] = args.transcript_rule
    records = parse_vcf(args.vcf, literature)
    if args.extra_vcf:
        records += parse_vcf(args.extra_vcf, literature, set(filter(None, args.ids.split(","))))
    constraint_map, constraint_record = None, None
    if args.gene_constraint and args.constraint_weight:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import gene_constraint
        table = gene_constraint.load(args.gene_constraint)
        constraint_map, constraint_record = constraint_bonus(
            records, table, args.constraint_weight, args.constraint_ceiling,
            args.constraint_metric)
        if args.constraint_report:
            args.constraint_report.write_text(json.dumps(constraint_record, indent=2) + "\n")
    bonus, phenotype_record = None, None
    if args.phenotype_weight and args.case_phenotype:
        if not args.hpo_obo or not args.genes_to_phenotype:
            parser.error("--phenotype-weight needs --hpo-obo and --genes-to-phenotype.")
        bonus, phenotype_record = phenotype_bonus(records, args.hpo_obo, args.genes_to_phenotype,
                                                  args.case_phenotype, args.phenotype_weight)
        if args.phenotype_report:
            args.phenotype_report.write_text(json.dumps(phenotype_record, indent=2) + "\n")
    rows = rows_for(hypotheses(records, merge_bonus(bonus, constraint_map)), args.proband,
                    args.max_rows, args.epcr_top, args.epcr_step, args.epcr_digits)
    report = None
    if args.reference:
        report = normalization_report(rows, args)
        rows = report.pop("normalized_rows")
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER)
        writer.writeheader()
        writer.writerows(rows)
    if args.norm_report and report is not None:
        args.norm_report.write_text(json.dumps(report, indent=2) + "\n")
    if args.summary:
        summary = {"candidates": len(records), "rows": len(rows), "profile": args.profile,
                   "literature_candidates": sum(r["literature"] for r in records)}
        if phenotype_record is not None:
            summary["phenotype"] = phenotype_record
        if constraint_record is not None:
            summary["gene_constraint"] = constraint_record
        if report is not None:
            summary["normalization_changed_variants"] = report["changed"]
            summary["normalization_backend"] = report["backend"]
        args.summary.write_text(json.dumps(summary) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
