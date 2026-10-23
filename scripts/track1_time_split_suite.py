#!/usr/bin/env python3
"""Build the held-out time split spike-in suite for Track 1. Data class A.

Defect R2-1 of the second jury round states that the suite `spikein_v1` draws 48 of its 50 scenarios
from ClinVar records that carry the same pathogenic significance term that the ranker reads as a
2.0 point feature. The ablation `tracking/receipts/t1_clinvar_ablation_v1.json` removes the feature.
It cannot remove the label from the allele, because every spiked allele is a ClinVar record of the
release that built the suite.

This script builds a suite that the label cannot solve. It uses a time split. The frozen annotation
resource is the Ensembl VEP 116.2 offline cache, whose header states `ClinVar="202509"`. The cutoff
release is therefore the last ClinVar weekly release of September 2025. An allele enters the
held-out pool when the current release calls it pathogenic and the cutoff release does not hold the
record at all. The frozen cache cannot carry a significance term for such an allele.

The design and the decision gates are in `docs/benchmarks/t1_time_split_v1.md`, recorded before any
result.

Every input is public. The background is GIAB HG002 and the script opens no patient file.

Usage:
    uv run python scripts/track1_time_split_suite.py \\
        --current tracking/cache/clinvar_20260913.vcf.gz \\
        --cutoff tracking/cache/t1_timesplit/clinvar_20250928.vcf.gz \\
        --out reports/benchmarks/timesplit_v1
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

import build_spikein_suite as base  # noqa: E402

SEED = 20260917
N_SCENARIOS = 50
CUTOFF_RELEASE = "clinvar_20250928"
CURRENT_RELEASE = "clinvar_20260913"


def release_index(path: Path) -> dict[str, str]:
    """Read every record of one ClinVar release as a key and its significance term.

    Args:
        path: The bgzipped ClinVar GRCh38 VCF.

    Returns:
        The key `chrom:pos:ref:alt` mapped to the `CLNSIG` value, which may be empty.
    """
    index: dict[str, str] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt, _, _, info = line.rstrip("\n").split("\t")[:8]
            if alt in (".", "") or "," in alt:
                continue
            index[f"{chrom}:{pos}:{ref}:{alt}"] = base.info_field(info, "CLNSIG")
    return index


def variant_key(variant: dict) -> str:
    """Build the release key of a candidate record.

    Args:
        variant: A candidate record with a `chr` prefixed chromosome.

    Returns:
        The key `chrom:pos:ref:alt` without the `chr` prefix.
    """
    return f"{variant['chrom'].removeprefix('chr')}:{variant['pos']}:{variant['ref']}:{variant['alt']}"


def split_pool(by_gene: dict[str, list[dict]], cutoff: dict[str, str]) -> dict[str, list[dict]]:
    """Split the eligible candidates of the current release against the cutoff release.

    Args:
        by_gene: The eligible candidates of the current release, by gene symbol.
        cutoff: The key to significance map of the cutoff release.

    Returns:
        Three lists of records under the keys `new`, `reclassified` and `pre_existing`.
    """
    groups: dict[str, list[dict]] = {"new": [], "reclassified": [], "pre_existing": []}
    for gene, variants in by_gene.items():
        for variant in variants:
            record = dict(variant) | {"gene": gene}
            key = variant_key(variant)
            if key not in cutoff:
                groups["new"].append(record)
            elif cutoff[key] in base.PLP:
                groups["pre_existing"].append(record)
            else:
                groups["reclassified"].append(record)
    return groups


def pool_by_gene(pool: list[dict]) -> dict[str, list[dict]]:
    """Group the held-out candidate records by gene symbol.

    Args:
        pool: The held-out candidate records, each with a `gene` field.

    Returns:
        The gene symbol mapped to its records.
    """
    grouped: dict[str, list[dict]] = defaultdict(list)
    for record in pool:
        grouped[record["gene"]].append(record)
    return grouped


def scenarios_from_pool(pool: list[dict], rng: random.Random, limit: int) -> list[dict]:
    """Build one compound heterozygous scenario per gene that holds two usable alleles.

    The suite mirrors the composition of `spikein_v1`. It takes up to 10 scenarios in the three MVA
    genes first, then fills the rest from the other autosomal recessive genes.

    Args:
        pool: The held-out candidate records, each with a `gene` field.
        rng: The seeded random generator.
        limit: The largest number of scenarios to build.

    Returns:
        The scenarios, each with a gene, a group and two variants.
    """
    grouped = pool_by_gene(pool)
    scenarios: list[dict] = []
    for gene in base.MVA_GENES:
        used: set[tuple[str, int]] = set()
        for _ in range(4):
            if sum(item["group"] == "heldout_mva" for item in scenarios) >= 10:
                break
            available = [item for item in grouped.get(gene, [])
                         if (item["chrom"], item["pos"]) not in used]
            pair = base.pick_pair(rng, available)
            if pair is None:
                break
            used.update((item["chrom"], item["pos"]) for item in pair)
            scenarios.append({"gene": gene, "group": "heldout_mva", "variants": pair})
    genes = sorted(gene for gene, items in grouped.items()
                   if len(items) >= 2 and gene not in base.MVA_GENES)
    rng.shuffle(genes)
    for gene in genes:
        if len(scenarios) >= limit:
            break
        pair = base.pick_pair(rng, grouped[gene])
        if pair:
            scenarios.append({"gene": gene, "group": "heldout_clinvar", "variants": pair})
    return scenarios


def literature_outside_clinvar(table: Path, current: dict[str, str]) -> list[dict]:
    """Find the resolved literature alleles that the current ClinVar release does not hold.

    Args:
        table: The resolved literature allele table.
        current: The key to significance map of the current release.

    Returns:
        The literature records that ClinVar does not hold, each with its gene.
    """
    outside = []
    for gene in ("BUB1B", "CEP57", "TRIP13", "BUB1", "BUB3", "CDC20", "PTPA"):
        for record in base.literature_alleles(table, gene):
            if variant_key(record) not in current:
                outside.append(dict(record) | {"gene": gene})
    return outside


def write_spikes(scenarios: list[dict], out_dir: Path) -> Path:
    """Write every spiked allele of the suite as a heterozygous VCF record.

    Args:
        scenarios: The suite scenarios.
        out_dir: The output directory.

    Returns:
        The path of the written VCF.
    """
    unique = {(v["chrom"], v["pos"], v["ref"], v["alt"]): v
              for scenario in scenarios for v in scenario["variants"]}
    spikes = sorted(unique.values(), key=base.genome_key)
    path = out_dir / "spikes.vcf"
    with path.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n##source=track1_time_split_suite.py\n")
        for chrom in range(1, 23):
            handle.write(f"##contig=<ID=chr{chrom}>\n")
        handle.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Synthetic">\n'
                     '##INFO=<ID=AQ,Number=A,Type=Integer,Description="Synthetic">\n'
                     '##INFO=<ID=FOUND_IN,Number=.,Type=String,Description="Source">\n'
                     '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                     '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
                     '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">\n'
                     '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
                     '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">\n'
                     '##FORMAT=<ID=RNC,Number=2,Type=Character,Description="No-call reason">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n")
        for variant in spikes:
            handle.write(base.spike_line(variant) + "\n")
    return path


def main() -> int:
    """Build the held-out suite and write its manifest.

    Returns:
        Zero on success, two when the pool supports no scenario.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--current", type=Path, required=True)
    parser.add_argument("--cutoff", type=Path, required=True)
    parser.add_argument("--hpoa", type=Path, default=ROOT / "tracking/cache/phenotype.hpoa")
    parser.add_argument("--genes-to-disease", type=Path,
                        default=ROOT / "tracking/cache/genes_to_disease.txt")
    parser.add_argument("--background", type=Path,
                        default=ROOT / "tracking/cache/hg002_s2_deepvariant.vcf.gz")
    parser.add_argument("--alleles", type=Path,
                        default=ROOT / "reports/known_alleles/known_alleles_grch38_2026_09_15.tsv")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    rng = random.Random(SEED)
    genes = base.recessive_genes(args.hpoa, args.genes_to_disease)
    background = base.background_positions(args.background)
    by_gene = base.clinvar_candidates(args.current, genes, background)
    cutoff = release_index(args.cutoff)
    current = release_index(args.current)
    groups = split_pool(by_gene, cutoff)
    scenarios = scenarios_from_pool(groups["new"], rng, N_SCENARIOS)
    for index, scenario in enumerate(scenarios, start=1):
        scenario["scenario_id"] = f"hld_{index:03d}"
    if not scenarios:
        print("The held-out pool supports no scenario.", file=sys.stderr)
        return 2

    literature = literature_outside_clinvar(args.alleles, current)
    literature_by_gene: dict[str, list[dict]] = defaultdict(list)
    for record in literature:
        literature_by_gene[record["gene"]].append(record)
    literature_scenarios = []
    for gene in sorted(literature_by_gene):
        pair = base.pick_pair(rng, literature_by_gene[gene])
        if pair:
            literature_scenarios.append({"gene": gene, "group": "heldout_literature",
                                         "variants": pair})
    for index, scenario in enumerate(literature_scenarios, start=1):
        scenario["scenario_id"] = f"lit_{index:03d}"

    args.out.mkdir(parents=True, exist_ok=True)
    spikes_path = write_spikes(scenarios + literature_scenarios, args.out)
    pool_genes = defaultdict(int)
    for record in groups["new"]:
        pool_genes[record["gene"]] += 1
    manifest = {
        "suite": "timesplit_v1",
        "created_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "data_class": "A",
        "seed": SEED,
        "gate_record": "docs/benchmarks/t1_time_split_v1.md",
        "time_split": {
            "cutoff_release": CUTOFF_RELEASE,
            "cutoff_reason": "The frozen VEP 116.2 offline cache header states ClinVar=202509. The "
                             "cutoff is the last ClinVar weekly release of September 2025.",
            "current_release": CURRENT_RELEASE,
            "cutoff_sha256": base.sha256(args.cutoff),
            "current_sha256": base.sha256(args.current),
            "records_in_the_cutoff_release": len(cutoff),
            "records_in_the_current_release": len(current),
        },
        "pool": {
            "eligible_genes_in_the_current_release": len(by_gene),
            "eligible_alleles_in_the_current_release": sum(len(v) for v in by_gene.values()),
            "new_after_the_cutoff": len(groups["new"]),
            "reclassified_to_pathogenic_after_the_cutoff": len(groups["reclassified"]),
            "pathogenic_before_the_cutoff": len(groups["pre_existing"]),
            "genes_with_two_or_more_new_alleles": sum(1 for count in pool_genes.values() if count >= 2),
        },
        "literature_arm": {
            "resolved_alleles_absent_from_the_current_clinvar_release": len(literature),
            "alleles_by_gene": {gene: len(items) for gene, items in sorted(literature_by_gene.items())},
            "pairs": len(literature_scenarios),
            "rule": "A pair needs two alleles of one gene more than 50 base pairs apart.",
        },
        "composition": {
            "heldout_mva": sum(1 for item in scenarios if item["group"] == "heldout_mva"),
            "heldout_clinvar": sum(1 for item in scenarios if item["group"] == "heldout_clinvar"),
            "new_alleles_in_the_mva_genes": {gene: len(pool_by_gene(groups["new"]).get(gene, []))
                                             for gene in base.MVA_GENES},
            "comparator_in_spikein_v1": "The group recessive_clinvar of spikein_v1 holds 38 "
                                        "scenarios in genes outside the MVA panel and the group "
                                        "mva_clinvar holds 10 scenarios inside it.",
        },
        "recessive_genes": len(genes),
        "n_scenarios": len(scenarios),
        "n_literature_scenarios": len(literature_scenarios),
        "spikes_vcf_sha256": base.sha256(spikes_path),
        "scenarios": scenarios + literature_scenarios,
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"pool": manifest["pool"], "time_split": {
        key: value for key, value in manifest["time_split"].items() if key != "cutoff_reason"},
        "literature_arm": manifest["literature_arm"],
        "n_scenarios": len(scenarios), "n_literature_scenarios": len(literature_scenarios)},
        indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
