"""Build the S4 spike-in benchmark suite on the HG002 S2 DeepVariant calls (ADR 0003). Data class A.

Each scenario adds one compound heterozygous pair of ClinVar pathogenic or likely pathogenic
variants in one autosomal recessive gene to the HG002 background, as heterozygous calls.

Sources:
    ClinVar GRCh38 VCF (tracking/cache/clinvar_<date>.vcf.gz)
    HPO phenotype.hpoa and genes_to_disease.txt (autosomal recessive diseases: HP:0000007)
    HG002 S2 calls (tracking/cache/hg002_s2_deepvariant.vcf.gz)
    Resolved literature alleles (reports/known_alleles/known_alleles_grch38_<date>.tsv)

Rules:
    ClinVar CLNSIG is Pathogenic, Likely_pathogenic or Pathogenic/Likely_pathogenic.
    CLNREVSTAT contains criteria_provided and no conflict.
    The variant is an SNV or an indel of 50 bp or less, with one gene in GENEINFO.
    The position is absent from the HG002 background.
    At least 10 scenarios use MVA genes (BUB1B, CEP57, TRIP13); two of them use literature non-coding alleles.

Outputs (reports/benchmarks/spikein_v1/):
    manifest.json      scenarios, truth pairs, sources, seed and hashes
    spikes.vcf         every spiked variant with a synthetic heterozygous DeepVariant-like genotype

Usage:
    uv run python scripts/build_spikein_suite.py --clinvar tracking/cache/clinvar_20260913.vcf.gz
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SEED = 20260915
# The benchmark MVA gene set of the frozen gene panel, one tuple per panel version. Panel version
# v1 held this list as an undeclared third copy of its three established genes, so
# harness/test_track1_gene_panel.py did not check it and the suite spikein_v1 carried no CENATAC
# and no MAD1L1 scenario. Panel version v2 declares this module as a consumer. The shipped default
# stays at v1, because the committed suite manifest of spikein_v1 reproduces from it.
PANEL_SETS = {
    "v1": ("BUB1B", "CEP57", "TRIP13"),
    "v2": ("BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1"),
}
SHIPPED_PANEL_VERSION = "v1"
MVA_GENES = PANEL_SETS[SHIPPED_PANEL_VERSION]
PLP = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
N_SCENARIOS = 50
N_MVA = 10
CHROMS = {str(number) for number in range(1, 23)}


def sha256(path: Path) -> str:
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


def position_key(variant: dict) -> int:
    """Sort key for variants on one chromosome.

    Args:
        variant: A variant record.

    Returns:
        The position.
    """
    return variant["pos"]


def genome_key(variant: dict) -> tuple[int, int]:
    """Sort key for autosomal variants across chromosomes.

    Args:
        variant: A variant record with a chr prefix.

    Returns:
        The chromosome number and the position.
    """
    return int(variant["chrom"].removeprefix("chr")), variant["pos"]


def recessive_genes(hpoa: Path, genes_to_disease: Path) -> set[str]:
    """Return gene symbols linked to at least one autosomal recessive disease.

    Args:
        hpoa: The HPO phenotype.hpoa file.
        genes_to_disease: The HPO genes_to_disease.txt file.

    Returns:
        The gene symbols.
    """
    diseases = set()
    with hpoa.open() as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("database_id"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 10 and fields[3] == "HP:0000007" and fields[10] == "I":
                diseases.add(fields[0])
    genes = set()
    with genes_to_disease.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["disease_id"] in diseases and row["association_type"] == "MENDELIAN":
                genes.add(row["gene_symbol"])
    return genes


def background_positions(vcf: Path) -> set[tuple[str, int]]:
    """Read the positions of the HG002 background calls.

    Args:
        vcf: The bgzipped VCF.

    Returns:
        The set of (chrom, pos) with a chr prefix.
    """
    positions = set()
    with gzip.open(vcf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos = line.split("\t", 2)[:2]
            positions.add((chrom, int(pos)))
    return positions


def info_field(info: str, key: str) -> str:
    """Read one INFO value.

    Args:
        info: The INFO column.
        key: The key.

    Returns:
        The value, or an empty string.
    """
    for item in info.split(";"):
        if item.startswith(key + "="):
            return item[len(key) + 1:]
    return ""


def clinvar_candidates(clinvar: Path, genes: set[str], background: set[tuple[str, int]]) -> dict[str, list[dict]]:
    """Collect eligible ClinVar P/LP variants by gene.

    Args:
        clinvar: The ClinVar GRCh38 VCF.
        genes: The allowed gene symbols.
        background: Positions present in HG002.

    Returns:
        Gene symbol to variant records.
    """
    by_gene: dict[str, list[dict]] = defaultdict(list)
    with gzip.open(clinvar, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos, vid, ref, alt, _, _, info = line.rstrip("\n").split("\t")[:8]
            if chrom not in CHROMS or alt in (".", "") or "," in alt:
                continue
            if info_field(info, "CLNSIG") not in PLP:
                continue
            status = info_field(info, "CLNREVSTAT")
            if "criteria_provided" not in status or "conflicting" in status:
                continue
            if len(ref) > 51 or len(alt) > 51 or not set(ref + alt) <= set("ACGT"):
                continue
            geneinfo = info_field(info, "GENEINFO")
            if not geneinfo or "|" in geneinfo:
                continue
            gene = geneinfo.split(":")[0]
            if gene not in genes and gene not in MVA_GENES:
                continue
            if ("chr" + chrom, int(pos)) in background:
                continue
            by_gene[gene].append({"chrom": "chr" + chrom, "pos": int(pos), "ref": ref, "alt": alt,
                                  "source": f"ClinVar:{vid}", "clnsig": info_field(info, "CLNSIG"),
                                  "mc": info_field(info, "MC").split("|")[-1] if info_field(info, "MC") else ""})
    return by_gene


def literature_alleles(table: Path, gene: str) -> list[dict]:
    """Read resolved literature alleles of one gene.

    Args:
        table: The S3 resolved allele table.
        gene: The gene symbol.

    Returns:
        Variant records.
    """
    rows = []
    with table.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["gene"] == gene and row["status"] == "resolved" and row["pos"]:
                rows.append({"chrom": "chr" + row["chrom"], "pos": int(row["pos"]), "ref": row["ref"], "alt": row["alt"],
                             "source": f"literature:{row['allele_id']}", "clnsig": "", "mc": ""})
    return rows


def pick_pair(rng: random.Random, variants: list[dict]) -> list[dict] | None:
    """Pick two variants at distinct positions more than 50 bp apart.

    Args:
        rng: The random generator.
        variants: Candidate variants of one gene.

    Returns:
        Two variants, or None.
    """
    pool = variants[:]
    rng.shuffle(pool)
    for first_index, first in enumerate(pool):
        for second in pool[first_index + 1:]:
            if abs(first["pos"] - second["pos"]) > 50:
                return sorted([first, second], key=position_key)
    return None


def spike_line(variant: dict) -> str:
    """Format one spiked variant as a DeepVariant-like heterozygous VCF line.

    Args:
        variant: The variant record.

    Returns:
        The VCF data line without a newline.
    """
    identifier = f"{variant['chrom']}_{variant['pos']}_{variant['ref']}_{variant['alt']}"
    return "\t".join([variant["chrom"], str(variant["pos"]), identifier, variant["ref"], variant["alt"], "50", "PASS",
                      "AF=0.5;AQ=50;FOUND_IN=spikein", "GT:DP:AD:GQ:PL:RNC", "0/1:36:18,18:50:50,0,50:.."])


def main() -> int:
    """Build the suite.

    Returns:
        Zero on success.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clinvar", type=Path, required=True)
    parser.add_argument("--hpoa", type=Path, default=ROOT / "tracking/cache/phenotype.hpoa")
    parser.add_argument("--genes-to-disease", type=Path, default=ROOT / "tracking/cache/genes_to_disease.txt")
    parser.add_argument("--background", type=Path, default=ROOT / "tracking/cache/hg002_s2_deepvariant.vcf.gz")
    parser.add_argument("--alleles", type=Path, default=ROOT / "reports/known_alleles/known_alleles_grch38_2026_09_15.tsv")
    parser.add_argument("--out", type=Path, default=ROOT / "reports/benchmarks/spikein_v1")
    args = parser.parse_args()
    rng = random.Random(SEED)
    genes = recessive_genes(args.hpoa, args.genes_to_disease)
    background = background_positions(args.background)
    by_gene = clinvar_candidates(args.clinvar, genes, background)
    scenarios = []
    # MVA scenarios: two with literature non-coding alleles, the rest from ClinVar pairs.
    bub1b_literature = {item["source"].split(":")[1]: item for item in literature_alleles(args.alleles, "BUB1B")}
    literature_pairs = [("bub1b_upstream_rs576524605", "bub1b_c2363_2364del"), ("bub1b_c1402_5A_G", "bub1b_c2386_11A_G")]
    for first, second in literature_pairs:
        pair = sorted([bub1b_literature[first], bub1b_literature[second]], key=position_key)
        if all((item["chrom"], item["pos"]) not in background for item in pair):
            scenarios.append({"gene": "BUB1B", "group": "mva_literature", "variants": pair})
    mva_quota = {"BUB1B": 3, "CEP57": 3, "TRIP13": 4}
    for gene, quota in mva_quota.items():
        used: set[tuple[str, int]] = set()
        for _ in range(quota):
            pool = [item for item in by_gene.get(gene, []) if (item["chrom"], item["pos"]) not in used]
            pair = pick_pair(rng, pool)
            if pair is None:
                break
            used.update((item["chrom"], item["pos"]) for item in pair)
            scenarios.append({"gene": gene, "group": "mva_clinvar", "variants": pair})
    other_genes = sorted(gene for gene, items in by_gene.items() if gene not in MVA_GENES and len(items) >= 2)
    rng.shuffle(other_genes)
    for gene in other_genes:
        if len(scenarios) >= N_SCENARIOS:
            break
        pair = pick_pair(rng, by_gene[gene])
        if pair:
            scenarios.append({"gene": gene, "group": "recessive_clinvar", "variants": pair})
    for index, scenario in enumerate(scenarios, start=1):
        scenario["scenario_id"] = f"spk_{index:03d}"
    n_mva = sum(scenario["group"].startswith("mva") for scenario in scenarios)
    if len(scenarios) != N_SCENARIOS or n_mva < N_MVA:
        print(f"The suite has {len(scenarios)} scenarios and {n_mva} MVA scenarios; rules not met.", file=sys.stderr)
        return 2
    args.out.mkdir(parents=True, exist_ok=True)
    spikes = sorted({(v["chrom"], v["pos"], v["ref"], v["alt"]): v for s in scenarios for v in s["variants"]}.values(),
                    key=genome_key)
    spike_path = args.out / "spikes.vcf"
    with spike_path.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n##source=build_spikein_suite.py\n")
        for chrom in range(1, 23):
            handle.write(f"##contig=<ID=chr{chrom}>\n")
        handle.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Synthetic">\n##INFO=<ID=AQ,Number=A,Type=Integer,Description="Synthetic">\n'
                     '##INFO=<ID=FOUND_IN,Number=.,Type=String,Description="Source">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                     '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">\n'
                     '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">\n'
                     '##FORMAT=<ID=RNC,Number=2,Type=Character,Description="No-call reason">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n")
        for variant in spikes:
            handle.write(spike_line(variant) + "\n")
    manifest = {"suite": "spikein_v1", "created_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"), "seed": SEED,
                "data_class": "A", "background": {"path": str(args.background),
                                                    "sha256": sha256(args.background)},
                "sources": {"clinvar": {"file": args.clinvar.name, "sha256": sha256(args.clinvar)},
                            "hpoa": {"file": args.hpoa.name, "sha256": sha256(args.hpoa)},
                            "genes_to_disease": {"file": args.genes_to_disease.name, "sha256": sha256(args.genes_to_disease)},
                            "alleles": {"file": str(args.alleles.relative_to(ROOT)), "sha256": sha256(args.alleles)}},
                "recessive_genes": len(genes), "eligible_genes": len(by_gene), "n_scenarios": len(scenarios), "n_mva": n_mva,
                "spikes_vcf_sha256": sha256(spike_path), "scenarios": scenarios}
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"{len(scenarios)} scenarios, {n_mva} MVA, {len(spikes)} spiked variants, {len(genes)} recessive genes, "
          f"{len(by_gene)} eligible genes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
