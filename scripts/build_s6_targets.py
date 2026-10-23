"""Build the S6 target sites and insertion scan intervals from public data (ADR 0003). Data class A.

Outputs (reports/known_alleles/s6_targets_<date>/):
    targets.vcf        sites-only VCF: resolved literature alleles and ClinVar P/LP SNVs and small indels
                       in BUB1B, CEP57, TRIP13 and SMC5, with a chr prefix
    intervals.bed      BUB1B MANE intron 2 and a 600 bp window around the intron 8 target site duplication
    manifest.json      counts, sources and hashes

Usage:
    uv run python scripts/build_s6_targets.py --clinvar tracking/cache/clinvar_20260913.vcf.gz
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# The forced genotyping target list of literature/track1_gene_panel_v1.yaml, which is the frozen
# gene panel of this work. SMC5 is a declared manual addition of that file on PMID 41374403, and
# the panel rule keeps it out of the structural annotation set because SMC5 is not a spindle
# assembly checkpoint component. harness/test_track1_gene_panel.py fails when this tuple no
# longer equals the operating set that the panel derives.
#
# PANEL_SETS holds one forced genotyping target list per panel version. Version v2 adds CENATAC and
# MAD1L1, which the established clause of v1 excluded through an enumerated disease list, and SLF2,
# TRIM37, MAD2L1BP and CEP192, which carry a published biallelic genotype with an aneuploidy or a
# hyperploidy phenotype. The shipped default stays at v1, because the run that produced the
# submission used the v1 list and queue row R30 blocks every patient rerun.
PANEL_SETS = {
    "v1": ("BUB1B", "CEP57", "TRIP13", "SMC5"),
    "v2": ("BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1", "SMC5", "SLF2", "TRIM37", "MAD2L1BP",
           "CEP192"),
}
SHIPPED_PANEL_VERSION = "v1"
GENES = PANEL_SETS[SHIPPED_PANEL_VERSION]
PLP = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}


def info_value(info: str, key: str) -> str:
    """Read one INFO value.

    Args:
        info: The INFO column.
        key: The key.

    Returns:
        The value or an empty string.
    """
    return next((item[len(key) + 1:] for item in info.split(";") if item.startswith(key + "=")), "")


def site_key(key: tuple) -> tuple[int, int]:
    """Sort key for (chrom, pos, ref, alt) sites.

    Args:
        key: The site tuple with a chr prefix.

    Returns:
        The chromosome number, 99 for other contigs, and the position.
    """
    name = key[0].removeprefix("chr")
    return (int(name) if name.isdigit() else 99), key[1]


def main() -> int:
    """Write targets, intervals and a manifest.

    Returns:
        Zero on success.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clinvar", type=Path, required=True)
    parser.add_argument("--alleles", type=Path, default=ROOT / "reports/known_alleles/known_alleles_grch38_2026_09_15.tsv")
    args = parser.parse_args()
    stamp = dt.datetime.now(dt.UTC).strftime("%Y_%m_%d")
    out = ROOT / f"reports/known_alleles/s6_targets_{stamp}"
    out.mkdir(parents=True, exist_ok=True)
    sites: dict[tuple, str] = {}
    intervals = []
    with args.alleles.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["status"] == "resolved":
                sites[("chr" + row["chrom"], int(row["pos"]), row["ref"], row["alt"])] = f"literature:{row['allele_id']}"
            elif row["status"] in ("interval_resolved", "intron_level"):
                start, end = int(row["interval_start"]), int(row["interval_end"])
                if row["status"] == "interval_resolved":
                    start, end = start - 300, end + 300
                intervals.append(("chr" + row["chrom"], start - 1, end, row["allele_id"]))
    n_literature = len(sites)
    with gzip.open(args.clinvar, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos, vid, ref, alt, _, _, info = line.rstrip("\n").split("\t")[:8]
            gene = info_value(info, "GENEINFO").split(":")[0]
            if gene not in GENES or info_value(info, "CLNSIG") not in PLP or "," in alt or alt in (".", ""):
                continue
            if len(ref) > 51 or len(alt) > 51 or not set(ref + alt) <= set("ACGT"):
                continue
            sites.setdefault(("chr" + chrom, int(pos), ref, alt), f"ClinVar:{vid}:{gene}")
    order = sorted(sites, key=site_key)
    targets = out / "targets.vcf"
    with targets.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        for chrom in sorted({key[0] for key in order}):
            handle.write(f"##contig=<ID={chrom}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for key in order:
            handle.write(f"{key[0]}\t{key[1]}\t{sites[key]}\t{key[2]}\t{key[3]}\t.\t.\t.\n")
    bed = out / "intervals.bed"
    bed.write_text("".join(f"{chrom}\t{start}\t{end}\t{name}\n" for chrom, start, end, name in intervals))
    manifest = {"created_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"), "data_class": "A",
                "genes": GENES, "literature_sites": n_literature, "total_sites": len(sites), "intervals": len(intervals),
                "clinvar": {"file": args.clinvar.name, "sha256": hashlib.sha256(args.clinvar.read_bytes()).hexdigest()},
                "alleles_sha256": hashlib.sha256(args.alleles.read_bytes()).hexdigest(),
                "targets_sha256": hashlib.sha256(targets.read_bytes()).hexdigest(),
                "intervals_sha256": hashlib.sha256(bed.read_bytes()).hexdigest()}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({key: manifest[key] for key in ("literature_sites", "total_sites", "intervals")}), out.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
