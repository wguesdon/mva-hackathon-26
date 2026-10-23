#!/usr/bin/env python3
"""Read the phase of every compound heterozygous pair of a Track 1 submission file.

The script takes a WhatsHap output and a submission file. For every submission row that names two
variants it reports whether the two variants share a phase block and, when they do, whether the two
alternate alleles lie on one haplotype or on opposite haplotypes. A recessive disease needs opposite
haplotypes.

The design and the decision gates are in `docs/benchmarks/t1_phase_v1.md`.

The script writes three files.

    out/class_b2_phase_rows.tsv   one row per submitted pair, with its phase result
    out/class_b2_phase_rows.md    the same content as a short note for owner review
    out/class_a_phase_gates.json  the process gates, the tool versions and the checksums

The class A file carries no result of the sample. It holds no count, no coordinate, no gene, no
genotype and no file size. The writer builds the class A file from a fixed key list.
Sample-specific results stay in the separate class B2 outputs.

Usage:
    phase_pairs.py --phased phased.vcf --submission slot1_b0.csv --out-dir out \\
        --run-id phase_pair_patient_v1 --accuracy-record docs/benchmarks/t1_phase_v1.md
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
from pathlib import Path

# The HG002 pair phase accuracy of each distance band, from tracking/receipts/t1_phase_control_v1.json.
BAND_ACCURACY = {
    "0_to_300": "0.9878 (0.9717 to 0.9948) over 409 phased pairs",
    "301_to_1000": "0.9862 (0.9749 to 0.9925) over 727 phased pairs",
    "1001_to_10000": "0.9752 (0.9639 to 0.9830) over 1049 phased pairs",
    "10001_to_100000": "no pair of this band phased in the HG002 validation",
    "above_100000": "no pair of this band phased in the HG002 validation",
}
BANDS = (("0_to_300", 0, 300), ("301_to_1000", 301, 1000), ("1001_to_10000", 1001, 10000),
         ("10001_to_100000", 10001, 100000), ("above_100000", 100001, 10 ** 12))
CLASS_A_KEYS = ("run_id", "created_utc", "data_class", "scope", "gates", "tools", "checksums",
                "verdict")


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


def open_text(path: Path):
    """Open a plain or a gzipped text file.

    Args:
        path: The file path.

    Returns:
        A text handle.
    """
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def band_of(distance: int) -> str:
    """Return the distance band of a pair.

    Args:
        distance: The distance in base pairs.

    Returns:
        The band name.
    """
    for name, low, high in BANDS:
        if low <= distance <= high:
            return name
    return BANDS[-1][0]


def haplotype_of_alt(genotype: str) -> int | None:
    """Return the haplotype index that carries the alternate allele of a phased heterozygote.

    Args:
        genotype: The GT value.

    Returns:
        0 or 1 for a phased heterozygote, and None otherwise.
    """
    if "|" not in genotype:
        return None
    alleles = genotype.split("|")
    if len(alleles) != 2 or set(alleles) != {"0", "1"}:
        return None
    return 0 if alleles[0] == "1" else 1


def read_phased(path: Path) -> dict[tuple, dict]:
    """Index the phased VCF by variant key.

    Args:
        path: The WhatsHap output.

    Returns:
        The key (chrom, pos, ref, alt) mapped to its phase set and its haplotype.
    """
    index = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            chrom, pos, _, ref, alt = columns[:5]
            if "," in alt or len(columns) < 10:
                continue
            sample = dict(zip(columns[8].split(":"), columns[9].split(":")))
            phase_set = sample.get("PS", ".")
            index[(chrom, int(pos), ref.upper(), alt.upper())] = {
                "ps": phase_set if phase_set not in (".", "") else None,
                "haplotype": haplotype_of_alt(sample.get("GT", "")),
                "gt": sample.get("GT", "")}
    return index


def submitted_pairs(path: Path) -> list[dict]:
    """Read every submission row that names two variants.

    Args:
        path: The submission CSV.

    Returns:
        The rows, each with its rank and its two variant keys.
    """
    pairs = []
    with path.open(newline="") as handle:
        for rank, row in enumerate(csv.DictReader(handle), start=1):
            if not row.get("chrom_2"):
                continue
            first = (row["chrom_1"], int(row["pos_1"]), row["ref_1"].upper(), row["alt_1"].upper())
            second = (row["chrom_2"], int(row["pos_2"]), row["ref_2"].upper(), row["alt_2"].upper())
            pairs.append({"rank": rank, "epcr": row.get("epcr", ""), "notes": row.get("notes", ""),
                          "first": first, "second": second})
    return pairs


def evaluate_pair(pair: dict, index: dict[tuple, dict]) -> dict:
    """Decide the phase result of one submitted pair.

    Args:
        pair: The submission row with its two variant keys.
        index: The phased VCF index.

    Returns:
        The row with its distance, its band, its phase state and its orientation.
    """
    first = index.get(pair["first"])
    second = index.get(pair["second"])
    distance = abs(pair["second"][1] - pair["first"][1])
    result = {"rank": pair["rank"], "epcr": pair["epcr"], "notes": pair["notes"],
              "chrom": pair["first"][0], "pos_1": pair["first"][1], "pos_2": pair["second"][1],
              "distance_bp": distance, "band": band_of(distance),
              "first_in_the_phasing_input": first is not None,
              "second_in_the_phasing_input": second is not None,
              "phase_set_1": first["ps"] if first else None,
              "phase_set_2": second["ps"] if second else None}
    if not first or not second:
        result |= {"phase_state": "NOT IN THE INPUT", "orientation": None,
                   "explains_a_recessive_disease": None}
        return result
    same_block = (first["ps"] is not None and first["ps"] == second["ps"]
                  and first["haplotype"] is not None and second["haplotype"] is not None)
    if not same_block:
        result |= {"phase_state": "NOT PHASED", "orientation": None,
                   "explains_a_recessive_disease": None}
        return result
    orientation = "cis" if first["haplotype"] == second["haplotype"] else "trans"
    result |= {"phase_state": "PHASED", "orientation": orientation,
               "explains_a_recessive_disease": orientation == "trans"}
    return result


def write_rows(results: list[dict], out_dir: Path) -> tuple[Path, Path]:
    """Write the class B2 table and the class B2 note.

    Args:
        results: The evaluated pairs.
        out_dir: The output directory.

    Returns:
        The path of the table and the path of the note.
    """
    fields = ["rank", "epcr", "chrom", "pos_1", "pos_2", "distance_bp", "band", "phase_state",
              "orientation", "explains_a_recessive_disease", "phase_set_1", "phase_set_2",
              "first_in_the_phasing_input", "second_in_the_phasing_input", "notes"]
    table = out_dir / "class_b2_phase_rows.tsv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    note = out_dir / "class_b2_phase_rows.md"
    lines = ["# Read backed phase of the submitted compound heterozygous pairs", "",
             "This file is class B2 content under ADR 0005. It stays in the llm_safe prefix and it",
             "goes to Amazon Bedrock only, after owner review.", "",
             "The method is WhatsHap 2.8 read backed phasing over the seven MVA target intervals.",
             "A pair explains a recessive disease only when the two alternate alleles lie on",
             "opposite haplotypes, which this table calls trans.", "",
             "| Rank | EPCR | Distance, bp | Band | Phase state | Orientation | Trans |",
             "|---|---|---|---|---|---|---|"]
    for item in results:
        lines.append(f"| {item['rank']} | {item['epcr']} | {item['distance_bp']} | "
                     f"{item['band']} | {item['phase_state']} | {item['orientation'] or 'none'} | "
                     f"{item['explains_a_recessive_disease']} |")
    lines += ["", "## The accuracy of the method in each band", "",
              "The numbers come from the GIAB HG002 validation, receipt",
              "`tracking/receipts/t1_phase_control_v1.json`.", "",
              "| Band | Pair phase accuracy on HG002 |", "|---|---|"]
    for name, value in BAND_ACCURACY.items():
        lines.append(f"| {name} | {value} |")
    lines += ["", "A pair with the state NOT PHASED carries no phase evidence. The reads do not",
              "connect the two positions, so the pair stays a hypothesis and its biallelic state is",
              "not demonstrated."]
    note.write_text("\n".join(lines) + "\n")
    return table, note


def build_gates(run_id: str, phased: Path, submission: Path, table: Path, note: Path,
                tools: dict) -> dict:
    """Build the class A gate record from a fixed key list.

    The record holds process gates, tool versions and checksums. It holds no result of the sample.

    Args:
        run_id: The run identifier.
        phased: The phased VCF.
        submission: The submission file.
        table: The class B2 table.
        note: The class B2 note.
        tools: The tool versions and image digests.

    Returns:
        The class A record.
    """
    record = {
        "run_id": run_id,
        "created_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "data_class": "A",
        "scope": "Process gates, tool versions and checksums. This file holds no result of the "
                 "sample, no count, no coordinate, no gene and no genotype. The phase result is "
                 "class B2 and it is in the llm_safe prefix.",
        "gates": {
            "phased_vcf_present": "PASS" if phased.exists() else "FAIL",
            "submission_present": "PASS" if submission.exists() else "FAIL",
            "pair_analysis_complete": "PASS" if table.exists() and note.exists() else "FAIL",
            "class_b2_written_to_llm_safe": "PASS" if table.exists() else "FAIL",
        },
        "tools": tools,
        "checksums": {
            "phased_vcf_sha256": sha256_of(phased),
            "submission_sha256": sha256_of(submission),
            "class_b2_phase_rows_tsv_sha256": sha256_of(table),
            "class_b2_phase_rows_md_sha256": sha256_of(note),
        },
    }
    record["verdict"] = ("COMPLETE" if all(value == "PASS" for value in record["gates"].values())
                         else "FAILED")
    return {key: record[key] for key in CLASS_A_KEYS}


def main() -> int:
    """Read the phase of every submitted pair and write the three files.

    Returns:
        Zero when every gate passes, two otherwise.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phased", required=True, type=Path)
    parser.add_argument("--submission", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--whatshap-version", default="2.8")
    parser.add_argument("--whatshap-image",
                        default="quay.io/biocontainers/whatshap:2.8--py312hf731ba3_0")
    parser.add_argument("--samtools-image",
                        default="quay.io/biocontainers/samtools:1.24--h9dcdb79_1")
    parser.add_argument("--bcftools-image",
                        default="quay.io/biocontainers/bcftools:1.24--h118bc1c_2")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    index = read_phased(args.phased)
    results = [evaluate_pair(pair, index) for pair in submitted_pairs(args.submission)]
    table, note = write_rows(results, args.out_dir)
    tools = {"whatshap": args.whatshap_version, "whatshap_image": args.whatshap_image,
             "samtools_image": args.samtools_image, "bcftools_image": args.bcftools_image,
             "method_record": "docs/benchmarks/t1_phase_v1.md",
             "validation_receipt": "tracking/receipts/t1_phase_control_v1.json"}
    gates = build_gates(args.run_id, args.phased, args.submission, table, note, tools)
    (args.out_dir / "class_a_phase_gates.json").write_text(json.dumps(gates, indent=2) + "\n")
    print(json.dumps({"verdict": gates["verdict"], "gates": gates["gates"],
                      "class_a_file": str(args.out_dir / "class_a_phase_gates.json")}, indent=2))
    return 0 if gates["verdict"] == "COMPLETE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
