#!/usr/bin/env python3
"""Summarize the S6 patient outputs into a class C table and a class B summary for owner review.

The S6 outputs hold patient genotypes, so they are class C. Run this script where class C data is
allowed: the owner machine, or an EC2 instance in the project account. Do not run it in an agent
session that is approved for class A only.

The script writes two files.

    class_c_genotypes.tsv   Every target site with the genotype, the depth and the allele depths.
                            This file stays out of the repository. Delete it with the data.
    class_b_summary.md      Only the sites that are not homozygous reference, the site counts, the
                            depth summary and the comparison of the BUB1B intron signal with the
                            HG002 control. Review this file before you give it to a model.

The HG002 control gives the baseline. HG002 is a public GIAB sample, so its numbers are class A.
A signal in the patient means little until it differs from the control.

Usage:
    uv run python pipelines/s6_summary/bin/summarize_s6_patient.py \\
        --forced proband_forced.vcf.gz \\
        --scan proband_insertion_scan.tsv \\
        --mei proband_mei_in_intervals.txt \\
        --targets reports/known_alleles/s6_targets_2026_09_15/targets.vcf \\
        --control tracking/receipts/s6_hg002_control_v2.json \\
        --out-dir tracking/private/s6_patient_v1
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import statistics
import sys
from pathlib import Path

MIN_DEPTH = 10
BIN_SIZE = 50


def target_names(path: Path) -> dict[tuple[str, int], str]:
    """Map each target position to its identifier.

    The forced genotyping step lost the ID column, so the name comes from the target VCF by position.

    Args:
        path: The target VCF.

    Returns:
        A map from chromosome and position to the target identifier.
    """
    names = {}
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        names[(fields[0], int(fields[1]))] = fields[2]
    return names


def open_maybe_gzip(path: Path):
    """Open a plain or a gzip compressed text file.

    Args:
        path: The file.

    Returns:
        A text mode file object.
    """
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def parse_forced(path: Path, names: dict[tuple[str, int], str]) -> list[dict]:
    """Read the forced genotyping VCF.

    Args:
        path: The forced VCF, bgzip compressed.
        names: The target identifiers by position.

    Returns:
        One record for each site with the genotype, the depth and the allele depths.
    """
    records = []
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            keys = fields[8].split(":")
            values = dict(zip(keys, fields[9].split(":")))
            position = (fields[0], int(fields[1]))
            depth = values.get("DP", ".")
            records.append({
                "target_id": names.get(position, "."),
                "chrom": fields[0], "pos": fields[1], "ref": fields[3], "alt": fields[4],
                "gt": values.get("GT", "./."), "dp": int(depth) if depth.isdigit() else 0,
                "ad": values.get("AD", "."), "qual": fields[5],
            })
    return records


def carrier_state(genotype: str) -> str:
    """Classify a genotype string.

    Args:
        genotype: The GT field.

    Returns:
        One of `hom_ref`, `het`, `hom_alt`, `no_call`.
    """
    alleles = genotype.replace("|", "/").split("/")
    if "." in alleles:
        return "no_call"
    if all(allele == "0" for allele in alleles):
        return "hom_ref"
    return "hom_alt" if all(allele != "0" for allele in alleles) else "het"


def parse_scan(path: Path) -> dict[str, list[dict]]:
    """Read the insertion evidence scan.

    Args:
        path: The scan table.

    Returns:
        The bins of each interval.
    """
    intervals: dict[str, list[dict]] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            numbers = {key: int(row[key]) for key in
                       ("bin_start", "bin_end", "reads", "softclip_20bp", "discordant_other_chrom", "supplementary_sa")}
            intervals.setdefault(row["interval"], []).append(numbers)
    return intervals


def interval_totals(bins: list[dict]) -> dict:
    """Total the evidence of one interval.

    Args:
        bins: The bins of the interval.

    Returns:
        The bin count and the sum of each evidence column.
    """
    totals = {"bins": len(bins)}
    for column in ("reads", "softclip_20bp", "discordant_other_chrom", "supplementary_sa"):
        totals[column] = sum(item[column] for item in bins)
    return totals


def peak_bins(bins: list[dict], column: str, floor: int = 5) -> list[str]:
    """Return the bins where one evidence column is high.

    A bin counts as a peak when the value is at least `floor` and at least three times the median of
    the interval. The median is the baseline and not the mean, because one large bin raises the mean
    enough to hide itself.

    Args:
        bins: The bins of the interval.
        column: The evidence column.
        floor: The smallest value that can be a peak.

    Returns:
        The position range of each peak bin.
    """
    values = [item[column] for item in bins]
    baseline = max(statistics.median(values), 1) if values else 1
    return [f"{item['bin_start']}-{item['bin_end']}" for item in bins
            if item[column] >= floor and item[column] >= 3 * baseline]


def control_interval(control: dict, name: str) -> dict:
    """Return the control numbers of one interval.

    Args:
        control: The HG002 control receipt.
        name: The interval name.

    Returns:
        The control record, or an empty record when the control has no such interval.
    """
    return control.get("insertion_scan", {}).get(name, {})


def summary_lines(records: list[dict], scan: dict[str, list[dict]], mei: list[str], control: dict) -> list[str]:
    """Build the class B summary.

    Args:
        records: The forced genotyping records.
        scan: The insertion scan bins by interval.
        mei: The lines of the mobile element file.
        control: The HG002 control receipt.

    Returns:
        The lines of the summary document.
    """
    states = [carrier_state(record["gt"]) for record in records]
    depths = [record["dp"] for record in records]
    low = [record for record in records if record["dp"] < MIN_DEPTH]
    called = [record for record, state in zip(records, states) if state in ("het", "hom_alt")]
    lines = [
        "# S6 patient summary, class B, owner review required",
        "",
        "This file holds a small number of named variants and QC numbers. It is class B under",
        "`docs/llm_data_privacy.md`. A model may read it only through Amazon Bedrock, and only after",
        "the owner approves it.",
        "",
        "## Forced genotyping at known allele sites",
        "",
        "| Field | Value |",
        "|---|---|",
        f"| Sites genotyped | {len(records)} |",
        f"| Homozygous reference | {states.count('hom_ref')} |",
        f"| Heterozygous | {states.count('het')} |",
        f"| Homozygous alternate | {states.count('hom_alt')} |",
        f"| No call | {states.count('no_call')} |",
        f"| Depth, minimum and median and maximum | {min(depths)}, {int(statistics.median(depths))}, {max(depths)} |",
        f"| Sites below depth {MIN_DEPTH} | {len(low)} |",
        "",
    ]
    if called:
        lines += ["### Sites with a non-reference genotype", "",
                  "| Target | Position | Change | Genotype | Depth | Allele depths |", "|---|---|---|---|---|---|"]
        lines += [f"| {r['target_id']} | {r['chrom']}:{r['pos']} | {r['ref']}>{r['alt']} | {r['gt']} | {r['dp']} | {r['ad']} |"
                  for r in called]
    else:
        lines.append("No site carries a non-reference genotype. The child carries no allele of the target list.")
    lines.append("")
    if low:
        lines += ["### Sites below the depth floor", "",
                  "A low depth means that the genotype is not reliable, not that the allele is absent.", "",
                  "| Target | Position | Depth |", "|---|---|---|"]
        lines += [f"| {r['target_id']} | {r['chrom']}:{r['pos']} | {r['dp']} |" for r in low]
        lines.append("")
    lines += ["## BUB1B intron insertion scan against the HG002 control", "",
              "| Interval | Measure | Patient | HG002 control |", "|---|---|---|---|"]
    for name in sorted(scan):
        patient = interval_totals(scan[name])
        baseline = control_interval(control, name)
        for column in ("bins", "reads", "softclip_20bp", "discordant_other_chrom", "supplementary_sa"):
            lines.append(f"| {name} | {column} | {patient[column]} | {baseline.get(column, 'not in control')} |")
    lines.append("")
    lines += ["### Peak bins in the patient", ""]
    for name in sorted(scan):
        for column in ("softclip_20bp", "discordant_other_chrom", "supplementary_sa"):
            peaks = peak_bins(scan[name], column)
            lines.append(f"- {name}, {column}: {', '.join(peaks) if peaks else 'no peak bin'}")
    lines += ["", "### Mobile element calls inside the intervals", ""]
    if mei:
        lines += [f"- {line}" for line in mei]
    else:
        lines.append("The mobile element caller placed no call inside the two intervals.")
    control_calls = [call for name in sorted(scan) for call in control_interval(control, name).get("retroseq_calls", [])]
    if control_calls:
        lines += ["", f"The HG002 control carries a call at {', '.join(control_calls)}. "
                      "A patient call at the same position is a common allele and not a finding."]
    return lines


def main() -> int:
    """Write the class C table and the class B summary.

    Returns:
        Zero on success.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--forced", type=Path, required=True)
    parser.add_argument("--scan", type=Path, required=True)
    parser.add_argument("--mei", type=Path, required=True)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--control", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    records = parse_forced(args.forced, target_names(args.targets))
    if not records:
        print("The forced VCF holds no record.", file=sys.stderr)
        return 1
    with (args.out_dir / "class_c_genotypes.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    mei = [line for line in args.mei.read_text().splitlines() if line.strip()]
    lines = summary_lines(records, parse_scan(args.scan), mei, json.loads(args.control.read_text()))
    (args.out_dir / "class_b_summary.md").write_text("\n".join(lines) + "\n")
    print(f"Wrote {args.out_dir}/class_c_genotypes.tsv and {args.out_dir}/class_b_summary.md")
    print("Read the summary yourself before you give it to a model.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
