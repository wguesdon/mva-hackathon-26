#!/usr/bin/env python3
"""Measure read backed phasing of variant pairs on GIAB HG002. Data class A.

Defect R2-4 of the second jury round states that the B0 ranker reports a compound heterozygous pair
with no phase. The disease of the case is recessive, so the two alleles must lie on opposite
haplotypes. This script measures whether read backed phasing answers that question on a public
sample, before the method touches the patient interval.

The quantity is a pair level quantity and not a switch error rate. For every pair of heterozygous
call set records inside one target interval the script asks two questions. Does the method place
both records in one phase block? When it does, does the relative orientation of the pair agree with
the published phased truth? The orientation is cis when the two alternate alleles lie on one
haplotype and trans when they lie on opposite haplotypes. A flip of a whole block does not change
the orientation of a pair, so the measure needs no block alignment step.

The design and the decision gates are in `docs/benchmarks/t1_phase_v1.md`, recorded before the
result.

The sample is GIAB HG002 and every input is public. The script opens no patient file.

Usage:
    uv run python scripts/track1_phase_control.py --phased work/hg002_phased.vcf \\
        --truth work/HG002_GRCh38_v5.0q_smvar.vcf.gz \\
        --truth-bed work/HG002_GRCh38_v5.0q_smvar.benchmark.bed \\
        --targets reports/r12_targets_2026_09_15/r12_targets.bed \\
        --receipt tracking/receipts/t1_phase_control_v1.json
"""

from __future__ import annotations

import argparse
import bisect
import datetime as dt
import gzip
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

import track1_interval_stats as stats  # noqa: E402

DISTANCE_BANDS = (("0_to_300", 0, 300), ("301_to_1000", 301, 1000),
                  ("1001_to_10000", 1001, 10000), ("10001_to_100000", 10001, 100000),
                  ("above_100000", 100001, 10 ** 12))
GATE_ACCURACY = 0.95
GATE_LOWER_BOUND = 0.90


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


def read_bed(path: Path) -> dict[str, list[tuple[int, int, str]]]:
    """Read a BED file as half-open intervals by contig.

    Args:
        path: The BED file.

    Returns:
        The contig mapped to its intervals, sorted by start.
    """
    intervals: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            name = fields[3] if len(fields) > 3 else f"{fields[0]}:{fields[1]}-{fields[2]}"
            intervals[fields[0]].append((int(fields[1]), int(fields[2]), name))
    return {contig: sorted(items) for contig, items in intervals.items()}


def locate(intervals: dict[str, list[tuple[int, int, str]]], chrom: str, pos: int) -> str | None:
    """Return the name of the interval that holds a one-based position.

    Args:
        intervals: The intervals by contig.
        chrom: The contig name.
        pos: The one-based position.

    Returns:
        The interval name, or None.
    """
    items = intervals.get(chrom)
    if not items:
        return None
    index = bisect.bisect_right(items, (pos - 1, float("inf"), "")) - 1
    while index >= 0:
        start, end, name = items[index]
        if start <= pos - 1 < end:
            return name
        if end <= pos - 1 - 10 ** 7:
            break
        index -= 1
    return None


def haplotype_of_alt(genotype: str) -> int | None:
    """Return the haplotype index that carries the alternate allele of a phased heterozygote.

    Args:
        genotype: The GT value.

    Returns:
        0 or 1 for a phased heterozygote, and None for any other value.
    """
    if "|" not in genotype:
        return None
    alleles = genotype.split("|")
    if len(alleles) != 2 or set(alleles) != {"0", "1"}:
        return None
    return 0 if alleles[0] == "1" else 1


def read_phased_calls(path: Path, targets: dict[str, list[tuple[int, int, str]]]) -> list[dict]:
    """Read the WhatsHap output over the target intervals.

    Args:
        path: The phased VCF.
        targets: The target intervals by contig.

    Returns:
        One record per call, with its phase set and its phased haplotype when it has one.
    """
    calls = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            chrom, pos, _, ref, alt = columns[:5]
            if "," in alt:
                continue
            interval = locate(targets, chrom, int(pos))
            if interval is None:
                continue
            sample = dict(zip(columns[8].split(":"), columns[9].split(":")))
            genotype = sample.get("GT", "")
            phase_set = sample.get("PS", ".")
            calls.append({"chrom": chrom, "pos": int(pos), "ref": ref.upper(), "alt": alt.upper(),
                          "interval": interval, "gt": genotype,
                          "ps": phase_set if phase_set not in (".", "") else None,
                          "haplotype": haplotype_of_alt(genotype)})
    return calls


def read_truth(path: Path, targets: dict[str, list[tuple[int, int, str]]],
               region: dict[str, list[tuple[int, int, str]]]) -> dict[tuple, int]:
    """Read the phased truth genotypes inside the benchmark region and the target intervals.

    Args:
        path: The truth VCF.
        targets: The target intervals by contig.
        region: The benchmark region by contig.

    Returns:
        The key (chrom, pos, ref, alt) mapped to the haplotype index of the alternate allele.
    """
    truth = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            chrom, pos, _, ref, alt = columns[:5]
            if "," in alt:
                continue
            position = int(pos)
            if locate(targets, chrom, position) is None:
                continue
            if locate(region, chrom, position) is None:
                continue
            sample = dict(zip(columns[8].split(":"), columns[9].split(":")))
            haplotype = haplotype_of_alt(sample.get("GT", ""))
            if haplotype is None:
                continue
            truth[(chrom, position, ref.upper(), alt.upper())] = haplotype
    return truth


def band_of(distance: int) -> str:
    """Return the distance band of a pair.

    Args:
        distance: The distance in base pairs.

    Returns:
        The band name.
    """
    for name, low, high in DISTANCE_BANDS:
        if low <= distance <= high:
            return name
    return DISTANCE_BANDS[-1][0]


def build_pairs(calls: list[dict], truth: dict[tuple, int]) -> list[dict]:
    """Enumerate every testable pair of heterozygous calls inside one target interval.

    A pair is testable when both calls carry a phased truth heterozygote inside the benchmark
    region. The truth orientation is trans when the two alternate alleles lie on opposite
    haplotypes.

    Args:
        calls: The phased call records.
        truth: The truth haplotype of every truth heterozygote.

    Returns:
        The pairs, each with its distance, its band, its truth orientation and its call orientation.
    """
    by_interval: dict[str, list[dict]] = defaultdict(list)
    for call in calls:
        key = (call["chrom"], call["pos"], call["ref"], call["alt"])
        if key in truth:
            by_interval[call["interval"]].append(call | {"truth_haplotype": truth[key]})
    pairs = []
    for interval, items in by_interval.items():
        ordered = sorted(items, key=lambda item: item["pos"])
        for first_index, first in enumerate(ordered):
            for second in ordered[first_index + 1:]:
                distance = second["pos"] - first["pos"]
                same_block = (first["ps"] is not None and first["ps"] == second["ps"]
                              and first["haplotype"] is not None
                              and second["haplotype"] is not None)
                called = None
                if same_block:
                    called = "cis" if first["haplotype"] == second["haplotype"] else "trans"
                expected = ("cis" if first["truth_haplotype"] == second["truth_haplotype"]
                            else "trans")
                pairs.append({"interval": interval, "distance": distance,
                              "band": band_of(distance), "phased": same_block,
                              "called_orientation": called, "truth_orientation": expected,
                              "agrees": called == expected if same_block else None})
    return pairs


def summarize(pairs: list[dict]) -> dict:
    """Summarize the three pair level rates with their Wilson intervals.

    Args:
        pairs: The pairs of one scope.

    Returns:
        The three rates and the counts behind them.
    """
    phased = [pair for pair in pairs if pair["phased"]]
    return {
        "pairs_tested": len(pairs),
        "pairs_phased": len(phased),
        "pair_phasing_rate": stats.wilson(len(phased), len(pairs)),
        "pair_phase_accuracy": stats.wilson(sum(1 for pair in phased if pair["agrees"]),
                                            len(phased)),
        "pair_trans_rate": stats.wilson(
            sum(1 for pair in phased if pair["called_orientation"] == "trans"), len(phased)),
        "truth_trans_rate": stats.wilson(
            sum(1 for pair in pairs if pair["truth_orientation"] == "trans"), len(pairs)),
    }


def band_gate(summary: dict) -> dict:
    """Evaluate gate P1 on one distance band.

    Args:
        summary: The summary of the band.

    Returns:
        The gate verdict of the band.
    """
    accuracy = summary["pair_phase_accuracy"]
    if accuracy["estimate"] is None:
        return {"outcome": "NO PAIR", "accuracy": None, "lower_bound": None}
    holds = accuracy["estimate"] >= GATE_ACCURACY and accuracy["low"] >= GATE_LOWER_BOUND
    return {"outcome": "GATE MET" if holds else "GATE NOT MET",
            "accuracy": accuracy["estimate"], "lower_bound": accuracy["low"],
            "required_accuracy": GATE_ACCURACY, "required_lower_bound": GATE_LOWER_BOUND}


def main() -> int:
    """Measure the pair level rates and write the receipt.

    Returns:
        Zero after the receipt is written.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phased", required=True, type=Path)
    parser.add_argument("--truth", required=True, type=Path)
    parser.add_argument("--truth-bed", required=True, type=Path)
    parser.add_argument("--targets", required=True, type=Path)
    parser.add_argument("--whatshap-version", default="2.8")
    parser.add_argument("--whatshap-image",
                        default="quay.io/biocontainers/whatshap:2.8--py312hf731ba3_0")
    parser.add_argument("--whatshap-digest",
                        default="sha256:e50ebc2cc57c796b41248139c384d8c9d6885aa7713e3685058ddeaa9c856d78")
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()

    targets = read_bed(args.targets)
    region = read_bed(args.truth_bed)
    calls = read_phased_calls(args.phased, targets)
    truth = read_truth(args.truth, targets, region)
    pairs = build_pairs(calls, truth)

    by_band = {name: [pair for pair in pairs if pair["band"] == name]
               for name, _, _ in DISTANCE_BANDS}
    by_interval: dict[str, list[dict]] = defaultdict(list)
    for pair in pairs:
        by_interval[pair["interval"]].append(pair)

    band_summaries = {name: summarize(items) for name, items in by_band.items()}
    receipt = {
        "receipt_id": "t1_phase_control_v1",
        "created_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "track": "t1_phase", "queue_item": "R22", "data_class": "A",
        "defect": "R2-4 of reports/jury/jury_review_track1_2026-09-17T154227Z.md",
        "gate_record": "docs/benchmarks/t1_phase_v1.md",
        "receipt_scope": "GIAB HG002 alone. The alignment, the call set and the truth set are "
                         "public. No patient file enters this run.",
        "tools": {"whatshap": args.whatshap_version, "whatshap_image": args.whatshap_image,
                  "whatshap_image_digest": args.whatshap_digest,
                  "samtools": "1.24, quay.io/biocontainers/samtools:1.24--h9dcdb79_1",
                  "bcftools": "1.24, quay.io/biocontainers/bcftools:1.24--h118bc1c_2",
                  "command": "whatshap phase --reference <GRCh38> --ignore-read-groups --indels"},
        "inputs": {
            "alignment": "Caller-supplied GIAB HG002 alignment; this program reads its phased output.",
            "call_set": "giab_hg002_snv.vcf.gz of run giab_hg002_calling_v1",
            "truth": "GIAB HG002 v5.0q small variant benchmark on GRCh38, "
                     "HG002_GRCh38_v5.0q_smvar.vcf.gz",
            "truth_sha256": sha256_of(args.truth),
            "truth_bed": args.truth_bed.name,
            "truth_bed_sha256": sha256_of(args.truth_bed),
            "targets": args.targets.name,
            "targets_sha256": sha256_of(args.targets),
            "phased_vcf_sha256": sha256_of(args.phased)},
        "counts": {
            "heterozygous_calls_in_the_target_intervals": len(calls),
            "calls_in_a_phase_block": sum(1 for call in calls if call["ps"] is not None),
            "truth_heterozygotes_inside_the_benchmark_region": len(truth),
            "calls_with_a_truth_heterozygote": sum(
                1 for call in calls
                if (call["chrom"], call["pos"], call["ref"], call["alt"]) in truth),
            "phase_blocks": len({call["ps"] for call in calls if call["ps"] is not None})},
        "overall": summarize(pairs),
        "by_distance_band": band_summaries,
        "by_interval": {name: summarize(items) for name, items in sorted(by_interval.items())},
        "gate_p1_by_band": {name: band_gate(summary)
                            for name, summary in band_summaries.items()},
        "gate_p1_overall": band_gate(summarize(pairs)),
        "verdict": "MEASURED",
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"counts": receipt["counts"], "overall": receipt["overall"],
                      "by_distance_band": receipt["by_distance_band"],
                      "gate_p1_by_band": receipt["gate_p1_by_band"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
