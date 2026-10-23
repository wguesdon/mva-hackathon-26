#!/usr/bin/env python3
"""Liftover control for the Track 1 structural variant benchmark. Class A.

The jury review of 2026-09-17 recorded defect T1-16. The R2 benchmark lifts the GIAB HG002
structural variant benchmark v0.6 from GRCh37 to GRCh38 with the UCSC hg19 to hg38 chain, and the
report asserts that the lift does not change the conclusion. The first repair pass removed the
claim that a similar CMRG precision validates the lift, but it ran no control at the same loci.
This script runs that control.

Three measurements run.

    Agreement of the lifted truth set with a native GRCh38 truth set
        The lifted v0.6 records and the native GIAB CMRG structural variant records v1.00 are
        compared inside the intersection of the two region files. A record of one set pairs with a
        record of the other set under the match rule of the benchmark, which is the same class, a
        breakpoint distance of 1000 base pairs or less and a size ratio of 0.70 or more. The
        measurement reports the agreement per class, the records that only one set holds and the
        breakpoint shift of each pair. The two sets come from different integration pipelines, so a
        disagreement here holds both the lift error and the difference between the two truth sets.

    The lift error alone
        The CMRG structural variant set v1.00 is published on GRCh37 and on GRCh38. The two files
        hold the same events from the same producer, and a sequence-resolved record carries the same
        REF and ALT allele on both builds. The script pairs the two files by that allele, lifts the
        GRCh37 record with the chain and the rule of the benchmark, and compares the lifted
        coordinate with the native GRCh38 coordinate of the same event. The difference is the lift
        error with no truth set difference inside it.

    The effect of the lift on a measured recall
        The prepared HG002 call set of the R2 run runs against the native GRCh38 CMRG truth and
        against the same CMRG truth lifted from GRCh37, with the same Truvari parameters. Both truth
        files are written as symbolic records by the same writer, so the two arms differ in the
        coordinate source and in nothing else. The paired truth records feed an exact McNemar test.

Every input is public and data class A. The script opens no patient file.

Usage:
    uv run --no-project --with "pyliftover==0.4.1" python scripts/track1_liftover_control.py \\
        --work-dir <dir> --receipt tracking/receipts/t1_liftover_control_v1.json
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import track1_benchmark_intervals as stats  # noqa: E402
import track1_sv_benchmark as bench  # noqa: E402
import track1_sv_common as common  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]

# The match rule of the R2 benchmark. Truvari runs with --refdist 1000 and --pctsize 0.7.
MATCH_REFDIST = 1000
MATCH_PCTSIZE = 0.70

# The classes that both truth sets carry.
TRUTH_CLASSES = ("DEL", "INS")

# The smallest structural variant that the benchmark scores.
MIN_SIZE = 50

# The public source of every input that this control adds to the R2 inputs.
SOURCE_URLS = {
    "cmrg_grch37_vcf": ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
                        "AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh37/StructuralVariant/"
                        "HG002_GRCh37_CMRG_SV_v1.00.vcf.gz"),
    "cmrg_grch37_bed": ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
                        "AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh37/StructuralVariant/"
                        "HG002_GRCh37_CMRG_SV_v1.00.bed"),
    "cmrg_grch38_vcf": ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
                        "AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/"
                        "HG002_GRCh38_CMRG_SV_v1.00.vcf.gz"),
    "cmrg_grch38_bed": ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
                        "AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/"
                        "HG002_GRCh38_CMRG_SV_v1.00.bed"),
    "chain": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
}

# The licence of every public input of this control.
LICENCES = {
    "cmrg": ("Public domain. NIST asks for a citation of doi:10.1038/s41587-021-01158-1."),
    "chain": "UCSC genome annotation data are free for all use, with attribution.",
}


def read_intervals(path: Path) -> list[tuple[str, int, int]]:
    """Read a BED file and merge its intervals.

    Args:
        path: The BED path.

    Returns:
        The merged intervals as contig, start and end.
    """
    return common.merge_intervals(common.read_bed(path))


def add_chr_prefix(intervals: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    """Give every contig name of an interval set the chr prefix.

    Args:
        intervals: The intervals as contig, start and end.

    Returns:
        The intervals with a chr prefix on every contig name.
    """
    return [(name if name.startswith("chr") else f"chr{name}", start, end)
            for name, start, end in intervals]


def intersect_interval_sets(first: list[tuple[str, int, int]],
                            second: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    """Intersect two merged interval sets.

    Args:
        first: The first merged interval set.
        second: The second merged interval set.

    Returns:
        The merged intervals that both sets cover.
    """
    index = common.index_intervals(second)
    shared: list[tuple[str, int, int]] = []
    for chrom, start, end in first:
        for other_start, other_end in index.get(chrom, []):
            if other_end <= start:
                continue
            if other_start >= end:
                break
            shared.append((chrom, max(start, other_start), min(end, other_end)))
    return common.merge_intervals(shared)


def prepare_truth_records(records: list[common.Variant], prefix: str) -> tuple[list[common.Variant], dict[str, int]]:
    """Select the truth records that the control can score and give each one an identifier.

    A record is kept when its FILTER is PASS or empty, its class is a deletion or an insertion and
    its size is 50 base pairs or larger. The identifier lets a record survive the lift and the
    Truvari output, so the same event can be followed through every arm.

    Args:
        records: The records of a truth VCF.
        prefix: The identifier prefix of this truth file.

    Returns:
        The kept records and the counts by reason.
    """
    kept: list[common.Variant] = []
    counts = {"input": len(records), "not_pass": 0, "wrong_type": 0, "below_50bp": 0, "kept": 0}
    for index, record in enumerate(records):
        if record.filter_field not in ("PASS", ".", ""):
            counts["not_pass"] += 1
            continue
        svtype = common.variant_type(record)
        if svtype not in TRUTH_CLASSES:
            counts["wrong_type"] += 1
            continue
        size = common.variant_size(record)
        if size is None or size < MIN_SIZE:
            counts["below_50bp"] += 1
            continue
        kept.append(common.Variant(
            chrom=record.chrom, pos=record.pos, identifier=f"{prefix}_{index}", ref=record.ref,
            alt=record.alt, qual=record.qual, filter_field="PASS", info=dict(record.info),
            tail=list(record.tail)))
        counts["kept"] += 1
    return kept, counts


def symbolic_record(record: common.Variant) -> common.Variant:
    """Write one sequence-resolved truth record as a symbolic record.

    The writer matches the writer of `track1_sv_benchmark.lift_truth_records`, so a native arm and
    a lifted arm carry the same representation and differ only in the coordinate source.

    Args:
        record: The sequence-resolved truth record.

    Returns:
        The symbolic record with SVTYPE, SVLEN and END.
    """
    svtype = common.variant_type(record)
    size = common.variant_size(record)
    end = record.pos + size if svtype == "DEL" else record.pos
    genotype = record.tail[1].split(":")[0] if len(record.tail) > 1 else "./."
    return common.Variant(
        chrom=record.chrom, pos=record.pos, identifier=record.identifier, ref="N",
        alt=f"<{svtype}>", qual=".", filter_field="PASS",
        info={"SVTYPE": svtype, "SVLEN": str(-size if svtype == "DEL" else size), "END": str(end)},
        tail=["GT", genotype])


def relabel(record: common.Variant, identifier: str) -> common.Variant:
    """Return a copy of one record under a new identifier.

    The lifted arm and the native arm must carry the same identifier for the same event, so that a
    paired test can follow one event through both Truvari runs.

    Args:
        record: The truth record.
        identifier: The identifier to write into the copy.

    Returns:
        The copy of the record.
    """
    return common.Variant(chrom=record.chrom, pos=record.pos, identifier=identifier,
                          ref=record.ref, alt=record.alt, qual=record.qual,
                          filter_field="PASS", info=dict(record.info), tail=list(record.tail))


def allele_key(record: common.Variant) -> tuple[str, str] | None:
    """Return the allele identity of a sequence-resolved record.

    Args:
        record: The truth record.

    Returns:
        The upper case REF and ALT pair, or None when the record is symbolic.
    """
    if record.alt.startswith("<") or record.alt in (".", "") or "," in record.alt:
        return None
    return record.ref.upper(), record.alt.upper()


def index_by_allele(records: list[common.Variant]) -> tuple[dict[tuple[str, str], common.Variant], int, int]:
    """Group records by their allele identity and drop every ambiguous key.

    Args:
        records: The truth records of one build.

    Returns:
        The unique key to record mapping, the number of symbolic records and the number of records
        that share a key with another record.
    """
    groups: dict[tuple[str, str], list[common.Variant]] = {}
    symbolic = 0
    for record in records:
        key = allele_key(record)
        if key is None:
            symbolic += 1
            continue
        groups.setdefault(key, []).append(record)
    unique = {key: items[0] for key, items in groups.items() if len(items) == 1}
    ambiguous = sum(len(items) for items in groups.values() if len(items) > 1)
    return unique, symbolic, ambiguous


def pair_builds_by_allele(source: list[common.Variant],
                          target: list[common.Variant]) -> dict:
    """Pair the records of two builds of the same truth set by their allele identity.

    Args:
        source: The GRCh37 records.
        target: The GRCh38 records.

    Returns:
        The pairs as source and target record tuples, and the counts of every unpaired group.
    """
    source_index, source_symbolic, source_ambiguous = index_by_allele(source)
    target_index, target_symbolic, target_ambiguous = index_by_allele(target)
    shared = sorted(set(source_index) & set(target_index))
    pairs = [(source_index[key], target_index[key]) for key in shared]
    return {
        "pairs": pairs,
        "counts": {
            "source_records": len(source),
            "target_records": len(target),
            "source_symbolic": source_symbolic,
            "target_symbolic": target_symbolic,
            "source_ambiguous_allele": source_ambiguous,
            "target_ambiguous_allele": target_ambiguous,
            "paired_by_allele": len(pairs),
            "source_only": len(source_index) - len(pairs),
            "target_only": len(target_index) - len(pairs),
        },
    }


def shift_statistics(shifts: list[int]) -> dict:
    """Summarise a set of signed coordinate shifts.

    Args:
        shifts: The signed shifts in base pairs.

    Returns:
        The count, the number of exact agreements, the median, the maximum and the number of shifts
        that exceed the Truvari breakpoint distance.
    """
    if not shifts:
        return {"pairs": 0, "exact": 0, "median_abs": None, "max_abs": None,
                "beyond_refdist": 0, "refdist": MATCH_REFDIST}
    absolute = sorted(abs(value) for value in shifts)
    return {
        "pairs": len(shifts),
        "exact": sum(1 for value in shifts if value == 0),
        "median_abs": statistics.median(absolute),
        "mean_abs": round(statistics.fmean(absolute), 4),
        "p95_abs": absolute[min(len(absolute) - 1, int(0.95 * len(absolute)))],
        "max_abs": absolute[-1],
        "beyond_refdist": sum(1 for value in absolute if value > MATCH_REFDIST),
        "refdist": MATCH_REFDIST,
    }


def records_in_regions(records: list[common.Variant],
                       index: dict[str, list[tuple[int, int]]]) -> list[common.Variant]:
    """Keep the records whose start position falls inside an interval set.

    Args:
        records: The truth records.
        index: The output of `common.index_intervals`.

    Returns:
        The records inside the intervals.
    """
    return [record for record in records
            if common.position_in_index(index, record.chrom, record.pos)]


def size_ratio(first: int, second: int) -> float:
    """Return the size ratio of two variants under the Truvari rule.

    Args:
        first: The size of the first variant in base pairs.
        second: The size of the second variant in base pairs.

    Returns:
        The smaller size divided by the larger size.
    """
    if max(first, second) == 0:
        return 0.0
    return min(first, second) / max(first, second)


def match_two_truth_sets(left: list[common.Variant], right: list[common.Variant]) -> dict:
    """Pair the records of two truth sets under the match rule of the benchmark.

    The pairing is greedy. Every candidate pair of the same class inside the breakpoint distance
    and the size ratio is sorted by the breakpoint distance, then by the size difference, and a
    record joins at most one pair.

    Args:
        left: The records of the first truth set.
        right: The records of the second truth set.

    Returns:
        The pairs, the unmatched records of each side and the signed shift of each pair.
    """
    candidates = []
    right_by_contig: dict[str, list[int]] = {}
    for index, record in enumerate(right):
        right_by_contig.setdefault(record.chrom, []).append(index)
    for left_index, left_record in enumerate(left):
        left_type = common.variant_type(left_record)
        left_size = common.variant_size(left_record) or 0
        for right_index in right_by_contig.get(left_record.chrom, []):
            right_record = right[right_index]
            if common.variant_type(right_record) != left_type:
                continue
            distance = abs(right_record.pos - left_record.pos)
            if distance > MATCH_REFDIST:
                continue
            right_size = common.variant_size(right_record) or 0
            if size_ratio(left_size, right_size) < MATCH_PCTSIZE:
                continue
            candidates.append((distance, abs(left_size - right_size), left_index, right_index))
    candidates.sort()
    used_left: set[int] = set()
    used_right: set[int] = set()
    pairs = []
    for distance, size_difference, left_index, right_index in candidates:
        if left_index in used_left or right_index in used_right:
            continue
        used_left.add(left_index)
        used_right.add(right_index)
        pairs.append((left_index, right_index, distance, size_difference))
    return {
        "pairs": pairs,
        "left_unmatched": [index for index in range(len(left)) if index not in used_left],
        "right_unmatched": [index for index in range(len(right)) if index not in used_right],
    }


def class_table(records: list[common.Variant], selected: set[int]) -> dict[str, dict[str, int]]:
    """Count records per class and split them into matched and unmatched.

    Args:
        records: The truth records.
        selected: The indexes of the matched records.

    Returns:
        A mapping from class to the total, the matched count and the unmatched count.
    """
    table: dict[str, dict[str, int]] = {}
    for index, record in enumerate(records):
        svtype = common.variant_type(record)
        row = table.setdefault(svtype, {"records": 0, "matched": 0, "unmatched": 0})
        row["records"] += 1
        if index in selected:
            row["matched"] += 1
        else:
            row["unmatched"] += 1
    return table


def describe_record(record: common.Variant) -> dict:
    """Return the fields of one truth record that the receipt reports.

    Args:
        record: The truth record.

    Returns:
        The contig, the position, the class, the size and the identifier.
    """
    return {"chrom": record.chrom, "pos": record.pos, "class": common.variant_type(record),
            "size": common.variant_size(record), "id": record.identifier}


def agreement_at_shared_loci(lifted: list[common.Variant], native: list[common.Variant],
                             shared: list[tuple[str, int, int]]) -> dict:
    """Compare the lifted v0.6 truth set with the native CMRG truth set at the shared loci.

    Args:
        lifted: The lifted GIAB v0.6 truth records on GRCh38.
        native: The native GRCh38 CMRG truth records.
        shared: The merged intervals that both region files cover.

    Returns:
        The region size, the record counts, the agreement per class, the shift statistics and the
        records that only one set holds.
    """
    index = common.index_intervals(shared)
    lifted_in = records_in_regions(lifted, index)
    native_in = records_in_regions(native, index)
    matched = match_two_truth_sets(lifted_in, native_in)
    lifted_matched = {pair[0] for pair in matched["pairs"]}
    native_matched = {pair[1] for pair in matched["pairs"]}
    shifts = [native_in[pair[1]].pos - lifted_in[pair[0]].pos for pair in matched["pairs"]]
    return {
        "shared_regions": {"intervals": len(shared), "bases": common.total_bases(shared)},
        "records": {"lifted_v06_in_shared_regions": len(lifted_in),
                    "native_cmrg_in_shared_regions": len(native_in),
                    "matched_pairs": len(matched["pairs"]),
                    "lifted_only": len(matched["left_unmatched"]),
                    "native_only": len(matched["right_unmatched"])},
        "by_class_lifted_v06": class_table(lifted_in, lifted_matched),
        "by_class_native_cmrg": class_table(native_in, native_matched),
        "breakpoint_shift_native_minus_lifted": shift_statistics(shifts),
        "lifted_only_records": [describe_record(lifted_in[index]) for index in matched["left_unmatched"]],
        "native_only_records": [describe_record(native_in[index]) for index in matched["right_unmatched"]],
        "note": ("The two sets come from different integration pipelines, so a record that only one "
                 "set holds is not by itself a lift error. The lift error alone is measured in the "
                 "next section."),
    }


def lift_error_on_the_same_set(source: list[common.Variant], target: list[common.Variant],
                               lifter) -> dict:
    """Measure the lift error on the CMRG set, which is published on both builds.

    Args:
        source: The prepared GRCh37 CMRG records.
        target: The prepared GRCh38 CMRG records.
        lifter: A pyliftover LiftOver instance for the hg19 to hg38 chain.

    Returns:
        The pairing counts, the lift loss, the position shift, the size change and every record
        whose lifted coordinate differs from its native coordinate.
    """
    paired = pair_builds_by_allele(source, target)
    source_records = [pair[0] for pair in paired["pairs"]]
    native_by_id = {pair[0].identifier: pair[1] for pair in paired["pairs"]}
    lifted, lift_counts = bench.lift_truth_records(source_records, lifter)
    shifts: list[int] = []
    size_changes: list[int] = []
    wrong_contig = 0
    disagreements = []
    for record in lifted:
        native = native_by_id[record.identifier]
        if record.chrom != native.chrom:
            wrong_contig += 1
            disagreements.append({"id": record.identifier, "reason": "contig differs",
                                  "lifted": describe_record(record),
                                  "native": describe_record(native)})
            continue
        shift = record.pos - native.pos
        lifted_size = common.variant_size(record) or 0
        native_size = common.variant_size(native) or 0
        shifts.append(shift)
        size_changes.append(lifted_size - native_size)
        if shift != 0 or lifted_size != native_size:
            disagreements.append({"id": record.identifier, "reason": "coordinate or size differs",
                                  "shift_bp": shift, "size_change_bp": lifted_size - native_size,
                                  "lifted": describe_record(record),
                                  "native": describe_record(native)})
    by_class: dict[str, dict[str, int]] = {}
    for record in lifted:
        svtype = common.variant_type(record)
        native = native_by_id[record.identifier]
        row = by_class.setdefault(svtype, {"lifted": 0, "exact_position": 0, "exact_size": 0})
        row["lifted"] += 1
        row["exact_position"] += int(record.chrom == native.chrom and record.pos == native.pos)
        row["exact_size"] += int((common.variant_size(record) or 0) == (common.variant_size(native) or 0))
    return {
        "pairing": paired["counts"],
        "lift_counts": lift_counts,
        "records_lost_by_the_lift": paired["counts"]["paired_by_allele"] - len(lifted),
        "lifted_to_a_different_contig": wrong_contig,
        "position_shift_lifted_minus_native": shift_statistics(shifts),
        "size_change_lifted_minus_native": shift_statistics(size_changes),
        "by_class": by_class,
        "disagreements": disagreements,
    }


def region_lift_error(lifted_regions: list[tuple[str, int, int]],
                      native_regions: list[tuple[str, int, int]]) -> dict:
    """Compare the lifted CMRG regions with the native CMRG regions.

    Args:
        lifted_regions: The GRCh37 CMRG regions lifted to GRCh38.
        native_regions: The native GRCh38 CMRG regions.

    Returns:
        The base counts of each set and of the intersection.
    """
    shared = intersect_interval_sets(lifted_regions, native_regions)
    lifted_bases = common.total_bases(lifted_regions)
    native_bases = common.total_bases(native_regions)
    shared_bases = common.total_bases(shared)
    return {
        "lifted_intervals": len(lifted_regions),
        "native_intervals": len(native_regions),
        "lifted_bases": lifted_bases,
        "native_bases": native_bases,
        "shared_bases": shared_bases,
        "lifted_only_bases": lifted_bases - shared_bases,
        "native_only_bases": native_bases - shared_bases,
        "shared_share_of_native": round(shared_bases / native_bases, 6) if native_bases else None,
    }


def bench_summary(bench_dir: Path) -> dict:
    """Read the summary of one Truvari output directory.

    Args:
        bench_dir: The Truvari output directory.

    Returns:
        The counts and the three measures of the run.
    """
    summary = json.loads((bench_dir / "summary.json").read_text())
    return bench.metrics(summary["TP-base"], summary["FN"], summary["TP-comp"], summary["FP"])


def matched_identifiers(bench_dir: Path) -> set[str]:
    """Read the identifiers of the truth records that one Truvari run recovered.

    Args:
        bench_dir: The Truvari output directory.

    Returns:
        The identifiers of the records in `tp-base.vcf.gz`.
    """
    _, records = common.read_vcf(bench_dir / "tp-base.vcf.gz")
    return {record.identifier for record in records}


def paired_recall_test(native_dir: Path, lifted_dir: Path, shared_ids: set[str]) -> dict:
    """Run the exact McNemar test on the truth records that both arms carry.

    Args:
        native_dir: The Truvari output of the native GRCh38 arm.
        lifted_dir: The Truvari output of the lifted arm.
        shared_ids: The identifiers that both arms scored.

    Returns:
        The two arm recalls, the discordant counts and the two-sided p value.
    """
    native_hits = matched_identifiers(native_dir) & shared_ids
    lifted_hits = matched_identifiers(lifted_dir) & shared_ids
    native_only = len(native_hits - lifted_hits)
    lifted_only = len(lifted_hits - native_hits)
    return {
        "paired_truth_records": len(shared_ids),
        "native_recall": stats.wilson_interval(len(native_hits), len(shared_ids)),
        "lifted_recall": stats.wilson_interval(len(lifted_hits), len(shared_ids)),
        "mcnemar_exact": stats.mcnemar_exact(native_only, lifted_only),
        "unit": "One truth record of the CMRG structural variant set v1.00.",
    }


def write_truth_file(work_dir: Path, name: str, header: list[str],
                     records: list[common.Variant]) -> str:
    """Write, sort, compress and index one symbolic truth VCF.

    Args:
        work_dir: The host work directory.
        name: The plain VCF file name.
        header: The header lines.
        records: The symbolic records.

    Returns:
        The file name of the bgzipped VCF.
    """
    common.write_vcf(work_dir / name, header, records)
    return bench.compress_and_index(name, work_dir)


def recall_effect(work_dir: Path, calls_name: str, header: list[str],
                  native_records: list[common.Variant], lifted_records: list[common.Variant],
                  native_regions: list[tuple[str, int, int]],
                  lifted_regions: list[tuple[str, int, int]]) -> dict:
    """Measure the recall of the call set against the native truth and against the lifted truth.

    Args:
        work_dir: The host work directory.
        calls_name: The file name of the prepared call set inside the work directory.
        header: The header lines for the symbolic truth files.
        native_records: The native GRCh38 CMRG records as symbolic records.
        lifted_records: The GRCh37 CMRG records lifted to GRCh38 as symbolic records.
        native_regions: The native GRCh38 CMRG regions.
        lifted_regions: The lifted CMRG regions.

    Returns:
        The two Truvari summaries, the third arm that swaps the region file, and the paired test.
    """
    common.write_bed(work_dir / "cmrg_native.bed", native_regions)
    common.write_bed(work_dir / "cmrg_lifted.bed", lifted_regions)
    native_file = write_truth_file(work_dir, "cmrg_native_symbolic.vcf", header, native_records)
    lifted_file = write_truth_file(work_dir, "cmrg_lifted_symbolic.vcf", header, lifted_records)
    native_dir = bench.run_bench(work_dir, native_file, calls_name, "cmrg_native.bed",
                                 "bench_cmrg_native", True)
    lifted_dir = bench.run_bench(work_dir, lifted_file, calls_name, "cmrg_lifted.bed",
                                 "bench_cmrg_lifted", True)
    lifted_same_bed = bench.run_bench(work_dir, lifted_file, calls_name, "cmrg_native.bed",
                                      "bench_cmrg_lifted_native_bed", True)
    native_index = common.index_intervals(native_regions)
    lifted_index = common.index_intervals(lifted_regions)
    native_ids = {record.identifier for record in native_records
                  if common.position_in_index(native_index, record.chrom, record.pos)}
    lifted_ids = {record.identifier for record in lifted_records
                  if common.position_in_index(lifted_index, record.chrom, record.pos)}
    return {
        "native_grch38_truth": bench_summary(native_dir),
        "lifted_from_grch37_truth": bench_summary(lifted_dir),
        "lifted_truth_with_the_native_region_file": bench_summary(lifted_same_bed),
        "paired_test": paired_recall_test(native_dir, lifted_dir, native_ids & lifted_ids),
        "truth_records_only_in_the_native_arm": len(native_ids - lifted_ids),
        "truth_records_only_in_the_lifted_arm": len(lifted_ids - native_ids),
    }


def build_receipt(agreement: dict, lift_error: dict, regions: dict, recall: dict,
                  inputs: dict, reproduction: dict) -> dict:
    """Assemble the receipt of the control.

    Args:
        agreement: The result of the shared loci comparison.
        lift_error: The result of the same set lift comparison.
        regions: The result of the region comparison.
        recall: The result of the recall comparison.
        inputs: The input files with their checksums.
        reproduction: The reproduction of the committed CMRG numbers.

    Returns:
        The receipt as a mapping.
    """
    return {
        "receipt_id": "t1_liftover_control_v1",
        "created_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "track": "t1_call",
        "data_class": "A",
        "queue_item": "R19",
        "defect": "T1-16 of reports/jury/defect_ledger_track1_2026-09-17.md",
        "question": ("Does the liftover of the GIAB v0.6 truth set from GRCh37 to GRCh38 change the "
                     "recall that the R2 benchmark reports?"),
        "acceptance_rule": (
            "The rule predates this run. tracking/receipts/t1_sv_mei_benchmark_v1.json of "
            "2026-09-15T21:16:11Z froze the match rule at --refdist 1000 and --pctsize 0.7. A lift "
            "that moves a breakpoint by 1000 base pairs or less and keeps the size ratio at 0.70 or "
            "more cannot change the match status of any call, so it cannot change a reported recall. "
            "No separate gate file was written before this run, because the run chooses no "
            "configuration and the thresholds it tests against are already committed."),
        "verdict": "MEASURED",
        "tools": {"liftover": "pyliftover 0.4.1", "truvari": bench.TRUVARI_IMAGE,
                  "bcftools": bench.BCFTOOLS_IMAGE},
        "truvari_arguments": list(bench.TRUVARI_ARGS) + ["--dup-to-ins"],
        "match_rule": {"refdist_bp": MATCH_REFDIST, "pctsize": MATCH_PCTSIZE,
                       "min_size_bp": MIN_SIZE, "classes": list(TRUTH_CLASSES)},
        "inputs": inputs,
        "reproduction_of_the_committed_cmrg_run": reproduction,
        "agreement_of_the_lifted_set_with_a_native_set_at_the_same_loci": agreement,
        "lift_error_measured_on_a_set_published_on_both_builds": lift_error,
        "region_lift": regions,
        "effect_on_recall": recall,
    }


def main() -> int:
    """Run the three measurements and write the receipt.

    Returns:
        Zero when every measurement completes.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--lifted-vcf", required=True, type=Path)
    parser.add_argument("--lifted-bed", required=True, type=Path)
    parser.add_argument("--cmrg38-vcf", required=True, type=Path)
    parser.add_argument("--cmrg38-bed", required=True, type=Path)
    parser.add_argument("--cmrg37-vcf", required=True, type=Path)
    parser.add_argument("--cmrg37-bed", required=True, type=Path)
    parser.add_argument("--calls-vcf", required=True, type=Path)
    parser.add_argument("--chain", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()

    from pyliftover import LiftOver

    work_dir = args.work_dir
    work_dir.mkdir(parents=True, exist_ok=True)
    lifter = LiftOver(str(args.chain))

    call_header, _ = common.read_vcf(args.calls_vcf)
    header = bench.contig_header(call_header, "HG002")
    calls_name = "calls_sv.vcf.gz"
    (work_dir / calls_name).write_bytes(args.calls_vcf.read_bytes())
    bench.run_container(bench.BCFTOOLS_IMAGE, ["bcftools", "index", "-f", "-t", calls_name], work_dir)

    _, lifted_v06 = common.read_vcf(args.lifted_vcf)
    _, cmrg38_raw = common.read_vcf(args.cmrg38_vcf)
    _, cmrg37_raw = common.read_vcf(args.cmrg37_vcf)
    cmrg38, cmrg38_counts = prepare_truth_records(cmrg38_raw, "cmrg38")
    cmrg37, cmrg37_counts = prepare_truth_records(cmrg37_raw, "cmrg37")

    lifted_bed = read_intervals(args.lifted_bed)
    cmrg38_bed = read_intervals(args.cmrg38_bed)
    cmrg37_bed = add_chr_prefix(read_intervals(args.cmrg37_bed))
    shared = intersect_interval_sets(lifted_bed, cmrg38_bed)

    agreement = agreement_at_shared_loci(lifted_v06, cmrg38, shared)
    lift_error = lift_error_on_the_same_set(cmrg37, cmrg38, lifter)
    lifted_regions, region_counts = bench.lift_bed(cmrg37_bed, lifter)
    regions = region_lift_error(lifted_regions, cmrg38_bed)
    regions["lift_counts"] = region_counts

    paired = pair_builds_by_allele(cmrg37, cmrg38)
    lifted_symbolic, _ = bench.lift_truth_records([pair[0] for pair in paired["pairs"]], lifter)
    lifted_ids = {record.identifier for record in lifted_symbolic}
    native_symbolic = [symbolic_record(relabel(target, source.identifier))
                       for source, target in paired["pairs"]
                       if source.identifier in lifted_ids]

    recall = recall_effect(work_dir, calls_name, header, native_symbolic, lifted_symbolic,
                           cmrg38_bed, lifted_regions)

    raw_native = write_truth_file(work_dir, "cmrg_native_raw.vcf",
                                  common.read_vcf(args.cmrg38_vcf)[0], cmrg38_raw)
    reproduction_dir = bench.run_bench(work_dir, raw_native, calls_name, "cmrg_native.bed",
                                       "bench_cmrg_reproduction", True)
    reproduction = {
        "arm": "The native GRCh38 CMRG truth file as the R2 benchmark ran it.",
        "committed": {"receipt": "tracking/receipts/t1_sv_mei_benchmark_v1.json",
                      "field": "results.cmrg_v1_00"},
        "measured": bench_summary(reproduction_dir),
    }

    inputs = {}
    for name, path in (("lifted_v06_grch38_vcf", args.lifted_vcf),
                       ("lifted_v06_grch38_bed", args.lifted_bed),
                       ("cmrg_grch38_vcf", args.cmrg38_vcf),
                       ("cmrg_grch38_bed", args.cmrg38_bed),
                       ("cmrg_grch37_vcf", args.cmrg37_vcf),
                       ("cmrg_grch37_bed", args.cmrg37_bed),
                       ("chain", args.chain),
                       ("calls_vcf", args.calls_vcf)):
        inputs[name] = {"path": path.name, "sha256": bench.sha256_of(path)}
        if name in SOURCE_URLS:
            inputs[name]["url"] = SOURCE_URLS[name]
            inputs[name]["licence"] = LICENCES["chain" if name == "chain" else "cmrg"]
    inputs["cmrg_grch38_preparation"] = cmrg38_counts
    inputs["cmrg_grch37_preparation"] = cmrg37_counts

    receipt = build_receipt(agreement, lift_error, regions, recall, inputs, reproduction)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"agreement": agreement["records"],
                      "shift": agreement["breakpoint_shift_native_minus_lifted"],
                      "lift_error": lift_error["position_shift_lifted_minus_native"],
                      "recall": {"native": recall["native_grch38_truth"],
                                 "lifted": recall["lifted_from_grch37_truth"],
                                 "mcnemar": recall["paired_test"]["mcnemar_exact"]},
                      "reproduction": reproduction["measured"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
