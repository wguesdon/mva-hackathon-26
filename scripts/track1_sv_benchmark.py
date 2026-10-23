"""Benchmark the Track 1 structural variant, copy number and mobile element callers on HG002.

The call sets come from the S2 run `giab_hg002_calling_v1`. That run used nf-core/raredisease
3.1.2 with the same caller parameters as the patient run S5, so the measured recall applies to
the patient call set. The truth set is the GIAB HG002 structural variant benchmark v0.6, which
is published on GRCh37. The script lifts the truth positions and the high confidence regions to
GRCh38 with a UCSC chain file, then compares the call set with Truvari inside a container.

The GIAB CMRG structural variant set v1.00 runs as a second comparison. That set is published on
GRCh38, so it needs no liftover and it checks the liftover step.

Every input is public and data class A. The script opens no patient file.

Usage:
    uv run python scripts/track1_sv_benchmark.py --work-dir <dir> --receipt <path>
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import track1_sv_common as common  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]

TRUVARI_IMAGE = "quay.io/biocontainers/truvari:4.3.1--pyhdfd78af_0"
BCFTOOLS_IMAGE = "quay.io/biocontainers/bcftools:1.22--h3a4d415_1"

# Truvari thresholds. The sequence check is off because TIDDIT and CNVnator report symbolic
# alleles with no inserted sequence. The GIAB v0.6 README recommends this setting for a call set
# that is not sequence resolved.
TRUVARI_ARGS = ("--refdist", "1000", "--pctseq", "0", "--pctsize", "0.7", "--pctovl", "0.0",
                "--sizemin", "50", "--sizefilt", "50", "--sizemax", "1000000", "--passonly")

# The classes that the call set reports. The truth set covers only DEL and INS.
CALL_CLASSES = ("DEL", "INS", "DUP", "INV", "BND")
TRUTH_CLASSES = ("DEL", "INS")

# A lifted deletion keeps its record when the lifted span differs from the source span by less
# than this fraction, or by less than 20 base pairs.
LIFT_LENGTH_TOLERANCE = 0.10
LIFT_LENGTH_FLOOR = 20

# A RetroSeq call and a truth insertion are one event when their positions are inside this window.
MEI_MATCH_WINDOW = 500

# The structural annotation set of the frozen gene panel, one tuple per panel version.
# harness/test_track1_gene_panel.py fails when a tuple no longer equals the operating set that the
# matching panel derives. The shipped default stays at v1, because the recorded structural
# benchmark ran over the v1 set.
PANEL_SETS = {
    "v1": ("BUB1B", "CEP57", "TRIP13", "BUB1", "BUB3", "CDC20", "PTPA"),
    "v2": ("BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1", "MAD2L1BP", "BUB1", "BUB3", "CDC20",
           "PTPA"),
}
SHIPPED_PANEL_VERSION = "v1"

# Ensembl release 116, GRCh38. The span of every gene of every panel version, read from the Ensembl
# REST lookup endpoint on 2026-09-17.
GENE_SPANS = {
    "BUB1B": ("chr15", 40160984, 40221137),
    "CEP57": ("chr11", 95789965, 95837070),
    "TRIP13": ("chr5", 892849, 919357),
    "BUB1": ("chr2", 110635468, 110678098),
    "BUB3": ("chr10", 123154365, 123313144),
    "CDC20": ("chr1", 43358907, 43363212),
    "PTPA": ("chr9", 129110550, 129148946),
    "CENATAC": ("chr11", 118998051, 119015811),
    "MAD1L1": ("chr7", 1815787, 2233243),
    "MAD2L1BP": ("chr6", 43629494, 43640960),
}

# The genes of the mosaic variegated aneuploidy panel that this benchmark reports, with their span.
MVA_GENES = {gene: GENE_SPANS[gene] for gene in PANEL_SETS[SHIPPED_PANEL_VERSION]}


# The call set of the S2 run. The callers and the parameters are the ones of the patient run S5.
CALLERS = {
    "workflow": "nf-core/raredisease v3.1.2-g83f2699 on Nextflow 26.04.6",
    "structural_variant_paired_and_split_read": "TIDDIT 3.9.5",
    "structural_variant_assembly": "Manta 1.6.0",
    "copy_number_read_depth": "CNVnator 0.4.1",
    "mobile_element": "RetroSeq 1.5",
    "merge": "SVDB 2.8.4, --merge --pass_only --same_order --priority tiddit,manta,cnvnator",
    "mobile_element_merge": "SVDB 2.8.4, --merge --bnd_distance 150 --overlap 0.5",
    "aligner": "bwa-mem2 2.2.1, duplicates marked with Picard 3.4.0",
    "skipped": "GATK GermlineCNVCaller, Gens, vcf2cytosure, MitoSalt, peddy, SMNCopyNumberCaller",
}

# Every input file is public. NIST publishes the GIAB files, and works of NIST carry no copyright
# in the United States. The consortium asks for a citation.
SOURCES = {
    "HG002_SVs_Tier1_v0.6.vcf.gz": {
        "set": "GIAB HG002 structural variant benchmark v0.6, Tier 1, GRCh37",
        "sample": "HG002 / NA24385, GIAB Ashkenazi son, Coriell NA24385",
        "url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
               "analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz",
        "mirror": "s3://giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/",
        "licence": "Public domain. NIST asks for a citation of doi:10.1101/664623.",
    },
    "HG002_SVs_Tier1_noVDJorXorY_v0.6.2.bed": {
        "set": "GIAB HG002 v0.6.2 high confidence regions, without the VDJ regions, chrX and chrY",
        "url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
               "analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_noVDJorXorY_v0.6.2.bed",
        "licence": "Public domain. NIST asks for a citation of doi:10.1101/664623.",
    },
    "HG002_GRCh38_CMRG_SV_v1.00.vcf.gz": {
        "set": "GIAB HG002 challenging medically relevant genes structural variant benchmark v1.00, GRCh38",
        "url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/"
               "HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz",
        "licence": "Public domain. NIST asks for a citation of doi:10.1038/s41587-021-01158-1.",
    },
    "HG002_GRCh38_CMRG_SV_v1.00.bed": {
        "set": "GIAB CMRG structural variant benchmark regions v1.00, GRCh38",
        "url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/"
               "HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.bed",
        "licence": "Public domain. NIST asks for a citation of doi:10.1038/s41587-021-01158-1.",
    },
    "hg19ToHg38.over.chain.gz": {
        "set": "UCSC chain file from hg19 to hg38",
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
        "licence": "UCSC genome annotation data are free for all use, with attribution.",
    },
    "giab_hg002_sv.vcf.gz": {
        "set": "Merged structural variant and copy number calls of the S2 run giab_hg002_calling_v1",
        "licence": "Derived from public GIAB reads. Data class A.",
    },
    "giab_hg002_mobile_elements.vcf.gz": {
        "set": "Mobile element calls of the S2 run giab_hg002_calling_v1",
        "licence": "Derived from public GIAB reads. Data class A.",
    },
}

# The limits that the GIAB README states for v0.6, and the limits that this comparison adds.
KNOWN_LIMITS = [
    "The v0.6 Tier 1 set holds insertions and deletions of 50 base pairs and larger only. It holds "
    "no duplication, no inversion and no break end, so those classes carry no truth value here.",
    "The v0.6 Tier 1 set labels no insertion as a mobile element insertion, so the set gives no "
    "recall by mobile element class.",
    "The GIAB README states that Tier 1 holds 50 percent or less of the structural variants of the "
    "genome, and that the set is biased towards the variants that are easier to call.",
    "The v0.6 set is published on GRCh37. This benchmark lifts the positions and the regions to "
    "GRCh38 with a UCSC chain file, and it reports the number of records that the lift loses.",
    "The comparison ignores the sequence of the allele, because TIDDIT and CNVnator report a "
    "symbolic allele. A call therefore matches on position, class and size only.",
    "chrX and chrY are outside the v0.6.2 regions, and the mitochondrion is outside every region.",
]


def sha256_of(path: Path) -> str:
    """Return the SHA-256 checksum of a file.

    Args:
        path: The file path.

    Returns:
        The checksum in lowercase hexadecimal.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def run_container(image: str, command: list[str], work_dir: Path) -> str:
    """Run one command inside a container with the work directory mounted.

    Args:
        image: The container image reference.
        command: The command and its arguments, with paths relative to the work directory.
        work_dir: The host directory to mount at /work.

    Returns:
        The standard output of the command.

    Raises:
        RuntimeError: The command exited with a non-zero status.
    """
    argv = ["podman", "run", "--rm", "-v", f"{work_dir}:/work:z", "-w", "/work", image]
    argv.extend(command)
    result = subprocess.run(argv, capture_output=True, text=True, timeout=7200)
    if result.returncode != 0:
        raise RuntimeError(f"{' '.join(command[:3])} failed: {result.stderr[-2000:]}")
    return result.stdout


def compress_and_index(name: str, work_dir: Path) -> str:
    """Sort, compress and index a plain VCF inside the work directory.

    Args:
        name: The file name of the plain VCF, relative to the work directory.
        work_dir: The host work directory.

    Returns:
        The file name of the bgzipped VCF.
    """
    output = name + ".gz"
    run_container(BCFTOOLS_IMAGE, ["bcftools", "sort", "-Oz", "-o", output, name], work_dir)
    run_container(BCFTOOLS_IMAGE, ["bcftools", "index", "-t", output], work_dir)
    return output


def contig_header(source_header: list[str], sample: str) -> list[str]:
    """Build a small GRCh38 VCF header from the contig lines of a call set.

    Args:
        source_header: The header lines of the GRCh38 call set.
        sample: The sample name for the column line.

    Returns:
        The header lines of the new VCF.
    """
    lines = ["##fileformat=VCFv4.2"]
    lines.extend(line for line in source_header if line.startswith("##contig="))
    lines.extend([
        '##ALT=<ID=DEL,Description="Deletion">',
        '##ALT=<ID=INS,Description="Insertion">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}",
    ])
    return lines


def lift_point(lifter, chrom: str, position: int) -> tuple[str, int] | None:
    """Lift one zero-based coordinate to GRCh38.

    Args:
        lifter: A pyliftover LiftOver instance for the hg19 to hg38 chain.
        chrom: The source contig name with the chr prefix.
        position: The zero-based source coordinate.

    Returns:
        The target contig and the zero-based target coordinate, or None when the point does not
        lift, lifts to the reverse strand, lifts to more than one place, or lands off a primary
        contig.
    """
    hits = lifter.convert_coordinate(chrom, position)
    if not hits or len(hits) != 1:
        return None
    target_chrom, target_pos, strand = hits[0][0], hits[0][1], hits[0][2]
    if strand != "+" or target_chrom not in common.PRIMARY_CONTIGS:
        return None
    return target_chrom, target_pos


def lift_truth_records(records: list[common.Variant], lifter) -> tuple[list[common.Variant], dict[str, int]]:
    """Lift the GIAB v0.6 truth records to GRCh38 and write them as symbolic records.

    Only PASS records of at least 50 base pairs are kept. A deletion keeps its record when both
    breakpoints lift to the same contig and the lifted span keeps its length. An insertion keeps
    its record when its position lifts, because the inserted sequence is a property of the sample
    and not of the reference.

    Args:
        records: The source records of the v0.6 VCF.
        lifter: A pyliftover LiftOver instance for the hg19 to hg38 chain.

    Returns:
        The lifted records and the counts of every reason for a loss.
    """
    kept: list[common.Variant] = []
    counts = {"input": len(records), "not_pass": 0, "wrong_type": 0, "below_50bp": 0,
              "pos_not_lifted": 0, "end_not_lifted": 0, "span_changed": 0, "kept": 0}
    for record in records:
        if record.filter_field != "PASS":
            counts["not_pass"] += 1
            continue
        svtype = common.variant_type(record)
        if svtype not in TRUTH_CLASSES:
            counts["wrong_type"] += 1
            continue
        size = common.variant_size(record)
        if size is None or size < 50:
            counts["below_50bp"] += 1
            continue
        source_chrom = record.chrom if record.chrom.startswith("chr") else f"chr{record.chrom}"
        lifted_start = lift_point(lifter, source_chrom, record.pos - 1)
        if lifted_start is None:
            counts["pos_not_lifted"] += 1
            continue
        target_chrom, target_pos = lifted_start[0], lifted_start[1] + 1
        end = target_pos
        if svtype == "DEL":
            lifted_end = lift_point(lifter, source_chrom, record.pos - 1 + size)
            if lifted_end is None or lifted_end[0] != target_chrom:
                counts["end_not_lifted"] += 1
                continue
            lifted_span = lifted_end[1] + 1 - target_pos
            tolerance = max(LIFT_LENGTH_FLOOR, size * LIFT_LENGTH_TOLERANCE)
            if lifted_span <= 0 or abs(lifted_span - size) > tolerance:
                counts["span_changed"] += 1
                continue
            end = target_pos + size
        genotype = record.tail[1].split(":")[0] if len(record.tail) > 1 else "./."
        kept.append(common.Variant(
            chrom=target_chrom, pos=target_pos, identifier=record.identifier, ref="N",
            alt=f"<{svtype}>", qual=".", filter_field="PASS",
            info={"SVTYPE": svtype, "SVLEN": str(-size if svtype == "DEL" else size), "END": str(end)},
            tail=["GT", genotype]))
        counts["kept"] += 1
    return kept, counts


def lift_bed(intervals: list[tuple[str, int, int]], lifter) -> tuple[list[tuple[str, int, int]], dict[str, int]]:
    """Lift the high confidence regions to GRCh38.

    An interval keeps its place when both edges lift to the same contig and the lifted length
    stays inside the tolerance. The result is sorted and merged.

    Args:
        intervals: The source intervals as contig, start and end.
        lifter: A pyliftover LiftOver instance for the hg19 to hg38 chain.

    Returns:
        The lifted intervals and the counts of the input, the loss and the bases.
    """
    lifted: list[tuple[str, int, int]] = []
    counts = {"input_intervals": len(intervals), "input_bases": common.total_bases(intervals),
              "edge_not_lifted": 0, "span_changed": 0}
    for chrom, start, end in intervals:
        source_chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
        lifted_start = lift_point(lifter, source_chrom, start)
        lifted_end = lift_point(lifter, source_chrom, end)
        if lifted_start is None or lifted_end is None or lifted_start[0] != lifted_end[0]:
            counts["edge_not_lifted"] += 1
            continue
        span = lifted_end[1] - lifted_start[1]
        source_span = end - start
        tolerance = max(LIFT_LENGTH_FLOOR, source_span * LIFT_LENGTH_TOLERANCE)
        if span <= 0 or abs(span - source_span) > tolerance:
            counts["span_changed"] += 1
            continue
        lifted.append((lifted_start[0], lifted_start[1], lifted_end[1]))
    merged = common.merge_intervals(lifted)
    counts["kept_intervals"] = len(merged)
    counts["kept_bases"] = common.total_bases(merged)
    return merged, counts


def prepare_call_set(records: list[common.Variant]) -> tuple[list[common.Variant], dict[str, int]]:
    """Select the call records that the comparison can score.

    PASS records on a primary contig are kept. An inversion or a break end record is dropped,
    because the truth set holds no such class and a kept record would only add a false positive
    that no truth can confirm. The counts report every dropped group.

    Args:
        records: The records of the merged structural variant VCF.

    Returns:
        The kept records and the counts by reason.
    """
    kept: list[common.Variant] = []
    counts = {"input": len(records), "not_pass": 0, "off_primary_contig": 0, "no_truth_class": 0,
              "inversion_calls": 0, "break_end_calls": 0, "insertion_of_unknown_length": 0,
              "unsized": 0, "below_50bp": 0, "kept": 0}
    for record in records:
        if record.filter_field not in ("PASS", "."):
            counts["not_pass"] += 1
            continue
        if record.chrom not in common.PRIMARY_CONTIGS:
            counts["off_primary_contig"] += 1
            continue
        svtype = common.variant_type(record)
        if svtype in ("INV", "BND", "NA"):
            counts["no_truth_class"] += 1
            counts["inversion_calls"] += int(svtype == "INV")
            counts["break_end_calls"] += int(svtype == "BND")
            continue
        size = common.variant_size(record)
        if size is None:
            counts["unsized"] += 1
            continue
        if size < 50:
            if svtype == "INS" and "SVLEN" not in record.info:
                counts["insertion_of_unknown_length"] += 1
            counts["below_50bp"] += 1
            continue
        kept.append(record)
        counts["kept"] += 1
    return kept, counts


def count_by_class_and_band(path: Path) -> dict[str, dict[str, int]]:
    """Count the records of a Truvari output VCF by class and by size band.

    Args:
        path: The Truvari output VCF.

    Returns:
        A mapping from class to a mapping from band to a count, with a "total" band.
    """
    _, records = common.read_vcf(path)
    table: dict[str, dict[str, int]] = {}
    for record in records:
        svtype = common.variant_type(record)
        band = common.size_band(common.variant_size(record))
        row = table.setdefault(svtype, {})
        row[band] = row.get(band, 0) + 1
        row["total"] = row.get("total", 0) + 1
    return table


def metrics(tp_base: int, false_negative: int, tp_comp: int, false_positive: int) -> dict[str, float | int]:
    """Compute recall, precision and F1 from the four Truvari counts.

    Recall comes from the truth side. Precision comes from the call side. Truvari writes the two
    sides in separate files, so the two true positive counts can differ.

    Args:
        tp_base: The number of matched truth records.
        false_negative: The number of missed truth records.
        tp_comp: The number of matched calls.
        false_positive: The number of unmatched calls.

    Returns:
        The counts and the three measures. A measure with an empty denominator is None.
    """
    recall = tp_base / (tp_base + false_negative) if (tp_base + false_negative) else None
    precision = tp_comp / (tp_comp + false_positive) if (tp_comp + false_positive) else None
    f1 = None
    if recall is not None and precision is not None and (recall + precision) > 0:
        f1 = 2 * recall * precision / (recall + precision)
    return {"tp_base": tp_base, "fn": false_negative, "tp_comp": tp_comp, "fp": false_positive,
            "truth_records": tp_base + false_negative,
            "recall": None if recall is None else round(recall, 4),
            "precision": None if precision is None else round(precision, 4),
            "f1": None if f1 is None else round(f1, 4)}


def stratify_bench(bench_dir: Path) -> dict:
    """Read one Truvari output directory and split the result by class and by size band.

    Recall uses the class and the size of the truth record. Precision uses the class and the
    size of the call. A class that the truth set does not hold gets no recall value.

    Args:
        bench_dir: The Truvari output directory.

    Returns:
        The summary, the per class table and the per band table.
    """
    tp_base = count_by_class_and_band(bench_dir / "tp-base.vcf.gz")
    false_negative = count_by_class_and_band(bench_dir / "fn.vcf.gz")
    tp_comp = count_by_class_and_band(bench_dir / "tp-comp.vcf.gz")
    false_positive = count_by_class_and_band(bench_dir / "fp.vcf.gz")
    summary = json.loads((bench_dir / "summary.json").read_text())

    by_class: dict[str, dict] = {}
    for name in CALL_CLASSES:
        row = metrics(tp_base.get(name, {}).get("total", 0),
                      false_negative.get(name, {}).get("total", 0),
                      tp_comp.get(name, {}).get("total", 0),
                      false_positive.get(name, {}).get("total", 0))
        if row["truth_records"] == 0:
            row["recall"] = None
            row["f1"] = None
            row["note"] = "The truth set holds no record of this class."
        by_class[name] = row

    bands = [f"{low}-{high}" if high else f"{low}+" for low, high in common.SIZE_BANDS]
    by_band: dict[str, dict] = {}
    for band in bands:
        by_band[band] = metrics(sum(tp_base.get(name, {}).get(band, 0) for name in TRUTH_CLASSES),
                                sum(false_negative.get(name, {}).get(band, 0) for name in TRUTH_CLASSES),
                                sum(tp_comp.get(name, {}).get(band, 0) for name in CALL_CLASSES),
                                sum(false_positive.get(name, {}).get(band, 0) for name in CALL_CLASSES))

    by_class_band: dict[str, dict[str, dict]] = {}
    for name in TRUTH_CLASSES:
        by_class_band[name] = {band: metrics(tp_base.get(name, {}).get(band, 0),
                                             false_negative.get(name, {}).get(band, 0),
                                             tp_comp.get(name, {}).get(band, 0),
                                             false_positive.get(name, {}).get(band, 0))
                               for band in bands}
    return {"truvari_summary": summary, "by_class": by_class, "by_size_band": by_band,
            "by_class_and_band": by_class_band}


def positional_reach(records: list[common.Variant], sites: dict[str, list[int]],
                     bed_index: dict[str, list[tuple[int, int]]]) -> tuple[int, int, set[tuple[str, int]]]:
    """Match calls to truth insertions by position alone.

    Args:
        records: The call records to test.
        sites: The sorted truth insertion positions by contig.
        bed_index: The indexed high confidence regions.

    Returns:
        The number of calls inside the regions, the number of those that match a truth insertion,
        and the set of truth insertions that the calls reach.
    """
    inside = 0
    matched = 0
    reached: set[tuple[str, int]] = set()
    for record in records:
        if not common.position_in_index(bed_index, record.chrom, record.pos):
            continue
        inside += 1
        hit = nearest_site(sites.get(record.chrom, []), record.pos)
        if hit is not None:
            matched += 1
            reached.add((record.chrom, hit))
    return inside, matched, reached


def mei_positional_check(mei_records: list[common.Variant], truth_records: list[common.Variant],
                         bed_index: dict[str, list[tuple[int, int]]], tp_base_path: Path,
                         sv_insertion_calls: list[common.Variant]) -> dict:
    """Compare the RetroSeq mobile element calls with the truth insertions by position.

    The v0.6 truth set holds no mobile element class, so it gives no mobile element recall. It
    does hold every insertion of at least 50 base pairs inside the high confidence regions. A
    mobile element insertion is such an insertion. A RetroSeq call inside those regions that
    matches no truth insertion is therefore a false positive.

    Args:
        mei_records: The RetroSeq call records.
        truth_records: The lifted truth records.
        bed_index: The indexed high confidence regions on GRCh38.
        tp_base_path: The Truvari tp-base VCF of the structural variant comparison.
        sv_insertion_calls: The insertion calls of the structural variant call set.

    Returns:
        The counts and the insertion reach of RetroSeq, of the structural variant callers and of
        the two sets together. Every value of this mapping comes from the positional rule, except
        the Truvari value that the last two keys name.
    """
    insertion_sites: dict[str, list[int]] = {}
    for record in truth_records:
        if common.variant_type(record) == "INS" and common.position_in_index(bed_index, record.chrom, record.pos):
            insertion_sites.setdefault(record.chrom, []).append(record.pos)
    for sites in insertion_sites.values():
        sites.sort()
    truth_total = sum(len(sites) for sites in insertion_sites.values())

    mei_inside, mei_matched, mei_reached = positional_reach(mei_records, insertion_sites, bed_index)
    sv_inside, sv_matched, sv_reached = positional_reach(sv_insertion_calls, insertion_sites, bed_index)
    union = mei_reached | sv_reached
    _, matched_by_truvari = common.read_vcf(tp_base_path)
    truvari_keys = {(record.chrom, record.pos) for record in matched_by_truvari
                    if common.variant_type(record) == "INS"}
    return {
        "match_window_bp": MEI_MATCH_WINDOW,
        "truth_insertions_in_regions": truth_total,
        "mei_calls_total": len(mei_records),
        "mei_calls_in_high_confidence_regions": mei_inside,
        "mei_calls_matching_a_truth_insertion": mei_matched,
        "mei_positional_match_rate": round(mei_matched / mei_inside, 4) if mei_inside else None,
        "sv_insertion_calls_in_high_confidence_regions": sv_inside,
        "sv_insertion_calls_matching_a_truth_insertion": sv_matched,
        "sv_positional_match_rate": round(sv_matched / sv_inside, 4) if sv_inside else None,
        "truth_insertions_reached_by_retroseq": len(mei_reached),
        "truth_insertions_reached_by_sv_callers": len(sv_reached),
        "truth_insertions_reached_by_either": len(union),
        "insertion_reach_retroseq": round(len(mei_reached) / truth_total, 4) if truth_total else None,
        "insertion_reach_sv_callers": round(len(sv_reached) / truth_total, 4) if truth_total else None,
        "insertion_reach_together": round(len(union) / truth_total, 4) if truth_total else None,
        "truth_insertions_matched_by_truvari": len(truvari_keys),
        "insertion_recall_truvari_rule": round(len(truvari_keys) / truth_total, 4) if truth_total else None,
    }


def nearest_site(sites: list[int], position: int) -> int | None:
    """Find a sorted site inside the match window of one position.

    Args:
        sites: The sorted truth positions of one contig.
        position: The call position.

    Returns:
        The matched site, or None when no site is inside the window.
    """
    if not sites:
        return None
    low, high = 0, len(sites) - 1
    best = None
    while low <= high:
        middle = (low + high) // 2
        if abs(sites[middle] - position) <= MEI_MATCH_WINDOW:
            best = sites[middle]
            break
        if sites[middle] < position:
            low = middle + 1
        else:
            high = middle - 1
    return best


def gene_report(truth_records: list[common.Variant], bed_index: dict[str, list[tuple[int, int]]],
                tp_base_path: Path) -> dict:
    """Measure the truth coverage and the recall inside the seven MVA genes.

    Args:
        truth_records: The lifted truth records.
        bed_index: The indexed high confidence regions on GRCh38.
        tp_base_path: The Truvari tp-base VCF, which holds the matched truth records.

    Returns:
        A mapping from gene name to the region coverage, the truth record count and the recall.
    """
    _, matched = common.read_vcf(tp_base_path)
    matched_keys = {(record.chrom, record.pos, common.variant_type(record)) for record in matched}
    report: dict[str, dict] = {}
    for gene, (chrom, start, end) in MVA_GENES.items():
        span = end - start
        covered = common.intersect_bases(bed_index, chrom, start, end)
        inside = [record for record in truth_records
                  if record.chrom == chrom and start < record.pos <= end
                  and common.position_in_index(bed_index, record.chrom, record.pos)]
        found = sum(1 for record in inside
                    if (record.chrom, record.pos, common.variant_type(record)) in matched_keys)
        report[gene] = {
            "region": f"{chrom}:{start}-{end}",
            "gene_bases": span,
            "high_confidence_bases": covered,
            "high_confidence_fraction": round(covered / span, 4) if span else 0.0,
            "truth_sv_records": len(inside),
            "truth_sv_recalled": found,
            "recall": round(found / len(inside), 4) if inside else None,
            "truth_sv_detail": [
                {"position": record.pos, "type": common.variant_type(record),
                 "size": common.variant_size(record),
                 "found": (record.chrom, record.pos, common.variant_type(record)) in matched_keys}
                for record in inside],
        }
    return report


def build_inputs(work_dir: Path, sv_vcf: Path, mei_vcf: Path, truth_vcf: Path, truth_bed: Path,
                 chain: Path) -> dict:
    """Lift the truth set, prepare the call set and write every file that Truvari needs.

    Args:
        work_dir: The directory that holds the prepared files.
        sv_vcf: The merged structural variant call set on GRCh38.
        mei_vcf: The mobile element call set on GRCh38.
        truth_vcf: The GIAB v0.6 truth VCF on GRCh37.
        truth_bed: The GIAB v0.6 high confidence BED on GRCh37.
        chain: The UCSC hg19 to hg38 chain file.

    Returns:
        The prepared file names and the counts of every preparation step.
    """
    from pyliftover import LiftOver

    lifter = LiftOver(str(chain))
    sv_header, sv_records = common.read_vcf(sv_vcf)
    _, mei_records = common.read_vcf(mei_vcf)
    _, truth_records = common.read_vcf(truth_vcf)

    lifted, truth_counts = lift_truth_records(truth_records, lifter)
    order = {name: index for index, name in enumerate(common.PRIMARY_CONTIGS)}
    lifted.sort(key=lambda record: (order[record.chrom], record.pos))
    common.write_vcf(work_dir / "truth_v06_grch38.vcf", contig_header(sv_header, "HG002"), lifted)

    intervals, bed_counts = lift_bed(common.read_bed(truth_bed), lifter)
    intervals.sort(key=lambda item: (order[item[0]], item[1]))
    common.write_bed(work_dir / "truth_v06_grch38.bed", intervals)

    calls, call_counts = prepare_call_set(sv_records)
    common.write_vcf(work_dir / "calls_sv.vcf", sv_header, calls)

    return {"lifted_truth": lifted, "intervals": intervals, "mei_records": mei_records,
            "calls": calls, "truth_counts": truth_counts, "bed_counts": bed_counts,
            "call_counts": call_counts}


def run_bench(work_dir: Path, base: str, comp: str, include_bed: str | None, output: str,
              dup_to_ins: bool) -> Path:
    """Run one Truvari comparison inside the container.

    Args:
        work_dir: The host work directory.
        base: The truth VCF file name.
        comp: The call VCF file name.
        include_bed: The high confidence BED file name, or None for the whole genome.
        output: The output directory name.
        dup_to_ins: Treat a duplication call as an insertion.

    Returns:
        The host path of the output directory.
    """
    command = ["truvari", "bench", "-b", base, "-c", comp, "-o", output]
    command.extend(TRUVARI_ARGS)
    if include_bed:
        command.extend(["--includebed", include_bed])
    if dup_to_ins:
        command.append("--dup-to-ins")
    run_container(TRUVARI_IMAGE, command, work_dir)
    return work_dir / output


def main() -> int:
    """Prepare the inputs, run every comparison and write the receipt.

    Returns:
        Zero when the benchmark completes.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--sv-vcf", required=True, type=Path)
    parser.add_argument("--mei-vcf", required=True, type=Path)
    parser.add_argument("--truth-vcf", required=True, type=Path)
    parser.add_argument("--truth-bed", required=True, type=Path)
    parser.add_argument("--cmrg-vcf", type=Path)
    parser.add_argument("--cmrg-bed", type=Path)
    parser.add_argument("--chain", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    parser.add_argument("--table", type=Path)
    args = parser.parse_args()

    args.work_dir.mkdir(parents=True, exist_ok=True)
    prepared = build_inputs(args.work_dir, args.sv_vcf, args.mei_vcf, args.truth_vcf,
                            args.truth_bed, args.chain)
    truth_gz = compress_and_index("truth_v06_grch38.vcf", args.work_dir)
    calls_gz = compress_and_index("calls_sv.vcf", args.work_dir)

    as_called = run_bench(args.work_dir, truth_gz, calls_gz, "truth_v06_grch38.bed",
                          "bench_v06_as_called", dup_to_ins=False)
    dup_as_ins = run_bench(args.work_dir, truth_gz, calls_gz, "truth_v06_grch38.bed",
                           "bench_v06_dup_as_ins", dup_to_ins=True)

    bed_index = common.index_intervals(prepared["intervals"])
    results = {
        "v06_as_called": stratify_bench(as_called),
        "v06_dup_as_insertion": stratify_bench(dup_as_ins),
        "mobile_elements": mei_positional_check(
            prepared["mei_records"], prepared["lifted_truth"], bed_index,
            dup_as_ins / "tp-base.vcf.gz",
            [record for record in prepared["calls"] if common.variant_type(record) == "INS"]),
        "mva_genes": gene_report(prepared["lifted_truth"], bed_index, dup_as_ins / "tp-base.vcf.gz"),
    }

    if args.cmrg_vcf and args.cmrg_bed:
        cmrg_vcf_name = args.cmrg_vcf.name
        cmrg_bed_name = args.cmrg_bed.name
        for source in (args.cmrg_vcf, args.cmrg_bed, Path(str(args.cmrg_vcf) + ".tbi")):
            if source.exists():
                (args.work_dir / source.name).write_bytes(source.read_bytes())
        cmrg = run_bench(args.work_dir, cmrg_vcf_name, calls_gz, cmrg_bed_name, "bench_cmrg",
                         dup_to_ins=True)
        results["cmrg_v1_00"] = json.loads((cmrg / "summary.json").read_text())

    inputs = {}
    for path in (args.sv_vcf, args.mei_vcf, args.truth_vcf, args.truth_bed, args.cmrg_vcf,
                 args.cmrg_bed, args.chain):
        if path is not None and path.exists():
            entry = dict(SOURCES.get(path.name, {}))
            entry["sha256"] = sha256_of(path)
            inputs[path.name] = entry

    receipt = {
        "receipt_id": args.receipt.stem,
        "created_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "track": "t1_call",
        "data_class": "A",
        "source_run": "giab_hg002_calling_v1 (S2). No caller ran again. The caller versions and the "
                      "caller parameters are the ones of the patient run patient_calling_v1 (S5).",
        "callers": CALLERS,
        "containers": {"truvari": TRUVARI_IMAGE, "truvari_licence": "MIT",
                       "bcftools": BCFTOOLS_IMAGE, "liftover": "pyliftover 0.4.1"},
        "truvari_arguments": list(TRUVARI_ARGS),
        "mobile_element_rule": f"A RetroSeq call and a truth insertion are one event when their "
                               f"positions are inside {MEI_MATCH_WINDOW} base pairs.",
        "gene_coordinates": "Ensembl release 116, GRCh38.",
        "known_limits": KNOWN_LIMITS,
        "inputs": inputs,
        "preparation": {"truth": prepared["truth_counts"], "high_confidence_bed": prepared["bed_counts"],
                        "call_set": prepared["call_counts"]},
        "results": results,
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n")
    if args.table:
        write_table(args.table, results)
    print(json.dumps(results, indent=2))
    return 0


def write_table(path: Path, results: dict) -> None:
    """Write the class and size band results as a tab separated table.

    Args:
        path: The output path.
        results: The result mapping of the benchmark.
    """
    header = ("comparison", "stratum", "truth_records", "tp_base", "fn", "tp_comp", "fp",
              "recall", "precision", "f1")
    rows: list[tuple] = []
    for key in ("v06_as_called", "v06_dup_as_insertion"):
        for name, row in results[key]["by_class"].items():
            rows.append((key, f"class:{name}", row["truth_records"], row["tp_base"], row["fn"],
                         row["tp_comp"], row["fp"], row["recall"], row["precision"], row["f1"]))
        for name, row in results[key]["by_size_band"].items():
            rows.append((key, f"band:{name}", row["truth_records"], row["tp_base"], row["fn"],
                         row["tp_comp"], row["fp"], row["recall"], row["precision"], row["f1"]))
        for svtype, bands in results[key]["by_class_and_band"].items():
            for name, row in bands.items():
                rows.append((key, f"class:{svtype}|band:{name}", row["truth_records"], row["tp_base"],
                             row["fn"], row["tp_comp"], row["fp"], row["recall"], row["precision"],
                             row["f1"]))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join("" if item is None else str(item) for item in row) + "\n")


if __name__ == "__main__":
    sys.exit(main())
