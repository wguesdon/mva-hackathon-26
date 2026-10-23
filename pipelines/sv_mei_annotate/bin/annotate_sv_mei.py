#!/usr/bin/env python3
"""Annotate, filter and rank the patient structural and mobile element calls (queue item R1).

The input calls are class C. Run this script where class C data is allowed: an AWS Batch task in the
project account, or the owner machine. Do not run it in an agent session that is approved for class A
only.

The script writes three groups of files.

    class_c_sv_mei_annotated.tsv      Every input call with its annotation. Class C. It stays in the
                                      run results prefix and it is deleted with the data.
    class_b2_sv_mei_ranked.tsv        The ranked candidate calls with the evidence fields. Class B2.
    class_b2_sv_mei_ranked.md         The same table as Markdown, for a Bedrock model. Class B2.
    class_b2_submission_rows.csv      A proposed row set in the organizer format. Class B2.
    class_b2_submission_rows.md       The reason for each proposed row. Class B2.
    class_a_gates.json                Checksums, tool versions, thresholds and gate verdicts. Class A.
                                      It holds no count of patient calls and no patient coordinate.

The class A file is the only output that a class A session may read. Every gate in it is a process
gate, so a verdict states that a step ran and not what the step found.

`--genes` names the disease genes. Its default is the seven MVA genes with BUB1B first, so a run that
omits the option behaves as the shipped run behaves. The option exists so that the causal ranking
benchmark `scripts/track1_sv_causal_benchmark.py` can score a public pathogenic event against the
gene of that event.

Usage:
    annotate_sv_mei.py --sv-vcf mva_case_sv.vcf.gz --mei-vcf mva_case_mobile_elements.vcf.gz \\
        --gtf gencode.v50.basic.annotation.gtf.gz --gnomad-sv gnomad.v4.1.sv.sites.vcf.gz \\
        --intervals intervals.bed --out-dir out
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

# The structural annotation set of literature/track1_gene_panel_v1.yaml, which is the frozen gene
# panel of this work. That file holds the inclusion rule, the ontology release it comes from and
# the reason for each manual addition. harness/test_track1_gene_panel.py fails when this tuple
# no longer equals the operating set that the panel derives. A Nextflow task runs this file in a
# container that holds no copy of literature/, so the tuple stays here and the test holds it.
# BUB1B, CEP57 and TRIP13 are the established disease genes of panel version v1. The other four are
# spindle assembly checkpoint partners that the queue item names. The first name is the primary gene.
# `--genes` replaces this tuple and its default is this tuple, so the shipped behaviour is unchanged.
#
# PANEL_SETS holds one structural annotation set per panel version. Version v2 adds CENATAC and
# MAD1L1, which the established clause of v1 excluded through an enumerated disease list, and
# MAD2L1BP, which silences the spindle assembly checkpoint. The shipped default stays at v1, because
# the run that produced the submission used the v1 set and queue row R30 blocks every patient rerun.
PANEL_SETS = {
    "v1": ("BUB1B", "CEP57", "TRIP13", "BUB1", "BUB3", "CDC20", "PTPA"),
    "v2": ("BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1", "MAD2L1BP", "BUB1", "BUB3", "CDC20",
           "PTPA"),
}
SHIPPED_PANEL_VERSION = "v1"
MVA_GENES = PANEL_SETS[SHIPPED_PANEL_VERSION]
PRIMARY_GENE = MVA_GENES[0]

# Filter thresholds. Every value is a class A constant and it is copied into the class A gate file.
MIN_SV_LENGTH = 50
MIN_GENOTYPE_QUALITY = 20
MIN_RETROSEQ_LEVEL = 6
MAX_GNOMAD_AF = 0.01
STRICT_GNOMAD_AF = 0.001
RECIPROCAL_OVERLAP = 0.5
BREAKPOINT_WINDOW = 100
REGULATORY_WINDOW = 50000
KNOWN_INTERVAL_PAD = 1000
CANDIDATE_SCORE_FLOOR = 25
SUBMISSION_SCORE_FLOOR = 60
MAX_CANDIDATE_ROWS = 25
# The ten-row output budget reserves at most three slots for small variants and at most seven
# slots for this writer. These are configuration limits, not counts from an individual run.
MAX_SUBMISSION_ROWS = 7
# Two eligible calls of one score map onto one EPCR, and the organizer scorer builds one threshold
# per distinct value. A tied row therefore falls by this step until its value is free. Seven rows
# need at most six steps, so no row leaves the interval (0, 1].
EPCR_TIE_STEP = 0.001
# The largest EPCR that this writer can print. score_call is a sum of five component terms and its
# largest value is 107, which maps onto 0.835, so the ceiling of 0.95 in the map is out of reach.
# The small-variant score schedule is configured independently. The merge check verifies its
# ordering relative to structural scores for each caller-supplied file.
EPCR_CEILING_OF_THIS_WRITER = 0.835
CONTEXT_MIN_LENGTH = 100000
CONTEXT_ROWS = 10

# Population labels of gnomAD SV v4.1 that carry an AF_<label> field.
GNOMAD_POPULATIONS = ("afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "sas")

# A tool name and a tool version are class B1, which ADR 0005 treats as class A. A command line can
# carry a private identifier, so read_tool_header drops it.
HEADER_LINE = re.compile(r"^##([A-Za-z0-9_]+)=(.*)$")
SYMBOLIC_ALT = re.compile(r"^<([A-Za-z0-9:_]+)>$")
BREAKEND_ALT = re.compile(r"[\[\]]")

CANDIDATE_FIELDS = ("rank", "score", "call_id", "source", "chrom", "pos", "end", "svtype", "svlen",
                    "genotype", "filter", "quality", "genotype_quality", "retroseq_level",
                    "caller_arm", "mva_gene", "gene_effect", "genes_overlapped", "gnomad_af",
                    "gnomad_popmax_af", "gnomad_match", "known_status", "score_components",
                    "alt_short", "ref_short")

CLASS_C_FIELDS = CANDIDATE_FIELDS + ("info", "format", "sample")

SUBMISSION_FIELDS = ("proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2",
                     "alt_2", "epcr", "finding_type", "notes")
PROBAND_ID = "PROBAND01"


def sha256_of(path: Path) -> str:
    """Compute the SHA-256 checksum of a file.

    Args:
        path: The file to read.

    Returns:
        The checksum in lower case hexadecimal.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path):
    """Open a plain or a gzip compressed text file.

    Args:
        path: The file to open.

    Returns:
        A file object in text mode.
    """
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def deny_tokens(raw: str) -> list[str]:
    """Split the private identifier list of the environment.

    Args:
        raw: The comma separated value of PRIVACY_DENY_TOKENS.

    Returns:
        One token for each non-empty entry.
    """
    return [token.strip() for token in raw.split(",") if token.strip()]


def redact(text: str, tokens: list[str]) -> str:
    """Replace each private identifier by a marker.

    Args:
        text: The text to clean.
        tokens: The private identifiers.

    Returns:
        The text with every identifier replaced.
    """
    for token in tokens:
        text = re.sub(re.escape(token), "[REDACTED]", text, flags=re.IGNORECASE)
    return text


def read_tool_header(path: Path, tokens: list[str]) -> list[str]:
    """Collect the tool name and the tool version lines of a VCF header.

    The function reads header lines only. It stops at the column header, so no data line is read. A
    command line entry is dropped, because it can carry a private identifier.

    Args:
        path: The VCF file.
        tokens: The private identifiers to remove.

    Returns:
        At most 40 header lines, in file order.
    """
    lines = []
    with open_text(path) as handle:
        for line in handle:
            if not line.startswith("##"):
                break
            match = HEADER_LINE.match(line.rstrip("\n"))
            if not match:
                continue
            key, value = match.group(1), match.group(2)
            lowered = key.lower()
            if lowered in ("cmdline", "commandline", "bcftools_viewcommand", "bcftools_mergecommand"):
                continue
            if lowered in ("source", "sourceversion", "fileformat", "reference") or "version" in lowered:
                lines.append(redact(f"##{key}={value}", tokens)[:300])
            if len(lines) >= 40:
                break
    return lines


def parse_info(field: str) -> dict[str, str]:
    """Split a VCF INFO column.

    Args:
        field: The INFO column.

    Returns:
        One entry for each key. A flag maps to an empty string.
    """
    info = {}
    if field in (".", ""):
        return info
    for item in field.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = ""
    return info


def to_int(value: str | None) -> int | None:
    """Convert a VCF field to an integer.

    Args:
        value: The raw field, or None.

    Returns:
        The integer, or None when the field is absent or not a number.
    """
    if value is None:
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def to_float(value: str | None) -> float | None:
    """Convert a VCF field to a float.

    Args:
        value: The raw field, or None.

    Returns:
        The float, or None when the field is absent or not a number.
    """
    if value is None:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def call_svtype(alt: str, info: dict[str, str]) -> str:
    """Decide the structural variant class of one record.

    Args:
        alt: The ALT column.
        info: The parsed INFO column.

    Returns:
        A class name such as DEL, DUP, INS, INV, CNV, BND or OTHER.
    """
    declared = info.get("SVTYPE", "").upper()
    if declared:
        return declared.split(":")[0]
    symbolic = SYMBOLIC_ALT.match(alt)
    if symbolic:
        return symbolic.group(1).upper().split(":")[0]
    if BREAKEND_ALT.search(alt):
        return "BND"
    return "OTHER"


def call_span(pos: int, ref: str, info: dict[str, str]) -> tuple[int, int | None]:
    """Decide the end coordinate and the length of one record.

    Args:
        pos: The POS column.
        ref: The REF column.
        info: The parsed INFO column.

    Returns:
        The end coordinate and the length, where the length is None for an unbounded record.
    """
    end = to_int(info.get("END"))
    length = to_int(info.get("SVLEN"))
    if length is not None:
        length = abs(length)
    if end is None and length is not None:
        end = pos + length
    if end is None:
        end = pos + max(len(ref) - 1, 0)
    if length is None:
        length = max(end - pos, 0)
    return end, length


def genotype_state(genotype: str) -> str:
    """Classify a GT field.

    Args:
        genotype: The GT field.

    Returns:
        One of hom_ref, het, hom_alt or no_call.
    """
    alleles = genotype.replace("|", "/").split("/")
    if not alleles or any(allele == "." for allele in alleles):
        return "no_call"
    if all(allele == "0" for allele in alleles):
        return "hom_ref"
    if all(allele != "0" for allele in alleles):
        return "hom_alt"
    return "het"


def read_calls(path: Path, source: str) -> list[dict]:
    """Read one patient VCF into records.

    Args:
        path: The VCF file, gzip compressed or plain.
        source: A label such as `sv` or `mei`.

    Returns:
        One record for each data line.
    """
    records = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            info = parse_info(fields[7])
            pos = int(fields[1])
            end, length = call_span(pos, fields[3], info)
            sample = {}
            if len(fields) >= 10:
                sample = dict(zip(fields[8].split(":"), fields[9].split(":")))
            records.append({
                "source": source,
                "call_id": fields[2] if fields[2] not in (".", "") else f"{source}:{fields[0]}:{pos}",
                "chrom": fields[0], "pos": pos, "end": end,
                # The full alleles reach a submission row, because the organizer scorer compares the
                # tuple of chromosome, position, reference allele and alternate allele and a
                # truncated allele can match nothing. The ranked table prints the short forms.
                "ref": fields[3], "alt": fields[4],
                "ref_short": fields[3][:40], "alt_short": fields[4][:60],
                "quality": fields[5], "filter": fields[6],
                "svtype": call_svtype(fields[4], info), "svlen": length,
                "info": fields[7][:400], "format": fields[8] if len(fields) >= 9 else "",
                "sample_fields": sample,
                "genotype": genotype_state(sample.get("GT", "./.")),
                "genotype_quality": to_int(sample.get("GQ")),
                "retroseq_level": to_int(sample.get("FL")),
                "caller_arm": info.get("set", "") or info.get("caller_set", ""),
            })
    return records


def gene_index(genes: dict[str, list[tuple]]) -> dict[str, tuple[list[tuple], list[int], int]]:
    """Build a search index over the gene intervals of each chromosome.

    Args:
        genes: The sorted gene intervals by chromosome.

    Returns:
        A map from chromosome to the interval list, the start coordinates and the longest gene.
    """
    index = {}
    for chrom, records in genes.items():
        starts = [record[0] for record in records]
        longest = max((record[1] - record[0] for record in records), default=0)
        index[chrom] = (records, starts, longest)
    return index


def gene_names(raw: str) -> tuple[str, ...]:
    """Split the disease gene option into gene names.

    Args:
        raw: The comma separated value of `--genes`.

    Returns:
        One name for each non-empty entry, in the order of the option. The first name is the primary
        gene. An empty option returns the shipped MVA gene tuple. The case of a name is kept, because
        `read_gene_models` matches the `gene_name` attribute of the GENCODE annotation exactly and a
        GENCODE name such as `C11orf65` is not upper case.
    """
    names = tuple(dict.fromkeys(item.strip() for item in raw.split(",") if item.strip()))
    return names or MVA_GENES


def read_gene_models(path: Path, genes_wanted: tuple[str, ...] = MVA_GENES
                     ) -> tuple[dict[str, list[tuple]], dict[str, dict]]:
    """Read gene intervals and the MANE exon structure of the disease genes from a GENCODE GTF.

    Args:
        path: The GENCODE annotation file, gzip compressed.
        genes_wanted: The disease gene names. The default is the shipped MVA gene tuple.

    Returns:
        A map from chromosome to a sorted list of gene intervals, and a map from disease gene name to
        its record with the exon list of its MANE Select transcript.
    """
    genes: dict[str, list[tuple]] = {}
    wanted = set(genes_wanted)
    mva: dict[str, dict] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("gene", "exon"):
                continue
            attributes = fields[8]
            name_match = re.search(r'gene_name "([^"]+)"', attributes)
            if not name_match:
                continue
            name = name_match.group(1)
            chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            if fields[2] == "gene":
                type_match = re.search(r'gene_type "([^"]+)"', attributes)
                genes.setdefault(chrom, []).append((start, end, name, type_match.group(1) if type_match else "."))
                if name in wanted:
                    mva[name] = {"chrom": chrom, "start": start, "end": end, "strand": strand,
                                 "gene_type": type_match.group(1) if type_match else ".", "exons": []}
            elif name in wanted and 'tag "MANE_Select"' in attributes:
                mva.setdefault(name, {"chrom": chrom, "start": start, "end": end, "strand": strand,
                                      "gene_type": ".", "exons": []})
                mva[name]["exons"].append((start, end))
    for chrom in genes:
        genes[chrom].sort()
    for record in mva.values():
        record["exons"].sort()
    return genes, mva


def overlapping_genes(index: dict[str, tuple[list[tuple], list[int], int]], chrom: str,
                      start: int, end: int) -> list[str]:
    """List the gene names that a call interval overlaps.

    Args:
        index: The gene search index.
        chrom: The chromosome of the call.
        start: The start coordinate.
        end: The end coordinate.

    Returns:
        The gene names, protein coding genes first, at most 12 names.
    """
    entry = index.get(chrom)
    if entry is None:
        return []
    records, starts, longest = entry
    first = bisect.bisect_left(starts, start - longest)
    found = []
    for gene_start, gene_end, name, gene_type in records[first:]:
        if gene_start > end:
            break
        if gene_end >= start:
            found.append((0 if gene_type == "protein_coding" else 1, name))
    found.sort()
    return [name for _, name in found[:12]]


def hits_exon(record: dict, start: int, end: int) -> bool:
    """Say whether a call interval touches a MANE exon of one gene.

    Args:
        record: The MVA gene record with its exon list.
        start: The start coordinate of the call.
        end: The end coordinate of the call.

    Returns:
        True when the call overlaps at least one exon.
    """
    for exon_start, exon_end in record["exons"]:
        if exon_start > end:
            return False
        if exon_end >= start:
            return True
    return False


def transcription_start(record: dict) -> int:
    """Return the transcription start coordinate of one gene.

    Args:
        record: The MVA gene record.

    Returns:
        The gene start on the plus strand and the gene end on the minus strand.
    """
    return record["start"] if record["strand"] == "+" else record["end"]


def gene_effect(record: dict, call: dict) -> str:
    """Classify the effect of one call on one MVA gene.

    Args:
        record: The MVA gene record.
        call: The patient call.

    Returns:
        An effect label such as whole_gene, exonic, intronic, upstream_regulatory, downstream_flank or
        breakend_in_gene.
    """
    start, end = call["pos"], call["end"]
    if call["svtype"] == "BND" and record["start"] <= start <= record["end"]:
        return "breakend_in_gene"
    if start <= record["start"] and end >= record["end"]:
        return "whole_gene"
    if start <= record["end"] and end >= record["start"]:
        return "exonic" if hits_exon(record, start, end) else "intronic"
    site = transcription_start(record)
    if record["strand"] == "+":
        if site - REGULATORY_WINDOW <= end < record["start"]:
            return "upstream_regulatory"
    elif record["end"] < start <= site + REGULATORY_WINDOW:
        return "upstream_regulatory"
    if record["start"] - REGULATORY_WINDOW <= end and start <= record["end"] + REGULATORY_WINDOW:
        return "downstream_flank"
    return "outside"


def read_gnomad_sv(path: Path) -> tuple[dict[str, list[tuple]], int]:
    """Read the gnomAD SV sites file into per chromosome frequency intervals.

    Args:
        path: The gnomAD SV sites VCF, gzip compressed.

    Returns:
        A map from chromosome to a sorted list of records, and the number of records read.
    """
    table: dict[str, list[tuple]] = {}
    total = 0
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 8)
            if len(fields) < 8:
                continue
            info = parse_info(fields[7])
            frequency = to_float(info.get("AF", "").split(",")[0]) or 0.0
            popmax = frequency
            for population in GNOMAD_POPULATIONS:
                value = to_float(info.get(f"AF_{population}", "").split(",")[0])
                if value is not None and value > popmax:
                    popmax = value
            pos = int(fields[1])
            end, _ = call_span(pos, fields[3], info)
            table.setdefault(fields[0], []).append(
                (pos, end, info.get("SVTYPE", "OTHER").upper(), frequency, popmax, fields[2]))
            total += 1
    for chrom in table:
        table[chrom].sort()
    return table, total


def compatible_type(patient: str, population: str) -> bool:
    """Say whether two structural variant classes may describe the same event.

    Args:
        patient: The class of the patient call.
        population: The class of the population record.

    Returns:
        True when the two classes belong to the same group.
    """
    gain = {"DUP", "CNV", "MCNV"}
    loss = {"DEL", "CNV", "MCNV"}
    insertion = {"INS", "MEI", "ALU", "LINE1", "SVA", "INS:ME"}
    if patient in gain and population in gain:
        return True
    if patient in loss and population in loss:
        return True
    if patient in insertion and population in insertion:
        return True
    return patient == population


def match_gnomad(table: dict[str, list[tuple]], call: dict) -> tuple[float | None, float | None, str]:
    """Find the population frequency of one patient call.

    A balanced or bounded event matches on reciprocal overlap. An insertion or a breakend matches on
    breakpoint distance, because its length is not comparable between call sets.

    Args:
        table: The gnomAD records by chromosome.
        call: The patient call.

    Returns:
        The allele frequency, the population maximum frequency and the identifier of the match. The
        three values are None, None and an empty string when no record matches.
    """
    records = table.get(call["chrom"])
    if not records:
        return None, None, ""
    start, end = call["pos"], call["end"]
    positional = call["svtype"] in ("INS", "BND", "MEI", "ALU", "LINE1", "SVA") or (end - start) < MIN_SV_LENGTH
    window = BREAKPOINT_WINDOW if positional else max(end - start, MIN_SV_LENGTH)
    # A record that reaches the reciprocal overlap floor is at most twice the length of the call, so a
    # back-off of twice that length finds every possible match.
    back = max(2 * (end - start), window) + 10
    left = bisect.bisect_left(records, (start - back,))
    best = (None, None, "")
    best_frequency = -1.0
    for record in records[left:]:
        if record[0] > end + window:
            break
        if not compatible_type(call["svtype"], record[2]):
            continue
        if positional:
            if abs(record[0] - start) > BREAKPOINT_WINDOW:
                continue
        else:
            shared = min(end, record[1]) - max(start, record[0])
            if shared <= 0:
                continue
            if shared < RECIPROCAL_OVERLAP * max(end - start, 1) or shared < RECIPROCAL_OVERLAP * max(record[1] - record[0], 1):
                continue
        if record[4] > best_frequency:
            best_frequency = record[4]
            best = (record[3], record[4], record[5])
    return best


def read_intervals(path: Path) -> list[tuple[str, int, int, str]]:
    """Read the known BUB1B intervals.

    Args:
        path: A BED file.

    Returns:
        One entry for each interval.
    """
    intervals = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        intervals.append((fields[0], int(fields[1]), int(fields[2]), fields[3] if len(fields) > 3 else "."))
    return intervals


def known_status(call: dict, intervals: list[tuple[str, int, int, str]]) -> str:
    """Say whether a call falls in an interval that the A5 work already covered.

    Args:
        call: The patient call.
        intervals: The known intervals.

    Returns:
        The interval name when the call touches one, otherwise the word `new`.
    """
    for chrom, start, end, name in intervals:
        if call["chrom"] != chrom:
            continue
        if call["pos"] <= end + KNOWN_INTERVAL_PAD and call["end"] >= start - KNOWN_INTERVAL_PAD:
            return f"known_interval:{name}"
    return "new"


def passes_quality(call: dict) -> bool:
    """Apply the caller quality filter to one call.

    Args:
        call: The patient call.

    Returns:
        True when the call meets every quality threshold.
    """
    if call["genotype"] in ("hom_ref", "no_call"):
        return False
    if call["svtype"] != "BND" and (call["svlen"] or 0) < MIN_SV_LENGTH and call["source"] == "sv":
        return False
    quality = call["genotype_quality"]
    if quality is not None and quality < MIN_GENOTYPE_QUALITY:
        return False
    level = call["retroseq_level"]
    if level is not None and level < MIN_RETROSEQ_LEVEL:
        return False
    return True


def score_call(call: dict, primary: str = PRIMARY_GENE) -> tuple[int, str]:
    """Score how well one call explains a recessive second hit in the primary gene.

    Args:
        call: The annotated patient call.
        primary: The primary disease gene. The default is the shipped primary gene BUB1B.

    Returns:
        The score and a text that names each component.
    """
    parts = []
    total = 0
    gene = call["mva_gene"]
    if gene == primary:
        total += 50
        parts.append(f"gene_{primary.lower()}=50")
    elif gene:
        total += 20
        parts.append(f"gene_{gene.lower()}=20")
    effect_points = {"whole_gene": 25, "exonic": 25, "breakend_in_gene": 20, "intronic": 12,
                     "upstream_regulatory": 15, "downstream_flank": 5}
    effect = call["gene_effect"]
    if effect == "intronic" and call["source"] == "mei":
        total += 20
        parts.append("effect_intronic_mei=20")
    elif effect in effect_points:
        total += effect_points[effect]
        parts.append(f"effect_{effect}={effect_points[effect]}")
    frequency = call["gnomad_popmax_af"]
    if frequency is None:
        total += 15
        parts.append("absent_from_gnomad=15")
    elif frequency <= STRICT_GNOMAD_AF:
        total += 12
        parts.append("gnomad_af_very_rare=12")
    elif frequency <= MAX_GNOMAD_AF:
        total += 6
        parts.append("gnomad_af_rare=6")
    if call["genotype"] == "hom_alt":
        total += 10
        parts.append("genotype_hom_alt=10")
    elif call["genotype"] == "het":
        total += 8
        parts.append("genotype_het=8")
    if call["filter"] in ("PASS", "."):
        total += 7
        parts.append("filter_pass=7")
    return total, ",".join(parts)


def annotate(calls: list[dict], index: dict[str, tuple[list[tuple], list[int], int]], mva: dict[str, dict],
             gnomad: dict[str, list[tuple]], intervals: list[tuple[str, int, int, str]],
             genes_wanted: tuple[str, ...] = MVA_GENES) -> list[dict]:
    """Add the gene, frequency, known and score fields to every call.

    Args:
        calls: The patient calls.
        index: The gene search index.
        mva: The disease gene records.
        gnomad: The gnomAD SV records by chromosome.
        intervals: The known BUB1B intervals.
        genes_wanted: The disease gene names. The first name is the primary gene. The default is the
            shipped MVA gene tuple.

    Returns:
        The same records with the annotation fields added.
    """
    primary = genes_wanted[0]
    for call in calls:
        call["genes_overlapped"] = ";".join(overlapping_genes(index, call["chrom"], call["pos"], call["end"]))
        best_gene, best_effect = "", "outside"
        for name in genes_wanted:
            record = mva.get(name)
            if record is None or record["chrom"] != call["chrom"]:
                continue
            effect = gene_effect(record, call)
            if effect == "outside":
                continue
            if not best_gene or name == primary:
                best_gene, best_effect = name, effect
                if name == primary:
                    break
        call["mva_gene"] = best_gene
        call["gene_effect"] = best_effect
        frequency, popmax, match = match_gnomad(gnomad, call)
        call["gnomad_af"] = frequency
        call["gnomad_popmax_af"] = popmax
        call["gnomad_match"] = match
        call["known_status"] = known_status(call, intervals)
        call["quality_pass"] = passes_quality(call)
        call["score"], call["score_components"] = score_call(call, primary)
    return calls


def select_candidates(calls: list[dict], limit: int | None = MAX_CANDIDATE_ROWS) -> list[dict]:
    """Apply the filters and order the surviving calls.

    Args:
        calls: The annotated calls.
        limit: The largest number of rows to return. None returns every row that passes the
            filters, which the R13 row budget count needs.

    Returns:
        The candidate calls, highest score first, at most `limit` entries.
    """
    kept = []
    for call in calls:
        if not call["quality_pass"]:
            continue
        if not call["mva_gene"]:
            continue
        frequency = call["gnomad_popmax_af"]
        if frequency is not None and frequency > MAX_GNOMAD_AF:
            continue
        if call["score"] < CANDIDATE_SCORE_FLOOR:
            continue
        kept.append(call)
    kept.sort(key=lambda call: (-call["score"], call["chrom"], call["pos"]))
    for position, call in enumerate(kept, start=1):
        call["rank"] = position
    return kept if limit is None else kept[:limit]


def select_context(calls: list[dict]) -> list[dict]:
    """Select the largest rare events outside the MVA genes.

    The ranked table answers one question about one gene. This second list protects the reading of
    that table, because a large rare event elsewhere changes what the child's genome says.

    Args:
        calls: The annotated calls.

    Returns:
        At most ten calls, longest first.
    """
    kept = []
    for call in calls:
        if call["mva_gene"] or not call["quality_pass"]:
            continue
        if call["filter"] not in ("PASS", "."):
            continue
        frequency = call["gnomad_popmax_af"]
        if frequency is not None and frequency > MAX_GNOMAD_AF:
            continue
        if (call["svlen"] or 0) < CONTEXT_MIN_LENGTH:
            continue
        kept.append(call)
    kept.sort(key=lambda call: -(call["svlen"] or 0))
    return kept[:CONTEXT_ROWS]


def row_for(call: dict, fields: tuple[str, ...]) -> dict:
    """Build one output row.

    Args:
        call: The annotated call.
        fields: The column names.

    Returns:
        One value for each column, with None written as an empty string.
    """
    row = {}
    for field in fields:
        if field == "sample":
            value = ":".join(f"{key}={item}" for key, item in call.get("sample_fields", {}).items())
        else:
            value = call.get(field, "")
        row[field] = "" if value is None else value
    return row


def write_table(path: Path, rows: list[dict], fields: tuple[str, ...]) -> None:
    """Write one tab separated table.

    Args:
        path: The output file.
        rows: The rows to write.
        fields: The column names.
    """
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def candidate_line(call: dict) -> str:
    """Render one candidate as a Markdown table row.

    Args:
        call: The annotated call.

    Returns:
        The table row.
    """
    frequency = "absent" if call["gnomad_popmax_af"] is None else f"{call['gnomad_popmax_af']:.5f}"
    return (f"| {call.get('rank', '-')} | {call['score']} | {call['call_id']} | {call['svtype']} | "
            f"{call['chrom']}:{call['pos']}-{call['end']} | {call['svlen']} | {call['genotype']} | "
            f"{call['filter']} | {call.get('caller_arm') or '-'} | {call['mva_gene'] or '-'} | "
            f"{call['gene_effect']} | {frequency} | "
            f"{call['known_status']} | {call['score_components'] or '-'} |")


def markdown_table(candidates: list[dict], context: list[dict], primary: str = PRIMARY_GENE) -> str:
    """Render the candidate table as Markdown for a Bedrock model.

    Args:
        candidates: The ranked candidates in the MVA genes.
        context: The largest rare events outside the MVA genes.
        primary: The primary disease gene. The default is the shipped primary gene BUB1B.

    Returns:
        The Markdown text.
    """
    columns = ("| rank | score | call | type | position | length | genotype | filter | arm | gene | "
               "effect | gnomAD popmax AF | known | components |")
    divider = "|" + "---|" * 14
    lines = ["# Ranked structural, copy number and mobile element candidates", "",
             "Data class B2. This table describes one child. Amazon Bedrock only.", "",
             "## 1. Candidates in the MVA genes", "",
             "The rank orders the calls by how well each one explains a recessive second hit in "
             f"{primary}. A call marked `known_interval` falls in one of the two {primary} "
             "intron intervals that the A5 work already covered, so it is not a new finding.", "",
             columns, divider]
    lines += [candidate_line(call) for call in candidates]
    if not candidates:
        lines.append("| - | - | no call passed the filters | - | - | - | - | - | - | - | - | - | - | - |")
    lines += ["", "## 2. The largest rare events outside the MVA genes", "",
              f"This list holds the longest calls above {CONTEXT_MIN_LENGTH} bp that pass the caller "
              "filter and are rare in the population. It is context for section 1 and it is not a "
              "ranking for the diagnosis.", "", columns, divider]
    lines += [candidate_line(call) for call in context]
    if not context:
        lines.append("| - | - | no call passed the filters | - | - | - | - | - | - | - | - | - | - | - |")
    lines += ["", "## 3. Filter thresholds", "",
              f"- The minimum structural variant length is {MIN_SV_LENGTH} bp.",
              f"- The minimum genotype quality is {MIN_GENOTYPE_QUALITY} where the caller reports one.",
              f"- The minimum RetroSeq call level is {MIN_RETROSEQ_LEVEL} where the caller reports one.",
              f"- The maximum gnomAD SV population frequency is {MAX_GNOMAD_AF}.",
              f"- A bounded call matches a gnomAD record at {RECIPROCAL_OVERLAP} reciprocal overlap.",
              f"- An insertion or a breakend matches a gnomAD record within {BREAKPOINT_WINDOW} bp.",
              f"- The regulatory window reaches {REGULATORY_WINDOW} bp upstream of the transcription "
              "start site.",
              f"- A call scores at least {CANDIDATE_SCORE_FLOOR} to enter section 1.",
              "- Section 1 keeps a call that the caller filter marks as failed, and the `filter` "
              "column shows the label. A caller filter must not hide the second allele, and the "
              "score gives the failed call no filter points.",
              "- A call with a homozygous reference genotype or with no call is dropped.", ""]
    return "\n".join(lines)


def epcr_of_score(score: int) -> float:
    """Map one call score onto the EPCR that this writer prints.

    The map is one to one on an integer score, and the largest score that `score_call` can return is
    107, so the largest EPCR that this writer can print is 0.835 and the ceiling of 0.95 is out of
    reach. `tracking/receipts/t1_epcr_ladder_v1.json` enumerates the reachable set.

    Args:
        score: The call score.

    Returns:
        The EPCR value, rounded to three digits.
    """
    return round(min(0.95, 0.30 + score / 200.0), 3)


def submission_rows(candidates: list[dict], limit: int | None = MAX_SUBMISSION_ROWS) -> list[dict]:
    """Build the proposed submission rows in the organizer format.

    A structural call has no small variant representation, so the row carries the VCF coordinate and
    the symbolic alternate allele. The owner decides whether the row belongs in the file.

    The organizer scorer builds one threshold per distinct EPCR, so two rows that share one value
    share one threshold. `tracking/receipts/s4_row_count_calibration_v1.json`, field
    `tie_counterfactual`, measures that cost: one shared value across ten rows drops mean F-max from
    0.960 to 0.238 over the 50 spike-in scenarios. Two eligible calls of one score map onto one
    EPCR, so this function lowers a tied row by one rounding unit of 0.001 until the value is free.
    The candidates arrive sorted by descending score, so the adjustment keeps the rank order and it
    never raises a row above the row before it.

    Args:
        candidates: The ranked candidates.
        limit: The largest number of rows to return. None returns every eligible row, which the
            R13 row budget count needs.

    Returns:
        One row for each call above the submission floor, with strictly decreasing distinct EPCR
        values.
    """
    rows = []
    used: set[float] = set()
    for call in candidates:
        if call["score"] < SUBMISSION_SCORE_FLOOR:
            continue
        if call["known_status"] != "new":
            continue
        epcr = epcr_of_score(call["score"])
        if rows:
            epcr = min(epcr, round(rows[-1]["epcr"] - EPCR_TIE_STEP, 3))
        while epcr in used:
            epcr = round(epcr - EPCR_TIE_STEP, 3)
        used.add(epcr)
        rows.append({
            "proband_id": PROBAND_ID,
            "chrom_1": call["chrom"], "pos_1": call["pos"], "ref_1": call["ref"], "alt_1": call["alt"],
            "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
            "epcr": epcr, "finding_type": "secondary",
            "notes": (f"{call['svtype']} {call['chrom']}:{call['pos']}-{call['end']}, "
                      f"{call['svlen']} bp, {call['genotype']}, {call['mva_gene']} "
                      f"{call['gene_effect']}, score {call['score']}"),
        })
        if limit is not None and len(rows) >= limit:
            break
    return rows


def submission_note(rows: list[dict], candidates: list[dict], primary: str = PRIMARY_GENE) -> str:
    """Write the reason for each proposed row.

    Args:
        rows: The proposed rows.
        candidates: The ranked candidates.
        primary: The primary disease gene. The default is the shipped primary gene BUB1B.

    Returns:
        The Markdown text.
    """
    lines = ["# Proposed Track 1 rows from the structural and mobile element analysis", "",
             "Data class B2. The owner reads this file and decides. Nothing here is uploaded, and "
             "`slot1_b0.csv` is unchanged.", "",
             f"A call enters this file when its score is at least {SUBMISSION_SCORE_FLOOR} and it does "
             f"not fall in a known {primary} interval. The file is written on every run, so its "
             "existence states nothing about the child.", ""]
    if not rows:
        lines += ["## Verdict", "",
                  "No call reached the submission floor. The current `slot1_b0.csv` needs no change "
                  "from this analysis.", ""]
        return "\n".join(lines)
    lines += ["## Rows", "", "| epcr | chrom | pos | ref | alt | reason |", "|---|---|---|---|---|---|"]
    by_id = {(call["chrom"], call["pos"]): call for call in candidates}
    for row in rows:
        call = by_id.get((row["chrom_1"], row["pos_1"]), {})
        reason = (f"{call.get('svtype', '')} in {call.get('mva_gene', '')}, "
                  f"{call.get('gene_effect', '')}, genotype {call.get('genotype', '')}, "
                  f"gnomAD {call.get('gnomad_popmax_af', 'absent')}, score {call.get('score', '')}")
        lines.append(f"| {row['epcr']} | {row['chrom_1']} | {row['pos_1']} | {row['ref_1']} | "
                     f"{row['alt_1']} | {reason} |")
    lines += ["", "## What the owner must check", "",
              "The organizer scorer compares the tuple of chromosome, position, reference allele and "
              "alternate allele. A structural call carries a symbolic alternate allele, so an exact "
              "match with a small variant answer key is unlikely. Decide whether the row earns a place "
              "against the row limit of 10 and against the F-max cost of an extra row.", ""]
    return "\n".join(lines)


def build_arguments() -> argparse.Namespace:
    """Read the command line.

    Returns:
        The parsed arguments.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sv-vcf", type=Path, required=True)
    parser.add_argument("--mei-vcf", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--gnomad-sv", type=Path, required=True)
    parser.add_argument("--intervals", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--genes", default=",".join(MVA_GENES),
                        help="Comma separated disease gene names. The first name is the primary gene "
                             "and it earns the gene points of the primary gene. The default is the "
                             "seven MVA genes with BUB1B first, so the default behaviour is the "
                             "shipped behaviour.")
    parser.add_argument("--deny-tokens", default="", help="Comma separated private identifiers.")
    parser.add_argument("--run-id", default="r1_sv_mei_v1")
    parser.add_argument("--min-frequency-records", type=int, default=100000,
                        help="The gate floor on the number of population records read.")
    return parser.parse_args()


def main() -> int:
    """Run the annotation, the filter, the ranking and the output writing.

    Returns:
        Zero when every gate passed.
    """
    args = build_arguments()
    tokens = deny_tokens(args.deny_tokens)
    disease_genes = gene_names(args.genes)
    primary = disease_genes[0]
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    genes, mva = read_gene_models(args.gtf, disease_genes)
    missing_genes = [name for name in disease_genes if name not in mva]
    gnomad, gnomad_records = read_gnomad_sv(args.gnomad_sv)
    intervals = read_intervals(args.intervals)

    calls = read_calls(args.sv_vcf, "sv") + read_calls(args.mei_vcf, "mei")
    calls = annotate(calls, gene_index(genes), mva, gnomad, intervals, disease_genes)
    candidates = select_candidates(calls)
    context = select_context(calls)

    class_c = out / "class_c_sv_mei_annotated.tsv"
    write_table(class_c, [row_for(call, CLASS_C_FIELDS) for call in calls], CLASS_C_FIELDS)
    ranked_tsv = out / "class_b2_sv_mei_ranked.tsv"
    write_table(ranked_tsv, [row_for(call, CANDIDATE_FIELDS) for call in candidates], CANDIDATE_FIELDS)
    ranked_md = out / "class_b2_sv_mei_ranked.md"
    ranked_md.write_text(markdown_table(candidates, context, primary) + "\n")
    rows = submission_rows(candidates)
    proposal_csv = out / "class_b2_submission_rows.csv"
    with proposal_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUBMISSION_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    proposal_md = out / "class_b2_submission_rows.md"
    proposal_md.write_text(submission_note(rows, candidates, primary) + "\n")

    gates = {
        "gene_models_loaded": "PASS" if genes else "FAIL",
        "mva_genes_found": "PASS" if not missing_genes else "FAIL",
        "frequency_source_loaded": "PASS" if gnomad_records >= args.min_frequency_records else "FAIL",
        "known_intervals_loaded": "PASS" if len(intervals) == 2 else "FAIL",
        "sv_vcf_parsed": "PASS" if args.sv_vcf.exists() else "FAIL",
        "mei_vcf_parsed": "PASS" if args.mei_vcf.exists() else "FAIL",
        "class_b2_written": "PASS" if ranked_md.exists() and proposal_md.exists() else "FAIL",
        "class_c_written": "PASS" if class_c.exists() else "FAIL",
        # Item B13 of reports/audit/submission_audit_plan.md. The verdict states that the proposed
        # rows carry distinct values. It carries no row count and no EPCR value.
        "submission_epcr_values_are_distinct": "PASS" if len(
            {row["epcr"] for row in rows}) == len(rows) else "FAIL",
        "submission_epcr_below_the_b0_ladder": "PASS" if all(
            row["epcr"] <= EPCR_CEILING_OF_THIS_WRITER for row in rows) else "FAIL",
    }
    report = {
        "run_id": args.run_id,
        "data_class_of_run": "C",
        "data_class_of_this_file": "A",
        "scope": ("Checksums, tool versions, reference counts, thresholds and process gates. No count "
                  "of patient calls, no patient coordinate and no patient genotype."),
        "inputs": {
            "sv_vcf": {"name": args.sv_vcf.name, "sha256": sha256_of(args.sv_vcf),
                       "header": read_tool_header(args.sv_vcf, tokens)},
            "mei_vcf": {"name": args.mei_vcf.name, "sha256": sha256_of(args.mei_vcf),
                        "header": read_tool_header(args.mei_vcf, tokens)},
        },
        "references": {
            "gene_models": {"name": args.gtf.name, "sha256": sha256_of(args.gtf),
                            "mva_genes_found": sorted(set(disease_genes) - set(missing_genes)),
                            "mva_genes_missing": missing_genes},
            "gnomad_sv": {"name": args.gnomad_sv.name, "sha256": sha256_of(args.gnomad_sv),
                          "records": gnomad_records},
            "known_intervals": {"name": args.intervals.name, "sha256": sha256_of(args.intervals),
                                "count": len(intervals)},
        },
        "thresholds": {
            "min_sv_length_bp": MIN_SV_LENGTH, "min_genotype_quality": MIN_GENOTYPE_QUALITY,
            "min_retroseq_level": MIN_RETROSEQ_LEVEL, "max_gnomad_af": MAX_GNOMAD_AF,
            "strict_gnomad_af": STRICT_GNOMAD_AF, "reciprocal_overlap": RECIPROCAL_OVERLAP,
            "breakpoint_window_bp": BREAKPOINT_WINDOW, "regulatory_window_bp": REGULATORY_WINDOW,
            "known_interval_pad_bp": KNOWN_INTERVAL_PAD, "candidate_score_floor": CANDIDATE_SCORE_FLOOR,
            "submission_score_floor": SUBMISSION_SCORE_FLOOR, "max_candidate_rows": MAX_CANDIDATE_ROWS,
            "max_submission_rows": MAX_SUBMISSION_ROWS, "mva_genes": list(disease_genes),
            "epcr_map": "round(min(0.95, 0.30 + score / 200), 3)",
            "epcr_tie_step": EPCR_TIE_STEP,
            "epcr_ceiling_of_this_writer": EPCR_CEILING_OF_THIS_WRITER,
            "context_min_length_bp": CONTEXT_MIN_LENGTH, "context_rows": CONTEXT_ROWS,
            "caller_filter_policy": ("A call in an MVA gene stays in the ranked table when the caller "
                                     "filter marks it as failed, and it earns no filter points. A "
                                     "context call must pass the caller filter."),
        },
        "outputs": {
            "class_c": [class_c.name],
            "class_b2": [ranked_tsv.name, ranked_md.name, proposal_csv.name, proposal_md.name],
            "note": ("Every class B2 file is written on every run, so the presence of a file states "
                     "nothing about the child."),
        },
        "gates": gates,
        "verdict": "PASS" if all(value == "PASS" for value in gates.values()) else "FAIL",
    }
    (out / "class_a_gates.json").write_text(json.dumps(report, indent=2) + "\n")
    for name, value in gates.items():
        print(f"{value}  {name}")
    print(f"verdict: {report['verdict']}")
    return 0 if report["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
