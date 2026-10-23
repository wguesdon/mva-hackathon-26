"""Shared readers and writers for the Track 1 structural variant benchmark.

The module holds the file operations that the benchmark driver needs: VCF record access,
INFO parsing, interval arithmetic and size bands. It reads no patient file. Every path that
it opens is a public GIAB file or a call set of the public HG002 run.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path

PRIMARY_CONTIGS = tuple(f"chr{name}" for name in list(range(1, 23)) + ["X", "Y"])

# The band edges follow the GIAB v0.6 report. The last band has no upper edge.
SIZE_BANDS = ((50, 99), (100, 299), (300, 999), (1000, 9999), (10000, None))


@dataclass
class Variant:
    """One VCF record with the fields that the benchmark uses.

    Attributes:
        chrom: The contig name.
        pos: The one-based position.
        identifier: The ID column.
        ref: The REF column.
        alt: The ALT column.
        qual: The QUAL column.
        filter_field: The FILTER column.
        info: The parsed INFO column.
        tail: The FORMAT column and the sample columns.
    """

    chrom: str
    pos: int
    identifier: str
    ref: str
    alt: str
    qual: str
    filter_field: str
    info: dict[str, str] = field(default_factory=dict)
    tail: list[str] = field(default_factory=list)


def open_text(path: Path):
    """Open a plain or a gzipped text file.

    Args:
        path: The file path.

    Returns:
        A text mode file handle.
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt")


def parse_info(text: str) -> dict[str, str]:
    """Split a VCF INFO column into a mapping.

    Args:
        text: The INFO column.

    Returns:
        The key to value mapping. A flag maps to an empty string.
    """
    values: dict[str, str] = {}
    for item in text.split(";"):
        if not item or item == ".":
            continue
        key, _, value = item.partition("=")
        values[key] = value
    return values


def format_info(values: dict[str, str]) -> str:
    """Join a mapping into a VCF INFO column.

    Args:
        values: The key to value mapping.

    Returns:
        The INFO column. A key with an empty value becomes a flag.
    """
    if not values:
        return "."
    parts = [key if value == "" else f"{key}={value}" for key, value in values.items()]
    return ";".join(parts)


def read_vcf(path: Path) -> tuple[list[str], list[Variant]]:
    """Read every record of a VCF.

    Args:
        path: The VCF path, plain or gzipped.

    Returns:
        The header lines and the records.
    """
    header: list[str] = []
    records: list[Variant] = []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header.append(line)
                continue
            fields = line.split("\t")
            records.append(Variant(chrom=fields[0], pos=int(fields[1]), identifier=fields[2],
                                   ref=fields[3], alt=fields[4], qual=fields[5],
                                   filter_field=fields[6], info=parse_info(fields[7]),
                                   tail=fields[8:]))
    return header, records


def variant_size(record: Variant) -> int | None:
    """Return the length of a structural variant in base pairs.

    SVLEN is used first. END minus POS is the fallback for a symbolic record without SVLEN.
    A sequence-resolved record without either field falls back to the allele lengths.

    Args:
        record: The VCF record.

    Returns:
        The absolute size, or None when no field gives a length.
    """
    raw = record.info.get("SVLEN")
    if raw not in (None, "", "."):
        return abs(int(float(raw.split(",")[0])))
    end = record.info.get("END")
    if end not in (None, "", "."):
        return abs(int(end) - record.pos)
    if not record.alt.startswith("<") and "," not in record.alt and record.alt != ".":
        return abs(len(record.alt) - len(record.ref))
    return None


def variant_type(record: Variant) -> str:
    """Return the structural variant class of a record.

    Args:
        record: The VCF record.

    Returns:
        DEL, INS, DUP, INV, BND or NA. A sub type such as DUP:TANDEM maps to DUP.
    """
    svtype = record.info.get("SVTYPE", "")
    if not svtype and not record.alt.startswith("<") and record.alt != ".":
        return "DEL" if len(record.ref) > len(record.alt) else "INS"
    return svtype.split(":")[0] or "NA"


def size_band(size: int | None) -> str:
    """Return the name of the size band of one variant.

    Args:
        size: The variant size in base pairs, or None.

    Returns:
        The band label, or "unsized" when the size is unknown.
    """
    if size is None:
        return "unsized"
    for low, high in SIZE_BANDS:
        if high is None:
            if size >= low:
                return f"{low}+"
        elif low <= size <= high:
            return f"{low}-{high}"
    return "below_50"


def read_bed(path: Path) -> list[tuple[str, int, int]]:
    """Read a BED file into intervals.

    Args:
        path: The BED path.

    Returns:
        The intervals as contig, start and end. The coordinates stay zero-based half open.
    """
    intervals: list[tuple[str, int, int]] = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            intervals.append((fields[0], int(fields[1]), int(fields[2])))
    return intervals


def write_bed(path: Path, intervals: list[tuple[str, int, int]]) -> None:
    """Write intervals to a BED file.

    Args:
        path: The output path.
        intervals: The intervals as contig, start and end.
    """
    with path.open("w") as handle:
        for chrom, start, end in intervals:
            handle.write(f"{chrom}\t{start}\t{end}\n")


def merge_intervals(intervals: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    """Sort intervals and join the ones that touch or overlap.

    Args:
        intervals: The intervals as contig, start and end.

    Returns:
        The sorted and merged intervals.
    """
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in sorted(intervals, key=lambda item: (item[0], item[1], item[2])):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            previous = merged[-1]
            merged[-1] = (chrom, previous[1], max(previous[2], end))
        else:
            merged.append((chrom, start, end))
    return merged


def total_bases(intervals: list[tuple[str, int, int]]) -> int:
    """Sum the length of a set of intervals.

    Args:
        intervals: The intervals as contig, start and end.

    Returns:
        The total number of bases.
    """
    return sum(end - start for _, start, end in intervals)


def index_intervals(intervals: list[tuple[str, int, int]]) -> dict[str, list[tuple[int, int]]]:
    """Group merged intervals by contig for a position lookup.

    Args:
        intervals: The merged intervals as contig, start and end.

    Returns:
        A mapping from contig to the sorted start and end pairs.
    """
    index: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end in intervals:
        index.setdefault(chrom, []).append((start, end))
    for spans in index.values():
        spans.sort()
    return index


def position_in_index(index: dict[str, list[tuple[int, int]]], chrom: str, pos: int) -> bool:
    """Test whether a one-based position falls inside an indexed interval set.

    Args:
        index: The output of index_intervals.
        chrom: The contig name.
        pos: The one-based position.

    Returns:
        True when the position is inside an interval.
    """
    spans = index.get(chrom)
    if not spans:
        return False
    low, high = 0, len(spans) - 1
    target = pos - 1
    while low <= high:
        middle = (low + high) // 2
        start, end = spans[middle]
        if target < start:
            high = middle - 1
        elif target >= end:
            low = middle + 1
        else:
            return True
    return False


def intersect_bases(index: dict[str, list[tuple[int, int]]], chrom: str, start: int, end: int) -> int:
    """Count the bases of one region that fall inside an indexed interval set.

    Args:
        index: The output of index_intervals.
        chrom: The contig name.
        start: The zero-based start of the region.
        end: The exclusive end of the region.

    Returns:
        The number of shared bases.
    """
    total = 0
    for span_start, span_end in index.get(chrom, []):
        if span_end <= start:
            continue
        if span_start >= end:
            break
        total += min(end, span_end) - max(start, span_start)
    return total


def write_vcf(path: Path, header: list[str], records: list[Variant]) -> None:
    """Write records to a plain text VCF.

    Args:
        path: The output path.
        header: The header lines, including the column line.
        records: The records to write.
    """
    with path.open("w") as handle:
        for line in header:
            handle.write(line + "\n")
        for record in records:
            fields = [record.chrom, str(record.pos), record.identifier, record.ref, record.alt,
                      record.qual, record.filter_field, format_info(record.info)]
            fields.extend(record.tail)
            handle.write("\t".join(fields) + "\n")
