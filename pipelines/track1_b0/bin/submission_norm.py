#!/usr/bin/env python3
"""Allele representation primitives and the normalization check of the Track 1 submission writer.

The Track 1 scorer matches a row only when the tuple (chrom, pos, REF, ALT) is equal to the tuple of
the answer key. One biological allele has several valid VCF representations, so a correct allele in
the wrong representation scores zero. This module holds the representation primitives and the check
that runs on a submission file before the file goes out.

Definitions:
    Parsimony      The representation carries the smallest number of bases. No shared base remains
                   at the start or at the end of REF and ALT, except the one padding base that VCF
                   4.2 section 1.6.1 requires for an indel.
    Left alignment The representation is at the smallest position that still describes the same
                   allele. Tan, Abecasis and Kang 2015, Bioinformatics 31(13):2202-2204, define the
                   pair of the two rules as normalization, and `vt normalize` and `bcftools norm`
                   implement it.

The module provides three normalizers, so that a caller can measure which rule recovers a match:
    left_align     Rule 1 alone. The allele shifts to the smallest position and keeps its length.
    parsimony      Rule 2 alone. The allele loses every redundant base and keeps its position.
    normalize      Both rules, which is the output of `bcftools norm -f <reference>`.

Usage:
    submission_norm.py --submission slot1_b0.csv --reference GRCh38.fasta \
        --out slot1_b0_normalized.csv --report norm_report.json [--backend bcftools]
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HEADER = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2", "alt_2",
          "epcr", "finding_type", "notes"]
BASES = set("ACGTN")
MAX_SHIFT = 1000


class Reference:
    """Random access to a FASTA file through its `.fai` index.

    Attributes:
        path: The FASTA path.
        index: Chromosome name to (length, offset, line_bases, line_width).
    """

    def __init__(self, path: Path, cache: bool = False) -> None:
        """Open the FASTA and read its index.

        Args:
            path: The FASTA path. A file `<path>.fai` must exist beside it.
            cache: When true, hold the sequence of the most recent chromosome in memory. A whole
                genome scan is much faster with the cache and it costs the size of one chromosome.

        Raises:
            FileNotFoundError: When the index is absent.
        """
        self.path = Path(path)
        self.cache = cache
        self._cached_name = ""
        self._cached_sequence = ""
        fai = Path(str(path) + ".fai")
        if not fai.exists():
            raise FileNotFoundError(f"No FASTA index at {fai}.")
        self.index: dict[str, tuple[int, int, int, int]] = {}
        with fai.open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                self.index[fields[0]] = (int(fields[1]), int(fields[2]), int(fields[3]), int(fields[4]))
        self._handle = self.path.open("rb")

    def close(self) -> None:
        """Close the FASTA handle."""
        self._handle.close()

    def _read_whole(self, chrom: str) -> str:
        """Read one whole chromosome into memory.

        Args:
            chrom: The chromosome name.

        Returns:
            The sequence in upper case.
        """
        size, offset, line_bases, line_width = self.index[chrom]
        lines = (size + line_bases - 1) // line_bases
        self._handle.seek(offset)
        block = self._handle.read(lines * line_width).decode("ascii")
        return block.replace("\n", "").replace("\r", "").upper()

    def length(self, chrom: str) -> int:
        """Return the length of one chromosome.

        Args:
            chrom: The chromosome name.

        Returns:
            The number of bases.
        """
        return self.index[chrom][0]

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Read reference bases on a 1-based inclusive interval.

        Args:
            chrom: The chromosome name.
            start: The first position, 1-based.
            end: The last position, 1-based and inclusive.

        Returns:
            The bases in upper case. The interval is clipped to the chromosome.
        """
        size, offset, line_bases, line_width = self.index[chrom]
        start = max(1, start)
        end = min(size, end)
        if end < start:
            return ""
        if self.cache:
            if self._cached_name != chrom:
                self._cached_name = chrom
                self._cached_sequence = self._read_whole(chrom)
            return self._cached_sequence[start - 1:end]
        first = offset + (start - 1) // line_bases * line_width + (start - 1) % line_bases
        last = offset + (end - 1) // line_bases * line_width + (end - 1) % line_bases
        self._handle.seek(first)
        block = self._handle.read(last - first + 1).decode("ascii")
        return "".join(character for character in block if character not in "\r\n").upper()


def variant_class(ref: str, alt: str) -> str:
    """Name the class of one allele.

    Args:
        ref: The REF allele.
        alt: The ALT allele.

    Returns:
        One of SNV, MNV, insertion, deletion and complex.
    """
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    shorter = min(len(ref), len(alt))
    if ref[:shorter] == alt[:shorter] or ref[-shorter:] == alt[-shorter:]:
        return "insertion" if len(alt) > len(ref) else "deletion"
    return "complex"


def parsimony(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Remove every redundant base without moving the allele leftwards.

    The function trims a shared suffix first and a shared prefix second, as Tan, Abecasis and Kang
    2015 define parsimony. It keeps one base in REF and one base in ALT.

    Args:
        pos: The 1-based position.
        ref: The REF allele.
        alt: The ALT allele.

    Returns:
        The position, REF and ALT of the parsimonious representation.
    """
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return pos, ref, alt


def left_align(reference: Reference, chrom: str, pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Move one allele to the smallest position that describes the same sequence change.

    The function applies rule 1 of Tan, Abecasis and Kang 2015 alone. While the last base of REF and
    the last base of ALT are equal, it drops that base and it adds the preceding reference base on
    the left. REF and ALT keep their lengths, so the result is left aligned and it is not always
    parsimonious.

    Args:
        reference: The reference genome.
        chrom: The chromosome name.
        pos: The 1-based position.
        ref: The REF allele.
        alt: The ALT allele.

    Returns:
        The position, REF and ALT of the left aligned representation.
    """
    shifts = 0
    while ref and alt and ref[-1] == alt[-1] and pos > 1 and shifts < MAX_SHIFT:
        base = reference.fetch(chrom, pos - 1, pos - 1)
        if not base:
            break
        pos -= 1
        ref, alt = base + ref[:-1], base + alt[:-1]
        shifts += 1
    return pos, ref, alt


def normalize(reference: Reference, chrom: str, pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Return the left aligned and parsimonious representation of one allele.

    This is the normalization of Tan, Abecasis and Kang 2015 and the output of
    `bcftools norm -f <reference>`.

    Args:
        reference: The reference genome.
        chrom: The chromosome name.
        pos: The 1-based position.
        ref: The REF allele.
        alt: The ALT allele.

    Returns:
        The position, REF and ALT of the normalized representation.
    """
    pos, ref, alt = left_align(reference, chrom, pos, ref, alt)
    return parsimony(pos, ref, alt)


def right_align(reference: Reference, chrom: str, pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Move one allele to the largest position that describes the same sequence change.

    A curator who follows the HGVS 3 prime rule writes an allele in this form. The VCF padding base
    stays on the left, so the result is a valid VCF record and not a left aligned one.

    Args:
        reference: The reference genome.
        chrom: The chromosome name.
        pos: The 1-based position.
        ref: The REF allele.
        alt: The ALT allele.

    Returns:
        The position, REF and ALT of the right aligned representation.
    """
    pos, ref, alt = normalize(reference, chrom, pos, ref, alt)
    if len(ref) != 1 and len(alt) != 1:
        return pos, ref, alt
    shifts = 0
    while shifts < MAX_SHIFT:
        unit = ref[1:] if len(alt) == 1 else alt[1:]
        if not unit:
            break
        following = reference.fetch(chrom, pos + len(ref), pos + len(ref))
        if not following or following != unit[0]:
            break
        pad = reference.fetch(chrom, pos + 1, pos + 1)
        unit = unit[1:] + following
        pos += 1
        ref = pad + unit if len(alt) == 1 else pad
        alt = pad if len(alt) == 1 else pad + unit
        shifts += 1
    return pos, ref, alt


def add_padding(reference: Reference, chrom: str, pos: int, ref: str, alt: str,
                left: int = 1, right: int = 1) -> tuple[int, str, str]:
    """Add shared reference bases on both sides of one allele.

    A caller that reports a wider haplotype writes an allele in this form. The representation is
    valid and it is not parsimonious.

    Args:
        reference: The reference genome.
        chrom: The chromosome name.
        pos: The 1-based position.
        ref: The REF allele.
        alt: The ALT allele.
        left: The number of bases to add on the left.
        right: The number of bases to add on the right.

    Returns:
        The position, REF and ALT of the padded representation.
    """
    prefix = reference.fetch(chrom, pos - left, pos - 1) if left else ""
    suffix = reference.fetch(chrom, pos + len(ref), pos + len(ref) + right - 1) if right else ""
    return pos - len(prefix), prefix + ref + suffix, prefix + alt + suffix


def spell_repeat_unit(reference: Reference, chrom: str, pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """Write one indel as the whole repeat tract that holds it.

    A caller that reassembles a tandem repeat reports the tract before the change as REF and the
    tract after the change as ALT. The representation is valid and it is not parsimonious.

    Args:
        reference: The reference genome.
        chrom: The chromosome name.
        pos: The 1-based position.
        ref: The REF allele.
        alt: The ALT allele.

    Returns:
        The position, REF and ALT of the spelled tract, or the input when no tract exists.
    """
    pos, ref, alt = normalize(reference, chrom, pos, ref, alt)
    if len(ref) != 1 and len(alt) != 1:
        return pos, ref, alt
    unit = ref[1:] if len(alt) == 1 else alt[1:]
    if not unit:
        return pos, ref, alt
    tract = reference.fetch(chrom, pos + 1, pos + len(unit) * 40)
    copies = 0
    while tract[copies * len(unit):(copies + 1) * len(unit)] == unit:
        copies += 1
    if copies < 2:
        return pos, ref, alt
    pad = reference.fetch(chrom, pos, pos)
    whole = unit * copies
    changed = unit * (copies - 1) if len(alt) == 1 else unit * (copies + 1)
    return pos, pad + whole, pad + changed


def homopolymer_run(reference: Reference, chrom: str, pos: int, ref: str) -> int:
    """Return the length of the longest homopolymer that touches one locus.

    Args:
        reference: The reference genome.
        chrom: The chromosome name.
        pos: The 1-based position.
        ref: The REF allele.

    Returns:
        The number of identical bases in the longest run that overlaps the locus.
    """
    first = max(1, pos - 30)
    window = reference.fetch(chrom, first, pos + len(ref) + 30)
    if not window:
        return 0
    start = pos - first
    best = 0
    for probe in range(max(0, start - 1), min(len(window), start + len(ref) + 1)):
        base = window[probe]
        low = probe
        while low > 0 and window[low - 1] == base:
            low -= 1
        high = probe
        while high + 1 < len(window) and window[high + 1] == base:
            high += 1
        best = max(best, high - low + 1)
    return best


def tandem_repeat_copies(reference: Reference, chrom: str, pos: int, span: int = 1) -> tuple[int, int]:
    """Return the unit length and the copy number of the tandem repeat that covers one locus.

    The function reads a window of 30 bases on each side of the locus. For every unit length from 2
    to 6 it grows the interval on which the window repeats with that period, starting at the locus,
    and it keeps the longest tract. A tract whose unit is one repeated base is a homopolymer and the
    function skips it.

    Args:
        reference: The reference genome.
        chrom: The chromosome name.
        pos: The 1-based position.
        span: The number of reference bases that the allele covers.

    Returns:
        The unit length and the number of copies. The pair is (0, 0) when no repeat is present.
    """
    first = max(1, pos - 30)
    window = reference.fetch(chrom, first, pos + span + 30)
    locus = pos - first
    best = (0, 0)
    for size in range(2, 7):
        seed = seed_offset(window, locus, span, size)
        if seed is None:
            continue
        low = seed
        while low > 0 and window[low - 1] == window[low - 1 + size]:
            low -= 1
        high = seed
        while high + size < len(window) and window[high] == window[high + size]:
            high += 1
        copies = (high - low + size) // size
        if len(set(window[low:low + size])) == 1:
            continue
        if copies >= 3 and high + size > locus and low <= locus + span + 1 and copies * size > best[0] * best[1]:
            best = (size, copies)
    return best


def seed_offset(window: str, locus: int, span: int, size: int) -> int | None:
    """Find a window offset near one locus where the sequence repeats with a given period.

    Args:
        window: The reference window.
        locus: The offset of the allele inside the window.
        span: The number of reference bases that the allele covers.
        size: The period.

    Returns:
        The first offset with a period match, or None.
    """
    for offset in range(max(0, locus - size), min(len(window) - size, locus + span + 2)):
        if window[offset] == window[offset + size]:
            return offset
    return None


def sequence_context(reference: Reference, chrom: str, pos: int, ref: str) -> str:
    """Name the sequence context of one locus.

    Args:
        reference: The reference genome.
        chrom: The chromosome name.
        pos: The 1-based position.
        ref: The REF allele.

    Returns:
        One of homopolymer_ge6, short_tandem_repeat and other.
    """
    if homopolymer_run(reference, chrom, pos, ref) >= 6:
        return "homopolymer_ge6"
    size, copies = tandem_repeat_copies(reference, chrom, pos, len(ref))
    if size and copies:
        return "short_tandem_repeat"
    return "other"


def read_submission(path: Path) -> list[dict]:
    """Read a Track 1 submission file.

    Args:
        path: The CSV path.

    Returns:
        One dictionary per row.

    Raises:
        ValueError: When the header is not the template header.
    """
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != HEADER:
            raise ValueError(f"The header is {reader.fieldnames} and not the template header.")
        return list(reader)


def row_variants(row: dict) -> list[tuple[int, str, int, str, str]]:
    """Return the variants of one submission row.

    Args:
        row: A submission row.

    Returns:
        One tuple per filled variant slot, holding the slot number, chromosome, position, REF and
        ALT.
    """
    found = []
    for slot in (1, 2):
        chrom = str(row[f"chrom_{slot}"]).strip()
        pos = str(row[f"pos_{slot}"]).strip()
        ref = str(row[f"ref_{slot}"]).strip().upper()
        alt = str(row[f"alt_{slot}"]).strip().upper()
        if not chrom or not pos or not ref or not alt:
            continue
        found.append((slot, chrom, int(pos), ref, alt))
    return found


def bcftools_command(bcftools: str, image: str | None, mounts: list[tuple[Path, str]]) -> list[str]:
    """Build the argument list that runs bcftools.

    Args:
        bcftools: The bcftools command on the host, or inside the container.
        image: A container image, or None to run on the host.
        mounts: Host directory and container directory pairs.

    Returns:
        The argument prefix.
    """
    if image is None:
        return [bcftools]
    argv = ["podman", "run", "--rm"]
    for host, guest in mounts:
        argv += ["-v", f"{host}:{guest}:z"]
    argv += [image, bcftools]
    return argv


def bcftools_version(bcftools: str, image: str | None) -> str:
    """Return the first line of the bcftools version banner.

    Args:
        bcftools: The bcftools command.
        image: A container image, or None.

    Returns:
        The version text.
    """
    argv = bcftools_command(bcftools, image, []) + ["--version"]
    result = subprocess.run(argv, capture_output=True, text=True, check=True)
    return result.stdout.splitlines()[0].strip()


def write_vcf(path: Path, reference: Reference, records: list[tuple[str, str, int, str, str]]) -> None:
    """Write the records of a check as a minimal VCF.

    Args:
        path: The output path.
        reference: The reference genome, which supplies the contig lines.
        records: Identifier, chromosome, position, REF and ALT tuples.
    """
    ordered = sorted(records, key=lambda item: (item[1], item[2], item[3], item[4]))
    with path.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        for name, values in reference.index.items():
            handle.write(f"##contig=<ID={name},length={values[0]}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for identifier, chrom, pos, ref, alt in ordered:
            handle.write(f"{chrom}\t{pos}\t{identifier}\t{ref}\t{alt}\t.\t.\t.\n")


def read_vcf(path: Path) -> dict[str, tuple[str, int, str, str]]:
    """Read the records of a normalized VCF by identifier.

    Args:
        path: The VCF path.

    Returns:
        The identifier mapped to chromosome, position, REF and ALT.
    """
    found = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos, identifier, ref, alt = line.split("\t")[:5]
            found[identifier] = (chrom, int(pos), ref.upper(), alt.upper())
    return found


def normalize_with_bcftools(records: list[tuple[str, str, int, str, str]], reference: Reference,
                            bcftools: str, image: str | None) -> tuple[dict, str, str]:
    """Normalize a set of records with `bcftools norm`.

    Args:
        records: Identifier, chromosome, position, REF and ALT tuples.
        reference: The reference genome.
        bcftools: The bcftools command.
        image: A container image, or None to run on the host.

    Returns:
        The normalized records by identifier, the bcftools version and the bcftools log.

    Raises:
        subprocess.CalledProcessError: When bcftools fails.
    """
    work = Path(tempfile.mkdtemp(prefix="submission_norm_"))
    write_vcf(work / "input.vcf", reference, records)
    mounts = [(work, "/work"), (reference.path.parent, "/ref")]
    argv = bcftools_command(bcftools, image, mounts)
    if image is None:
        source, target, fasta = str(work / "input.vcf"), str(work / "output.vcf"), str(reference.path)
    else:
        source, target, fasta = "/work/input.vcf", "/work/output.vcf", f"/ref/{reference.path.name}"
    argv += ["norm", "--no-version", "--check-ref", "w", "--fasta-ref", fasta, "--output", target,
             "--output-type", "v", source]
    result = subprocess.run(argv, capture_output=True, text=True, check=True)
    return read_vcf(work / "output.vcf"), bcftools_version(bcftools, image), result.stderr.strip()


def check_submission(rows: list[dict], reference: Reference, bcftools: str | None,
                     image: str | None) -> dict:
    """Compare every submission variant with its normalized form.

    The function does not decide anything. It returns the normalized rows and the record of every
    variant that changed, so the owner reads the record and decides.

    Args:
        rows: The submission rows.
        reference: The reference genome.
        bcftools: The bcftools command, or None to use the internal normalizer only.
        image: A container image, or None to run bcftools on the host.

    Returns:
        A report with the normalized rows, the changed variants and the backend agreement.
    """
    records = []
    for index, row in enumerate(rows):
        for slot, chrom, pos, ref, alt in row_variants(row):
            records.append((f"r{index}s{slot}", chrom, pos, ref, alt))
    internal = {}
    for identifier, chrom, pos, ref, alt in records:
        internal[identifier] = (chrom,) + normalize(reference, chrom, pos, ref, alt)
    external: dict = {}
    version = "not run"
    log = ""
    if bcftools is not None and records:
        external, version, log = normalize_with_bcftools(records, reference, bcftools, image)
    disagreements = [identifier for identifier in external if external[identifier] != internal[identifier]]
    chosen = external if external else internal
    changes = []
    alternatives = []
    normalized_rows = [dict(row) for row in rows]
    for identifier, chrom, pos, ref, alt in records:
        index = int(identifier[1:identifier.index("s")])
        slot = int(identifier[identifier.index("s") + 1:])
        new_chrom, new_pos, new_ref, new_alt = chosen[identifier]
        reference_base = reference.fetch(chrom, pos, pos + len(ref) - 1)
        row = normalized_rows[index]
        row[f"chrom_{slot}"], row[f"pos_{slot}"] = new_chrom, str(new_pos)
        row[f"ref_{slot}"], row[f"alt_{slot}"] = new_ref, new_alt
        right = right_align(reference, new_chrom, new_pos, new_ref, new_alt)
        if right != (new_pos, new_ref, new_alt):
            alternatives.append({"row": index + 1, "slot": slot,
                                 "left_aligned": f"{new_chrom}:{new_pos}:{new_ref}:{new_alt}",
                                 "right_aligned": f"{new_chrom}:{right[0]}:{right[1]}:{right[2]}",
                                 "shift_bp": right[0] - new_pos,
                                 "variant_class": variant_class(new_ref, new_alt),
                                 "sequence_context": sequence_context(reference, new_chrom, new_pos, new_ref)})
        if (chrom, pos, ref, alt) == (new_chrom, new_pos, new_ref, new_alt):
            continue
        changes.append({"row": index + 1, "slot": slot,
                        "submitted": f"{chrom}:{pos}:{ref}:{alt}",
                        "normalized": f"{new_chrom}:{new_pos}:{new_ref}:{new_alt}",
                        "variant_class": variant_class(ref, alt),
                        "sequence_context": sequence_context(reference, chrom, pos, ref),
                        "reference_matches_ref": reference_base == ref})
    return {"reference": str(reference.path), "bcftools": version, "bcftools_log": log,
            "backend": "bcftools" if external else "internal",
            "variants": len(records), "changed": len(changes),
            "with_a_second_representation": len(alternatives),
            "internal_bcftools_disagreements": disagreements,
            "changes": changes, "alternatives": alternatives, "normalized_rows": normalized_rows}


def write_rows(path: Path, rows: list[dict]) -> None:
    """Write submission rows to a CSV file.

    Args:
        path: The output path.
        rows: The rows.
    """
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    """Run the normalization check on one submission file.

    Returns:
        Zero when no variant changed, and 3 when at least one variant changed.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--submission", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--backend", choices=("bcftools", "internal"), default="bcftools")
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--container-image")
    args = parser.parse_args()
    reference = Reference(args.reference)
    rows = read_submission(args.submission)
    report = check_submission(rows, reference, args.bcftools if args.backend == "bcftools" else None,
                              args.container_image)
    normalized = report.pop("normalized_rows")
    if args.out:
        write_rows(args.out, normalized)
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(f"{report['variants']} variants, {report['changed']} changed under normalization, "
          f"{report['with_a_second_representation']} have a right aligned form at another position, "
          f"backend {report['backend']}, {report['bcftools']}")
    for change in report["changes"]:
        print(f"  row {change['row']} slot {change['slot']}: {change['submitted']} -> "
              f"{change['normalized']} ({change['variant_class']}, {change['sequence_context']})")
    if report["internal_bcftools_disagreements"]:
        print(f"The two normalizers disagree on {report['internal_bcftools_disagreements']}.", file=sys.stderr)
        return 4
    reference.close()
    return 3 if report["changed"] else 0


if __name__ == "__main__":
    sys.exit(main())
