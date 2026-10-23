"""Measure how often the Track 1 allele representation differs from the public source form. Class A.

The Track 1 scorer matches a row only when the tuple (chrom, pos, REF, ALT) is equal to the tuple of
the answer key. The organizers do not publish how the key normalizes an allele, so a correct allele
in another representation scores zero. This benchmark measures the size of that risk on public data
only. It reads ClinVar, the S3 literature allele table, the S4 spike-in manifest and the GIAB HG002
call set that forms the suite background. It reads no patient file.

Arms:
    source     The spike-in suite alleles, which come from ClinVar and from the literature list.
               It compares the raw source tuple with the left aligned and parsimonious tuple.
    caller     The ClinVar records that the HG002 DeepVariant call set also holds. It joins the two
               sets on the normalized tuple and compares the two raw tuples. This is the direct
               measurement of the pipeline form against the public source form.
    clinvar    Every eligible ClinVar record. The suite holds no MNV and no complex allele, so this
               arm supplies the counts of those two classes.
    recovery   Generated representation classes of every suite indel. For each class it records
               whether left alignment alone, parsimony alone, both rules or no rule recovers the
               normalized tuple.

Usage:
    uv run python scripts/track1_allele_representation.py \
        --reference work/ref/Homo_sapiens_assembly38.fasta \
        --clinvar tracking/cache/clinvar_20260913.vcf.gz \
        --background tracking/cache/hg002_s2_deepvariant.vcf.gz \
        --out reports/benchmarks/allele_representation_v1
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipelines/track1_b0/bin"))

import submission_norm as norm  # noqa: E402

MAX_ALLELE = 51
CHROMS = tuple(f"chr{number}" for number in range(1, 23))
CLASSES = ("SNV", "MNV", "insertion", "deletion", "complex")
CONTEXTS = ("homopolymer_ge6", "short_tandem_repeat", "other")
REPRESENTATIONS = ("left_parsimonious", "right_aligned", "padded_both_sides", "repeat_unit_spelled",
                   "left_aligned_unparsimonious", "mnv_merged", "complex_merged")


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


def open_text(path: Path):
    """Open a plain or gzipped text file.

    Args:
        path: The file path.

    Returns:
        A text handle.
    """
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def suite_alleles(manifest: Path) -> list[dict]:
    """Read the unique alleles of the S4 spike-in suite.

    Args:
        manifest: The suite manifest.

    Returns:
        One record per allele with its chromosome, position, REF, ALT and source.
    """
    found: dict[tuple, dict] = {}
    data = json.loads(manifest.read_text())
    for scenario in data["scenarios"]:
        for variant in scenario["variants"]:
            key = (variant["chrom"], variant["pos"], variant["ref"].upper(), variant["alt"].upper())
            found.setdefault(key, {"chrom": key[0], "pos": key[1], "ref": key[2], "alt": key[3],
                                   "source": variant["source"], "gene": scenario["gene"],
                                   "group": scenario["group"]})
    return sorted(found.values(), key=lambda item: (item["chrom"], item["pos"]))


def eligible(chrom: str, ref: str, alt: str) -> bool:
    """Say whether one VCF record can enter the benchmark.

    Args:
        chrom: The chromosome with a chr prefix.
        ref: The REF allele.
        alt: The ALT allele.

    Returns:
        True when the record is an autosomal single-allele record of plain bases.
    """
    if chrom not in CHROMS or "," in alt or alt in (".", ""):
        return False
    if len(ref) > MAX_ALLELE or len(alt) > MAX_ALLELE:
        return False
    return set(ref) <= set("ACGT") and set(alt) <= set("ACGT")


def read_clinvar(path: Path) -> list[tuple[str, int, str, str, str]]:
    """Read every eligible ClinVar record.

    Args:
        path: The ClinVar GRCh38 VCF.

    Returns:
        Chromosome, position, REF, ALT and variation identifier tuples.
    """
    found = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = "chr" + fields[0].removeprefix("chr")
            ref, alt = fields[3].upper(), fields[4].upper()
            if not eligible(chrom, ref, alt):
                continue
            found.append((chrom, int(fields[1]), ref, alt, fields[2]))
    return found


def read_calls(path: Path) -> list[tuple[str, int, str, str, str]]:
    """Read every eligible PASS record of a single-sample call set.

    Args:
        path: The bgzipped VCF.

    Returns:
        Chromosome, position, REF, ALT and an empty identifier.
    """
    found = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[6] not in ("PASS", "."):
                continue
            chrom = "chr" + fields[0].removeprefix("chr")
            ref, alt = fields[3].upper(), fields[4].upper()
            if not eligible(chrom, ref, alt):
                continue
            found.append((chrom, int(fields[1]), ref, alt, ""))
    return found


def normalized_index(reference: norm.Reference,
                     records: list[tuple[str, int, str, str, str]]) -> dict[tuple, tuple]:
    """Map the normalized tuple of every record to its raw tuple.

    Args:
        reference: The reference genome.
        records: Chromosome, position, REF, ALT and identifier tuples.

    Returns:
        The normalized tuple mapped to the raw tuple. A normalized tuple that two raw records share
        is dropped, because the join would be ambiguous.
    """
    index: dict[tuple, tuple] = {}
    ambiguous = set()
    for chrom, pos, ref, alt, _ in sorted(records):
        key = (chrom,) + norm.normalize(reference, chrom, pos, ref, alt)
        raw = (chrom, pos, ref, alt)
        if key in index and index[key] != raw:
            ambiguous.add(key)
        index[key] = raw
    for key in ambiguous:
        index.pop(key, None)
    return index


def measure_raw_against_normalized(reference: norm.Reference,
                                   records: list[tuple[str, int, str, str, str]]) -> dict:
    """Count the records whose raw tuple differs from its normalized tuple.

    Args:
        reference: The reference genome.
        records: Chromosome, position, REF, ALT and identifier tuples.

    The variant class and the sequence context come from the normalized allele, so that one allele
    keeps one class whatever representation the source used.

    Returns:
        Totals and differing counts by variant class and by sequence context, with the first 500
        differing records listed.
    """
    total: Counter = Counter()
    differ: Counter = Counter()
    by_context_total: Counter = Counter()
    by_context_differ: Counter = Counter()
    rules: Counter = Counter()
    listed = []
    for chrom, pos, ref, alt, identifier in sorted(records):
        canonical = norm.normalize(reference, chrom, pos, ref, alt)
        new_pos, new_ref, new_alt = canonical
        variant = norm.variant_class(new_ref, new_alt)
        context = norm.sequence_context(reference, chrom, new_pos, new_ref)
        total[variant] += 1
        by_context_total[(variant, context)] += 1
        if (pos, ref, alt) != (new_pos, new_ref, new_alt):
            rule = recovery_of(reference, chrom, canonical, (pos, ref, alt))
            rules[f"{variant}|{rule}"] += 1
            differ[variant] += 1
            by_context_differ[(variant, context)] += 1
            if len(listed) < 500:
                listed.append({"raw": f"{chrom}:{pos}:{ref}:{alt}",
                               "normalized": f"{chrom}:{new_pos}:{new_ref}:{new_alt}",
                               "variant_class": variant, "sequence_context": context,
                               "recovery": rule, "id": identifier})
    return {"total": dict(total), "differ": dict(differ),
            "by_context_total": {f"{key[0]}|{key[1]}": value for key, value in by_context_total.items()},
            "by_context_differ": {f"{key[0]}|{key[1]}": value for key, value in by_context_differ.items()},
            "recovery_rules": dict(rules), "records": listed}


def join_on_normalized(reference: norm.Reference, source: dict[tuple, tuple],
                       caller: dict[tuple, tuple]) -> dict:
    """Compare the raw tuples of two sets that describe the same normalized allele.

    Args:
        reference: The reference genome.
        source: The normalized to raw map of the public source.
        caller: The normalized to raw map of the call set.

    Returns:
        Totals, disagreement counts by class and by context, and the disagreeing records.
    """
    total: Counter = Counter()
    differ: Counter = Counter()
    by_context_total: Counter = Counter()
    by_context_differ: Counter = Counter()
    rules: Counter = Counter()
    listed = []
    for key in sorted(set(source) & set(caller)):
        chrom, pos, ref, alt = key
        variant = norm.variant_class(ref, alt)
        context = norm.sequence_context(reference, chrom, pos, ref)
        total[variant] += 1
        by_context_total[(variant, context)] += 1
        if source[key] != caller[key]:
            rule = recovery_of(reference, chrom, (pos, ref, alt), source[key][1:])
            rules[f"{variant}|{rule}"] += 1
            differ[variant] += 1
            by_context_differ[(variant, context)] += 1
            if len(listed) < 200:
                listed.append({"normalized": f"{chrom}:{pos}:{ref}:{alt}",
                               "source": "{}:{}:{}:{}".format(*source[key]),
                               "caller": "{}:{}:{}:{}".format(*caller[key]),
                               "variant_class": variant, "sequence_context": context,
                               "recovery": rule})
    return {"total": dict(total), "differ": dict(differ),
            "by_context_total": {f"{key[0]}|{key[1]}": value for key, value in by_context_total.items()},
            "by_context_differ": {f"{key[0]}|{key[1]}": value for key, value in by_context_differ.items()},
            "recovery_rules": dict(rules), "records": listed,
            "records_total": sum(differ.values())}


def measure_ambiguity(reference: norm.Reference,
                      records: list[tuple[str, int, str, str, str]]) -> dict:
    """Count the alleles that have a second valid representation at another position.

    An allele whose right aligned form differs from its left aligned form is exposed to the
    representation risk, because a curator who follows the HGVS 3 prime rule writes the right
    aligned form. An allele whose two forms are equal carries no risk of that kind.

    Args:
        reference: The reference genome.
        records: Chromosome, position, REF, ALT and identifier tuples.

    Returns:
        Totals and exposed counts by variant class and by sequence context, with the shift sizes.
    """
    total: Counter = Counter()
    exposed: Counter = Counter()
    by_context_total: Counter = Counter()
    by_context_exposed: Counter = Counter()
    shifts = []
    for chrom, pos, ref, alt, _ in records:
        canonical = norm.normalize(reference, chrom, pos, ref, alt)
        variant = norm.variant_class(canonical[1], canonical[2])
        context = norm.sequence_context(reference, chrom, canonical[0], canonical[1])
        total[variant] += 1
        by_context_total[(variant, context)] += 1
        right = norm.right_align(reference, chrom, *canonical)
        if right != canonical:
            exposed[variant] += 1
            by_context_exposed[(variant, context)] += 1
            shifts.append(right[0] - canonical[0])
    shifts.sort()
    return {"total": dict(total), "differ": dict(exposed),
            "by_context_total": {f"{key[0]}|{key[1]}": value for key, value in by_context_total.items()},
            "by_context_differ": {f"{key[0]}|{key[1]}": value for key, value in by_context_exposed.items()},
            "shift_bp": {"n": len(shifts), "min": shifts[0] if shifts else 0,
                         "median": shifts[len(shifts) // 2] if shifts else 0,
                         "max": shifts[-1] if shifts else 0}}


def measure_atomization(reference: norm.Reference, clinvar: list[tuple[str, int, str, str, str]],
                        called: set[tuple]) -> dict:
    """Count the merged ClinVar records whose single base changes the caller reports one by one.

    A curator writes two adjacent substitutions as one MNV record. A short read caller writes them
    as two records. No normalization rule joins the two records, so this class of allele cannot be
    recovered after the fact.

    Args:
        reference: The reference genome.
        clinvar: Chromosome, position, REF, ALT and identifier tuples of the public source.
        called: The normalized tuples of the call set.

    Returns:
        The counts of examined, atomized and merged records for the MNV class and the complex class.
    """
    result = {"MNV": {"records": 0, "merged_called": 0, "all_atoms_called": 0},
              "complex": {"records": 0, "merged_called": 0}}
    for chrom, pos, ref, alt, _ in clinvar:
        canonical = norm.normalize(reference, chrom, pos, ref, alt)
        variant = norm.variant_class(canonical[1], canonical[2])
        if variant not in ("MNV", "complex"):
            continue
        result[variant]["records"] += 1
        if (chrom,) + canonical in called:
            result[variant]["merged_called"] += 1
        if variant == "MNV":
            atoms = [(chrom, canonical[0] + offset, canonical[1][offset], canonical[2][offset])
                     for offset in range(len(canonical[1])) if canonical[1][offset] != canonical[2][offset]]
            if atoms and all(atom in called for atom in atoms):
                result["MNV"]["all_atoms_called"] += 1
    return result


def representations(reference: norm.Reference, chrom: str, pos: int, ref: str,
                    alt: str) -> dict[str, tuple[int, str, str]]:
    """Build one representation of each class for a single allele.

    Args:
        reference: The reference genome.
        chrom: The chromosome name.
        pos: The 1-based position.
        ref: The REF allele.
        alt: The ALT allele.

    Returns:
        The class name mapped to its position, REF and ALT. A class that does not apply to the
        allele is absent.
    """
    canonical = norm.normalize(reference, chrom, pos, ref, alt)
    found = {"left_parsimonious": canonical}
    right = norm.right_align(reference, chrom, *canonical)
    if right != canonical:
        found["right_aligned"] = right
    padded = norm.add_padding(reference, chrom, *canonical)
    if padded != canonical:
        found["padded_both_sides"] = padded
    spelled = norm.spell_repeat_unit(reference, chrom, *canonical)
    if spelled != canonical:
        found["repeat_unit_spelled"] = spelled
    unparsimonious = norm.left_align(reference, chrom, *padded)
    if unparsimonious != canonical:
        found["left_aligned_unparsimonious"] = unparsimonious
    merged = merge_with_neighbour(reference, chrom, *canonical)
    if merged is not None and merged != canonical:
        found["mnv_merged" if norm.variant_class(*merged[1:]) == "MNV" else "complex_merged"] = merged
    return found


def merge_with_neighbour(reference: norm.Reference, chrom: str, pos: int, ref: str,
                         alt: str) -> tuple[int, str, str] | None:
    """Merge one allele with a substitution of the next reference base.

    A caller that reports a haplotype writes two changes that sit side by side as one record. The
    record is an MNV when both changes are substitutions, and a complex allele otherwise.

    Args:
        reference: The reference genome.
        chrom: The chromosome name.
        pos: The 1-based position.
        ref: The REF allele.
        alt: The ALT allele.

    Returns:
        The merged position, REF and ALT, or None when the next base cannot be read.
    """
    following = reference.fetch(chrom, pos + len(ref), pos + len(ref))
    if not following:
        return None
    changed = "A" if following != "A" else "C"
    return pos, ref + following, alt + changed


def recovery_of(reference: norm.Reference, chrom: str, canonical: tuple[int, str, str],
                candidate: tuple[int, str, str]) -> str:
    """Name the rule that recovers the canonical tuple from one representation.

    Args:
        reference: The reference genome.
        chrom: The chromosome name.
        canonical: The left aligned and parsimonious position, REF and ALT.
        candidate: The position, REF and ALT of the representation under test.

    Returns:
        One of identical, left_alignment_alone, parsimony_alone, both_rules and none.
    """
    if candidate == canonical:
        return "identical"
    aligned = norm.left_align(reference, chrom, *candidate)
    trimmed = norm.parsimony(*candidate)
    if aligned == canonical:
        return "left_alignment_alone"
    if trimmed == canonical:
        return "parsimony_alone"
    if norm.normalize(reference, chrom, *candidate) == canonical:
        return "both_rules"
    return "none"


def measure_recovery(reference: norm.Reference, alleles: list[dict]) -> dict:
    """Build every representation class of every allele and record its recovery rule.

    Args:
        reference: The reference genome.
        alleles: Allele records with chrom, pos, ref and alt.

    Returns:
        The recovery counts per representation class and the alleles that no rule recovers.
    """
    counts: dict[str, Counter] = defaultdict(Counter)
    unrecovered = []
    for allele in alleles:
        chrom, pos, ref, alt = allele["chrom"], allele["pos"], allele["ref"], allele["alt"]
        canonical = norm.normalize(reference, chrom, pos, ref, alt)
        for name, candidate in representations(reference, chrom, pos, ref, alt).items():
            rule = recovery_of(reference, chrom, canonical, candidate)
            counts[name][rule] += 1
            if rule == "none":
                unrecovered.append({"allele": f"{chrom}:{canonical[0]}:{canonical[1]}:{canonical[2]}",
                                    "representation": name,
                                    "written_as": f"{chrom}:{candidate[0]}:{candidate[1]}:{candidate[2]}",
                                    "variant_class": norm.variant_class(ref, alt),
                                    "sequence_context": norm.sequence_context(reference, chrom, canonical[0], canonical[1]),
                                    "source": allele.get("source", "")})
    return {"counts": {name: dict(value) for name, value in counts.items()},
            "unrecovered": unrecovered}


def cross_check(reference: norm.Reference, records: list[tuple[str, int, str, str, str]],
                bcftools: str, image: str | None) -> dict:
    """Compare the internal normalizer with `bcftools norm` on a record set.

    Args:
        reference: The reference genome.
        records: Chromosome, position, REF, ALT and identifier tuples.
        bcftools: The bcftools command.
        image: A container image, or None to run on the host.

    Returns:
        The record count, the bcftools version and the disagreeing tuples.
    """
    numbered = [(f"v{index}", item[0], item[1], item[2], item[3]) for index, item in enumerate(records)]
    external, version, log = norm.normalize_with_bcftools(numbered, reference, bcftools, image)
    disagreements = []
    for identifier, chrom, pos, ref, alt in numbered:
        internal = (chrom,) + norm.normalize(reference, chrom, pos, ref, alt)
        if external.get(identifier) != internal:
            disagreements.append({"raw": f"{chrom}:{pos}:{ref}:{alt}",
                                  "internal": "{}:{}:{}:{}".format(*internal),
                                  "bcftools": "{}:{}:{}:{}".format(*external.get(identifier, ("", 0, "", "")))})
    return {"records": len(numbered), "bcftools": version, "disagreements": disagreements,
            "log": log[:2000]}


def write_table(path: Path, rows: list[dict], header: list[str]) -> None:
    """Write one result table as a comma separated file.

    Args:
        path: The output path.
        rows: The rows.
        header: The column names.
    """
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)


def rate_rows(result: dict) -> list[dict]:
    """Turn a measurement into one row per variant class and sequence context.

    Args:
        result: A measurement from measure_raw_against_normalized or join_on_normalized.

    Returns:
        One row per class and context with the counts and the rate.
    """
    rows = []
    for variant in CLASSES:
        total = result["total"].get(variant, 0)
        differ = result["differ"].get(variant, 0)
        if total:
            rows.append({"variant_class": variant, "sequence_context": "all", "n": total,
                         "differing": differ, "rate": round(differ / total, 4)})
        for context in CONTEXTS:
            key = f"{variant}|{context}"
            sub_total = result["by_context_total"].get(key, 0)
            if not sub_total:
                continue
            sub_differ = result["by_context_differ"].get(key, 0)
            rows.append({"variant_class": variant, "sequence_context": context, "n": sub_total,
                         "differing": sub_differ, "rate": round(sub_differ / sub_total, 4)})
    return rows


def main() -> int:
    """Run the benchmark and write the result files.

    Returns:
        Zero on success.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--clinvar", type=Path, required=True)
    parser.add_argument("--background", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=ROOT / "reports/benchmarks/spikein_v1/manifest.json")
    parser.add_argument("--out", type=Path, default=ROOT / "reports/benchmarks/allele_representation_v1")
    parser.add_argument("--bcftools", default="bcftools")
    parser.add_argument("--container-image")
    parser.add_argument("--cross-check-sample", type=int, default=20000,
                        help="The number of ClinVar indel, MNV and complex records that bcftools "
                             "normalizes beside the suite alleles, as a check of the internal "
                             "normalizer.")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    reference = norm.Reference(args.reference, cache=True)

    alleles = suite_alleles(args.manifest)
    suite_records = [(item["chrom"], item["pos"], item["ref"], item["alt"], item["source"]) for item in alleles]
    source_arm = measure_raw_against_normalized(reference, suite_records)
    print(f"source arm: {len(suite_records)} suite alleles, "
          f"{sum(source_arm['differ'].values())} differ from the normalized form")

    clinvar = read_clinvar(args.clinvar)
    clinvar_arm = measure_raw_against_normalized(reference, clinvar)
    print(f"clinvar arm: {len(clinvar)} records, {sum(clinvar_arm['differ'].values())} differ")

    calls = read_calls(args.background)
    caller_arm_self = measure_raw_against_normalized(reference, calls)
    print(f"caller arm: {len(calls)} HG002 calls, {sum(caller_arm_self['differ'].values())} differ")

    called_index = normalized_index(reference, calls)
    joined = join_on_normalized(reference, normalized_index(reference, clinvar), called_index)
    print(f"join arm: {sum(joined['total'].values())} shared alleles, "
          f"{sum(joined['differ'].values())} written differently")

    suite_ambiguity = measure_ambiguity(reference, suite_records)
    clinvar_ambiguity = measure_ambiguity(reference, clinvar)
    print(f"ambiguity arm: {sum(suite_ambiguity['differ'].values())} of {len(suite_records)} suite "
          f"alleles and {sum(clinvar_ambiguity['differ'].values())} of {len(clinvar)} ClinVar "
          f"records have a right aligned form at another position")

    atomization = measure_atomization(reference, clinvar, set(called_index))
    print(f"atomization arm: {atomization}")

    recovery = measure_recovery(reference, alleles)
    print(f"recovery arm: {sum(sum(value.values()) for value in recovery['counts'].values())} representations")

    sample = [item for item in clinvar if len(item[2]) != 1 or len(item[3]) != 1]
    step = max(1, len(sample) // args.cross_check_sample)
    check = cross_check(reference, suite_records + sample[::step], args.bcftools, args.container_image)
    print(f"cross check: {check['records']} alleles, {len(check['disagreements'])} disagreements, {check['bcftools']}")

    write_table(args.out / "source_arm_rates.csv", rate_rows(source_arm),
                ["variant_class", "sequence_context", "n", "differing", "rate"])
    write_table(args.out / "clinvar_arm_rates.csv", rate_rows(clinvar_arm),
                ["variant_class", "sequence_context", "n", "differing", "rate"])
    write_table(args.out / "caller_arm_rates.csv", rate_rows(caller_arm_self),
                ["variant_class", "sequence_context", "n", "differing", "rate"])
    write_table(args.out / "join_arm_rates.csv", rate_rows(joined),
                ["variant_class", "sequence_context", "n", "differing", "rate"])
    write_table(args.out / "suite_ambiguity_rates.csv", rate_rows(suite_ambiguity),
                ["variant_class", "sequence_context", "n", "differing", "rate"])
    write_table(args.out / "clinvar_ambiguity_rates.csv", rate_rows(clinvar_ambiguity),
                ["variant_class", "sequence_context", "n", "differing", "rate"])

    summary = {"created_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
               "data_class": "A",
               "reference": {"file": args.reference.name, "accession": "GCA_000001405.15 GRCh38, GATK bundle"},
               "inputs": {"clinvar": {"file": args.clinvar.name, "sha256": sha256(args.clinvar)},
                          "background": {"file": args.background.name, "sha256": sha256(args.background)},
                          "manifest": {"file": str(args.manifest.relative_to(ROOT)), "sha256": sha256(args.manifest)}},
               "source_arm": {key: value for key, value in source_arm.items() if key != "records"},
               "source_arm_records": source_arm["records"],
               "clinvar_arm": {key: value for key, value in clinvar_arm.items() if key != "records"},
               "caller_arm": {key: value for key, value in caller_arm_self.items() if key != "records"},
               "join_arm": joined,
               "suite_ambiguity": suite_ambiguity,
               "clinvar_ambiguity": clinvar_ambiguity,
               "atomization": atomization,
               "recovery": recovery,
               "cross_check": check}
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    reference.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
