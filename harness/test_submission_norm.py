"""Verify the normalization check of the Track 1 submission writer on a synthetic reference.

The fixture reference holds a homopolymer of 8 A bases, a period 3 tandem repeat, a period 2 tandem
repeat, a run of 6 G bases and ordinary sequence. Each test writes one submission row in one
representation class and asserts that the check recovers the left aligned and parsimonious tuple and
records the change. The test uses the internal normalizer, so it needs no container and no network.
"""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipelines/track1_b0/bin"))
sys.path.insert(0, str(ROOT / "scripts"))
import submission_norm as norm  # noqa: E402
import track1_allele_representation as bench  # noqa: E402

# chrT holds 1-40 ACGT units, 41-48 a homopolymer of 8 A bases, 49-63 a period 3 repeat of CGT,
# 64-75 a period 2 repeat of CA, 76-79 a spacer, 80-85 a run of 6 G bases and 86-133 a sequence with
# no homopolymer of 6 bases and no tandem repeat of 3 copies.
SEQUENCE = ("ACGT" * 10 + "A" * 8 + "CGT" * 5 + "CA" * 6 + "TTAC" + "G" * 6
            + "TGCGTCTGACATCGCATTTGAAACCCAGTAGGTACTGCCTTAGTTGCA")
LINE = 60


def write_reference(directory: Path) -> Path:
    """Write the fixture FASTA and its index.

    Args:
        directory: The directory that receives the files.

    Returns:
        The FASTA path.
    """
    path = directory / "chrT.fa"
    lines = [SEQUENCE[start:start + LINE] for start in range(0, len(SEQUENCE), LINE)]
    path.write_text(">chrT\n" + "\n".join(lines) + "\n")
    header = len(">chrT\n")
    (directory / "chrT.fa.fai").write_text(f"chrT\t{len(SEQUENCE)}\t{header}\t{LINE}\t{LINE + 1}\n")
    return path


def write_submission(path: Path, variants: list[tuple[str, int, str, str]]) -> None:
    """Write one submission row per variant.

    Args:
        path: The CSV path.
        variants: Chromosome, position, REF and ALT tuples.
    """
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=norm.HEADER)
        writer.writeheader()
        for index, (chrom, pos, ref, alt) in enumerate(variants):
            writer.writerow({"proband_id": "PROBAND01", "chrom_1": chrom, "pos_1": pos,
                             "ref_1": ref, "alt_1": alt, "chrom_2": "", "pos_2": "", "ref_2": "",
                             "alt_2": "", "epcr": f"{0.95 - 0.09 * index:.2f}",
                             "finding_type": "primary", "notes": "fixture"})


class NormalizerTest(unittest.TestCase):
    """The three normalizers behave as Tan, Abecasis and Kang 2015 define them."""

    def setUp(self) -> None:
        """Write the fixture reference."""
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.reference = norm.Reference(write_reference(self.root))
        self.addCleanup(self.reference.close)

    def test_fetch_reads_the_sequence(self) -> None:
        """The index reader returns the fixture bases."""
        self.assertEqual(self.reference.fetch("chrT", 1, len(SEQUENCE)), SEQUENCE)
        self.assertEqual(self.reference.fetch("chrT", 41, 48), "A" * 8)

    def test_variant_class_names_every_class(self) -> None:
        """Each representation class carries the expected name."""
        self.assertEqual(norm.variant_class("A", "G"), "SNV")
        self.assertEqual(norm.variant_class("TA", "GC"), "MNV")
        self.assertEqual(norm.variant_class("T", "TA"), "insertion")
        self.assertEqual(norm.variant_class("TA", "T"), "deletion")
        self.assertEqual(norm.variant_class("CG", "A"), "complex")

    def test_homopolymer_deletion_left_aligns(self) -> None:
        """A deletion inside the A run moves to the first base of the run."""
        for pos, ref, alt in [(47, "AA", "A"), (44, "AA", "A"), (40, "TA", "T")]:
            self.assertEqual(norm.normalize(self.reference, "chrT", pos, ref, alt), (40, "TA", "T"))

    def test_homopolymer_insertion_left_aligns(self) -> None:
        """An insertion inside the A run moves to the first base of the run."""
        for pos, ref, alt in [(48, "A", "AA"), (44, "A", "AA"), (40, "T", "TA")]:
            self.assertEqual(norm.normalize(self.reference, "chrT", pos, ref, alt), (40, "T", "TA"))

    def test_tandem_repeat_deletion_left_aligns(self) -> None:
        """A CA deletion inside the period 2 repeat moves to the base before the tract."""
        for pos, ref, alt in [(73, "ACA", "A"), (69, "ACA", "A"), (63, "TCA", "T")]:
            self.assertEqual(norm.normalize(self.reference, "chrT", pos, ref, alt), (63, "TCA", "T"))

    def test_tandem_repeat_insertion_left_aligns(self) -> None:
        """A CA insertion inside the period 2 repeat moves to the base before the tract."""
        for pos, ref, alt in [(75, "A", "ACA"), (69, "A", "ACA"), (63, "T", "TCA")]:
            self.assertEqual(norm.normalize(self.reference, "chrT", pos, ref, alt), (63, "T", "TCA"))

    def test_parsimony_alone_keeps_the_position(self) -> None:
        """Parsimony removes a redundant base and does not shift the allele."""
        self.assertEqual(norm.parsimony(90, "TCT", "AGT"), (90, "TC", "AG"))
        self.assertEqual(norm.parsimony(47, "AA", "A"), (47, "AA", "A"))

    def test_left_alignment_alone_keeps_the_length(self) -> None:
        """Left alignment shifts an allele and leaves a redundant leading base in place."""
        self.assertEqual(norm.left_align(self.reference, "chrT", 39, "GTAA", "GTA"), (38, "CGTA", "CGT"))

    def test_right_alignment_is_the_mirror(self) -> None:
        """The right aligned form sits at the last position of the tract."""
        self.assertEqual(norm.right_align(self.reference, "chrT", 40, "TA", "T"), (47, "AA", "A"))
        self.assertEqual(norm.right_align(self.reference, "chrT", 63, "TCA", "T"), (73, "ACA", "A"))

    def test_sequence_context_names_the_locus(self) -> None:
        """The context reader separates a homopolymer, a tandem repeat and ordinary sequence."""
        self.assertEqual(norm.sequence_context(self.reference, "chrT", 40, "TA"), "homopolymer_ge6")
        self.assertEqual(norm.sequence_context(self.reference, "chrT", 63, "TCA"), "short_tandem_repeat")
        self.assertEqual(norm.sequence_context(self.reference, "chrT", 90, "TC"), "other")

    def test_complex_allele_trims_on_both_sides(self) -> None:
        """A complex allele keeps its two changed bases and loses its shared bases."""
        self.assertEqual(norm.normalize(self.reference, "chrT", 60, "TCGT", "TAT"), (61, "CG", "A"))
        self.assertEqual(norm.variant_class("CG", "A"), "complex")


class CheckTest(unittest.TestCase):
    """The writer check normalizes a submission file and records every row that changed."""

    def setUp(self) -> None:
        """Write the fixture reference."""
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.reference = norm.Reference(write_reference(self.root))
        self.addCleanup(self.reference.close)

    def run_check(self, variants: list[tuple[str, int, str, str]]) -> dict:
        """Write a submission and run the check with the internal normalizer.

        Args:
            variants: Chromosome, position, REF and ALT tuples, one per row.

        Returns:
            The check report.
        """
        path = self.root / "submission.csv"
        write_submission(path, variants)
        return norm.check_submission(norm.read_submission(path), self.reference, None, None)

    def test_each_representation_class_is_recorded(self) -> None:
        """Every non-canonical representation changes and the record names the class."""
        cases = [(("chrT", 47, "AA", "A"), "deletion", "homopolymer_ge6"),
                 (("chrT", 48, "A", "AA"), "insertion", "homopolymer_ge6"),
                 (("chrT", 73, "ACA", "A"), "deletion", "short_tandem_repeat"),
                 (("chrT", 75, "A", "ACA"), "insertion", "short_tandem_repeat"),
                 (("chrT", 39, "GTAA", "GTA"), "deletion", "homopolymer_ge6"),
                 (("chrT", 90, "TCT", "AGT"), "MNV", "other"),
                 (("chrT", 60, "TCGT", "TAT"), "complex", "short_tandem_repeat")]
        report = self.run_check([case[0] for case in cases])
        self.assertEqual(report["variants"], len(cases))
        self.assertEqual(report["changed"], len(cases))
        for index, (_, variant_class, context) in enumerate(cases):
            change = report["changes"][index]
            self.assertEqual(change["row"], index + 1)
            self.assertEqual(change["variant_class"], variant_class)
            self.assertEqual(change["sequence_context"], context)
            self.assertTrue(change["reference_matches_ref"])

    def test_a_canonical_file_reports_no_change(self) -> None:
        """A file that already holds the normalized form reports zero changes."""
        report = self.run_check([("chrT", 40, "TA", "T"), ("chrT", 63, "T", "TCA"),
                                 ("chrT", 90, "TC", "AG"), ("chrT", 5, "A", "G")])
        self.assertEqual(report["changed"], 0)
        self.assertEqual(report["backend"], "internal")

    def test_the_check_emits_the_normalized_rows(self) -> None:
        """The normalized rows carry the left aligned and parsimonious tuple."""
        report = self.run_check([("chrT", 47, "AA", "A")])
        row = report["normalized_rows"][0]
        self.assertEqual((row["chrom_1"], row["pos_1"], row["ref_1"], row["alt_1"]),
                         ("chrT", "40", "TA", "T"))
        self.assertEqual(row["epcr"], "0.95")
        self.assertEqual(row["notes"], "fixture")

    def test_both_variant_slots_are_checked(self) -> None:
        """A compound heterozygous row carries two variants and the check reads both."""
        path = self.root / "pair.csv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=norm.HEADER)
            writer.writeheader()
            writer.writerow({"proband_id": "PROBAND01", "chrom_1": "chrT", "pos_1": 47,
                             "ref_1": "AA", "alt_1": "A", "chrom_2": "chrT", "pos_2": 73,
                             "ref_2": "ACA", "alt_2": "A", "epcr": "0.95",
                             "finding_type": "primary", "notes": "pair"})
        report = norm.check_submission(norm.read_submission(path), self.reference, None, None)
        self.assertEqual(report["variants"], 2)
        self.assertEqual(report["changed"], 2)
        self.assertEqual([change["slot"] for change in report["changes"]], [1, 2])

    def test_a_second_representation_is_reported(self) -> None:
        """A row whose right aligned form sits elsewhere is named with its shift."""
        report = self.run_check([("chrT", 40, "TA", "T"), ("chrT", 90, "TC", "AG")])
        self.assertEqual(report["with_a_second_representation"], 1)
        alternative = report["alternatives"][0]
        self.assertEqual(alternative["left_aligned"], "chrT:40:TA:T")
        self.assertEqual(alternative["right_aligned"], "chrT:47:AA:A")
        self.assertEqual(alternative["shift_bp"], 7)

    def test_a_wrong_reference_base_is_reported(self) -> None:
        """A row whose REF does not match the reference is flagged and not hidden."""
        report = self.run_check([("chrT", 47, "GG", "G")])
        self.assertEqual(report["changed"], 1)
        self.assertFalse(report["changes"][0]["reference_matches_ref"])


class WriterTest(unittest.TestCase):
    """The B0 writer runs the check when a reference is given."""

    def setUp(self) -> None:
        """Write the fixture reference and an annotated VCF with one right aligned deletion."""
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.reference_path = write_reference(self.root)
        csq = ("##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from "
               "Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|PICK|gnomADg_AF|"
               "gnomADe_AF|CLIN_SIG\">\n")
        lines = [csq, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"]
        for pos, ref, alt in [(47, "AA", "A"), (73, "ACA", "A")]:
            info = f"CSQ={alt}|frameshift|HIGH|GENEA|G1|1|||pathogenic"
            lines.append("\t".join(["chrT", str(pos), ".", ref, alt, "50", "PASS", info,
                                    "GT:GQ", "0/1:50"]) + "\n")
        self.vcf = self.root / "case.vcf"
        self.vcf.write_text("".join(lines))

    def test_the_writer_emits_the_normalized_form(self) -> None:
        """The submission carries the normalized tuple and the record names the two changes."""
        sys.argv = ["b0_rank.py", "--vcf", str(self.vcf), "--out", str(self.root / "out.csv"),
                    "--reference", str(self.reference_path), "--norm-backend", "internal",
                    "--norm-report", str(self.root / "report.json"),
                    "--summary", str(self.root / "summary.json")]
        import b0_rank
        self.assertEqual(b0_rank.main(), 0)
        with (self.root / "out.csv").open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual((rows[0]["chrom_1"], rows[0]["pos_1"], rows[0]["ref_1"], rows[0]["alt_1"]),
                         ("chrT", "40", "TA", "T"))
        self.assertEqual((rows[0]["chrom_2"], rows[0]["pos_2"], rows[0]["ref_2"], rows[0]["alt_2"]),
                         ("chrT", "63", "TCA", "T"))
        report = json.loads((self.root / "report.json").read_text())
        self.assertEqual(report["changed"], 2)
        summary = json.loads((self.root / "summary.json").read_text())
        self.assertEqual(summary["normalization_changed_variants"], 2)


class RecoveryTest(unittest.TestCase):
    """The benchmark names the rule that recovers each representation class."""

    def setUp(self) -> None:
        """Write the fixture reference."""
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.reference = norm.Reference(write_reference(self.root))
        self.addCleanup(self.reference.close)

    def test_right_alignment_needs_left_alignment(self) -> None:
        """Parsimony alone cannot recover a right aligned indel."""
        canonical = (40, "TA", "T")
        self.assertEqual(bench.recovery_of(self.reference, "chrT", canonical, (47, "AA", "A")),
                         "left_alignment_alone")

    def test_padding_needs_parsimony(self) -> None:
        """A padded allele at the canonical position is recovered by parsimony alone."""
        canonical = (90, "TC", "AG")
        self.assertEqual(bench.recovery_of(self.reference, "chrT", canonical, (90, "TCT", "AGT")),
                         "parsimony_alone")

    def test_a_padded_and_shifted_allele_needs_both_rules(self) -> None:
        """An allele that is padded and not left aligned needs the two rules together."""
        canonical = (40, "TA", "T")
        self.assertEqual(bench.recovery_of(self.reference, "chrT", canonical, (46, "AAA", "AA")),
                         "both_rules")

    def test_a_merged_neighbour_is_never_recovered(self) -> None:
        """An allele merged with a neighbouring substitution is a different allele."""
        canonical = (40, "TA", "T")
        merged = bench.merge_with_neighbour(self.reference, "chrT", *canonical)
        self.assertEqual(bench.recovery_of(self.reference, "chrT", canonical, merged), "none")

    def test_the_generator_covers_the_representation_classes(self) -> None:
        """Every generated class is a valid representation and its name is known."""
        found = bench.representations(self.reference, "chrT", 63, "TCA", "T")
        self.assertIn("left_parsimonious", found)
        self.assertIn("right_aligned", found)
        self.assertIn("padded_both_sides", found)
        self.assertIn("repeat_unit_spelled", found)
        for name in found:
            self.assertIn(name, bench.REPRESENTATIONS)


if __name__ == "__main__":
    unittest.main(verbosity=2)
