"""Verify the B0 ranker and its submission file against the vendored official scorer, without network."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "pipelines/track1_b0/bin"
sys.path.insert(0, str(BIN))
import b0_rank  # noqa: E402
import hackathon_evaluation as evaluation  # noqa: E402

CSQ = ("##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. "
       "Format: Allele|Consequence|IMPACT|SYMBOL|Gene|PICK|gnomADg_AF|gnomADe_AF|CLIN_SIG\">\n")


def vcf_line(chrom: str, pos: int, ref: str, alt: str, gt: str, impact: str, gene: str, af: str, clin: str) -> str:
    """Build one annotated VCF line.

    Args:
        chrom: Chromosome with chr prefix.
        pos: Position.
        ref: Reference allele.
        alt: Alternate allele.
        gt: Genotype.
        impact: VEP impact.
        gene: Gene symbol.
        af: gnomAD genome AF text.
        clin: CLIN_SIG text.

    Returns:
        The VCF line.
    """
    csq = f"{alt}|x|{impact}|{gene}|G|1|{af}||{clin}"
    return "\t".join([chrom, str(pos), f"{chrom}_{pos}", ref, alt, "50", "PASS", f"CSQ={csq}",
                      "GT:GQ", f"{gt}:50"]) + "\n"


class B0Test(unittest.TestCase):
    """Rank a synthetic case with a pathogenic pair, a common variant and a benign gene."""

    def setUp(self) -> None:
        """Write a small annotated VCF."""
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.vcf = self.root / "case.vcf"
        lines = ["##fileformat=VCFv4.2\n", CSQ, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                 vcf_line("chr15", 40185238, "TG", "T", "0/1", "HIGH", "BUB1B", "", "pathogenic"),
                 vcf_line("chr15", 40212554, "G", "A", "0/1", "MODERATE", "BUB1B", "0.00001", "pathogenic"),
                 vcf_line("chr1", 1000, "A", "G", "0/1", "HIGH", "COMMON", "0.2", ""),
                 vcf_line("chr2", 2000, "C", "T", "0/1", "MODERATE", "LONE", "0.0005", ""),
                 vcf_line("chr3", 3000, "G", "C", "0/1", "LOW", "WEAK1", "0.005", ""),
                 vcf_line("chr3", 3100, "G", "T", "0/1", "LOW", "WEAK1", "0.005", "")]
        self.vcf.write_text("".join(lines))

    def test_pair_ranks_first_and_scores_100(self) -> None:
        """The pathogenic BUB1B pair is rank 1 and the official scorer gives 100 rank points."""
        out = self.root / "submission.csv"
        records = b0_rank.parse_vcf(self.vcf, {})
        rows = b0_rank.rows_for(b0_rank.hypotheses(records), "PROBAND01")
        with out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=b0_rank.HEADER)
            writer.writeheader()
            writer.writerows(rows)
        submission = evaluation.load_submission(str(out))
        truth = frozenset([("chr15", 40185238, "TG", "T"), ("chr15", 40212554, "G", "A")])
        result = evaluation.score_proband("PROBAND01", submission["PROBAND01"], truth)
        self.assertEqual(result.full_match_rank, 1)
        self.assertEqual(result.rank_points, 100.0)
        self.assertLessEqual(len(rows), 10)

    def test_common_variant_is_dropped(self) -> None:
        """A variant with gnomAD AF above 0.01 is not a candidate."""
        records = b0_rank.parse_vcf(self.vcf, {})
        self.assertNotIn("COMMON", {record["gene"] for record in records})

    def test_literature_allele_counts_as_high(self) -> None:
        """A MODIFIER literature allele becomes a candidate with impact 3."""
        lit_vcf = self.root / "lit.vcf"
        lit_vcf.write_text("".join(["##fileformat=VCFv4.2\n", CSQ, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                                    vcf_line("chr15", 40117088, "G", "A", "0/1", "MODIFIER", "BUB1B", "", "")]))
        self.assertEqual(b0_rank.parse_vcf(lit_vcf, {}), [])
        records = b0_rank.parse_vcf(lit_vcf, {("chr15", 40117088, "G", "A"): "BUB1B"})
        self.assertEqual(records[0]["impact"], 3.0)

    def test_intergenic_literature_allele_joins_its_gene(self) -> None:
        """An allele with no VEP gene still pairs with a variant of its literature gene."""
        vcf = self.root / "upstream.vcf"
        upstream = "\t".join(["chr15", "40117088", "x", "G", "A", "50", "PASS", "CSQ=A|intergenic_variant|MODIFIER||||||",
                              "GT:GQ", "0/1:50"]) + "\n"
        vcf.write_text("".join(["##fileformat=VCFv4.2\n", CSQ, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                                upstream, vcf_line("chr15", 40210185, "ACT", "A", "0/1", "HIGH", "BUB1B", "", "")]))
        records = b0_rank.parse_vcf(vcf, {("chr15", 40117088, "G", "A"): "BUB1B"})
        pairs = [item for item in b0_rank.hypotheses(records) if len(item["variants"]) == 2]
        self.assertEqual(pairs[0]["gene"], "BUB1B")


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
