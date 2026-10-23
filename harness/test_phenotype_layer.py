"""Verify the Resnik and best match average phenotype layer on a small synthetic ontology."""

from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "pipelines/track1_b0/bin"
sys.path.insert(0, str(BIN))
import b0_rank  # noqa: E402
import phenotype_layer  # noqa: E402

# A four level chain and one side branch.
#   HP:0000001 root
#     HP:0000118 phenotypic abnormality
#       HP:0000100 organ abnormality        -> genes A, B, C, D
#         HP:0000200 specific finding       -> genes A, B
#           HP:0000300 very specific        -> gene A
#         HP:0000400 other finding          -> genes C, D
OBO = """format-version: 1.2
data-version: hp/releases/2026-09-01

[Term]
id: HP:0000001
name: All

[Term]
id: HP:0000118
name: Phenotypic abnormality
is_a: HP:0000001 ! All

[Term]
id: HP:0000100
name: Organ abnormality
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:0000200
name: Specific finding
is_a: HP:0000100 ! Organ abnormality

[Term]
id: HP:0000300
name: Very specific finding
is_a: HP:0000200 ! Specific finding

[Term]
id: HP:0000400
name: Other finding
is_a: HP:0000100 ! Organ abnormality

[Term]
id: HP:0000999
name: Retired finding
is_obsolete: true
replaced_by: HP:0000400
"""

GENES = """ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id
1\tGENEA\tHP:0000300\tVery specific finding\t-\tOMIM:1
2\tGENEB\tHP:0000200\tSpecific finding\t-\tOMIM:2
3\tGENEC\tHP:0000400\tOther finding\t-\tOMIM:3
4\tGENED\tHP:0000400\tOther finding\t-\tOMIM:4
"""

CSQ = ("##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. "
       "Format: Allele|Consequence|IMPACT|SYMBOL|Gene|PICK|gnomADg_AF|gnomADe_AF|CLIN_SIG\">\n")


def vcf_line(chrom: str, pos: int, ref: str, alt: str, gene: str, impact: str) -> str:
    """Build one annotated VCF line with no population frequency.

    Args:
        chrom: Chromosome with a chr prefix.
        pos: Position.
        ref: Reference allele.
        alt: Alternate allele.
        gene: Gene symbol.
        impact: VEP impact.

    Returns:
        The VCF line.
    """
    csq = f"{alt}|x|{impact}|{gene}|G|1|||"
    return "\t".join([chrom, str(pos), f"{chrom}_{pos}", ref, alt, "50", "PASS", f"CSQ={csq}",
                      "GT:GQ", "0/1:50"]) + "\n"


class PhenotypeLayerTest(unittest.TestCase):
    """Check the ontology parse, the information content and the score of the layer."""

    def setUp(self) -> None:
        """Write the synthetic ontology and the synthetic gene annotation."""
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.obo = self.root / "hp.obo"
        self.obo.write_text(OBO)
        self.genes = self.root / "genes_to_phenotype.txt"
        self.genes.write_text(GENES)
        self.model = phenotype_layer.build_model(self.obo, self.genes)

    def test_obsolete_term_maps_to_its_replacement(self) -> None:
        """A retired identifier resolves to the live term."""
        self.assertEqual(self.model.normalize(["HP:0000999"]), ["HP:0000400"])

    def test_information_content_follows_the_annotation_counts(self) -> None:
        """A term that every gene carries has no information content."""
        self.assertAlmostEqual(self.model.content["HP:0000100"], 0.0)
        self.assertAlmostEqual(self.model.content["HP:0000200"], -math.log(2 / 4))
        self.assertAlmostEqual(self.model.content["HP:0000300"], -math.log(1 / 4))

    def test_the_exact_match_scores_one(self) -> None:
        """A gene whose profile equals the case profile reaches the largest score."""
        score = self.model.score_gene("GENEA", ["HP:0000300"])
        self.assertAlmostEqual(score, 1.0)

    def test_a_related_gene_scores_between_the_bounds(self) -> None:
        """A gene on the same branch scores above a gene on another branch and below the exact match."""
        exact = self.model.score_gene("GENEA", ["HP:0000300"])
        related = self.model.score_gene("GENEB", ["HP:0000300"])
        unrelated = self.model.score_gene("GENEC", ["HP:0000300"])
        self.assertLess(unrelated, related)
        self.assertLess(related, exact)
        self.assertAlmostEqual(unrelated, 0.0)

    def test_the_score_is_symmetric_in_the_best_match_average(self) -> None:
        """The two directions of the best match average carry the same weight."""
        forward = self.model.score_gene("GENEB", ["HP:0000300"])
        backward = self.model.score_gene("GENEA", ["HP:0000200"])
        self.assertAlmostEqual(forward, backward)

    def test_an_unannotated_gene_scores_zero(self) -> None:
        """A gene with no phenotype annotation contributes nothing."""
        self.assertEqual(self.model.score_gene("GENEZ", ["HP:0000300"]), 0.0)


class RankerIntegrationTest(unittest.TestCase):
    """Check that the layer is off by default and that it orders two equal genes when it is on."""

    def setUp(self) -> None:
        """Write a VCF with one heterozygous pair in each of two genes."""
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.obo = self.root / "hp.obo"
        self.obo.write_text(OBO)
        self.genes = self.root / "genes_to_phenotype.txt"
        self.genes.write_text(GENES)
        self.vcf = self.root / "case.vcf"
        lines = ["##fileformat=VCFv4.2\n", CSQ,
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                 vcf_line("chr1", 1000, "A", "G", "GENEC", "HIGH"),
                 vcf_line("chr1", 2000, "A", "G", "GENEC", "HIGH"),
                 vcf_line("chr2", 1000, "A", "G", "GENEA", "HIGH"),
                 vcf_line("chr2", 2000, "A", "G", "GENEA", "HIGH")]
        self.vcf.write_text("".join(lines))
        self.records = b0_rank.parse_vcf(self.vcf, {})

    def test_the_default_run_ignores_the_phenotype(self) -> None:
        """Without a bonus map the two genes keep the same score and the position decides."""
        found = b0_rank.hypotheses(self.records)
        pairs = [item for item in found if len(item["variants"]) == 2]
        self.assertEqual(len(pairs), 2)
        self.assertEqual(pairs[0]["score"], pairs[1]["score"])
        self.assertEqual(pairs[0]["gene"], "GENEC")

    def test_the_layer_orders_two_genes_with_the_same_variant_score(self) -> None:
        """With the layer on, the gene that matches the case profile comes first."""
        model = phenotype_layer.build_model(self.obo, self.genes)
        bonus = {gene: 2.0 * model.score_gene(gene, ["HP:0000300"])
                 for gene in {record["gene"] for record in self.records}}
        found = b0_rank.hypotheses(self.records, bonus)
        pairs = [item for item in found if len(item["variants"]) == 2]
        self.assertEqual(pairs[0]["gene"], "GENEA")
        self.assertEqual(pairs[0]["variant_score"], pairs[1]["variant_score"])
        self.assertGreater(pairs[0]["phenotype_points"], pairs[1]["phenotype_points"])


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
