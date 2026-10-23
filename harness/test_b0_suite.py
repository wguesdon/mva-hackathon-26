"""Run b0_suite.py end to end on a tiny synthetic suite, without network or VEP."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "pipelines/track1_b0/bin"
CSQ = ("##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. "
       "Format: Allele|Consequence|IMPACT|SYMBOL|Gene|PICK|gnomADg_AF|gnomADe_AF|CLIN_SIG\">\n")
HEAD = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"


def line(chrom: str, pos: int, ref: str, alt: str, impact: str, gene: str, clin: str = "") -> str:
    """Build one annotated heterozygous VCF line with the suite ID form.

    Args:
        chrom: Chromosome.
        pos: Position.
        ref: Reference allele.
        alt: Alternate allele.
        impact: VEP impact.
        gene: Gene symbol.
        clin: CLIN_SIG.

    Returns:
        The line.
    """
    return "\t".join([chrom, str(pos), f"{chrom}_{pos}_{ref}_{alt}", ref, alt, "50", "PASS",
                      f"CSQ={alt}|x|{impact}|{gene}|G|1|||{clin}", "GT:GQ", "0/1:50"]) + "\n"


class SuiteTest(unittest.TestCase):
    """Two scenarios: a pathogenic pair above the noise (score 12 against 6) and a weak pair behind the noise (score 5 against 6)."""

    def test_suite_receipt(self) -> None:
        """The receipt holds both arms, the checks and per-scenario ranks from the official scorer."""
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            background = root / "bg.vcf"
            background.write_text("##fileformat=VCFv4.2\n" + CSQ + HEAD
                                  + line("chr1", 100, "A", "T", "MODERATE", "NOISE")
                                  + line("chr1", 200, "C", "G", "MODERATE", "NOISE"))
            spikes = root / "spikes.vcf"
            spikes.write_text("##fileformat=VCFv4.2\n" + CSQ + HEAD
                              + line("chr2", 10, "G", "A", "HIGH", "GENEA", "pathogenic")
                              + line("chr2", 900, "T", "C", "HIGH", "GENEA", "pathogenic")
                              + line("chr3", 10, "G", "A", "MODERATE", "GENEB")
                              + line("chr3", 900, "T", "C", "LOW", "GENEB"))
            manifest = root / "manifest.json"
            manifest.write_text(json.dumps({"suite": "tiny", "scenarios": [
                {"scenario_id": "spk_001", "group": "recessive_clinvar", "gene": "GENEA",
                 "variants": [{"chrom": "chr2", "pos": 10, "ref": "G", "alt": "A"}, {"chrom": "chr2", "pos": 900, "ref": "T", "alt": "C"}]},
                {"scenario_id": "spk_002", "group": "recessive_clinvar", "gene": "GENEB",
                 "variants": [{"chrom": "chr3", "pos": 10, "ref": "G", "alt": "A"}, {"chrom": "chr3", "pos": 900, "ref": "T", "alt": "C"}]}]}))
            literature = root / "lit.tsv"
            literature.write_text("allele_id\tgene\tstatus\tchrom\tpos\tref\talt\n")
            out = root / "out"
            result = subprocess.run([sys.executable, str(BIN / "b0_suite.py"), "--background", str(background),
                                     "--spikes", str(spikes), "--manifest", str(manifest), "--literature", str(literature),
                                     "--min-mean-rank-points", "25", "--min-pair-top10", "0.5", "--out-dir", str(out)],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            receipt = json.loads((out / "s4_b0_receipt.json").read_text())
            ranks = {row["scenario_id"]: row for row in receipt["per_scenario"]["b0"]}
            self.assertEqual(ranks["spk_001"]["full_match_rank"], 1)
            self.assertEqual(ranks["spk_002"]["full_match_rank"], 2)
            self.assertEqual(receipt["arms"]["b0"]["metrics"]["mean_rank_points"], 75.0)
            self.assertEqual(receipt["verdict"], "PASS")


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
