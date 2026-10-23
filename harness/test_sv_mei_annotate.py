"""Verify pipelines/sv_mei_annotate/bin/annotate_sv_mei.py on synthetic class A input.

The test builds a small GENCODE extract, a small gnomAD SV extract, a structural variant VCF and a
mobile element VCF in a temporary directory. Every coordinate is invented. No patient data enters the
test.
"""

from __future__ import annotations

import csv
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipelines" / "sv_mei_annotate" / "bin"))
import annotate_sv_mei as annotate  # noqa: E402

# A synthetic BUB1B on chr15 with three exons, a synthetic CEP57 on chr1 and one filler gene.
GTF = "\n".join([
    '##description: synthetic test extract',
    'chr15\ttest\tgene\t40160000\t40220000\t.\t+\t.\tgene_id "G1"; gene_type "protein_coding"; gene_name "BUB1B";',
    'chr15\ttest\texon\t40160000\t40160500\t.\t+\t.\tgene_id "G1"; gene_name "BUB1B"; tag "MANE_Select";',
    'chr15\ttest\texon\t40185000\t40185400\t.\t+\t.\tgene_id "G1"; gene_name "BUB1B"; tag "MANE_Select";',
    'chr15\ttest\texon\t40219000\t40220000\t.\t+\t.\tgene_id "G1"; gene_name "BUB1B"; tag "MANE_Select";',
    'chr1\ttest\tgene\t95000000\t95040000\t.\t-\t.\tgene_id "G2"; gene_type "protein_coding"; gene_name "CEP57";',
    'chr1\ttest\texon\t95039000\t95040000\t.\t-\t.\tgene_id "G2"; gene_name "CEP57"; tag "MANE_Select";',
    'chr5\ttest\tgene\t900000\t950000\t.\t+\t.\tgene_id "G3"; gene_type "protein_coding"; gene_name "TRIP13";',
    'chr5\ttest\texon\t900000\t901000\t.\t+\t.\tgene_id "G3"; gene_name "TRIP13"; tag "MANE_Select";',
    'chr8\ttest\tgene\t1000000\t1400000\t.\t+\t.\tgene_id "G4"; gene_type "protein_coding"; gene_name "FILLER";',
    'chr2\ttest\tgene\t10000\t20000\t.\t+\t.\tgene_id "G5"; gene_type "protein_coding"; gene_name "BUB1";',
    'chr10\ttest\tgene\t10000\t20000\t.\t+\t.\tgene_id "G6"; gene_type "protein_coding"; gene_name "BUB3";',
    'chr1\ttest\tgene\t43000\t53000\t.\t+\t.\tgene_id "G7"; gene_type "protein_coding"; gene_name "CDC20";',
    'chr9\ttest\tgene\t10000\t20000\t.\t+\t.\tgene_id "G8"; gene_type "protein_coding"; gene_name "PTPA";',
]) + "\n"

# A common deletion that matches call 4, and filler records so that the frequency gate passes.
GNOMAD_HEAD = "\n".join([
    '##fileformat=VCFv4.2',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
]) + "\n"
GNOMAD_COMMON = ("chr15\t40186000\tgnomad_common_del\tN\t<DEL>\t.\tPASS\t"
                 "END=40187000;SVTYPE=DEL;SVLEN=1000;AF=0.22;AF_nfe=0.25\n")

# Call 1 is a mobile element insertion in BUB1B intron 2, inside the known interval.
# Call 2 is a mobile element insertion in a BUB1B intron outside the known intervals.
# Call 3 is a deletion that removes a BUB1B exon.
# Call 4 is a deletion in a BUB1B intron that is common in gnomAD.
# Call 5 is a large rare deletion outside every MVA gene.
# Call 6 is a homozygous reference call that the genotype gate must drop.
SV_VCF = "\n".join([
    '##fileformat=VCFv4.2',
    '##source=GenerateSVCandidates 1.6.0',
    '##cmdline=/opt/manta/bin/configManta.py --bam proband_sorted_md.bam',
    '##reference=file:///references/Homo_sapiens_assembly38.fasta',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tproband',
    'chr15\t40185100\tsv_exon_del\tN\t<DEL>\t320\tPASS\tSVTYPE=DEL;END=40185300;SVLEN=-200\tGT:GQ\t0/1:60',
    'chr15\t40186100\tsv_common_del\tN\t<DEL>\t300\tPASS\tSVTYPE=DEL;END=40187100;SVLEN=-1000\tGT:GQ\t0/1:55',
    'chr8\t1100000\tsv_big_del\tN\t<DEL>\t400\tPASS\tSVTYPE=DEL;END=1300000;SVLEN=-200000\tGT:GQ\t0/1:70',
    'chr15\t40190000\tsv_homref\tN\t<DEL>\t10\tPASS\tSVTYPE=DEL;END=40190400;SVLEN=-400\tGT:GQ\t0/0:30',
    'chr15\t40170000\tsv_low_gq\tN\t<DUP>\t10\tPASS\tSVTYPE=DUP;END=40171000;SVLEN=1000\tGT:GQ\t0/1:5',
]) + "\n"

MEI_VCF = "\n".join([
    '##fileformat=VCFv4.2',
    '##source=RetroSeq v1.5',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tproband',
    'chr15\t40166377\tmei_known\tN\t<INS:ME:ALU>\t.\tPASS\tSVTYPE=INS;MEINFO=ALU,40166377,40166677,+\tGT:GQ:FL\t0/1:40:7',
    'chr15\t40190500\tmei_new\tN\t<INS:ME:SVA>\t.\tPASS\tSVTYPE=INS;MEINFO=SVA,40190500,40190800,+\tGT:GQ:FL\t0/1:38:7',
    'chr15\t40195000\tmei_lowlevel\tN\t<INS:ME:ALU>\t.\tPASS\tSVTYPE=INS;MEINFO=ALU,40195000,40195300,+\tGT:GQ:FL\t0/1:35:2',
]) + "\n"

INTERVALS = ("chr15\t40165196\t40170061\tbub1b_intron2_alu_sva_insertion\n"
             "chr15\t40196226\t40196842\tbub1b_intron8_alu_insertion\n")


def build_inputs(folder: Path) -> dict[str, Path]:
    """Write every synthetic input file.

    Args:
        folder: A temporary directory.

    Returns:
        The path of each input.
    """
    paths = {"gtf": folder / "genes.gtf.gz", "gnomad": folder / "gnomad.vcf.gz",
             "sv": folder / "sv.vcf.gz", "mei": folder / "mei.vcf.gz",
             "intervals": folder / "intervals.bed", "out": folder / "out"}
    with gzip.open(paths["gtf"], "wt") as handle:
        handle.write(GTF)
    filler = "".join(f"chr21\t{10000 + index * 300}\tfill{index}\tN\t<DEL>\t.\tPASS\t"
                     f"END={10200 + index * 300};SVTYPE=DEL;SVLEN=200;AF=0.30\n"
                     for index in range(50))
    with gzip.open(paths["gnomad"], "wt") as handle:
        handle.write(GNOMAD_HEAD + GNOMAD_COMMON + filler)
    with gzip.open(paths["sv"], "wt") as handle:
        handle.write(SV_VCF)
    with gzip.open(paths["mei"], "wt") as handle:
        handle.write(MEI_VCF)
    paths["intervals"].write_text(INTERVALS)
    return paths


def run_script(paths: dict[str, Path]) -> dict:
    """Run the annotation script on the synthetic input.

    Args:
        paths: The input paths.

    Returns:
        The parsed class A gate file.
    """
    argv = ["annotate_sv_mei.py", "--sv-vcf", str(paths["sv"]), "--mei-vcf", str(paths["mei"]),
            "--gtf", str(paths["gtf"]), "--gnomad-sv", str(paths["gnomad"]),
            "--intervals", str(paths["intervals"]), "--out-dir", str(paths["out"]),
            "--deny-tokens", "proband", "--run-id", "test_run", "--min-frequency-records", "10"]
    saved = sys.argv
    sys.argv = argv
    try:
        code = annotate.main()
    finally:
        sys.argv = saved
    if code != 0:
        raise AssertionError("The script returned a FAIL verdict on the synthetic input.")
    return json.loads((paths["out"] / "class_a_gates.json").read_text())


def read_candidates(folder: Path) -> list[dict]:
    """Read the ranked candidate table.

    Args:
        folder: The output directory.

    Returns:
        One row for each candidate.
    """
    with (folder / "class_b2_sv_mei_ranked.tsv").open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class AnnotateTest(unittest.TestCase):
    """Check the parser, the filters, the ranking and the class A boundary."""

    @classmethod
    def setUpClass(cls) -> None:
        """Build the synthetic input once and run the script."""
        cls.folder = Path(tempfile.mkdtemp())
        cls.paths = build_inputs(cls.folder)
        cls.gates = run_script(cls.paths)
        cls.candidates = read_candidates(cls.paths["out"])
        cls.by_id = {row["call_id"]: row for row in cls.candidates}

    def test_every_gate_passes(self) -> None:
        """Every process gate reads PASS on valid input."""
        self.assertEqual(self.gates["verdict"], "PASS")
        self.assertTrue(all(value == "PASS" for value in self.gates["gates"].values()))

    def test_gate_file_holds_no_patient_count(self) -> None:
        """The class A gate file holds no count of patient calls and no patient coordinate."""
        text = json.dumps(self.gates)
        self.assertNotIn("40185100", text)
        self.assertNotIn("sv_exon_del", text)
        self.assertNotIn("candidates", text)
        self.assertNotIn("proband", text)

    def test_exon_deletion_outranks_intron_insertion(self) -> None:
        """A deletion that removes a BUB1B exon scores above an intronic insertion."""
        self.assertIn("sv_exon_del", self.by_id)
        self.assertIn("mei_new", self.by_id)
        self.assertGreater(int(self.by_id["sv_exon_del"]["score"]), int(self.by_id["mei_new"]["score"]))
        self.assertEqual(self.by_id["sv_exon_del"]["gene_effect"], "exonic")
        self.assertEqual(self.by_id["mei_new"]["gene_effect"], "intronic")

    def test_known_interval_is_marked(self) -> None:
        """A call inside a known BUB1B interval carries the known label."""
        self.assertTrue(self.by_id["mei_known"]["known_status"].startswith("known_interval:"))
        self.assertEqual(self.by_id["mei_new"]["known_status"], "new")

    def test_common_call_is_removed(self) -> None:
        """A call that matches a common gnomAD record does not reach the candidate table."""
        self.assertNotIn("sv_common_del", self.by_id)

    def test_low_quality_calls_are_removed(self) -> None:
        """A homozygous reference call, a low genotype quality call and a low RetroSeq level call go."""
        for call_id in ("sv_homref", "sv_low_gq", "mei_lowlevel"):
            self.assertNotIn(call_id, self.by_id)

    def test_context_list_holds_the_large_event(self) -> None:
        """The large rare deletion outside the MVA genes appears in the context section."""
        text = (self.paths["out"] / "class_b2_sv_mei_ranked.md").read_text()
        self.assertIn("sv_big_del", text)
        self.assertIn("## 2. The largest rare events outside the MVA genes", text)

    def test_submission_proposal_is_always_written(self) -> None:
        """The proposal files exist on every run, so their presence states nothing about the child."""
        self.assertTrue((self.paths["out"] / "class_b2_submission_rows.csv").exists())
        self.assertTrue((self.paths["out"] / "class_b2_submission_rows.md").exists())

    def test_header_extract_drops_the_command_line(self) -> None:
        """The tool header extract keeps the source line and drops the command line."""
        header = self.gates["inputs"]["sv_vcf"]["header"]
        self.assertTrue(any("GenerateSVCandidates" in line for line in header))
        self.assertFalse(any("configManta" in line for line in header))

    def test_checksums_are_recorded(self) -> None:
        """Every input and every reference carries a SHA-256 checksum."""
        self.assertEqual(len(self.gates["inputs"]["sv_vcf"]["sha256"]), 64)
        self.assertEqual(len(self.gates["inputs"]["mei_vcf"]["sha256"]), 64)
        self.assertEqual(len(self.gates["references"]["gnomad_sv"]["sha256"]), 64)

    def test_all_mva_genes_are_found(self) -> None:
        """The gene model gate finds every one of the seven MVA genes."""
        self.assertEqual(self.gates["references"]["gene_models"]["mva_genes_missing"], [])
        self.assertEqual(len(self.gates["references"]["gene_models"]["mva_genes_found"]), 7)

    def test_two_calls_of_one_score_do_not_share_one_epcr(self) -> None:
        """Item B13. A tied row falls by one rounding unit, so the values stay distinct."""
        tied = [{"score": 94, "known_status": "new", "chrom": "chr15", "pos": 40200000 + step,
                 "end": 40200500 + step, "ref": "N", "alt": "<DEL>", "svtype": "DEL",
                 "svlen": -500, "genotype": "het", "mva_gene": "BUB1B",
                 "gene_effect": "exonic"} for step in (0, 1000, 2000)]
        rows = annotate.submission_rows(tied)
        values = [row["epcr"] for row in rows]
        self.assertEqual(len(set(values)), len(values))
        self.assertEqual(values, [0.77, 0.769, 0.768])
        self.assertEqual(values, sorted(values, reverse=True))

    def test_the_reachable_epcr_ceiling_is_below_the_shipped_b0_ladder(self) -> None:
        """The largest score of 107 maps onto 0.835, which sits below the B0 value of 0.97."""
        best = {"mva_gene": "BUB1B", "gene_effect": "exonic", "source": "sv",
                "gnomad_popmax_af": None, "genotype": "hom_alt", "filter": "PASS"}
        score, _ = annotate.score_call(best)
        self.assertEqual(score, 107)
        self.assertEqual(annotate.epcr_of_score(score), annotate.EPCR_CEILING_OF_THIS_WRITER)
        self.assertEqual(annotate.epcr_of_score(score), 0.835)
        self.assertLess(annotate.epcr_of_score(score), 0.97)

    def test_the_two_collisions_with_the_b0_suite_ladder(self) -> None:
        """The suite ladder values 0.77 and 0.68 are reachable, and 0.95 and 0.86 are not."""
        self.assertEqual(annotate.epcr_of_score(94), 0.77)
        self.assertEqual(annotate.epcr_of_score(76), 0.68)
        ladder = [round(0.95 - 0.09 * rank, 2) for rank in range(10)]
        reachable = {annotate.epcr_of_score(score)
                     for score in range(annotate.SUBMISSION_SCORE_FLOOR, 108)}
        self.assertEqual(sorted(value for value in ladder if value in reachable), [0.68, 0.77])
        self.assertNotIn(0.95, reachable)
        self.assertNotIn(0.86, reachable)


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=0).result
    print(f"{result.testsRun - len(result.failures) - len(result.errors)} of {result.testsRun} checks passed.")
    sys.exit(0 if result.wasSuccessful() else 1)
