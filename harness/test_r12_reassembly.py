"""Verify the R12 targeted local reassembly on a synthetic class A fixture.

The test builds a small reference FASTA with its index, a target BED file and a PAF file with `cs`
tags that carry a known insertion and a known deletion. Every coordinate is invented. No patient
file and no real alignment enters the test.

The test also verifies the target list builder on a small GENCODE extract, because the strand rule
of the 50 kb upstream pad and the 10 kb downstream pad decides which reads the pipeline reads.
"""

from __future__ import annotations

import gzip
import json
import random
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipelines" / "local_reassembly" / "bin"))
sys.path.insert(0, str(ROOT / "scripts"))

import build_r12_targets as targets_builder  # noqa: E402
import r12_reassembly as reassembly  # noqa: E402

CONTIG = "chrT"
CONTIG_LENGTH = 6000
LINE_WIDTH = 60

INSERTED = "ACGT" * 50           # 200 base pairs
DELETED = "TTGCA" * 60           # 300 base pairs

# A synthetic GENCODE extract with one gene on each strand.
GTF_LINES = [
    "##description: synthetic test extract",
    'chr15\ttest\tgene\t400000\t460000\t.\t+\t.\tgene_id "G1"; gene_type "protein_coding"; gene_name "BUB1B";',
    'chr11\ttest\tgene\t300000\t340000\t.\t+\t.\tgene_id "G2"; gene_type "protein_coding"; gene_name "CEP57";',
    'chr5\ttest\tgene\t200000\t220000\t.\t+\t.\tgene_id "G3"; gene_type "protein_coding"; gene_name "TRIP13";',
    'chr2\ttest\tgene\t500000\t540000\t.\t-\t.\tgene_id "G4"; gene_type "protein_coding"; gene_name "BUB1";',
    'chr10\ttest\tgene\t700000\t760000\t.\t+\t.\tgene_id "G5"; gene_type "protein_coding"; gene_name "BUB3";',
    'chr1\ttest\tgene\t800000\t804000\t.\t+\t.\tgene_id "G6"; gene_type "protein_coding"; gene_name "CDC20";',
    'chr9\ttest\tgene\t900000\t940000\t.\t-\t.\tgene_id "G7"; gene_type "protein_coding"; gene_name "PTPA";',
    'chr15_alt_random\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "G8"; gene_type "protein_coding"; gene_name "BUB1B";',
]

KNOWN_BED = "chr15\t410000\t411000\tbub1b_known_one\nchr15\t420000\t420600\tbub1b_known_two\n"


def write_reference(path: Path) -> str:
    """Write a synthetic reference FASTA with its FASTA index.

    Args:
        path: The FASTA path. The index goes to `<path>.fai`.

    Returns:
        The sequence of the single contig.
    """
    generator = random.Random(12345)
    sequence = "".join(generator.choice("ACGT") for _ in range(CONTIG_LENGTH))
    lines = [sequence[start:start + LINE_WIDTH] for start in range(0, CONTIG_LENGTH, LINE_WIDTH)]
    header = f">{CONTIG}\n"
    path.write_text(header + "\n".join(lines) + "\n")
    Path(str(path) + ".fai").write_text(
        f"{CONTIG}\t{CONTIG_LENGTH}\t{len(header)}\t{LINE_WIDTH}\t{LINE_WIDTH + 1}\n")
    return sequence


def paf_line(query: str, target_start: int, target_end: int, cs_value: str, mapq: int = 60,
             alignment_type: str = "P") -> str:
    """Build one PAF line with a `cs` tag.

    Args:
        query: The contig name.
        target_start: The 0-based reference start.
        target_end: The 0-based exclusive reference end.
        cs_value: The value of the `cs` tag.
        mapq: The mapping quality.
        alignment_type: The value of the `tp` tag.

    Returns:
        The PAF line without a line break.
    """
    length = target_end - target_start
    return "\t".join([query, str(length), "0", str(length), "+", CONTIG, str(CONTIG_LENGTH),
                      str(target_start), str(target_end), str(length), str(length), str(mapq),
                      f"tp:A:{alignment_type}", f"cs:Z:{cs_value}"])


class TestCsParsing(unittest.TestCase):
    """The `cs` tag reader must return every operation in order."""

    def test_operations_split_in_order(self) -> None:
        """The reader splits a mixed `cs` value into its operations."""
        operations = reassembly.parse_cs(":100*ac+GGG-TTTT:50")
        self.assertEqual(operations, [(":", "100"), ("*", "ac"), ("+", "GGG"),
                                      ("-", "TTTT"), (":", "50")])

    def test_a_substitution_moves_one_base(self) -> None:
        """A substitution advances the reference position by one base."""
        alignment = reassembly.parse_paf_line(paf_line("c1", 0, 2000, f":500*ac:499+{INSERTED}:1000"))
        events = reassembly.events_of_alignment(alignment, 50)
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0]["position"], 1000)


class TestEventReading(unittest.TestCase):
    """The event reader must find the insertion and the deletion with the right size."""

    def setUp(self) -> None:
        """Build one alignment that holds one insertion and one deletion."""
        cs_value = f":1000+{INSERTED}:1000-{DELETED}:1000"
        self.alignment = reassembly.parse_paf_line(paf_line("c1", 0, 3300, cs_value))
        self.events = reassembly.events_of_alignment(self.alignment, 50)

    def test_both_events_are_found(self) -> None:
        """The reader finds two events."""
        self.assertEqual(len(self.events), 2)

    def test_the_insertion_carries_its_sequence(self) -> None:
        """The insertion event holds the inserted sequence and its length."""
        insertion = next(item for item in self.events if item["svtype"] == "INS")
        self.assertEqual(insertion["position"], 1000)
        self.assertEqual(insertion["size"], len(INSERTED))
        self.assertEqual(insertion["sequence"], INSERTED)

    def test_the_deletion_advances_the_reference(self) -> None:
        """The deletion event sits after the insertion and moves the reference position."""
        deletion = next(item for item in self.events if item["svtype"] == "DEL")
        self.assertEqual(deletion["position"], 2000)
        self.assertEqual(deletion["size"], len(DELETED))

    def test_a_small_event_is_dropped(self) -> None:
        """An event below the size limit does not reach the output."""
        alignment = reassembly.parse_paf_line(paf_line("c1", 0, 2000, ":1000+ACGTACGT:1000"))
        self.assertEqual(reassembly.events_of_alignment(alignment, 50), [])

    def test_an_event_at_the_edge_is_dropped(self) -> None:
        """An event inside the edge margin does not reach the output."""
        alignment = reassembly.parse_paf_line(paf_line("c1", 0, 1050, f":20+{INSERTED}:1030"))
        self.assertEqual(reassembly.events_of_alignment(alignment, 50), [])

    def test_a_secondary_alignment_is_dropped(self) -> None:
        """A secondary alignment and a low mapping quality do not reach the reader."""
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "one.paf"
            path.write_text("\n".join([
                paf_line("c1", 0, 3000, f":1000+{INSERTED}:2000", 60, "S"),
                paf_line("c2", 0, 3000, f":1000+{INSERTED}:2000", 3, "P"),
            ]) + "\n")
            self.assertEqual(reassembly.read_alignments(path), [])


class TestMerging(unittest.TestCase):
    """Two redundant contigs that report one allele must produce one record."""

    def test_two_contigs_of_one_allele_merge(self) -> None:
        """Two events inside the merge window and the size ratio become one record."""
        cs_one = f":1000+{INSERTED}:1000"
        cs_two = f":1040+{INSERTED[:190]}:1000"
        alignments = [reassembly.parse_paf_line(paf_line("c1", 0, 2200, cs_one)),
                      reassembly.parse_paf_line(paf_line("c2", 0, 2230, cs_two))]
        events = []
        for alignment in alignments:
            events.extend(reassembly.events_of_alignment(alignment, 50))
        merged = reassembly.merge_events(events)
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0]["support"], 2)
        self.assertEqual(merged[0]["size"], len(INSERTED))

    def test_two_distinct_alleles_stay_apart(self) -> None:
        """Two events far apart stay as two records."""
        cs_value = f":1000+{INSERTED}:1000-{DELETED}:1000"
        alignment = reassembly.parse_paf_line(paf_line("c1", 0, 3300, cs_value))
        merged = reassembly.merge_events(reassembly.events_of_alignment(alignment, 50))
        self.assertEqual(len(merged), 2)


class TestVcfOutput(unittest.TestCase):
    """The call step must write a valid, sequence-resolved VCF inside the target list."""

    def setUp(self) -> None:
        """Build the fixture and run the call step and the merge step."""
        self.folder = tempfile.TemporaryDirectory()
        base = Path(self.folder.name)
        self.reference = base / "ref.fa"
        self.sequence = write_reference(self.reference)
        self.bed = base / "targets.bed"
        self.bed.write_text(f"{CONTIG}\t500\t4000\tTESTGENE\tgene\n")
        self.paf = base / "TESTGENE.paf"
        self.paf.write_text(paf_line("c1", 0, 3300, f":1000+{INSERTED}:1000-{DELETED}:1000") + "\n")
        self.calls = base / "calls"
        self.out_vcf = self.calls / "TESTGENE.vcf"
        reassembly.call_interval(self.paf, self.reference, self.bed, "TESTGENE", 50,
                                 self.out_vcf, "TESTSAMPLE")
        self.records = [line for line in self.out_vcf.read_text().splitlines()
                        if not line.startswith("#")]
        self.merged_vcf = base / "out" / "r12_calls.vcf"
        self.gates_path = base / "out" / "class_a_gates.json"
        reassembly.merge_vcf_files(self.calls, self.reference, self.bed, "TESTSAMPLE", self.merged_vcf)

    def tearDown(self) -> None:
        """Remove the temporary directory."""
        self.folder.cleanup()

    def test_two_records_are_written(self) -> None:
        """The call step writes one record per event."""
        self.assertEqual(len(self.records), 2)

    def test_the_insertion_record_holds_the_sequence(self) -> None:
        """The insertion record carries the anchor base and the inserted sequence in ALT."""
        record = next(line for line in self.records if "SVTYPE=INS" in line)
        fields = record.split("\t")
        self.assertEqual(int(fields[1]), 1000)
        self.assertEqual(fields[3], self.sequence[999])
        self.assertEqual(fields[4], self.sequence[999] + INSERTED)
        self.assertIn(f"SVLEN={len(INSERTED)}", fields[7])

    def test_the_deletion_record_holds_the_reference_span(self) -> None:
        """The deletion record carries the deleted reference sequence in REF."""
        record = next(line for line in self.records if "SVTYPE=DEL" in line)
        fields = record.split("\t")
        self.assertEqual(int(fields[1]), 2000)
        self.assertEqual(fields[3], self.sequence[1999:2000 + len(DELETED)])
        self.assertEqual(fields[4], self.sequence[1999])
        self.assertIn(f"SVLEN=-{len(DELETED)}", fields[7])

    def test_a_call_outside_the_target_list_is_dropped(self) -> None:
        """An event outside every target interval does not reach the output."""
        base = Path(self.folder.name)
        narrow = base / "narrow.bed"
        narrow.write_text(f"{CONTIG}\t500\t1500\tNARROW\tgene\n")
        out = base / "narrow.vcf"
        reassembly.call_interval(self.paf, self.reference, narrow, "NARROW", 50, out, "TESTSAMPLE")
        kept = [line for line in out.read_text().splitlines() if not line.startswith("#")]
        self.assertEqual(len(kept), 1)
        self.assertIn("SVTYPE=INS", kept[0])

    def test_the_genotype_reads_homozygous_with_one_contig(self) -> None:
        """One spanning contig that carries the event gives a homozygous genotype."""
        record = next(line for line in self.records if "SVTYPE=INS" in line)
        self.assertEqual(record.split("\t")[9].split(":")[0], "1/1")

    def test_a_second_reference_contig_gives_a_heterozygous_genotype(self) -> None:
        """A second spanning contig without the event gives a heterozygous genotype."""
        base = Path(self.folder.name)
        paf = base / "HET.paf"
        paf.write_text("\n".join([
            paf_line("c1", 0, 3300, f":1000+{INSERTED}:1000-{DELETED}:1000"),
            paf_line("c2", 0, 3000, ":3000"),
        ]) + "\n")
        out = base / "het.vcf"
        reassembly.call_interval(paf, self.reference, self.bed, "HET", 50, out, "TESTSAMPLE")
        record = next(line for line in out.read_text().splitlines() if "SVTYPE=INS" in line)
        self.assertEqual(record.split("\t")[9].split(":")[0], "0/1")

    def test_the_merge_gate_reads_pass(self) -> None:
        """The merge step writes a gate file with the PASS verdict."""
        self.assertTrue(self.merged_vcf.exists())
        merged = [line for line in self.merged_vcf.read_text().splitlines()
                  if not line.startswith("#")]
        self.assertEqual(len(merged), 2)


class TestTargetList(unittest.TestCase):
    """The target list builder must pad each gene in the direction of its strand."""

    def setUp(self) -> None:
        """Build the synthetic gene model and the target list."""
        self.folder = tempfile.TemporaryDirectory()
        base = Path(self.folder.name)
        self.gtf = base / "gencode.gtf.gz"
        with gzip.open(self.gtf, "wt") as handle:
            handle.write("\n".join(GTF_LINES) + "\n")
        self.known = base / "known.bed"
        self.known.write_text(KNOWN_BED)
        self.merged, self.gene_table = targets_builder.build_targets(self.gtf, self.known)
        self.manifest = targets_builder.write_outputs(
            self.merged, self.gene_table, self.gtf, self.known,
            base / "targets.bed", base / "manifest.json")

    def tearDown(self) -> None:
        """Remove the temporary directory."""
        self.folder.cleanup()

    def test_every_mva_gene_is_present(self) -> None:
        """The builder finds all seven genes."""
        self.assertEqual(len(self.gene_table), 7)

    def test_a_forward_gene_pads_upstream_of_the_start(self) -> None:
        """A gene on the plus strand takes 50 kb before the start and 10 kb after the end."""
        self.assertEqual(self.gene_table["TRIP13"]["target"], "chr5:150000-230000")

    def test_a_reverse_gene_pads_upstream_of_the_end(self) -> None:
        """A gene on the minus strand takes 50 kb after the end and 10 kb before the start."""
        self.assertEqual(self.gene_table["BUB1"]["target"], "chr2:490000-590000")

    def test_a_scaffold_record_is_skipped(self) -> None:
        """A gene record on a scaffold does not replace the primary record."""
        self.assertTrue(self.gene_table["BUB1B"]["gene_body"].startswith("chr15:"))

    def test_the_known_intervals_merge_into_the_bub1b_target(self) -> None:
        """The two known intervals sit inside the BUB1B target and carry the known label."""
        record = next(item for item in self.merged if item["contig"] == "chr15")
        self.assertTrue(record["known"])
        self.assertIn("bub1b_known_one", record["name"])
        self.assertIn("bub1b_known_two", record["name"])

    def test_the_manifest_records_the_gene_model_checksum(self) -> None:
        """The manifest carries a SHA-256 of the gene model and the total size."""
        self.assertEqual(len(self.manifest["gene_model"]["sha256"]), 64)
        self.assertEqual(self.manifest["interval_count"], 7)
        self.assertGreater(self.manifest["total_bp"], 0)


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=0).result
    print(f"{result.testsRun - len(result.failures) - len(result.errors)} of {result.testsRun} checks passed.")
    sys.exit(0 if result.wasSuccessful() else 1)
