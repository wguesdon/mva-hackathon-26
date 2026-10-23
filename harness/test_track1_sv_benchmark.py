"""Verify the Track 1 structural variant benchmark on a synthetic fixture.

The fixture holds a small GRCh37 truth set, a small high confidence BED, a small GRCh38 call set
and a small Truvari output directory. No file of the test comes from a real sample, so the test
runs on any machine and reads no patient data.
"""

from __future__ import annotations

import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import track1_sv_benchmark as bench  # noqa: E402
import track1_sv_common as common  # noqa: E402

# The stub shifts every GRCh37 coordinate of chr1 by 1000 and refuses chr2 entirely.
SHIFT = 1000


class StubLifter:
    """A chain file stand-in that shifts chr1 and drops every other contig.

    Attributes:
        broken: A source coordinate that the stub refuses, to test the loss counters.
    """

    def __init__(self, broken: int | None = None) -> None:
        """Store the coordinate that must fail to lift.

        Args:
            broken: The zero-based source coordinate that returns no hit.
        """
        self.broken = broken

    def convert_coordinate(self, chrom: str, position: int) -> list[tuple[str, int, str, int]]:
        """Return the shifted coordinate of a chr1 position.

        Args:
            chrom: The source contig name.
            position: The zero-based source coordinate.

        Returns:
            One hit for chr1, and an empty list otherwise.
        """
        if chrom != "chr1" or position == self.broken:
            return []
        return [("chr1", position + SHIFT, "+", 0)]


def write_lines(path: Path, lines: list[str]) -> None:
    """Write text lines to a file.

    Args:
        path: The output path.
        lines: The lines without a line break.
    """
    path.write_text("\n".join(lines) + "\n")


def write_gzip_vcf(path: Path, body: list[str]) -> None:
    """Write a gzipped VCF with a minimal header.

    Args:
        path: The output path.
        body: The data lines without a line break.
    """
    header = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002"]
    with gzip.open(path, "wt") as handle:
        handle.write("\n".join(header + body) + "\n")


def truth_line(pos: int, svtype: str, size: int, filter_field: str = "PASS") -> str:
    """Build one sequence-resolved GRCh37 truth record.

    Args:
        pos: The one-based position.
        svtype: DEL or INS.
        size: The variant size in base pairs.
        filter_field: The FILTER column.

    Returns:
        The tab separated VCF line.
    """
    if svtype == "DEL":
        ref, alt, svlen = "A" * (size + 1), "A", -size
    else:
        ref, alt, svlen = "A", "A" + "T" * size, size
    info = f"SVTYPE={svtype};SVLEN={svlen}"
    return f"1\t{pos}\tid{pos}\t{ref}\t{alt}\t20\t{filter_field}\t{info}\tGT:AD\t0/1:10,10"


class CommonTest(unittest.TestCase):
    """Check the VCF and interval helpers."""

    def test_info_round_trip(self) -> None:
        """A parsed INFO column rebuilds to the same string."""
        text = "SVTYPE=DEL;SVLEN=-120;IMPRECISE;END=500"
        self.assertEqual(common.format_info(common.parse_info(text)), text)

    def test_variant_size_sources(self) -> None:
        """SVLEN wins, END is the fallback and the alleles are the last resort."""
        with_svlen = common.Variant("chr1", 10, ".", "N", "<DEL>", ".", "PASS", {"SVLEN": "-300", "END": "9999"})
        with_end = common.Variant("chr1", 10, ".", "N", "<DUP>", ".", "PASS", {"END": "310"})
        resolved = common.Variant("chr1", 10, ".", "A", "A" + "T" * 60, ".", "PASS", {})
        symbolic = common.Variant("chr1", 10, ".", "N", "<INS:ME>", ".", "PASS", {})
        self.assertEqual(common.variant_size(with_svlen), 300)
        self.assertEqual(common.variant_size(with_end), 300)
        self.assertEqual(common.variant_size(resolved), 60)
        self.assertIsNone(common.variant_size(symbolic))

    def test_variant_type_drops_the_sub_type(self) -> None:
        """A tandem duplication reports the DUP class."""
        record = common.Variant("chr1", 10, ".", "N", "<DUP:TANDEM>", ".", "PASS", {"SVTYPE": "DUP:TANDEM"})
        self.assertEqual(common.variant_type(record), "DUP")

    def test_size_bands(self) -> None:
        """Each size falls in one band and an unknown size is separate."""
        self.assertEqual(common.size_band(50), "50-99")
        self.assertEqual(common.size_band(299), "100-299")
        self.assertEqual(common.size_band(1000), "1000-9999")
        self.assertEqual(common.size_band(50000), "10000+")
        self.assertEqual(common.size_band(None), "unsized")

    def test_intervals(self) -> None:
        """Merging joins the overlap and the index answers a position question."""
        merged = common.merge_intervals([("chr1", 10, 20), ("chr1", 15, 30), ("chr2", 0, 5)])
        self.assertEqual(merged, [("chr1", 10, 30), ("chr2", 0, 5)])
        self.assertEqual(common.total_bases(merged), 25)
        index = common.index_intervals(merged)
        self.assertTrue(common.position_in_index(index, "chr1", 20))
        self.assertFalse(common.position_in_index(index, "chr1", 40))
        self.assertEqual(common.intersect_bases(index, "chr1", 0, 15), 5)


class LiftoverTest(unittest.TestCase):
    """Check the truth and region liftover rules."""

    def test_truth_liftover_keeps_only_scorable_records(self) -> None:
        """A non PASS record, a small record and a record that fails to lift are dropped."""
        records = [
            common.Variant("1", 100, "a", "A" * 201, "A", ".", "PASS", {"SVTYPE": "DEL", "SVLEN": "-200"}, ["GT", "0/1"]),
            common.Variant("1", 500, "b", "A", "A" + "T" * 300, ".", "PASS", {"SVTYPE": "INS", "SVLEN": "300"}, ["GT", "1/1"]),
            common.Variant("1", 900, "c", "A" * 61, "A", ".", "LongReadHomRef", {"SVTYPE": "DEL", "SVLEN": "-60"}, ["GT", "0/1"]),
            common.Variant("1", 1200, "d", "A", "AT", ".", "PASS", {"SVTYPE": "INS", "SVLEN": "1"}, ["GT", "0/1"]),
            common.Variant("2", 100, "e", "A" * 201, "A", ".", "PASS", {"SVTYPE": "DEL", "SVLEN": "-200"}, ["GT", "0/1"]),
        ]
        lifted, counts = bench.lift_truth_records(records, StubLifter())
        self.assertEqual(counts["kept"], 2)
        self.assertEqual(counts["not_pass"], 1)
        self.assertEqual(counts["below_50bp"], 1)
        self.assertEqual(counts["pos_not_lifted"], 1)
        self.assertEqual([record.alt for record in lifted], ["<DEL>", "<INS>"])
        self.assertEqual(lifted[0].pos, 100 + SHIFT)
        self.assertEqual(lifted[0].info["END"], str(100 + SHIFT + 200))
        self.assertEqual(lifted[1].tail, ["GT", "1/1"])

    def test_truth_liftover_drops_a_deletion_with_a_broken_end(self) -> None:
        """A deletion whose far breakpoint does not lift is counted and dropped."""
        record = common.Variant("1", 100, "a", "A" * 201, "A", ".", "PASS",
                                {"SVTYPE": "DEL", "SVLEN": "-200"}, ["GT", "0/1"])
        lifted, counts = bench.lift_truth_records([record], StubLifter(broken=299))
        self.assertEqual(lifted, [])
        self.assertEqual(counts["end_not_lifted"], 1)

    def test_bed_liftover_reports_the_loss(self) -> None:
        """An interval on a contig with no chain is dropped and the bases are recounted."""
        intervals = [("1", 0, 100), ("1", 50, 200), ("2", 0, 100)]
        lifted, counts = bench.lift_bed(intervals, StubLifter())
        self.assertEqual(lifted, [("chr1", SHIFT, SHIFT + 200)])
        self.assertEqual(counts["edge_not_lifted"], 1)
        self.assertEqual(counts["input_bases"], 350)
        self.assertEqual(counts["kept_bases"], 200)


class CallSetTest(unittest.TestCase):
    """Check the call set filter and the measures."""

    def test_prepare_call_set_counts_every_drop(self) -> None:
        """An inversion, a break end, an unsized insertion and a small call are dropped."""
        records = [
            common.Variant("chr1", 10, ".", "N", "<DEL>", ".", "PASS", {"SVTYPE": "DEL", "SVLEN": "-200"}),
            common.Variant("chr1", 20, ".", "N", "<INV>", ".", "PASS", {"SVTYPE": "INV", "SVLEN": "900"}),
            common.Variant("chr1", 30, ".", "N", "N[chr2:100[", ".", "PASS", {"SVTYPE": "BND"}),
            common.Variant("chr1", 40, ".", "N", "<INS>", ".", "PASS", {"SVTYPE": "INS", "END": "40"}),
            common.Variant("chr1", 50, ".", "N", "<DEL>", ".", "PASS", {"SVTYPE": "DEL", "SVLEN": "-10"}),
            common.Variant("chrUn_KI270442v1", 60, ".", "N", "<DEL>", ".", "PASS", {"SVTYPE": "DEL", "SVLEN": "-200"}),
        ]
        kept, counts = bench.prepare_call_set(records)
        self.assertEqual(counts["kept"], 1)
        self.assertEqual(len(kept), 1)
        self.assertEqual(counts["inversion_calls"], 1)
        self.assertEqual(counts["break_end_calls"], 1)
        self.assertEqual(counts["insertion_of_unknown_length"], 1)
        self.assertEqual(counts["off_primary_contig"], 1)
        self.assertEqual(counts["below_50bp"], 2)

    def test_metrics(self) -> None:
        """Recall reads the truth side and precision reads the call side."""
        row = bench.metrics(tp_base=8, false_negative=2, tp_comp=8, false_positive=2)
        self.assertEqual(row["recall"], 0.8)
        self.assertEqual(row["precision"], 0.8)
        self.assertEqual(row["f1"], 0.8)
        self.assertEqual(row["truth_records"], 10)
        empty = bench.metrics(tp_base=0, false_negative=0, tp_comp=0, false_positive=5)
        self.assertIsNone(empty["recall"])
        self.assertIsNone(empty["f1"])
        self.assertEqual(empty["precision"], 0.0)

    def test_positional_reach(self) -> None:
        """Two calls near one truth insertion reach one site, and a far call reaches none."""
        index = common.index_intervals([("chr1", 0, 100000)])
        sites = {"chr1": [5000, 20000]}
        records = [common.Variant("chr1", 5100, ".", "N", "<INS:ME>", ".", "PASS", {}),
                   common.Variant("chr1", 5400, ".", "N", "<INS:ME>", ".", "PASS", {}),
                   common.Variant("chr1", 60000, ".", "N", "<INS:ME>", ".", "PASS", {})]
        inside, matched, reached = bench.positional_reach(records, sites, index)
        self.assertEqual((inside, matched), (3, 2))
        self.assertEqual(reached, {("chr1", 5000)})


class StratifyTest(unittest.TestCase):
    """Check that the stratification reads a Truvari output directory correctly."""

    def setUp(self) -> None:
        """Write a synthetic Truvari output directory."""
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name)
        write_gzip_vcf(self.path / "tp-base.vcf.gz", [
            "chr1\t100\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-200\tGT\t0/1",
            "chr1\t900\t.\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=60\tGT\t0/1"])
        write_gzip_vcf(self.path / "fn.vcf.gz", [
            "chr1\t300\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-200\tGT\t0/1",
            "chr1\t500\t.\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=5000\tGT\t0/1"])
        write_gzip_vcf(self.path / "tp-comp.vcf.gz", [
            "chr1\t100\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-200\tGT\t0/1",
            "chr1\t900\t.\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=60\tGT\t0/1"])
        write_gzip_vcf(self.path / "fp.vcf.gz", [
            "chr1\t700\t.\tN\t<DUP:TANDEM>\t.\tPASS\tSVTYPE=DUP:TANDEM;SVLEN=400\tGT\t0/1",
            "chr1\t800\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-200\tGT\t0/1"])
        (self.path / "summary.json").write_text(json.dumps({"recall": 0.5, "precision": 0.5}))

    def test_stratify(self) -> None:
        """The class table and the band table follow the synthetic counts."""
        result = bench.stratify_bench(self.path)
        deletions = result["by_class"]["DEL"]
        self.assertEqual(deletions["recall"], 0.5)
        self.assertEqual(deletions["precision"], 0.5)
        insertions = result["by_class"]["INS"]
        self.assertEqual(insertions["recall"], 0.5)
        self.assertEqual(insertions["precision"], 1.0)
        duplications = result["by_class"]["DUP"]
        self.assertIsNone(duplications["recall"])
        self.assertEqual(duplications["precision"], 0.0)
        self.assertIn("note", duplications)
        self.assertEqual(result["by_size_band"]["100-299"]["truth_records"], 2)
        self.assertEqual(result["by_class_and_band"]["INS"]["1000-9999"]["recall"], 0.0)


if __name__ == "__main__":
    outcome = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if outcome.wasSuccessful() else 1)
