"""Verify pipelines/s6_summary/bin/summarize_s6_patient.py on synthetic class A input.

The test builds a small forced VCF, a scan table and a control receipt in a temporary directory.
No patient data enters the test.
"""

from __future__ import annotations

import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipelines" / "s6_summary" / "bin"))
import summarize_s6_patient as summary  # noqa: E402

TARGETS = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr15\t40117088\tlit:bub1b_upstream_rs576524605\tG\tA\t.\t.\t.
chr15\t40160000\tClinVar:1:BUB1B\tC\tT\t.\t.\t.
chr5\t900524\tClinVar:2:TRIP13\tT\tTA\t.\t.\t.
"""

FORCED = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tproband
chr15\t40117088\t.\tG\tA\t45.2\t.\t.\tGT:DP:AD\t0/1:42:21,21
chr15\t40160000\t.\tC\tT\t0.5\t.\t.\tGT:DP:AD\t0/0:33:33,0
chr5\t900524\t.\tT\tTA\t0.3\t.\t.\tGT:DP:AD\t0/0:4:4,0
"""

SCAN = """interval\tbin_start\tbin_end\treads\tsoftclip_20bp\tdiscordant_other_chrom\tsupplementary_sa
bub1b_intron2_alu_sva_insertion\t40165196\t40165246\t50\t0\t0\t0
bub1b_intron2_alu_sva_insertion\t40165246\t40165296\t48\t1\t0\t1
bub1b_intron2_alu_sva_insertion\t40165296\t40165346\t52\t30\t12\t9
"""

CONTROL = {"insertion_scan": {"bub1b_intron2_alu_sva_insertion": {
    "bins": 98, "reads": 1234, "softclip_20bp": 45, "discordant_other_chrom": 23,
    "supplementary_sa": 12, "retroseq_calls": ["chr15:40166000"]}}}


def write_inputs(folder: Path) -> dict[str, Path]:
    """Write the synthetic input files.

    Args:
        folder: A temporary directory.

    Returns:
        The path of each input.
    """
    paths = {"targets": folder / "targets.vcf", "forced": folder / "proband_forced.vcf.gz",
             "scan": folder / "scan.tsv", "mei": folder / "mei.txt", "control": folder / "control.json"}
    paths["targets"].write_text(TARGETS)
    with gzip.open(paths["forced"], "wt") as handle:
        handle.write(FORCED)
    paths["scan"].write_text(SCAN)
    paths["mei"].write_text("")
    paths["control"].write_text(json.dumps(CONTROL))
    return paths


class SummaryTest(unittest.TestCase):
    """Check the genotype classes, the target names, the peak bins and the class B content."""

    def test_carrier_state(self) -> None:
        """Each genotype string maps to one state."""
        self.assertEqual(summary.carrier_state("0/0"), "hom_ref")
        self.assertEqual(summary.carrier_state("0|1"), "het")
        self.assertEqual(summary.carrier_state("1/1"), "hom_alt")
        self.assertEqual(summary.carrier_state("./."), "no_call")

    def test_peak_bins_finds_the_high_bin_only(self) -> None:
        """A bin at three times the median and above the floor is a peak."""
        bins = summary.parse_scan(self.paths["scan"])["bub1b_intron2_alu_sva_insertion"]
        self.assertEqual(summary.peak_bins(bins, "softclip_20bp"), ["40165296-40165346"])
        self.assertEqual(summary.peak_bins(bins, "supplementary_sa"), ["40165296-40165346"])

    def test_target_name_comes_from_position(self) -> None:
        """The forced VCF lost the ID column, so the name comes from the target VCF."""
        records = summary.parse_forced(self.paths["forced"], summary.target_names(self.paths["targets"]))
        self.assertEqual(records[0]["target_id"], "lit:bub1b_upstream_rs576524605")
        self.assertEqual(records[0]["gt"], "0/1")
        self.assertEqual(records[2]["dp"], 4)

    def test_summary_reports_the_carrier_the_low_depth_site_and_the_control(self) -> None:
        """The class B summary names the het site, flags the low depth site and quotes the control call."""
        records = summary.parse_forced(self.paths["forced"], summary.target_names(self.paths["targets"]))
        text = "\n".join(summary.summary_lines(records, summary.parse_scan(self.paths["scan"]), [],
                                               json.loads(self.paths["control"].read_text())))
        self.assertIn("lit:bub1b_upstream_rs576524605 | chr15:40117088 | G>A | 0/1 | 42 | 21,21", text)
        self.assertIn("Sites below the depth floor", text)
        self.assertIn("ClinVar:2:TRIP13 | chr5:900524 | 4", text)
        self.assertIn("placed no call inside the two intervals", text)
        self.assertIn("chr15:40166000", text)
        self.assertNotIn("| ClinVar:1:BUB1B | chr15:40160000 | C>T |", text)

    def setUp(self) -> None:
        """Write the synthetic inputs for one test."""
        self.folder = tempfile.TemporaryDirectory()
        self.paths = write_inputs(Path(self.folder.name))

    def tearDown(self) -> None:
        """Remove the temporary directory."""
        self.folder.cleanup()


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
