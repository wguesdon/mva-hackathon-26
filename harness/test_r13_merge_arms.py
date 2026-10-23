"""Verify the R13 merge module on a synthetic class A fixture.

The test covers the tag handling and the command that the module builds. It does not run SVDB,
because a unit test must not need a container. Every coordinate is invented and no patient file
enters the test.
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipelines" / "local_reassembly" / "bin"))

import r13_merge_arms as merge  # noqa: E402

HEADER = [
    "##fileformat=VCFv4.1",
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
]

PIPELINE_RECORDS = [
    "chrT\t1000\tcall_one\tN\t<DEL>\t.\tPASS\tEND=1500;SVTYPE=DEL;SVLEN=-500;set=tiddit;"
    "svdb_origin=tiddit\tGT\t0/1",
    "chrT\t4000\tcall_two\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=300;set=tiddit-manta;"
    "svdb_origin=tiddit|manta\tGT\t1/1",
]

REASSEMBLY_RECORDS = [
    "chrT\t4010\tr12_one\tN\tNACGT\t.\tPASS\tSVTYPE=INS;SVLEN=300\tGT\t1/1",
]

MERGED_RECORDS = [
    "chrT\t1000\tcall_one\tN\t<DEL>\t.\tPASS\tEND=1500;SVTYPE=DEL;caller_set=tiddit;set=pipeline\tGT\t0/1",
    "chrT\t4000\tcall_two\tN\t<INS>\t.\tPASS\tSVTYPE=INS;caller_set=tiddit-manta;"
    "set=pipeline|reassembly\tGT\t1/1",
]


def write_vcf(path: Path, records: list[str]) -> None:
    """Write a small VCF fixture.

    Args:
        path: The output path.
        records: The data lines.
    """
    path.write_text("\n".join(HEADER + records) + "\n")


class MergeArmsTest(unittest.TestCase):
    """Check the tag handling, the command and the composition count."""

    def setUp(self) -> None:
        """Write the three fixture files."""
        self.folder = Path(tempfile.mkdtemp())
        self.pipeline = self.folder / "pipeline.vcf"
        self.reassembly = self.folder / "reassembly.vcf"
        self.merged = self.folder / "merged.vcf"
        write_vcf(self.pipeline, PIPELINE_RECORDS)
        write_vcf(self.reassembly, REASSEMBLY_RECORDS)
        write_vcf(self.merged, MERGED_RECORDS)

    def test_the_svdb_tags_are_renamed(self) -> None:
        """The `set` and `svdb_origin` entries become `caller_set` and `caller_origin`."""
        line = merge.retag_info("END=1500;SVTYPE=DEL;set=tiddit;svdb_origin=tiddit", "pipeline")
        self.assertIn("caller_set=tiddit", line)
        self.assertIn("caller_origin=tiddit", line)
        self.assertNotIn(";set=", ";" + line)
        self.assertNotIn("svdb_origin=", line)

    def test_the_arm_is_named_on_every_record(self) -> None:
        """A record with no SVDB tag takes the arm name as its caller."""
        line = merge.retag_info("SVTYPE=INS;SVLEN=300", "reassembly")
        self.assertIn("caller_arm=reassembly", line)
        self.assertIn("caller_set=reassembly", line)

    def test_a_flag_entry_survives_the_rename(self) -> None:
        """An INFO flag with no value stays a flag."""
        line = merge.retag_info("IMPRECISE;SVTYPE=DEL;set=manta", "pipeline")
        self.assertIn("IMPRECISE", line)
        self.assertNotIn("IMPRECISE=", line)

    def test_the_retagged_file_keeps_every_record(self) -> None:
        """The rewrite keeps the record count and adds the INFO definitions."""
        target = self.folder / "arm_pipeline.vcf"
        count = merge.retag_file(self.pipeline, target, merge.ARM_PIPELINE)
        text = target.read_text()
        self.assertEqual(count, len(PIPELINE_RECORDS))
        self.assertIn("ID=caller_set", text)
        self.assertIn("ID=caller_arm", text)
        self.assertEqual(text.count("caller_arm=pipeline"), len(PIPELINE_RECORDS))

    def test_the_merge_command_carries_the_pinned_parameters(self) -> None:
        """The command names the overlap, the breakpoint window, the priority and --no_intra."""
        command = merge.merge_command("a.vcf.gz", "b.vcf.gz", merge.MERGE_OVERLAP,
                                      merge.MERGE_BND_DISTANCE, merge.MERGE_INS_DISTANCE)
        self.assertEqual(command[0], "svdb")
        self.assertIn("--pass_only", command)
        self.assertIn("--same_order", command)
        self.assertIn("--no_intra", command)
        self.assertEqual(command[command.index("--overlap") + 1], "0.5")
        self.assertEqual(command[command.index("--bnd_distance") + 1], "150")
        self.assertEqual(command[command.index("--ins_distance") + 1], "50")
        self.assertEqual(command[command.index("--priority") + 1], "pipeline,reassembly")
        self.assertIn("a.vcf.gz:pipeline", command)
        self.assertIn("b.vcf.gz:reassembly", command)

    def test_the_gate_file_carries_no_count(self) -> None:
        """The class A parameter section holds the thresholds and no record count."""
        report = {"tool": "SVDB 2.8.4", "arguments": ["--merge"], "reciprocal_overlap": 0.5,
                  "breakpoint_window_bp": 150, "insertion_distance_bp": 50,
                  "source_caller_kept_in": "caller_set", "input_records": {"pipeline": 7},
                  "merged_records": 7, "records_removed_by_the_merge": 0,
                  "composition": {"pipeline": 7}}
        section = merge.parameter_section(report)
        self.assertEqual(section["reciprocal_overlap"], 0.5)
        for key in ("input_records", "merged_records", "records_removed_by_the_merge",
                    "composition"):
            self.assertNotIn(key, section)

    def test_the_composition_counts_each_arm(self) -> None:
        """The composition reads the `set` entry that SVDB writes on a merged record."""
        counts = merge.composition(self.merged)
        self.assertEqual(counts, {"pipeline": 1, "pipeline|reassembly": 1})


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=0).result
    print(f"{result.testsRun - len(result.failures) - len(result.errors)} of {result.testsRun} checks passed.")
    sys.exit(0 if result.wasSuccessful() else 1)
