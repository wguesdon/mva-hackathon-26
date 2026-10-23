"""Gate G1 of ADR 0008. Every class A writer of a class B2 pipeline passes the publish gate.

The test builds the class A record of each writer on synthetic inputs, then runs
`privacy_guard.class_a_publish_problems` over the bytes that the writer would publish. A writer that
adds a count, a size or any key outside its declared list fails this test.

Three writers expose the record builder as a function and the test calls it. Two writers build the
record inside a long main, and the test uses a fixture whose key set is quoted from the writer with
the file and the reason.
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import privacy_guard

ROOT = Path(__file__).resolve().parent.parent
RULES = privacy_guard.load_publish_rules(ROOT / "literature" / "class_a_publish_rules.json")

sys.path.insert(0, str(ROOT / "pipelines" / "phase_pair" / "bin"))
sys.path.insert(0, str(ROOT / "pipelines" / "local_reassembly" / "bin"))

import phase_pairs  # noqa: E402  the pipeline bin path is set above
import r12_reassembly  # noqa: E402
import r13_merge_arms  # noqa: E402


def gate(name: str, record: dict) -> list[str]:
    """Run the class A publish gate over a record that a writer would publish.

    Args:
        name: The published file name.
        record: The class A record of the writer.

    Returns:
        The problems that the gate reports.
    """
    return privacy_guard.class_a_publish_problems(Path(name), json.dumps(record, indent=2), RULES)


class WriterTest(unittest.TestCase):
    """Each writer that exposes its record builder passes the gate."""

    def test_phase_pairs_writer_passes(self) -> None:
        """`pipelines/phase_pair/bin/phase_pairs.py` builds a record that passes."""
        with tempfile.TemporaryDirectory() as name:
            work = Path(name)
            for stem in ("phased.vcf", "slot1.csv", "class_b2_phase_rows.tsv",
                         "class_b2_phase_rows.md"):
                (work / stem).write_text("synthetic\n")
            tools = {"whatshap": "2.8", "whatshap_image": "quay.io/biocontainers/whatshap:2.8",
                     "method_record": "docs/benchmarks/t1_phase_v1.md"}
            record = phase_pairs.build_gates(
                "phase_pair_patient_v1", work / "phased.vcf", work / "slot1.csv",
                work / "class_b2_phase_rows.tsv", work / "class_b2_phase_rows.md", tools)
            self.assertEqual(gate("class_a_phase_gates.json", record), [])
            self.assertEqual(record["verdict"], "COMPLETE")

    def test_r13_merge_writer_passes(self) -> None:
        """`r13_merge_arms.build_class_a_gates` keeps every count out of the class A record."""
        report = {"tool": "SVDB 2.8.4", "arguments": "--no_intra", "reciprocal_overlap": 0.5,
                  "breakpoint_window_bp": 150, "insertion_distance_bp": 50,
                  "source_caller_kept_in": "caller_set",
                  "input_records": {r13_merge_arms.ARM_PIPELINE: 41,
                                    r13_merge_arms.ARM_REASSEMBLY: 37},
                  "merged_records": 63,
                  "records_removed_by_the_merge": 15,
                  "composition": {"arm_pipeline": 20}}
        record = r13_merge_arms.build_class_a_gates(report, "r13_merged_sv_v1")
        self.assertEqual(gate("class_a_merge_gates.json", record), [])
        self.assertNotIn("input_records", json.dumps(record))
        self.assertNotIn("merged_records", json.dumps(record))

    def test_r12_reassembly_writer_passes(self) -> None:
        """`r12_reassembly.build_class_a_gates` carries public target quantities only."""
        targets = [{"start": 40161000, "end": 40162000}, {"start": 40200000, "end": 40200500}]
        record = r12_reassembly.build_class_a_gates(targets, 2, True, "r12_local_reassembly_v1")
        self.assertEqual(gate("class_a_reassembly_gates.json", record), [])

    def test_sv_mei_writer_schema_passes(self) -> None:
        """The key set of `annotate_sv_mei.py` passes.

        The record is built inside the main of that script. The fixture below quotes its key set
        from `pipelines/sv_mei_annotate/bin/annotate_sv_mei.py`, where the report dictionary names
        run_id, data_class_of_run, data_class_of_this_file, scope, inputs, references, thresholds,
        outputs, gates and verdict.
        """
        record = {
            "run_id": "r13_merged_sv_v1", "data_class_of_run": "C",
            "data_class_of_this_file": "A",
            "scope": "Checksums, tool versions, reference counts, thresholds and process gates.",
            "inputs": {"sv_vcf": {"name": "class_c_merged_sv.vcf", "sha256": "0" * 64,
                                  "header": ["##fileformat=VCFv4.2"]}},
            "references": {"gene_models": {"name": "gencode.gtf", "sha256": "1" * 64,
                                           "mva_genes_found": ["BUB1B"], "mva_genes_missing": []},
                           "gnomad_sv": {"name": "gnomad_sv.vcf", "sha256": "2" * 64,
                                         "records": 1199117},
                           "known_intervals": {"name": "intervals.bed", "sha256": "3" * 64,
                                               "count": 2}},
            "thresholds": {"min_sv_length_bp": 50, "min_genotype_quality": 20,
                           "max_gnomad_af": 0.01, "max_submission_rows": 10,
                           "mva_genes": ["BUB1B"]},
            "outputs": {"class_c": ["class_c_sv_mei_annotated.tsv"],
                        "class_b2": ["class_b2_sv_mei_ranked.tsv"],
                        "note": "Every class B2 file is written on every run."},
            "gates": {"gene_models_loaded": "PASS", "class_b2_written": "PASS"},
            "verdict": "PASS",
        }
        self.assertEqual(gate("class_a_gates.json", record), [])

    def test_b0_format_check_schema_passes(self) -> None:
        """The repaired key set of the format check of `pipelines/track1_b0/main.nf` passes."""
        record = {"file": "slot1/slot1_b0_phenotype.csv", "sha256": "0" * 64, "data_class": "A",
                  "scope": "The shape of the file, not its content.", "row_budget": 3,
                  "format_check": "PASS", "proband_ids": ["PROBAND01"],
                  "rows_present": "PASS", "rows_within_budget": "PASS",
                  "epcr_values_are_distinct": True,
                  "epcr_minimum_above_the_structural_ceiling": True}
        self.assertEqual(gate("format_check.json", record), [])

    def test_undeclared_b0_format_fields_fail(self) -> None:
        """Reject undeclared row counts and EPCR values in a synthetic gate record."""
        # All counts and values below are invented for this fixture.
        record = {"file": "slot1/slot1_b0_phenotype.csv", "sha256": "0" * 64, "data_class": "A",
                  "scope": "The shape of the file, not its content.",
                  "format_check": "PASS", "proband_ids": ["PROBAND01"], "rows": 7,
                  "epcr_values": ["0.71", "0.70", "0.69"],
                  "epcr_values_are_distinct": True}
        problems = gate("format_check.json", record)
        self.assertTrue(any("rows" in problem for problem in problems))
        self.assertTrue(any("epcr_values" in problem for problem in problems))

    def test_every_pipeline_process_runs_the_gate(self) -> None:
        """Each workflow that publishes a class A gate file calls the gate in the same task."""
        expected = {
            "pipelines/phase_pair/main.nf": ["out/class_a_phase_gates.json"],
            "pipelines/local_reassembly/main.nf": ["out/class_a_reassembly_gates.json",
                                                   "out/class_a_merge_gates.json",
                                                   "out/class_a_gates.json"],
            "pipelines/sv_mei_annotate/main.nf": ["out/class_a_gates.json"],
            "pipelines/track1_b0/main.nf": ["slot1/format_check.json",
                                            "phenotype/phenotype_extraction_gate.json"],
        }
        for workflow, targets in expected.items():
            flat = " ".join((ROOT / workflow).read_text().replace("\\\\\n", " ").split())
            for target in targets:
                with self.subTest(workflow=workflow, target=target):
                    call = f"--class-a-publish --rules ${{publish_rules}} {target}"
                    self.assertIn(call, flat, f"{workflow} does not gate {target}.")


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
