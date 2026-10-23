"""Test the content and listing rules of the conservative public dispatcher.

The content gate checks fields exported to a class A path. The listing gate checks
whether a restricted-path listing includes object sizes. Every fixture is synthetic.
The tests cover invented record counts, monitor lines, and prefix listings.
Writer schema tests are in harness/test_class_a_writers.py.
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import privacy_guard

ROOT = Path(__file__).resolve().parent.parent
RULES = privacy_guard.load_publish_rules(ROOT / "literature" / "class_a_publish_rules.json")

# The allowed class A gate file of the merge step. It carries parameters, process gates and a
# verdict, and no count.
MERGE_GATES_PASS = {
    "run_id": "r13_merged_sv_v1",
    "data_class_of_run": "C",
    "data_class_of_this_file": "A",
    "scope": "The merge parameters and the process gates.",
    "merge": {"tool": "SVDB 2.8.4", "arguments": "--no_intra", "reciprocal_overlap": 0.5,
              "breakpoint_window_bp": 150, "insertion_distance_bp": 50,
              "source_caller_kept_in": "caller_set"},
    "gates": {"pipeline_vcf_parsed": "PASS", "reassembly_vcf_read": "PASS",
              "merged_vcf_written": "PASS", "source_caller_kept": "PASS"},
    "verdict": "PASS",
}

# This synthetic fixture includes undeclared record counts. All counts are invented.
MERGE_GATES_WITH_COUNTS = dict(MERGE_GATES_PASS)
MERGE_GATES_WITH_COUNTS["merge"] = dict(MERGE_GATES_PASS["merge"])
MERGE_GATES_WITH_COUNTS["merge"]["input_records"] = {"arm_pipeline": 41, "arm_reassembly": 37}
MERGE_GATES_WITH_COUNTS["merge"]["merged_records"] = 63

# A class A gate file that carries a byte size of a class B2 output.
PHASE_GATES_WITH_SIZE = {
    "run_id": "phase_pair_patient_v1",
    "created_utc": "2026-09-17T16:40:00Z",
    "data_class": "A",
    "scope": "Process gates, tool versions and checksums.",
    "gates": {"phased_vcf_present": "PASS", "submission_present": "PASS",
              "pair_analysis_complete": "PASS", "class_b2_written_to_llm_safe": "PASS"},
    "tools": {"whatshap": "2.8"},
    "checksums": {"phased_vcf_sha256": "0" * 64,
                  "class_b2_phase_rows_tsv_size_bytes": 812},
    "verdict": "COMPLETE",
}

PHASE_GATES_PASS = json.loads(json.dumps(PHASE_GATES_WITH_SIZE))
PHASE_GATES_PASS["checksums"] = {"phased_vcf_sha256": "0" * 64,
                                 "class_b2_phase_rows_tsv_sha256": "1" * 64}

# Synthetic output of `aws s3 ls --recursive` over a class B2 prefix. Every key and every size
# below is invented for this fixture.
S3_LISTING_WITH_SIZES = (
    "2026-09-17 19:14:58        731 derived/llm_safe/b0_phenotype_patient_v2/out/class_b2_shipped_rows.csv\n"
    "2026-09-17 19:14:58       2184 derived/llm_safe/b0_phenotype_patient_v2/out/class_b2_shipped_rows.md\n"
)

# The output of the approved listing command. Keys alone, no size.
S3_LISTING_SAFE = (
    "derived/llm_safe/b0_phenotype_patient_v2/out/class_b2_shipped_rows.csv\n"
    "derived/llm_safe/b0_phenotype_patient_v2/out/class_b2_shipped_rows.md\n"
    "2 keys\n"
)

S3_LISTING_HUMAN = (
    "2026-09-17 19:14:58    0.7 KiB derived/llm_safe/b0_phenotype_patient_v2/out/class_b2_shipped_rows.csv\n"
)

S3_LISTING_SUMMARIZE = S3_LISTING_SAFE + "\nTotal Objects: 2\n   Total Size: 2915\n"

MONITOR_LINE_WITH_SIZE = "derived/runs/s6_patient_v1/results/out/proband_mei_in_intervals.txt size 0 bytes\n"

PUBLIC_LISTING_WITH_SIZE = (
    "2026-09-13 08:00:00      12345 reference/GRCh38/GRCh38.fa.fai\n"
)


def gate(name: str, payload: dict) -> list[str]:
    """Run the class A publish gate over one synthetic file.

    Args:
        name: The published file name.
        payload: The content of the file.

    Returns:
        The problems that the gate reports.
    """
    return privacy_guard.class_a_publish_problems(Path(name), json.dumps(payload, indent=2), RULES)


def listing(text: str) -> list[str]:
    """Run the listing gate over one synthetic listing.

    Args:
        text: The listing output.

    Returns:
        The problems that the gate reports.
    """
    import re
    return privacy_guard.listing_problems(
        text, RULES["class_b2_prefixes"], re.compile(RULES["size_key_pattern"], re.IGNORECASE))


class ClassAPublishGateTest(unittest.TestCase):
    """Control 7 on synthetic class A gate files."""

    def test_allowed_merge_gate_passes(self) -> None:
        """The merge gate file with parameters, gates and a verdict passes."""
        self.assertEqual(gate("class_a_merge_gates.json", MERGE_GATES_PASS), [])

    def test_record_count_fails(self) -> None:
        """The record count of each arm is an undeclared quantity."""
        problems = gate("class_a_merge_gates.json", MERGE_GATES_WITH_COUNTS)
        self.assertTrue(problems)
        self.assertTrue(any("input_records" in problem for problem in problems))
        self.assertTrue(any("merged_records" in problem for problem in problems))

    def test_byte_size_of_a_class_b2_output_fails(self) -> None:
        """A byte size of a class B2 file in a class A gate file fails."""
        problems = gate("class_a_phase_gates.json", PHASE_GATES_WITH_SIZE)
        self.assertTrue(problems)
        self.assertEqual(gate("class_a_phase_gates.json", PHASE_GATES_PASS), [])

    def test_an_undeclared_file_name_fails_closed(self) -> None:
        """A file name with no entry in the rules file fails."""
        problems = gate("class_a_new_gates.json", {"verdict": "PASS"})
        self.assertTrue(any("no entry" in problem for problem in problems))

    def test_a_key_outside_the_fixed_list_fails(self) -> None:
        """A key that the rules file does not allow fails, even with no numeric value."""
        payload = dict(MERGE_GATES_PASS)
        payload["top_candidate_gene"] = "BUB1B"
        problems = gate("class_a_merge_gates.json", payload)
        self.assertTrue(any("top_candidate_gene" in problem for problem in problems))

    def test_a_declared_public_quantity_passes(self) -> None:
        """The lexicon label count of the public ontology passes, because the rules declare it."""
        payload = {"gate_id": "case_phenotype_extraction_v1", "data_class": "A",
                   "scope": "A verdict and the checksums of the inputs.",
                   "document_sha256": "0" * 64, "hp_obo_sha256": "1" * 64,
                   "lexicon_labels": 41595, "method": "A dictionary lookup.",
                   "verdict": "PASS", "condition": "At least one HPO term.",
                   "out_terms": "phenotype/case_hpo_terms.tsv", "out_terms_sha256": "2" * 64}
        self.assertEqual(gate("phenotype_extraction_gate.json", payload), [])

    def test_a_term_count_of_the_child_fails(self) -> None:
        """A count of the extracted terms of the child fails, because no rule declares it."""
        payload = {"gate_id": "case_phenotype_extraction_v1", "data_class": "A",
                   "scope": "A verdict and the checksums of the inputs.",
                   "document_sha256": "0" * 64, "hp_obo_sha256": "1" * 64,
                   "lexicon_labels": 41595, "method": "A dictionary lookup.",
                   "verdict": "PASS", "condition": "At least one HPO term.",
                   "out_terms": "phenotype/case_hpo_terms.tsv", "out_terms_sha256": "2" * 64,
                   "matched_terms": 11}
        self.assertTrue(gate("phenotype_extraction_gate.json", payload))

    def test_a_genomic_record_in_a_class_a_file_fails(self) -> None:
        """A genotype field in a class A gate file fails the record scan of the same gate."""
        payload = dict(MERGE_GATES_PASS)
        payload["scope"] = "The merge parameters. " + "0/1" + ":12,10:22"
        self.assertTrue(gate("class_a_merge_gates.json", payload))

    def test_a_declared_file_that_is_not_json_fails(self) -> None:
        """A declared file whose bytes are not JSON fails, because its key list cannot be read."""
        problems = privacy_guard.class_a_publish_problems(
            Path("class_a_merge_gates.json"), "verdict: PASS\n", RULES)
        self.assertTrue(any("not\nvalid JSON" in problem or "valid JSON" in problem
                            for problem in problems))


class ListingGateTest(unittest.TestCase):
    """Control 8 on synthetic listings."""

    def test_object_sizes_of_a_class_b2_prefix_fail(self) -> None:
        """The recursive listing prints a byte count for each class B2 object."""
        problems = listing(S3_LISTING_WITH_SIZES)
        self.assertTrue(problems)
        self.assertEqual(len(problems), 2)

    def test_the_safe_listing_passes(self) -> None:
        """A listing of keys alone passes."""
        self.assertEqual(listing(S3_LISTING_SAFE), [])

    def test_a_human_readable_size_fails(self) -> None:
        """A `--human-readable` listing of a class B2 prefix fails."""
        self.assertTrue(listing(S3_LISTING_HUMAN))

    def test_a_total_size_fails(self) -> None:
        """A `--summarize` footer over a class B2 prefix fails."""
        self.assertTrue(listing(S3_LISTING_SUMMARIZE))

    def test_a_monitor_line_with_a_size_fails(self) -> None:
        """A monitor line that carries a size beside a class B2 path fails."""
        self.assertTrue(listing(MONITOR_LINE_WITH_SIZE))

    def test_a_public_prefix_is_not_flagged(self) -> None:
        """A listing of a public reference prefix passes, because the rule covers class B2 only."""
        self.assertEqual(listing(PUBLIC_LISTING_WITH_SIZE), [])


class CommandLineTest(unittest.TestCase):
    """The gate exits non-zero and prints no value of a refused file."""

    def test_the_gate_exits_two_and_prints_no_value(self) -> None:
        """The command line gate fails with exit code 2 and prints the key and not the count."""
        with tempfile.TemporaryDirectory(prefix="synthetic_path_41_63_") as name:
            path = Path(name) / "class_a_merge_gates.json"
            path.write_text(json.dumps(MERGE_GATES_WITH_COUNTS, indent=2))
            result = subprocess.run(
                [sys.executable, str(ROOT / "harness" / "privacy_guard.py"),
                 "--class-a-publish", "--rules",
                 str(ROOT / "literature" / "class_a_publish_rules.json"), str(path)],
                capture_output=True, text=True, cwd=str(ROOT))
            self.assertEqual(result.returncode, 2)
            self.assertEqual(result.stderr, "")
            lines = result.stdout.splitlines()
            self.assertEqual(len(lines), 1, result.stdout)
            status, reported_path, diagnostic = lines[0].split("\t", 2)
            self.assertEqual(status, "FAIL")
            self.assertEqual(reported_path, str(path))
            self.assertIn("input_records", diagnostic)
            for value in ("41", "37", "63"):
                with self.subTest(value=value):
                    self.assertNotIn(value, diagnostic)

    def test_the_gate_passes_the_allowed_file(self) -> None:
        """The command line gate exits zero on the allowed merge gate file."""
        with tempfile.TemporaryDirectory() as name:
            path = Path(name) / "class_a_merge_gates.json"
            path.write_text(json.dumps(MERGE_GATES_PASS, indent=2))
            result = subprocess.run(
                [sys.executable, str(ROOT / "harness" / "privacy_guard.py"),
                 "--class-a-publish", "--rules",
                 str(ROOT / "literature" / "class_a_publish_rules.json"), str(path)],
                capture_output=True, text=True, cwd=str(ROOT))
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_the_listing_gate_reads_standard_input(self) -> None:
        """The listing gate reads a dash as standard input and refuses the synthetic size listing."""
        result = subprocess.run(
            [sys.executable, str(ROOT / "harness" / "privacy_guard.py"), "--listing", "--rules",
             str(ROOT / "literature" / "class_a_publish_rules.json"), "-"],
            input=S3_LISTING_WITH_SIZES, capture_output=True, text=True, cwd=str(ROOT))
        self.assertEqual(result.returncode, 2)


class CommittedWriterTest(unittest.TestCase):
    """Gate G1 of ADR 0008. Every declared file name has a committed writer."""

    def test_every_declared_file_names_its_writer(self) -> None:
        """Each entry of the rules file names the writer that produces the file."""
        for name, entry in RULES["files"].items():
            with self.subTest(name=name):
                self.assertTrue(entry.get("writer"), f"{name} declares no writer.")
                self.assertTrue(entry.get("allowed_key_paths"), f"{name} declares no key list.")

    def test_every_declared_quantity_carries_a_reason(self) -> None:
        """Each declared quantity states why the value is not a property of the child."""
        for name, entry in RULES["files"].items():
            for key, reason in entry.get("declared_quantities", {}).items():
                with self.subTest(name=name, key=key):
                    self.assertGreater(len(reason), 20, f"{name}:{key} has no reason.")


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
