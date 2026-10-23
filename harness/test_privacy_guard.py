"""Verify that the privacy guard refuses raw genomic records and arms above their ceiling.

The tests cover direct file reads, shell commands and wildcard expansion against protected paths.
Every record and private path in these fixtures is synthetic.
"""

from __future__ import annotations

import csv
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import class_b2_command_guard
import privacy_guard

VCF_TEXT = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr15\t40161100\t.\tG\tA\t50\tPASS\t.\tGT:AD:DP\t0/1:12,10:22\n"
VCF_ROW_ONLY = "Look at this row:\n15\t40161100\trs1\tG\tA\t50\tPASS\tDP=22\n"
GENOTYPE_TABLE = "gene\tgenotype\nBUB1B\t0/1:12,10:22\n"
FASTQ_TEXT = "@A00123:8:HGW:1:1101:1000:1000 1:N:0:ACGT\nACGTACGTACGTACGTACGTACGTACGTACGTACGT\n+\nFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF\n"
SAM_TEXT = "read1\t99\tchr15\t40161000\t60\t151M\t=\t40161200\t350\tACGT\tFFFF\n"
LITERATURE = ("Rio Frio et al. reported a homozygous BUB1B intronic variant c.2386-11A>G that creates a "
              "de novo splice site. Ochiai et al. found ss802470619 44 kb upstream of BUB1B. "
              "Compound heterozygous p.Arg814His and a truncating allele reduce BUBR1.")


class ScanTest(unittest.TestCase):
    """Check the pattern scan on synthetic records and on normal literature text."""

    def test_refuses_raw_records(self) -> None:
        """Each raw genomic format fails the scan."""
        for text in (VCF_TEXT, VCF_ROW_ONLY, GENOTYPE_TABLE, FASTQ_TEXT, SAM_TEXT):
            with self.subTest(text=text[:20]):
                self.assertTrue(privacy_guard.scan(text, tokens=[]))

    def test_refuses_json_escaped_records(self) -> None:
        """A VCF record inside a JSON string fails the scan."""
        import json
        self.assertTrue(privacy_guard.scan(json.dumps({"body": VCF_ROW_ONLY}), tokens=[]))

    def test_passes_literature(self) -> None:
        """Variant nomenclature in literature text passes the scan."""
        self.assertEqual(privacy_guard.scan(LITERATURE, tokens=[]), [])

    def test_refuses_private_identifier(self) -> None:
        """A deny token fails the scan without regard to case."""
        self.assertTrue(privacy_guard.scan("The file sample_zz9999 is large.", tokens=["SAMPLE_ZZ9999"]))


class CeilingTest(unittest.TestCase):
    """Check the class ceiling of each arm."""

    def test_subscription_arms_take_class_a_only(self) -> None:
        """The exported dispatcher assigns class A ceilings to subscription backends."""
        with mock.patch.object(privacy_guard, "deny_tokens", return_value=[]):
            privacy_guard.check_call("claude", "opus", "A", LITERATURE)
            for backend in ("claude", "codex"):
                with self.assertRaises(ValueError):
                    privacy_guard.check_call(backend, "model", "B", LITERATURE)

    def test_bedrock_takes_class_b_but_not_c(self) -> None:
        """OpenCode on Bedrock accepts class B and refuses class C."""
        with mock.patch.object(privacy_guard, "deny_tokens", return_value=[]):
            model = "amazon-bedrock/us.anthropic.claude-opus-5"
            self.assertEqual(privacy_guard.check_call("opencode", model, "B", LITERATURE)["ceiling"], "B")
            with self.assertRaises(ValueError):
                privacy_guard.check_call("opencode", model, "C", LITERATURE)
            with self.assertRaises(ValueError):
                privacy_guard.check_call("opencode", "openrouter/some-model", "B", LITERATURE)

    def test_scan_applies_to_class_a(self) -> None:
        """A class A call with a VCF record fails."""
        with mock.patch.object(privacy_guard, "deny_tokens", return_value=[]):
            with self.assertRaises(ValueError):
                privacy_guard.check_call("codex", "gpt", "A", VCF_TEXT)


class ModelCallLogTest(unittest.TestCase):
    """Check that sent and refused calls both reach the model call log."""

    def test_records_sent_and_refused(self) -> None:
        """One SENT row and one REFUSED row, with a prompt hash and no prompt text."""
        with tempfile.TemporaryDirectory() as name, mock.patch.object(privacy_guard, "deny_tokens", return_value=[]):
            root = Path(name)
            privacy_guard.guard_call(root, "claude", "opus", "A", LITERATURE, "t1", "run_a")
            with self.assertRaises(ValueError):
                privacy_guard.guard_call(root, "claude", "opus", "A", FASTQ_TEXT, "t1", "run_b")
            with (root / privacy_guard.CALL_LOG).open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["verdict"] for row in rows], ["SENT", "REFUSED"])
            self.assertNotIn("ACGTACGT", (root / privacy_guard.CALL_LOG).read_text())


def build_class_b2_tree(root: Path) -> None:
    """Write a synthetic report tree that holds one class B2 include.

    Args:
        root: The directory that receives the tree.
    """
    report = root / "writeups" / "track1_methods"
    report.mkdir(parents=True)
    (report / "track1_methods.qmd").write_text("The report text of the fixture.\n")
    (report / "_class_b2_results.qmd").write_text("A synthetic line with no value of any person.\n")
    (report / "_results_public_stub.qmd").write_text("The public stub of the fixture.\n")
    prefix = root / "derived" / "llm_safe" / "b0_fixture" / "out"
    prefix.mkdir(parents=True)
    (prefix / "rows.csv").write_text("column\nvalue\n")
    other = root / "writeups" / "track2_repurposing"
    other.mkdir(parents=True)
    (other / "track2_repurposing.qmd").write_text("The second report text of the fixture.\n")


class ClassB2CommandTest(unittest.TestCase):
    """Check that the command guard refuses the reader forms of control 9."""

    def setUp(self) -> None:
        """Build the synthetic tree that each fixture reads."""
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        build_class_b2_tree(self.root)

    def tearDown(self) -> None:
        """Remove the synthetic tree."""
        self.directory.cleanup()

    def test_refuses_restricted_include_grep(self) -> None:
        """The conservative command guard refuses grep over a restricted include."""
        command = "grep -c -- '-' writeups/track1_methods/_class_b2_results.qmd"
        verdict, reason = class_b2_command_guard.decide(command, self.root)
        self.assertEqual(verdict, class_b2_command_guard.REFUSE)
        self.assertIn("_class_b2_results.qmd", reason)

    def test_allows_the_same_command_on_the_public_stub(self) -> None:
        """The same reader over the public stub passes."""
        command = "grep -c -- '-' writeups/track1_methods/_results_public_stub.qmd"
        verdict, _reason = class_b2_command_guard.decide(command, self.root)
        self.assertEqual(verdict, class_b2_command_guard.ALLOW)

    def test_allows_the_one_approved_reader(self) -> None:
        """The privacy guard reads a class B2 path, because it prints a verdict and no content."""
        command = "uv run python harness/privacy_guard.py writeups/track1_methods/_class_b2_results.qmd"
        verdict, _reason = class_b2_command_guard.decide(command, self.root)
        self.assertEqual(verdict, class_b2_command_guard.ALLOW)

    def test_refuses_a_recursive_reader_on_the_parent_directory(self) -> None:
        """A recursive grep reaches the include without naming it, so the directory is refused."""
        command = "grep -rn 'text' writeups/track1_methods"
        verdict, reason = class_b2_command_guard.decide(command, self.root)
        self.assertEqual(verdict, class_b2_command_guard.REFUSE)
        self.assertIn("_class_b2_results.qmd", reason)

    def test_refuses_a_reader_of_the_class_b2_prefix(self) -> None:
        """A command that measures a class B2 output is refused."""
        for command in ("wc -c derived/llm_safe/b0_fixture/out/rows.csv",
                        "aws s3 ls --human-readable s3://bucket/derived/llm_safe/b0_fixture/out/"):
            with self.subTest(command=command):
                verdict, _reason = class_b2_command_guard.decide(command, self.root)
                self.assertEqual(verdict, class_b2_command_guard.REFUSE)

    def test_allows_a_pattern_that_names_the_class(self) -> None:
        """A search for the string of the class over a class A file passes."""
        command = "grep -n '_class_b2_' writeups/track1_methods/track1_methods.qmd"
        verdict, _reason = class_b2_command_guard.decide(command, self.root)
        self.assertEqual(verdict, class_b2_command_guard.ALLOW)

    def test_refuses_a_class_b2_path_in_a_second_segment(self) -> None:
        """A refused reader after a pipe or a separator is still refused."""
        command = "echo start && cat writeups/track1_methods/_class_b2_results.qmd"
        verdict, _reason = class_b2_command_guard.decide(command, self.root)
        self.assertEqual(verdict, class_b2_command_guard.REFUSE)

    def test_refuses_glob_matching_restricted_include(self) -> None:
        """The command guard refuses a word count when its glob reaches a restricted include."""
        command = "wc -lw writeups/track1_methods/*.qmd writeups/track2_repurposing/*.qmd"
        verdict, reason = class_b2_command_guard.decide(command, self.root)
        self.assertEqual(verdict, class_b2_command_guard.REFUSE)
        self.assertIn("_class_b2_results.qmd", reason)

    def test_allows_a_glob_that_reaches_no_class_b2_file(self) -> None:
        """A glob whose expansion holds no class B2 path passes."""
        command = "wc -lw writeups/track1_methods/track1_*.qmd"
        verdict, _reason = class_b2_command_guard.decide(command, self.root)
        self.assertEqual(verdict, class_b2_command_guard.ALLOW)

    def test_refuses_a_recursive_reader_on_a_glob_of_directories(self) -> None:
        """A recursive reader over a glob of directories reaches the include, so it is refused."""
        command = "grep -rn 'text' writeups/*"
        verdict, reason = class_b2_command_guard.decide(command, self.root)
        self.assertEqual(verdict, class_b2_command_guard.REFUSE)
        self.assertIn("_class_b2_results.qmd", reason)

    def test_the_hook_payload_returns_the_refusal_code(self) -> None:
        """The hook exits 2 when the synthetic command targets a restricted include."""
        payload = ('{"tool_input": {"command": '
                   '"grep -c x writeups/track1_methods/_class_b2_results.qmd"}}')
        with mock.patch.object(sys, "stdin", io.StringIO(payload)):
            code = class_b2_command_guard.main(["--hook", "--root", str(self.root)])
        self.assertEqual(code, 2)


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
