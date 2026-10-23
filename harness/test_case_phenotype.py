"""Verify the clinical document to HPO term mapping on a synthetic public document."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "pipelines/track1_b0/bin"
sys.path.insert(0, str(BIN))
import extract_case_phenotype as extract  # noqa: E402

OBO = """format-version: 1.2

[Term]
id: HP:0000001
name: All

[Term]
id: HP:0000118
name: Phenotypic abnormality
is_a: HP:0000001 ! All

[Term]
id: HP:0001249
name: Intellectual disability
is_a: HP:0000118 ! Phenotypic abnormality
synonym: "Mental retardation" EXACT []

[Term]
id: HP:0000252
name: Microcephaly
is_a: HP:0000118 ! Phenotypic abnormality
synonym: "Small head circumference" EXACT []

[Term]
id: HP:0001250
name: Seizure
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:0002373
name: Febrile seizure
is_a: HP:0001250 ! Seizure

[Term]
id: HP:0004322
name: Short stature
is_a: HP:0000118 ! Phenotypic abnormality
"""

# The text below is invented for this test. It describes no real person.
DOCUMENT = (
    "Clinical summary of a fictional test case.\n"
    "The child shows microcephaly and short stature.\n"
    "A febrile seizure was recorded at the age of two years.\n"
    "No cardiac finding was reported.\n"
)


def write_docx(path: Path, paragraphs: list[str]) -> None:
    """Write a minimal Word document with the given paragraphs.

    Args:
        path: The output `.docx` file.
        paragraphs: The paragraph text.
    """
    namespace = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
    body = "".join(f"<w:p><w:r><w:t>{line}</w:t></w:r></w:p>" for line in paragraphs)
    document = f'<?xml version="1.0"?><w:document xmlns:w="{namespace}"><w:body>{body}</w:body></w:document>'
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("word/document.xml", document)


class ExtractionTest(unittest.TestCase):
    """Check the lexicon, the matching rule and the gate file."""

    def setUp(self) -> None:
        """Write the synthetic ontology and the synthetic document."""
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.obo = self.root / "hp.obo"
        self.obo.write_text(OBO)
        self.text_document = self.root / "case.txt"
        self.text_document.write_text(DOCUMENT)
        self.docx = self.root / "case.docx"
        write_docx(self.docx, DOCUMENT.splitlines())

    def run_extraction(self, document: Path) -> tuple[list[str], dict]:
        """Run the extraction program on one document.

        Args:
            document: The document file.

        Returns:
            The matched term identifiers and the gate record.
        """
        terms = self.root / f"{document.stem}_terms.tsv"
        gate = self.root / f"{document.stem}_gate.json"
        subprocess.run([sys.executable, str(BIN / "extract_case_phenotype.py"),
                        "--document", str(document), "--obo", str(self.obo),
                        "--out-terms", str(terms), "--out-gate", str(gate)],
                       check=True, capture_output=True)
        rows = terms.read_text().splitlines()[1:]
        return [row.split("\t")[0] for row in rows], json.loads(gate.read_text())

    def test_the_text_document_gives_the_specific_terms(self) -> None:
        """The document maps to the three findings that it names."""
        found, _ = self.run_extraction(self.text_document)
        self.assertEqual(found, ["HP:0000252", "HP:0002373", "HP:0004322"])

    def test_the_general_term_is_dropped_when_a_child_matches(self) -> None:
        """Seizure is not reported beside febrile seizure."""
        found, _ = self.run_extraction(self.text_document)
        self.assertNotIn("HP:0001250", found)

    def test_a_term_that_the_document_does_not_name_is_absent(self) -> None:
        """Intellectual disability does not appear in the document, so it does not appear here."""
        found, _ = self.run_extraction(self.text_document)
        self.assertNotIn("HP:0001249", found)

    def test_the_word_document_gives_the_same_terms(self) -> None:
        """The Word reader and the text reader agree."""
        from_text, _ = self.run_extraction(self.text_document)
        from_docx, _ = self.run_extraction(self.docx)
        self.assertEqual(from_text, from_docx)

    def test_the_gate_file_carries_no_count_of_the_terms(self) -> None:
        """The class A gate file holds a verdict and checksums and no count of the matched terms."""
        found, gate = self.run_extraction(self.text_document)
        self.assertEqual(gate["verdict"], "PASS")
        self.assertEqual(gate["data_class"], "A")
        text = json.dumps(gate)
        for term in found:
            self.assertNotIn(term, text)
        self.assertNotIn("terms_matched", gate)
        self.assertNotIn("matched", " ".join(key for key in gate))

    def test_normalization_drops_the_punctuation(self) -> None:
        """The normalizer lowers the case and keeps word boundaries."""
        self.assertEqual(extract.normalize("Short-stature, severe!"), " short stature severe ")


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
