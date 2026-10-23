"""Verify the Track 2 one command runner, without a network call.

`scripts/track2_run_disease.py` is the single entry point of the Track 2 method. It must fail closed
when a prerequisite is missing, and it must resolve a disease identifier without a hand written
field. The tests below check three groups of behaviour.

The first group checks that each preflight check fails when its resource is missing, changed,
unwritable or too small, and that each failure carries a repair instruction. The second group checks
the identifier classifier, the disease name matcher and the hit selector of the resolver on fixed
records, so no autocomplete call runs. The third group checks the cost model fit and the cost
projection on a fixture whose answer is known by construction.

Every fixture holds public class A content only: MeSH identifiers, gene symbols, OMIM identifiers
and counts of a public graph. No patient data enters this test and no network call runs.
"""

from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import track2_common as common  # noqa: E402
import track2_cost_model_fit as cost_model  # noqa: E402
import track2_disease_resolver as resolver  # noqa: E402
import track2_peak_memory as peak_memory  # noqa: E402
import track2_run_disease as runner  # noqa: E402

# Two rows of the Human Phenotype Ontology release 2026-09-01, copied field for field.
GENES_TO_DISEASE = (
    "ncbi_gene_id\tgene_symbol\tassociation_type\tdisease_id\tsource\n"
    "NCBIGene:157570\tESCO2\tMENDELIAN\tOMIM:268300\t"
    "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/mim2gene_medgen\n"
    "NCBIGene:4683\tNBN\tMENDELIAN\tOMIM:251260\t"
    "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/mim2gene_medgen\n"
    "NCBIGene:672\tBRCA1\tPOLYGENIC\tOMIM:114480\t"
    "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/mim2gene_medgen\n"
)
PHENOTYPE_HPOA = (
    "#version: 2026-09-02\n"
    "database_id\tdisease_name\tqualifier\thpo_id\n"
    "OMIM:268300\tRoberts-SC phocomelia syndrome\t\tHP:0000175\n"
    "OMIM:251260\tNijmegen breakage syndrome\t\tHP:0000252\n"
)

# The autocomplete record of `@DISEASE_Roberts_Syndrome`, copied field for field.
ROBERTS_HIT = {"_id": "@DISEASE_Roberts_Syndrome", "db": "ncbi_mesh", "db_id": "C535687",
               "name": "Roberts Syndrome", "biotype": "disease"}
LISSENCEPHALY_HIT = {"_id": "@DISEASE_Norman_Roberts_lissencephaly_syndrome", "db": "ncbi_mesh",
                     "db_id": "C537848", "name": "Norman Roberts lissencephaly syndrome",
                     "biotype": "disease"}


class FakeResponse:
    """A response that carries a status code alone, for the throttle tests."""

    def __init__(self, status_code: int) -> None:
        """Store the status code.

        Args:
            status_code: The HTTP status code.
        """
        self.status_code = status_code
        self.headers: dict[str, str] = {}


class FakeStatusError(Exception):
    """An error that carries a response, as `httpx.HTTPStatusError` does."""

    def __init__(self, status_code: int) -> None:
        """Store the response.

        Args:
            status_code: The HTTP status code of the response.
        """
        super().__init__(f"HTTP {status_code}")
        self.response = FakeResponse(status_code)


def write_release(folder: str) -> Path:
    """Write the two release files of the fixture.

    Args:
        folder: A temporary directory.

    Returns:
        The directory that holds the two files.
    """
    directory = Path(folder) / "hpo"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "genes_to_disease.txt").write_text(GENES_TO_DISEASE, encoding="utf-8")
    (directory / "phenotype.hpoa").write_text(PHENOTYPE_HPOA, encoding="utf-8")
    return directory


class PreflightFailsClosed(unittest.TestCase):
    """Every preflight check refuses a missing or changed resource."""

    def test_a_missing_release_file_fails_with_a_download_command(self) -> None:
        """A release file that is absent fails and names the download command."""
        entry = {"id": "hpo", "path": "tracking/cache/does_not_exist.txt", "release": "2026-09-01",
                 "sha256": "0" * 64, "download": "https://example.invalid/genes_to_disease.txt"}
        check = runner.check_release_file(entry)
        self.assertFalse(check["ok"])
        self.assertIn("is missing", check["detail"])
        self.assertIn("curl", check["fix"])

    def test_a_changed_release_file_fails_on_its_digest(self) -> None:
        """A release file whose SHA-256 differs from the pinned value fails."""
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "genes_to_disease.txt"
            path.write_text(GENES_TO_DISEASE, encoding="utf-8")
            check = runner.check_release_file(
                {"id": "hpo", "path": str(path), "release": "2026-09-01", "sha256": "0" * 64,
                 "download": "https://example.invalid/genes_to_disease.txt"})
        self.assertFalse(check["ok"])
        self.assertIn("has the SHA-256", check["detail"])

    def test_an_unwritable_output_directory_fails(self) -> None:
        """An output directory that the account cannot create fails."""
        check = runner.check_writable(Path("/proc/track2_output_probe"))
        self.assertFalse(check["ok"])
        self.assertIn("permission", check["fix"])

    def test_a_disk_below_the_floor_fails(self) -> None:
        """A free disk floor that no machine reaches fails."""
        with tempfile.TemporaryDirectory() as folder:
            check = runner.check_disk(Path(folder), 10 ** 9)
        self.assertFalse(check["ok"])

    def test_a_missing_cost_model_fails_with_the_fit_command(self) -> None:
        """A missing cost model fails and names the command that fits one."""
        check = runner.check_cost_model(Path("/nonexistent/track2_cost_model.yaml"))
        self.assertFalse(check["ok"])
        self.assertIn("track2_cost_model_fit.py", check["fix"])

    def test_a_missing_tool_fails(self) -> None:
        """A command line tool that is not installed fails."""
        check = runner.check_tool("nosuchtool", ["nosuchtool_that_does_not_exist", "--version"])
        self.assertFalse(check["ok"])

    def test_the_python_check_reads_the_running_interpreter(self) -> None:
        """The Python check passes on a floor below the running version and fails above it."""
        self.assertTrue(runner.check_python("3.9")["ok"])
        self.assertFalse(runner.check_python("99.0")["ok"])


class ResolverReadsTheSources(unittest.TestCase):
    """The resolver reads every field from a source and writes none by hand."""

    def test_the_classifier_reads_the_three_identifier_forms(self) -> None:
        """An OMIM identifier, a MeSH identifier and a name are told apart."""
        self.assertEqual(resolver.classify("OMIM:268300"), ("omim", "OMIM:268300"))
        self.assertEqual(resolver.classify("268300"), ("omim", "OMIM:268300"))
        self.assertEqual(resolver.classify("MESH:C535687"), ("mesh", "C535687"))
        self.assertEqual(resolver.classify("D003550"), ("mesh", "D003550"))
        self.assertEqual(resolver.classify("Roberts syndrome"), ("name", "Roberts syndrome"))

    def test_the_gene_reader_keeps_a_mendelian_row_alone(self) -> None:
        """A causal gene row is kept and a polygenic row is not."""
        with tempfile.TemporaryDirectory() as folder:
            directory = write_release(folder)
            rows = resolver.read_gene_rows(directory / "genes_to_disease.txt", "OMIM:268300")
            polygenic = resolver.read_gene_rows(directory / "genes_to_disease.txt", "OMIM:114480")
        self.assertEqual([row["symbol"] for row in rows], ["ESCO2"])
        self.assertEqual(rows[0]["ncbi_gene"], "NCBIGene:157570")
        self.assertEqual(polygenic, [])

    def test_a_missing_release_file_raises_with_the_repair(self) -> None:
        """A missing release file raises an error that names the download target."""
        with self.assertRaises(resolver.ResolutionError) as caught:
            resolver.read_gene_rows(Path("/nonexistent/genes_to_disease.txt"), "OMIM:268300")
        self.assertIn("make track2-resources", str(caught.exception))

    def test_the_name_matcher_reads_case_and_punctuation_as_noise(self) -> None:
        """The MeSH term "Roberts-Sc Phocomelia Syndrome" matches the HPO name."""
        with tempfile.TemporaryDirectory() as folder:
            directory = write_release(folder)
            found = resolver.match_omim_by_name(directory / "phenotype.hpoa", "Roberts Syndrome",
                                                "Roberts-Sc Phocomelia Syndrome")
        self.assertEqual(found, "OMIM:268300")

    def test_an_ambiguous_autocomplete_answer_refuses_to_guess(self) -> None:
        """Two hits with no exact name raise an error that asks for the MeSH identifier."""
        with self.assertRaises(resolver.ResolutionError) as caught:
            resolver.pick_disease_hit([ROBERTS_HIT, LISSENCEPHALY_HIT], "Roberts")
        self.assertIn("MESH:", str(caught.exception))

    def test_an_exact_name_selects_one_hit_of_two(self) -> None:
        """A query that equals one hit name selects that hit."""
        hit = resolver.pick_disease_hit([ROBERTS_HIT, LISSENCEPHALY_HIT], "Roberts Syndrome")
        self.assertEqual(hit["db_id"], "C535687")

    def test_an_empty_autocomplete_answer_names_the_cause(self) -> None:
        """A disease that PubTator3 does not recognise raises an error that says so."""
        with self.assertRaises(resolver.ResolutionError) as caught:
            resolver.pick_disease_hit([], "Warsaw breakage syndrome")
        self.assertIn("PubTator3 recognises no entity", str(caught.exception))

    def test_the_mesh_unique_identifier_follows_the_record_type(self) -> None:
        """A descriptor takes the prefix 68 and a supplementary record takes 67."""
        self.assertEqual(resolver.mesh_uid_of("D003550"), "68003550")
        self.assertEqual(resolver.mesh_uid_of("C535687"), "67535687")
        self.assertEqual(resolver.mesh_uid_of("X1"), "")


class CostModelFits(unittest.TestCase):
    """The cost model reproduces a relation that the fixture builds by construction."""

    def test_the_line_recovers_a_known_slope_and_intercept(self) -> None:
        """Points on a line return that line and a coefficient of determination of one."""
        xs = [1.0, 2.0, 3.0, 4.0]
        ys = [2.0 * x + 0.5 for x in xs]
        fit = cost_model.ols(xs, ys)
        self.assertAlmostEqual(fit["slope"], 2.0, places=3)
        self.assertAlmostEqual(fit["intercept"], 0.5, places=3)
        self.assertAlmostEqual(fit["r_squared"], 1.0, places=3)

    def test_one_point_is_refused(self) -> None:
        """A fit with fewer than two points raises."""
        with self.assertRaises(ValueError):
            cost_model.ols([1.0], [1.0])

    def test_the_origin_fit_recovers_a_known_rate(self) -> None:
        """Points through the origin return that rate."""
        fit = cost_model.origin_slope([10.0, 20.0, 30.0], [27.0, 54.0, 81.0])
        self.assertAlmostEqual(fit["slope"], 2.7, places=3)

    def test_the_prediction_inverts_the_power_law(self) -> None:
        """A graph of ten thousand edges returns the value that the coefficients state."""
        block = {"slope": 2.0, "intercept": -6.0}
        self.assertAlmostEqual(cost_model.predict_seconds_per_ranking(block, 10_000),
                               10.0 ** (-6.0 + 2.0 * math.log10(10_000)), places=6)

    def test_the_projection_multiplies_the_work_by_the_rate(self) -> None:
        """The projected seconds equal the rankings times the predicted seconds per ranking."""
        model = {"version": 1, "fitted_on": ["a"],
                 "ranking_seconds": {"slope": 2.0, "intercept": -6.0, "valid_edges_from": 1_000,
                                     "valid_edges_to": 100_000, "largest_absolute_ratio": 1.1}}
        manifest = {"edge_count": 10_000, "node_count": 2_000, "diseases_with_an_indication": 50}
        values = runner.projection(model, manifest, 5)
        self.assertEqual(values["rankings"], 300)
        self.assertAlmostEqual(values["predicted_wall_seconds"],
                               round(values["predicted_seconds_per_ranking"] * 300, 1), places=1)
        self.assertTrue(values["edges_inside_the_measured_range"])

    def test_a_graph_outside_the_fitted_range_is_marked(self) -> None:
        """A graph larger than every fitted run is marked as an extrapolation."""
        model = {"version": 1, "fitted_on": ["a"],
                 "ranking_seconds": {"slope": 2.0, "intercept": -6.0, "valid_edges_from": 1_000,
                                     "valid_edges_to": 5_000, "largest_absolute_ratio": 1.1}}
        manifest = {"edge_count": 10_000, "node_count": 2_000, "diseases_with_an_indication": 50}
        self.assertFalse(runner.projection(model, manifest, 5)["edges_inside_the_measured_range"])


class ThrottlePacing(unittest.TestCase):
    """A server that refuses a burst slows the process and keeps its attempt budget."""

    def setUp(self) -> None:
        """Start each test from a process with no extra pause."""
        common.note_success()
        while common.extra_pause() > 0.0:
            common.note_success()

    def test_a_throttle_status_is_told_from_a_dropped_connection(self) -> None:
        """HTTP 429 and HTTP 503 read as a throttle and HTTP 404 does not."""
        self.assertTrue(common.is_throttle(FakeStatusError(429)))
        self.assertTrue(common.is_throttle(FakeStatusError(503)))
        self.assertFalse(common.is_throttle(FakeStatusError(404)))
        self.assertFalse(common.is_throttle(ValueError("no response")))

    def test_the_pause_rises_with_each_refusal_and_stops_at_the_ceiling(self) -> None:
        """Repeated refusals raise the process pause to the recorded ceiling and no further."""
        for _ in range(100):
            common.note_throttle()
        self.assertEqual(common.extra_pause(), common.THROTTLE_PAUSE_MAX_S)

    def test_the_pause_falls_after_answered_calls(self) -> None:
        """A run of answered calls lowers the pause back to zero."""
        common.note_throttle()
        raised = common.extra_pause()
        self.assertGreater(raised, 0.0)
        for _ in range(500):
            common.note_success()
        self.assertEqual(common.extra_pause(), 0.0)

    def test_a_throttled_call_carries_the_larger_attempt_budget(self) -> None:
        """A call that only met a throttle status keeps trying past the ordinary budget."""
        self.assertTrue(common.attempt_budget(common.MAX_ATTEMPTS, common.MAX_ATTEMPTS))
        self.assertFalse(common.attempt_budget(common.MAX_THROTTLE_ATTEMPTS,
                                               common.MAX_THROTTLE_ATTEMPTS))
        self.assertFalse(common.attempt_budget(common.MAX_ATTEMPTS, 0))
        self.assertFalse(common.attempt_budget(common.MAX_ATTEMPTS, common.MAX_ATTEMPTS - 1))


class PeakMemorySampler(unittest.TestCase):
    """The memory sampler reads a process tree or reports that it cannot."""

    def test_the_sampler_reads_its_own_process(self) -> None:
        """The sampler returns a positive size for the process that runs the test."""
        if not peak_memory.available():
            self.skipTest("this system has no /proc interface.")
        self.assertGreater(peak_memory.tree_rss_kb(__import__("os").getpid()), 0)

    def test_an_exited_process_reads_as_zero(self) -> None:
        """A process identifier that no process holds reads as zero and raises nothing."""
        self.assertEqual(peak_memory.read_rss_kb(2 ** 22), 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
