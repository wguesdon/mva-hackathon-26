"""Verify scripts/track1_sv_causal_benchmark.py on synthetic fixtures. Data class A.

The test builds a synthetic GENCODE extract, a synthetic gnomAD SV extract, a synthetic background
call set and a synthetic pathogenic record in a temporary directory. Every coordinate is invented.
The test needs no network and it opens no patient file.

The three scenarios that the test runs through the real scorer are:

    rank_one   The spike removes an exon and every background call in the gene is intronic.
    rank_two   The spike is intronic and one background call removes an exon.
    budget     Seven background calls tie with the spike and they sit at a lower position, so the
               spike holds row eight of a seven row budget.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import struct
import sys
import tempfile
import unittest
import zlib
from contextlib import redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "pipelines" / "sv_mei_annotate" / "bin"))

import annotate_sv_mei as scorer  # noqa: E402
import track1_sv_causal_benchmark as causal  # noqa: E402

# A synthetic gene TESTG1 on chr1 with three MANE exons, and a second gene TESTG2 that holds TESTG1
# to test the containment rule. A third gene TESTG3 sits on chr2.
GTF = "\n".join([
    "##description: synthetic test extract",
    'chr1\ttest\tgene\t100000\t200000\t.\t+\t.\tgene_id "T1"; gene_type "protein_coding"; gene_name "TESTG1";',
    'chr1\ttest\texon\t100000\t100500\t.\t+\t.\tgene_id "T1"; gene_name "TESTG1"; tag "MANE_Select";',
    'chr1\ttest\texon\t150000\t150400\t.\t+\t.\tgene_id "T1"; gene_name "TESTG1"; tag "MANE_Select";',
    'chr1\ttest\texon\t199000\t200000\t.\t+\t.\tgene_id "T1"; gene_name "TESTG1"; tag "MANE_Select";',
    'chr1\ttest\tgene\t90000\t260000\t.\t+\t.\tgene_id "T2"; gene_type "protein_coding"; gene_name "TESTG2";',
    'chr1\ttest\texon\t90000\t90500\t.\t+\t.\tgene_id "T2"; gene_name "TESTG2"; tag "MANE_Select";',
    'chr2\ttest\tgene\t10000\t20000\t.\t+\t.\tgene_id "T3"; gene_type "protein_coding"; gene_name "TESTG3";',
    'chr2\ttest\tgene\t30000\t40000\t.\t+\t.\tgene_id "T4"; gene_type "lncRNA"; gene_name "TESTL1";',
    # The seven MVA genes with synthetic coordinates. They let the comparison run with the shipped
    # default gene set. None of them shares a contig with TESTG1, TESTG2 or TESTG3.
    'chr15\ttest\tgene\t40160000\t40220000\t.\t+\t.\tgene_id "M1"; gene_type "protein_coding"; gene_name "BUB1B";',
    'chr15\ttest\texon\t40185000\t40185400\t.\t+\t.\tgene_id "M1"; gene_name "BUB1B"; tag "MANE_Select";',
    'chr16\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id "M2"; gene_type "protein_coding"; gene_name "CEP57";',
    'chr17\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id "M3"; gene_type "protein_coding"; gene_name "TRIP13";',
    'chr18\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id "M4"; gene_type "protein_coding"; gene_name "BUB1";',
    'chr19\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id "M5"; gene_type "protein_coding"; gene_name "BUB3";',
    'chr20\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id "M6"; gene_type "protein_coding"; gene_name "CDC20";',
    'chr21\ttest\tgene\t1000\t2000\t.\t+\t.\tgene_id "M7"; gene_type "protein_coding"; gene_name "PTPA";',
]) + "\n"

GNOMAD = "\n".join([
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t160000\tcommon_del\tN\t<DEL>\t.\tPASS\tEND=161000;SVTYPE=DEL;SVLEN=1000;AF=0.22;AF_nfe=0.25",
]) + "\n"

# The two known intervals. The scorer gate needs exactly two.
INTERVALS = ("chr9\t1000\t2000\tsynthetic_interval_one\n"
             "chr9\t5000\t6000\tsynthetic_interval_two\n")

BACKGROUND_INTRONIC = [
    "chr1\t120000\tbg_intronic_1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-300;END=120300\tGT:GQ\t0/1:60",
    "chr1\t130000\tbg_intronic_2\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-400;END=130400\tGT:GQ\t0/1:60",
    "chr2\t12000\tbg_other_gene\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-500;END=12500\tGT:GQ\t0/1:60",
]

# One call inside a MANE exon of the synthetic BUB1B, so the comparison of the default option
# against no option writes a candidate row and a submission row.
BACKGROUND_DEFAULT_ARM = BACKGROUND_INTRONIC + [
    "chr15\t40185100\tbg_bub1b_exon\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-200;END=40185300\tGT:GQ\t0/1:60",
]

BACKGROUND_EXONIC = BACKGROUND_INTRONIC + [
    "chr1\t150100\tbg_exonic\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-200;END=150300\tGT:GQ\t0/1:60",
]

SPIKE_EXONIC = {"spike_id": "nssv_test_exonic", "chrom": "chr1", "pos": 199100, "end": 199400,
                "svlen": 300, "svtype": "DEL", "gene": "TESTG1", "gene_count": 1,
                "clnsig": "pathogenic", "clinvar_accession": "VCV000000001.1", "imprecise": False}

# The early intronic spike sits before every intronic background call, so only the exonic background
# call can outrank it. The late intronic spike sits after the seven tied background calls of the
# budget scenario, so the tie breaker of the scorer puts it in row eight.
SPIKE_INTRONIC = {"spike_id": "nssv_test_intronic", "chrom": "chr1", "pos": 105000, "end": 105500,
                  "svlen": 500, "svtype": "DEL", "gene": "TESTG1", "gene_count": 1,
                  "clnsig": "pathogenic", "clinvar_accession": "VCV000000002.1", "imprecise": True}

SPIKE_INTRONIC_LATE = dict(SPIKE_INTRONIC, spike_id="nssv_test_intronic_late", pos=170000,
                           end=170500)


def bgzf_block(payload: bytes) -> bytes:
    """Compress one payload into a single BGZF block.

    Args:
        payload: The bytes to compress.

    Returns:
        The BGZF block, with the gzip header, the BC extra field, the deflate body and the trailer.
    """
    engine = zlib.compressobj(6, zlib.DEFLATED, -15)
    body = engine.compress(payload) + engine.flush()
    total = len(body) + 26
    header = (b"\x1f\x8b\x08\x04" + b"\x00\x00\x00\x00" + b"\x00" + b"\xff"
              + struct.pack("<H", 6) + b"BC" + struct.pack("<H", 2) + struct.pack("<H", total - 1))
    trailer = struct.pack("<I", zlib.crc32(payload)) + struct.pack("<I", len(payload) & 0xFFFFFFFF)
    return header + body + trailer


def budget_background() -> list[str]:
    """Build seven intronic background calls that tie with the intronic spike.

    Every call sits at a position below the spike, so the tie breaker of the scorer puts each one
    above the spike and the spike holds row eight of the seven row budget.

    Returns:
        The background data lines.
    """
    lines = []
    for index in range(7):
        start = 110000 + index * 1000
        lines.append(f"chr1\t{start}\tbg_tie_{index}\tN\t<DEL>\t.\tPASS\t"
                     f"SVTYPE=DEL;SVLEN=-500;END={start + 500}\tGT:GQ\t0/1:60")
    return lines


def build_inputs(folder: Path) -> dict[str, Path]:
    """Write every synthetic input file.

    Args:
        folder: A temporary directory.

    Returns:
        The path of each input.
    """
    paths = {"gtf": folder / "genes.gtf.gz", "gnomad": folder / "gnomad.vcf.gz",
             "intervals": folder / "intervals.bed", "mei": folder / "empty_mei.vcf"}
    with gzip.open(paths["gtf"], "wt") as handle:
        handle.write(GTF)
    with gzip.open(paths["gnomad"], "wt") as handle:
        handle.write(GNOMAD)
    paths["intervals"].write_text(INTERVALS)
    paths["mei"].write_text("\n".join(causal.COMBINED_HEADER) + "\n")
    return paths


def run_scenario(folder: Path, paths: dict[str, Path], name: str, background: list[str],
                 spike: dict) -> dict:
    """Place one spike, run the real scorer and read the result row.

    Args:
        folder: The working directory.
        paths: The synthetic input paths.
        name: The scenario name.
        background: The background data lines.
        spike: The pathogenic record to add.

    Returns:
        The result row of `measure_spike`.
    """
    sv_vcf = folder / f"{name}.vcf"
    sv_vcf.write_text("\n".join(causal.combined_lines(background, spike)) + "\n")
    out_dir = folder / f"out_{name}"
    outputs = causal.run_scorer(sv_vcf, paths["mei"], paths["gtf"], paths["gnomad"],
                                paths["intervals"], out_dir, spike["gene"], f"test_{name}")
    return causal.measure_spike(spike, outputs)


def output_checksums(folder: Path, paths: dict[str, Path], name: str,
                     extra: list[str]) -> dict[str, str]:
    """Run the scorer on the shared background file and return the checksum of every output.

    Both runs read one input file, because the class A gate file records the input name and the
    input checksum. Two input files would differ in that record for a reason that is not the option.

    Args:
        folder: The working directory.
        paths: The synthetic input paths, including the shared `background` VCF.
        name: The run name.
        extra: The extra command line arguments.

    Returns:
        The SHA-256 of each output file, keyed by file name.

    Raises:
        AssertionError: When the scorer returns a FAIL verdict.
    """
    sv_vcf = paths["background"]
    out_dir = folder / f"cmp_{name}"
    out_dir.mkdir(parents=True, exist_ok=True)
    argv = ["annotate_sv_mei.py", "--sv-vcf", str(sv_vcf), "--mei-vcf", str(paths["mei"]),
            "--gtf", str(paths["gtf"]), "--gnomad-sv", str(paths["gnomad"]),
            "--intervals", str(paths["intervals"]), "--out-dir", str(out_dir),
            "--run-id", "cmp", "--min-frequency-records", "1"] + extra
    saved = sys.argv
    sys.argv = argv
    try:
        with redirect_stdout(io.StringIO()):
            code = scorer.main()
    finally:
        sys.argv = saved
    if code != 0:
        raise AssertionError(f"The scorer returned a FAIL verdict for {name}")
    return {item.name: hashlib.sha256(item.read_bytes()).hexdigest()
            for item in sorted(out_dir.iterdir())}


class DefaultBehaviourTest(unittest.TestCase):
    """Check that the new option changes nothing when a run omits it."""

    @classmethod
    def setUpClass(cls) -> None:
        """Build the synthetic inputs once."""
        cls.folder = Path(tempfile.mkdtemp())
        cls.paths = build_inputs(cls.folder)
        cls.paths["background"] = cls.folder / "background.vcf"
        cls.paths["background"].write_text(
            "\n".join(causal.combined_lines(BACKGROUND_DEFAULT_ARM, None)) + "\n")

    def test_the_default_option_reproduces_the_run_that_omits_it(self) -> None:
        """A run with no `--genes` writes the same six files as a run with the default value."""
        without = output_checksums(self.folder, self.paths, "without", [])
        with_default = output_checksums(self.folder, self.paths, "with_default",
                                        ["--genes", ",".join(scorer.MVA_GENES)])
        self.assertEqual(len(without), 6)
        self.assertEqual(without, with_default)

    def test_the_comparison_run_is_not_empty(self) -> None:
        """The comparison writes at least one candidate row, so it compares real content."""
        checksums = output_checksums(self.folder, self.paths, "content", [])
        table = (self.folder / "cmp_content" / "class_b2_sv_mei_ranked.tsv").read_text()
        self.assertIn("bg_bub1b_exon", table)
        self.assertIn("class_b2_submission_rows.csv", checksums)


class SpikePlacementTest(unittest.TestCase):
    """Check that a spike enters the call set once, in coordinate order and with the right fields."""

    def test_spike_line_of_a_deletion(self) -> None:
        """A deletion carries a negative SVLEN, its end, a PASS filter and the chosen genotype."""
        line = causal.spike_line(SPIKE_EXONIC)
        fields = line.split("\t")
        self.assertEqual(fields[0], "chr1")
        self.assertEqual(fields[1], "199100")
        self.assertEqual(fields[2], "nssv_test_exonic")
        self.assertEqual(fields[4], "<DEL>")
        self.assertEqual(fields[6], "PASS")
        self.assertIn("SVLEN=-300", fields[7])
        self.assertIn("END=199400", fields[7])
        self.assertEqual(fields[9], f"{causal.SPIKE_GENOTYPE}:{causal.SPIKE_GENOTYPE_QUALITY}")

    def test_spike_line_of_a_duplication_keeps_a_positive_length(self) -> None:
        """A duplication carries a positive SVLEN."""
        spike = dict(SPIKE_EXONIC, svtype="DUP")
        self.assertIn("SVLEN=300", causal.spike_line(spike).split("\t")[7])

    def test_spike_enters_in_coordinate_order(self) -> None:
        """The combined call set holds every background line once and the spike in position order."""
        lines = causal.combined_lines(BACKGROUND_INTRONIC, SPIKE_INTRONIC_LATE)
        data = [line for line in lines if not line.startswith("#")]
        self.assertEqual(len(data), len(BACKGROUND_INTRONIC) + 1)
        self.assertEqual([line.split("\t")[2] for line in data],
                         ["bg_intronic_1", "bg_intronic_2", "nssv_test_intronic_late",
                          "bg_other_gene"])

    def test_spike_after_every_background_call_is_appended(self) -> None:
        """A spike that sits after every background call is the last data line."""
        spike = dict(SPIKE_INTRONIC, chrom="chrY", pos=5_000_000)
        data = [line for line in causal.combined_lines(BACKGROUND_INTRONIC, spike)
                if not line.startswith("#")]
        self.assertEqual(data[-1].split("\t")[2], "nssv_test_intronic")

    def test_no_spike_returns_the_background_alone(self) -> None:
        """The negative control call set holds the background lines and no added line."""
        data = [line for line in causal.combined_lines(BACKGROUND_INTRONIC, None)
                if not line.startswith("#")]
        self.assertEqual(data, BACKGROUND_INTRONIC)


class RankAndBudgetTest(unittest.TestCase):
    """Check the rank computation and the budget check on the real scorer."""

    @classmethod
    def setUpClass(cls) -> None:
        """Build the synthetic inputs once and run the three scenarios."""
        cls.folder = Path(tempfile.mkdtemp())
        cls.paths = build_inputs(cls.folder)
        cls.rank_one = run_scenario(cls.folder, cls.paths, "rank_one", BACKGROUND_INTRONIC,
                                    SPIKE_EXONIC)
        cls.rank_two = run_scenario(cls.folder, cls.paths, "rank_two", BACKGROUND_EXONIC,
                                    SPIKE_INTRONIC)
        cls.budget = run_scenario(cls.folder, cls.paths, "budget", budget_background(),
                                  SPIKE_INTRONIC_LATE)

    def test_exonic_spike_reaches_rank_one(self) -> None:
        """An exonic spike outranks every intronic background call and holds a submission row."""
        self.assertEqual(self.rank_one["rank"], 1)
        self.assertEqual(self.rank_one["gene_effect"], "exonic")
        self.assertEqual(self.rank_one["score"], 105)
        self.assertEqual(self.rank_one["clears_candidate_floor"], 1)
        self.assertEqual(self.rank_one["clears_submission_floor"], 1)
        self.assertEqual(self.rank_one["in_row_budget"], 1)

    def test_intronic_spike_falls_behind_an_exonic_background_call(self) -> None:
        """An intronic spike takes rank 2 when a background call removes an exon."""
        self.assertEqual(self.rank_two["rank"], 2)
        self.assertEqual(self.rank_two["score"], 92)
        self.assertEqual(self.rank_two["top_background_score"], 105)
        self.assertEqual(self.rank_two["in_row_budget"], 1)

    def test_row_eight_falls_outside_the_seven_row_budget(self) -> None:
        """A spike at rank 8 clears the submission floor and holds no row of the budget."""
        self.assertEqual(self.budget["rank"], scorer.MAX_SUBMISSION_ROWS + 1)
        self.assertEqual(self.budget["clears_submission_floor"], 1)
        self.assertEqual(self.budget["in_row_budget"], 0)
        self.assertEqual(self.budget["submission_rows"], scorer.MAX_SUBMISSION_ROWS)

    def test_calls_outside_the_disease_gene_stay_out_of_the_table(self) -> None:
        """A background call in another gene enters no candidate row of this run."""
        self.assertEqual(self.rank_one["candidate_rows"], len(BACKGROUND_INTRONIC))
        self.assertEqual(self.rank_one["background_calls_in_gene"], 2)

    def test_rank_of_reads_the_table(self) -> None:
        """The rank lookup returns the rank of a present call and None for an absent call."""
        rows = [{"call_id": "a", "rank": "1"}, {"call_id": "b", "rank": "2"}]
        self.assertEqual(causal.rank_of(rows, "b"), 2)
        self.assertIsNone(causal.rank_of(rows, "c"))

    def test_budget_check_matches_on_the_coordinate(self) -> None:
        """The budget check finds a row by its contig and its position and nothing else."""
        rows = [{"chrom_1": "chr1", "pos_1": "199100"}]
        self.assertTrue(causal.in_row_budget(rows, "chr1", 199100))
        self.assertFalse(causal.in_row_budget(rows, "chr1", 199101))
        self.assertFalse(causal.in_row_budget(rows, "chr2", 199100))


class DropReasonTest(unittest.TestCase):
    """Check that the filter attribution names the filter that the scorer applied."""

    def base_row(self, **changes) -> dict:
        """Build one annotated row.

        Args:
            **changes: The fields to replace.

        Returns:
            The row.
        """
        row = {"genotype": "het", "svtype": "DEL", "svlen": "300", "source": "sv",
               "genotype_quality": "60", "retroseq_level": "", "mva_gene": "TESTG1",
               "gnomad_popmax_af": "", "score": "105"}
        row.update(changes)
        return row

    def test_no_filter_removed_the_call(self) -> None:
        """A call that passes every filter carries no reason."""
        self.assertEqual(causal.drop_reason(self.base_row()), "")

    def test_caller_quality_filter(self) -> None:
        """A homozygous reference genotype names the caller quality filter."""
        self.assertEqual(causal.drop_reason(self.base_row(genotype="hom_ref")),
                         "caller_quality_filter")

    def test_no_disease_gene_overlap(self) -> None:
        """A call outside the disease gene names the gene filter."""
        self.assertEqual(causal.drop_reason(self.base_row(mva_gene="")),
                         "no_disease_gene_overlap")

    def test_frequency_filter(self) -> None:
        """A call above the population frequency ceiling names the frequency filter."""
        self.assertEqual(causal.drop_reason(self.base_row(gnomad_popmax_af="0.05")),
                         "gnomad_frequency_above_0.01")

    def test_score_floor(self) -> None:
        """A call below the candidate score floor names the score floor."""
        self.assertEqual(causal.drop_reason(self.base_row(score="20")),
                         "below_candidate_score_floor_25")


class SelectionTest(unittest.TestCase):
    """Check the dbVar selection rules and the draw."""

    def test_clinical_terms_unescape_and_split(self) -> None:
        """A quoted and URL escaped CLNSIG value becomes lower case terms."""
        self.assertEqual(causal.clinical_terms('"Likely%20pathogenic","Pathogenic"'),
                         ["likely pathogenic", "pathogenic"])

    def test_record_span_of_a_deletion_and_an_insertion(self) -> None:
        """A deletion takes its end from END. An insertion keeps its position as its end."""
        self.assertEqual(causal.record_span(1000, {"END": "1500", "SVLEN": "-500"}, "DEL"),
                         (1500, 500))
        self.assertEqual(causal.record_span(1000, {"END": "1001", "SVLEN": "300"}, "INS"),
                         (1000, 300))
        self.assertEqual(causal.record_span(1000, {}, "DEL"), (None, None))

    def test_containment_rule(self) -> None:
        """A record enters the benchmark only when one gene span holds it whole."""
        genes = causal.read_protein_coding_genes(self.gtf)
        index = causal.gene_span_index(genes)
        self.assertEqual(causal.containing_genes(index, "chr1", 150000, 150200),
                         ["TESTG1", "TESTG2"])
        self.assertEqual(causal.containing_genes(index, "chr1", 95000, 250000), ["TESTG2"])
        self.assertEqual(causal.containing_genes(index, "chr1", 95000, 300000), [])
        self.assertEqual(causal.containing_genes(index, "chr3", 100, 200), [])

    def test_a_non_coding_gene_is_not_a_container(self) -> None:
        """A long non-coding RNA gene holds no spike."""
        genes = causal.read_protein_coding_genes(self.gtf)
        self.assertEqual(causal.containing_genes(causal.gene_span_index(genes), "chr2", 31000,
                                                 31100), [])

    def test_draw_keeps_one_locus_for_each_gene(self) -> None:
        """The draw keeps one record for each locus and one locus for each gene."""
        records = [
            {"spike_id": "nssv3", "chrom": "chr1", "pos": 100, "end": 200, "svtype": "DEL",
             "gene": "A"},
            {"spike_id": "nssv1", "chrom": "chr1", "pos": 100, "end": 200, "svtype": "DEL",
             "gene": "A"},
            {"spike_id": "nssv2", "chrom": "chr1", "pos": 500, "end": 900, "svtype": "DEL",
             "gene": "A"},
            {"spike_id": "nssv4", "chrom": "chr2", "pos": 100, "end": 200, "svtype": "DUP",
             "gene": "B"},
        ]
        drawn = causal.draw_spikes(records, count=10, seed=1)
        self.assertEqual(sorted(item["gene"] for item in drawn), ["A", "B"])
        self.assertEqual(drawn, causal.draw_spikes(records, count=10, seed=1))

    def test_draw_respects_the_count(self) -> None:
        """The draw returns at most the requested number of spikes."""
        records = [{"spike_id": f"nssv{index}", "chrom": "chr1", "pos": index * 1000,
                    "end": index * 1000 + 300, "svtype": "DEL", "gene": f"G{index}"}
                   for index in range(1, 21)]
        self.assertEqual(len(causal.draw_spikes(records, count=5, seed=7)), 5)

    @classmethod
    def setUpClass(cls) -> None:
        """Write the synthetic annotation once."""
        cls.folder = Path(tempfile.mkdtemp())
        cls.gtf = cls.folder / "genes.gtf.gz"
        with gzip.open(cls.gtf, "wt") as handle:
            handle.write(GTF)


class RemoteReaderTest(unittest.TestCase):
    """Check the BGZF reader and the tabix region arithmetic with no network."""

    def test_two_blocks_decode_in_order(self) -> None:
        """The reader joins every complete block and reports the bytes it consumed."""
        payload = b"chr1\t100\ta\tN\t<DEL>\n"
        second = b"chr1\t200\tb\tN\t<DEL>\n"
        buffer = bgzf_block(payload) + bgzf_block(second)
        text, used = causal.decode_bgzf_blocks(buffer)
        self.assertEqual(text, payload + second)
        self.assertEqual(used, len(buffer))

    def test_a_truncated_final_block_is_left_alone(self) -> None:
        """A block that the fetched range cuts short yields no partial record."""
        first = bgzf_block(b"one\n")
        buffer = first + bgzf_block(b"two\n")[:-4]
        text, used = causal.decode_bgzf_blocks(buffer)
        self.assertEqual(text, b"one\n")
        self.assertEqual(used, len(first))

    def test_region_bins_cover_every_level(self) -> None:
        """A one base region occupies one bin of each of the six levels."""
        self.assertEqual(causal.region_bins(0, 1), [0, 1, 9, 73, 585, 4681])

    def test_merge_ranges_joins_a_gap_inside_the_slack(self) -> None:
        """Two ranges join when the gap is inside the slack and they stay apart otherwise."""
        self.assertEqual(causal.merge_ranges([(0, 10), (5, 20), (1000, 1100)], 0),
                         [(0, 20), (1000, 1100)])
        self.assertEqual(causal.merge_ranges([(0, 10), (1000, 1100)], 2000), [(0, 1100)])
        self.assertEqual(causal.merge_ranges([], 0), [])

    def test_region_filter_keeps_only_an_overlapping_record(self) -> None:
        """The region filter keeps a record that overlaps a wanted span and drops the rest."""
        text = (b"#comment\n"
                b"chr1\t1000\tin\tN\t<DEL>\t.\tPASS\tEND=1500;SVTYPE=DEL\n"
                b"chr1\t9000\tout\tN\t<DEL>\t.\tPASS\tEND=9100;SVTYPE=DEL\n"
                b"chr7\t1000\telsewhere\tN\t<DEL>\t.\tPASS\tEND=1500;SVTYPE=DEL\n")
        sink: dict = {}
        causal.collect_region_lines(text, {"chr1": [(900, 2000)]}, sink)
        self.assertEqual([key[2] for key in sink], ["in"])

    def test_the_trimmed_info_keeps_the_keys_that_the_scorer_reads(self) -> None:
        """The gnomAD subset keeps the class, the span and the frequency keys and drops the rest."""
        line = ("chr1\t1000\tg1\tN\t<DEL>\t10\tPASS\t"
                "ALGORITHMS=wham;SVTYPE=DEL;END=1500;SVLEN=500;AF=0.01;AF_nfe=0.02;NCR=0.5")
        trimmed = causal.trim_gnomad_info(line)
        self.assertEqual(len(trimmed.split("\t")), 8)
        info = scorer.parse_info(trimmed.split("\t")[7])
        self.assertEqual(set(info), {"SVTYPE", "END", "SVLEN", "AF", "AF_nfe"})


class GateTest(unittest.TestCase):
    """Check the gate verdicts and the default gene set of the scorer."""

    def test_the_default_gene_set_is_the_shipped_tuple(self) -> None:
        """The `--genes` default of the scorer is the shipped MVA gene tuple with BUB1B first."""
        self.assertEqual(causal.scorer_default_genes(), list(scorer.MVA_GENES))
        self.assertEqual(causal.scorer_default_genes()[0], "BUB1B")

    def test_a_gene_name_keeps_its_case(self) -> None:
        """The option parser keeps the case of a GENCODE gene name such as C11orf65."""
        self.assertEqual(scorer.gene_names("C11orf65"), ("C11orf65",))
        self.assertEqual(scorer.gene_names(""), scorer.MVA_GENES)

    def test_a_failing_gate_reads_fail(self) -> None:
        """A negative control that proposes a row fails gate six and the verdict stays in the map."""
        summary = {"spikes": 200,
                   "rank_1_fraction": {"estimate": 1.0, "low": 1.0, "high": 1.0},
                   "in_row_budget_fraction": {"estimate": 1.0, "low": 1.0, "high": 1.0},
                   "dropped_before_scoring_fraction": {"estimate": 0.0, "low": 0.0, "high": 0.1}}
        negative = {"default_gene_set": {"submission_rows": 2}}
        verdicts = causal.gate_verdicts(summary, negative, list(scorer.MVA_GENES))
        self.assertEqual(verdicts["g6_negative_control_clean"]["verdict"], "FAIL")
        self.assertEqual(verdicts["g3_rank_one_fraction"]["verdict"], "PASS")
        negative["default_gene_set"]["submission_rows"] = 0
        self.assertEqual(causal.gate_verdicts(summary, negative, list(scorer.MVA_GENES))
                         ["g6_negative_control_clean"]["verdict"], "PASS")

    def test_a_low_rank_one_fraction_fails_gate_three(self) -> None:
        """A rank 1 fraction below the floor fails gate three."""
        summary = {"spikes": 200,
                   "rank_1_fraction": {"estimate": 0.80, "low": 0.74, "high": 0.85},
                   "in_row_budget_fraction": {"estimate": 1.0, "low": 1.0, "high": 1.0},
                   "dropped_before_scoring_fraction": {"estimate": 0.0, "low": 0.0, "high": 0.1}}
        verdicts = causal.gate_verdicts(summary, {"default_gene_set": {"submission_rows": 0}},
                                        list(scorer.MVA_GENES))
        self.assertEqual(verdicts["g3_rank_one_fraction"]["verdict"], "FAIL")


class SummaryTest(unittest.TestCase):
    """Check the reported fractions."""

    def test_the_fractions_count_the_right_rows(self) -> None:
        """The three reported fractions count rank 1, the budget and the pre-scoring drop."""
        results = [
            {"rank": 1, "in_row_budget": 1, "reached_candidate_table": 1, "score": 105,
             "clears_candidate_floor": 1, "clears_submission_floor": 1, "candidate_rows": 1,
             "drop_reason": "", "svtype": "DEL", "gene_effect": "exonic", "clnsig": "pathogenic",
             "known_status": "new", "imprecise": 0, "background_calls_in_gene": 0,
             "top_background_score": ""},
            {"rank": 2, "in_row_budget": 1, "reached_candidate_table": 1, "score": 92,
             "clears_candidate_floor": 1, "clears_submission_floor": 1, "candidate_rows": 2,
             "drop_reason": "", "svtype": "DEL", "gene_effect": "intronic", "clnsig": "pathogenic",
             "known_status": "new", "imprecise": 1, "background_calls_in_gene": 1,
             "top_background_score": 105},
            {"rank": "", "in_row_budget": 0, "reached_candidate_table": 0, "score": 90,
             "clears_candidate_floor": 1, "clears_submission_floor": 1, "candidate_rows": 0,
             "drop_reason": "gnomad_frequency_above_0.01", "svtype": "DEL", "gene_effect": "exonic",
             "clnsig": "pathogenic", "known_status": "new", "imprecise": 0,
             "background_calls_in_gene": 0, "top_background_score": ""},
        ]
        summary = causal.summarise(results)
        self.assertEqual(summary["spikes"], 3)
        self.assertEqual(summary["rank_1_fraction"]["successes"], 1)
        self.assertEqual(summary["in_row_budget_fraction"]["successes"], 2)
        self.assertEqual(summary["dropped_before_scoring_fraction"]["successes"], 1)
        self.assertEqual(summary["drop_reasons"], {"gnomad_frequency_above_0.01": 1})
        self.assertEqual(summary["score_margin_over_top_background_call"]["minimum"], -13)


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=0).result
    passed = result.testsRun - len(result.failures) - len(result.errors)
    print(f"{passed} of {result.testsRun} checks passed.")
    sys.exit(0 if result.wasSuccessful() else 1)
