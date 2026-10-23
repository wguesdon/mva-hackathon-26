#!/usr/bin/env python3
"""Causal ranking benchmark for the Track 1 structural and mobile element score. Data class A.

The benchmark answers defect R2-2 of `reports/jury/defect_ledger_track1_round2_2026-09-17.md`. That
defect states that the eight point terms and the two score floors of
`pipelines/sv_mei_annotate/bin/annotate_sv_mei.py` carry no causal ranking measurement. This script
supplies one. It spikes a known pathogenic structural variant into a public background call set, runs
the shipped scorer once for each spike, and measures where the causal call lands.

Every input is public. The script opens no patient file and it writes no patient quantity.

The five steps are:

    1. Read the ClinVar structural variants of dbVar study nstd102 on GRCh38. Keep a record whose
       clinical significance is pathogenic or likely pathogenic, whose length is at least 50 base
       pairs and whose whole span falls inside one protein coding gene of GENCODE 50.
    2. Lift the GIAB HG002 structural variant benchmark v0.6 Tier 1 call set from GRCh37 to GRCh38
       with the liftover rule of `scripts/track1_sv_benchmark.py`. That lifted set is the background.
    3. Add one pathogenic record to the background, once for each spike, and run the shipped scorer
       with `--genes` set to the gene of that record.
    4. Read the rank, the score, the two score floors and the seven row submission budget.
    5. Run the scorer on the background with no spike, as the negative control.

Usage:
    uv run --with "pyliftover==0.4.1" python scripts/track1_sv_causal_benchmark.py all
    uv run python scripts/track1_sv_causal_benchmark.py fetch
    uv run python scripts/track1_sv_causal_benchmark.py select
    uv run python scripts/track1_sv_causal_benchmark.py run
"""

from __future__ import annotations

import argparse
import bisect
import collections
import concurrent.futures
import csv
import datetime as dt
import gzip
import hashlib
import json
import platform
import random
import re
import io
import statistics
import shutil
import struct
import sys
import threading
import urllib.parse
import zlib
from contextlib import redirect_stdout
from pathlib import Path

import httpx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "pipelines" / "sv_mei_annotate" / "bin"))

import annotate_sv_mei as scorer  # noqa: E402
import track1_interval_stats as stats  # noqa: E402
import track1_sv_benchmark as bench  # noqa: E402
import track1_sv_common as common  # noqa: E402

RECEIPT_ID = "t1_sv_causal_benchmark_v1"
CACHE = ROOT / "tracking" / "cache" / "t1_sv_causal"
RESULTS_TABLE = ROOT / "tracking" / "t1_sv_causal_benchmark_2026_09_17.tsv"
RECEIPT = ROOT / "tracking" / "receipts" / f"{RECEIPT_ID}.json"
KNOWN_INTERVALS = ROOT / "reports" / "known_alleles" / "s6_targets_2026_09_15" / "intervals.bed"

# Public sources. Every entry carries the URL, the release name and the licence.
SOURCES = {
    "dbvar_nstd102": {
        "url": ("https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/"
                "nstd102.GRCh38.variant_call.vcf.gz"),
        "name": "nstd102.GRCh38.variant_call.vcf.gz",
        "release": "dbVar study nstd102, the ClinVar structural variants, GRCh38 variant calls",
        "published": "The file header records fileDate 20260616. The server reports Last-Modified 2026-07-09.",
        "licence": "Public domain. NCBI asks for a citation of dbVar, PMID 23193291.",
    },
    "dbvar_nstd102_index": {
        "url": ("https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/"
                "nstd102.GRCh38.variant_call.vcf.gz.tbi"),
        "name": "nstd102.GRCh38.variant_call.vcf.gz.tbi",
        "release": ("The tabix index of dbVar study nstd102 on GRCh38. The benchmark streams the "
                    "whole 5,706,478 byte study file and it reads no range of it, so this index is "
                    "recorded for provenance and it is not read."),
        "published": "The server reports Last-Modified 2026-07-09.",
        "licence": "Public domain.",
    },
    "gencode": {
        "url": ("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/"
                "gencode.v50.basic.annotation.gtf.gz"),
        "name": "gencode.v50.basic.annotation.gtf.gz",
        "release": "GENCODE 50 basic annotation, GRCh38",
        "published": "GENCODE release 50. The server reports Last-Modified 2026-06-09.",
        "licence": "Free for all use. GENCODE asks for a citation of PMID 30357393.",
    },
    "giab_v06": {
        "url": ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
                "analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz"),
        "name": "HG002_SVs_Tier1_v0.6.vcf.gz",
        "release": "GIAB HG002 structural variant benchmark v0.6, Tier 1, GRCh37",
        "published": "The server reports Last-Modified 2018-06-12.",
        "licence": "Public domain. NIST asks for a citation of doi:10.1101/664623.",
    },
    "chain": {
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
        "name": "hg19ToHg38.over.chain.gz",
        "release": "UCSC liftOver chain file, hg19 to hg38",
        "published": "The server reports Last-Modified 2014-01-01.",
        "licence": "UCSC genome annotation data are free for all use, with attribution.",
    },
    "gnomad_sv_index": {
        "url": ("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/"
                "gnomad.v4.1.sv.sites.vcf.gz.tbi"),
        "name": "gnomad.v4.1.sv.sites.vcf.gz.tbi",
        "release": "The tabix index of gnomAD SV v4.1 sites, GRCh38",
        "published": "gnomAD release 4.1. The server reports Last-Modified 2024-04-17.",
        "licence": "Free for all use under CC0. gnomAD asks for a citation of PMID 38057664.",
    },
}

GNOMAD_SV_URL = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/"
                 "gnomad.v4.1.sv.sites.vcf.gz")

# The selection parameters. Every value is fixed before the run.
SPIKE_COUNT = 200
SEED = 20260917
SPIKE_GENOTYPE = "0/1"
SPIKE_GENOTYPE_QUALITY = 99
MIN_SPIKE_LENGTH = 50
PATHOGENIC_TERMS = ("pathogenic", "likely pathogenic")
PRIMARY_CONTIGS = tuple(str(name) for name in list(range(1, 23)) + ["X", "Y"])
# The regulatory window of the scorer. A frequency lookup can reach this far outside a gene.
GNOMAD_PAD = scorer.REGULATORY_WINDOW
# The gnomAD INFO keys that `annotate_sv_mei.read_gnomad_sv` reads. The subset keeps these and drops
# every other key, which changes no computed value.
GNOMAD_KEPT_KEYS = ("SVTYPE", "END", "SVLEN", "AF") + tuple(
    f"AF_{population}" for population in scorer.GNOMAD_POPULATIONS)

# The gate conditions. They are fixed before the result is known and they stay in the record when
# they fail.
GATES = {
    "g1_benchmark_has_power": {
        "condition": "At least 100 pathogenic structural variants enter the benchmark as spikes.",
        "floor": 100,
    },
    "g2_default_gene_set_unchanged": {
        "condition": ("The `--genes` default of annotate_sv_mei.py equals the shipped MVA_GENES "
                      "tuple, the first default name is BUB1B, harness/test_sv_mei_annotate.py "
                      "passes with no change, and a run that omits the option writes the same six "
                      "output files as a run that passes the default value. The last condition is "
                      "the check DefaultBehaviourTest of "
                      "harness/test_track1_sv_causal_benchmark.py, which compares the SHA-256 of "
                      "every output file of the two runs."),
        "expected": list(scorer.MVA_GENES),
    },
    "g3_rank_one_fraction": {
        "condition": ("The fraction of spikes that reach rank 1 of the candidate table is at least "
                      "0.90 and the lower limit of its Wilson interval is at least 0.85."),
        "floor": 0.90,
        "interval_floor": 0.85,
    },
    "g4_budget_fraction": {
        "condition": ("The fraction of spikes that land inside the seven row submission budget is at "
                      "least 0.90 and the lower limit of its Wilson interval is at least 0.85."),
        "floor": 0.90,
        "interval_floor": 0.85,
    },
    "g5_pre_scoring_drop": {
        "condition": ("The fraction of spikes that a filter drops before the ranking is at most "
                      "0.05."),
        "ceiling": 0.05,
    },
    "g6_negative_control_clean": {
        "condition": ("The background call set with no spike proposes no submission row under the "
                      "shipped default gene set."),
        "expected": 0,
    },
}

BGZF_MAGIC = b"\x1f\x8b\x08\x04"
BGZF_MAX_BLOCK = 65536
RANGE_MERGE_SLACK = 65536
RANGE_WORKERS = 8

COMBINED_HEADER = (
    "##fileformat=VCFv4.2",
    "##reference=GRCh38",
    '##ALT=<ID=DEL,Description="Deletion">',
    '##ALT=<ID=DUP,Description="Duplication">',
    '##ALT=<ID=INS,Description="Insertion">',
    '##ALT=<ID=INV,Description="Inversion">',
    '##ALT=<ID=CNV,Description="Copy number variable region">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002",
)

RESULT_FIELDS = ("spike_id", "gene", "chrom", "pos", "end", "svlen", "svtype", "clnsig",
                 "clinvar_accession", "imprecise", "reached_candidate_table", "rank",
                 "candidate_rows", "score", "score_components", "gene_effect",
                 "gnomad_popmax_af", "clears_candidate_floor", "clears_submission_floor",
                 "in_row_budget", "submission_rows", "drop_reason", "known_status",
                 "background_calls_in_gene", "top_background_score")


def utc_now() -> str:
    """Return the current time in UTC.

    Returns:
        The time as an ISO 8601 string with a Z suffix.
    """
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_of(path: Path) -> str:
    """Compute the SHA-256 checksum of a file.

    Args:
        path: The file to read.

    Returns:
        The checksum in lower case hexadecimal.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


# --- Step 1. Fetch the public files. -----------------------------------------------------------


def download(url: str, target: Path, client: httpx.Client) -> None:
    """Download one public file to the cache.

    Args:
        url: The source URL.
        target: The cache path.
        client: An HTTP client.
    """
    temporary = target.with_suffix(target.suffix + ".part")
    with client.stream("GET", url) as reply:
        reply.raise_for_status()
        with temporary.open("wb") as handle:
            for block in reply.iter_bytes(1 << 20):
                handle.write(block)
    temporary.replace(target)


def fetch_sources() -> dict:
    """Download every public source file that the cache does not hold.

    Returns:
        One manifest entry for each source with the URL, the release, the licence, the checksum and
        the fetch date.
    """
    CACHE.mkdir(parents=True, exist_ok=True)
    manifest = {}
    client = httpx.Client(timeout=900.0, follow_redirects=True)
    try:
        for key, source in SOURCES.items():
            target = CACHE / source["name"]
            if not target.exists():
                print(f"fetch {source['name']}")
                download(source["url"], target, client)
            manifest[key] = {"url": source["url"], "name": source["name"],
                             "release": source["release"], "licence": source["licence"],
                             "published": source["published"],
                             "sha256": sha256_of(target), "fetched_utc": utc_now()}
    finally:
        client.close()
    manifest["gnomad_sv_sites"] = {
        "url": GNOMAD_SV_URL,
        "name": "gnomad.v4.1.sv.sites.vcf.gz",
        "release": "gnomAD SV v4.1 sites, GRCh38",
        "published": "gnomAD release 4.1. The server reports Last-Modified 2024-04-17.",
        "licence": "Free for all use under CC0. gnomAD asks for a citation of PMID 38057664.",
        "access": ("The whole file is 1,743,544,245 bytes and it is not downloaded. The subset is "
                   "read by tabix range queries over https with HTTP range requests, using the "
                   "published tabix index."),
        "fetched_utc": utc_now(),
    }
    (CACHE / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


# --- A tabix range query over https, with no htslib. -------------------------------------------


def decode_bgzf_blocks(buffer: bytes) -> tuple[bytes, int]:
    """Decompress every complete BGZF block of a byte buffer.

    The buffer must start on a BGZF block boundary. A truncated final block is left alone, so a
    caller that over-fetches by one block length loses no record.

    Args:
        buffer: The raw bytes of one or more BGZF blocks.

    Returns:
        The decompressed bytes and the number of input bytes consumed.
    """
    parts = []
    offset = 0
    size = len(buffer)
    while offset + 18 <= size:
        if buffer[offset:offset + 4] != BGZF_MAGIC:
            break
        extra_length = int.from_bytes(buffer[offset + 10:offset + 12], "little")
        extra = buffer[offset + 12:offset + 12 + extra_length]
        block_size = None
        cursor = 0
        while cursor + 4 <= len(extra):
            field_length = int.from_bytes(extra[cursor + 2:cursor + 4], "little")
            if extra[cursor] == 0x42 and extra[cursor + 1] == 0x43 and field_length == 2:
                block_size = int.from_bytes(extra[cursor + 4:cursor + 6], "little") + 1
                break
            cursor += 4 + field_length
        if block_size is None or offset + block_size > size:
            break
        parts.append(zlib.decompress(buffer[offset:offset + block_size], 31))
        offset += block_size
    return b"".join(parts), offset


def read_tabix_index(path: Path) -> dict:
    """Read a tabix index file.

    Args:
        path: The `.tbi` file.

    Returns:
        A map with the column layout and, for each contig, the bin to chunk map and the linear
        index.

    Raises:
        ValueError: When the file is not a tabix index.
    """
    raw = gzip.decompress(path.read_bytes())
    if raw[:4] != b"TBI\x01":
        raise ValueError(f"{path} is not a tabix index")
    file_format, column_sequence, column_begin, column_end, meta, skip, name_length = \
        struct.unpack_from("<7i", raw, 8)
    names = [item.decode() for item in raw[36:36 + name_length].split(b"\x00") if item]
    cursor = 36 + name_length
    refs: dict[str, dict] = {}
    for name in names:
        (bin_count,) = struct.unpack_from("<i", raw, cursor)
        cursor += 4
        bins: dict[int, list[tuple[int, int]]] = {}
        for _ in range(bin_count):
            bin_id, chunk_count = struct.unpack_from("<Ii", raw, cursor)
            cursor += 8
            chunks = []
            for _ in range(chunk_count):
                begin, end = struct.unpack_from("<QQ", raw, cursor)
                cursor += 16
                chunks.append((begin, end))
            bins[bin_id] = chunks
        (interval_count,) = struct.unpack_from("<i", raw, cursor)
        cursor += 4
        linear = list(struct.unpack_from(f"<{interval_count}Q", raw, cursor)) if interval_count else []
        cursor += 8 * interval_count
        refs[name] = {"bins": bins, "linear": linear}
    return {"format": file_format, "column_sequence": column_sequence,
            "column_begin": column_begin, "column_end": column_end,
            "meta": meta, "skip": skip, "refs": refs}


def region_bins(begin: int, end: int) -> list[int]:
    """List the tabix bin numbers that a region can occupy.

    Args:
        begin: The zero based start of the region.
        end: The zero based end of the region, exclusive.

    Returns:
        Every bin number of the six level binning scheme of tabix.
    """
    result = [0]
    end -= 1
    for shift, offset in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        result.extend(range(offset + (begin >> shift), offset + (end >> shift) + 1))
    return result


def region_chunks(index: dict, chrom: str, begin: int, end: int) -> list[tuple[int, int]]:
    """List the index chunks that can hold a record of one region.

    Args:
        index: The parsed tabix index.
        chrom: The contig name as the index spells it.
        begin: The zero based start of the region.
        end: The zero based end of the region, exclusive.

    Returns:
        The virtual offset pairs of every candidate chunk.
    """
    ref = index["refs"].get(chrom)
    if ref is None:
        return []
    linear = ref["linear"]
    floor = linear[min(begin >> 14, len(linear) - 1)] if linear else 0
    found = []
    for bin_id in region_bins(begin, end):
        for chunk in ref["bins"].get(bin_id, ()):
            if chunk[1] > floor:
                found.append(chunk)
    return found


def merge_ranges(ranges: list[tuple[int, int]], slack: int = RANGE_MERGE_SLACK) -> list[tuple[int, int]]:
    """Merge overlapping or nearby byte ranges.

    Args:
        ranges: The byte ranges as a start and an end.
        slack: Two ranges join when the gap between them is at most this many bytes.

    Returns:
        The merged ranges, in order.
    """
    if not ranges:
        return []
    ordered = sorted(ranges)
    merged = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start <= merged[-1][1] + slack:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


CLIENT_STORE = threading.local()


def thread_client() -> httpx.Client:
    """Return the HTTP client of the calling thread.

    One client for each thread keeps the connection open across range requests.

    Returns:
        The client of this thread.
    """
    client = getattr(CLIENT_STORE, "client", None)
    if client is None:
        client = httpx.Client(timeout=300.0, follow_redirects=True)
        CLIENT_STORE.client = client
    return client


def fetch_range(job: tuple[str, int, int]) -> bytes:
    """Fetch one byte range of a remote file.

    Args:
        job: The URL, the first byte and the last byte.

    Returns:
        The raw bytes of the range.
    """
    url, start, end = job
    reply = thread_client().get(url, headers={"Range": f"bytes={start}-{end}"})
    reply.raise_for_status()
    return reply.content


def tabix_query(url: str, index: dict, regions: list[tuple[str, int, int]],
                workers: int = RANGE_WORKERS) -> tuple[list[str], dict]:
    """Read every record of a remote bgzip file that falls inside a set of regions.

    Args:
        url: The URL of the bgzip compressed file.
        index: The parsed tabix index of that file.
        regions: The regions as a contig, a one based start and a one based end.
        workers: The number of parallel range requests.

    Returns:
        The matching data lines and the counts of the transfer.
    """
    wanted: dict[str, list[tuple[int, int]]] = {}
    ranges: list[tuple[int, int]] = []
    chunk_total = 0
    for chrom, start, end in regions:
        wanted.setdefault(chrom, []).append((start, end))
        chunks = region_chunks(index, chrom, max(0, start - 1), end)
        chunk_total += len(chunks)
        for begin, finish in chunks:
            ranges.append((begin >> 16, (finish >> 16) + BGZF_MAX_BLOCK))
    for spans in wanted.values():
        spans.sort()
    merged = merge_ranges(ranges)
    jobs = [(url, start, end) for start, end in merged]
    lines: dict[tuple[str, int, str], str] = {}
    fetched = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for position, payload in enumerate(pool.map(fetch_range, jobs), start=1):
            fetched += len(payload)
            text, _ = decode_bgzf_blocks(payload)
            collect_region_lines(text, wanted, lines)
            if position % 200 == 0:
                print(f"  range {position} of {len(jobs)}, {fetched / 1e6:.0f} MB")
    counts = {"regions": len(regions), "index_chunks": chunk_total, "range_requests": len(jobs),
              "bytes_fetched": fetched, "records_kept": len(lines)}
    return [lines[key] for key in sorted(lines, key=lambda item: (item[0], item[1]))], counts


def collect_region_lines(text: bytes, wanted: dict[str, list[tuple[int, int]]],
                         sink: dict[tuple[str, int, str], str]) -> None:
    """Keep the VCF data lines of a decompressed buffer that fall inside the wanted regions.

    Args:
        text: The decompressed bytes of one fetched range.
        wanted: The sorted one based regions by contig.
        sink: The record store, keyed by contig, position and identifier.
    """
    for line in text.decode("utf-8", errors="replace").split("\n"):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t", 8)
        if len(fields) < 8:
            continue
        spans = wanted.get(fields[0])
        if not spans:
            continue
        try:
            pos = int(fields[1])
        except ValueError:
            continue
        info = scorer.parse_info(fields[7])
        end, _ = scorer.call_span(pos, fields[3], info)
        if not spans_overlap(spans, pos, end):
            continue
        sink[(fields[0], pos, fields[2])] = line


def spans_overlap(spans: list[tuple[int, int]], start: int, end: int) -> bool:
    """Say whether a record overlaps any span of a span list.

    The list holds the benchmark regions of one contig, so it is short and a full scan is cheap.

    Args:
        spans: The one based spans.
        start: The record start.
        end: The record end.

    Returns:
        True when one span overlaps the record.
    """
    for span_start, span_end in spans:
        if span_start <= end and start <= span_end:
            return True
    return False


# --- Step 2. Select the pathogenic structural variants. ----------------------------------------


def clinical_terms(raw: str) -> list[str]:
    """Split a dbVar CLNSIG value into lower case terms.

    Args:
        raw: The CLNSIG value, URL escaped and quoted.

    Returns:
        One lower case term for each entry.
    """
    text = urllib.parse.unquote(raw).replace('"', "")
    return [item.strip().lower() for item in text.split(",") if item.strip()]


def record_length(pos: int, info: dict[str, str]) -> int | None:
    """Return the length of a dbVar structural variant record.

    Args:
        pos: The POS column.
        info: The parsed INFO column.

    Returns:
        The length in base pairs, or None when neither SVLEN nor END gives one.
    """
    raw = info.get("SVLEN", "")
    if raw not in ("", "."):
        try:
            return abs(int(float(raw.split(",")[0])))
        except ValueError:
            pass
    end = info.get("END", "")
    if end not in ("", "."):
        try:
            return abs(int(end) - pos)
        except ValueError:
            return None
    return None


def record_span(pos: int, info: dict[str, str], svtype: str) -> tuple[int | None, int | None]:
    """Return the reference end and the length of one dbVar record.

    An insertion carries no reference span, so its end is its position and its length is the
    inserted length. Every other class takes its end from the END key when the record holds one.

    Args:
        pos: The POS column.
        info: The parsed INFO column.
        svtype: The structural variant class.

    Returns:
        The one based end and the length in base pairs. Both are None when neither key gives a
        length.
    """
    length = record_length(pos, info)
    if svtype == "INS":
        return (pos, length) if length is not None else (None, None)
    raw = info.get("END", "")
    if raw not in ("", "."):
        try:
            end = int(raw)
        except ValueError:
            return None, None
    elif length is not None:
        end = pos + length
    else:
        return None, None
    return end, max(end - pos, 0)


def read_protein_coding_genes(path: Path) -> dict[str, list[tuple[int, int, str]]]:
    """Read the span of every protein coding gene of a GENCODE GTF.

    Args:
        path: The GENCODE annotation file, gzip compressed.

    Returns:
        A map from contig to a sorted list of the start, the end and the gene name.
    """
    genes: dict[str, list[tuple[int, int, str]]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            if len(fields) < 9 or fields[2] != "gene":
                continue
            if 'gene_type "protein_coding"' not in fields[8]:
                continue
            match = re.search(r'gene_name "([^"]+)"', fields[8])
            if not match:
                continue
            genes.setdefault(fields[0], []).append((int(fields[3]), int(fields[4]), match.group(1)))
    for spans in genes.values():
        spans.sort()
    return genes


def gene_span_index(genes: dict[str, list[tuple[int, int, str]]]
                    ) -> dict[str, tuple[list[tuple[int, int, str]], list[int], int]]:
    """Build a search index over the gene spans of each contig.

    Args:
        genes: The sorted gene spans by contig.

    Returns:
        A map from contig to the span list, the start coordinates and the longest span.
    """
    index = {}
    for chrom, spans in genes.items():
        index[chrom] = (spans, [start for start, _, _ in spans],
                        max((end - start for start, end, _ in spans), default=0))
    return index


def unique_gene_names(genes: dict[str, list[tuple[int, int, str]]]) -> tuple[set[str], int]:
    """Return the gene names that carry exactly one gene line in the annotation.

    `annotate_sv_mei.read_gene_models` keys its disease gene record by name, so a name that the
    annotation gives twice is ambiguous and it cannot enter the benchmark.

    Args:
        genes: The gene spans by contig.

    Returns:
        The unambiguous names and the number of names that the annotation gives more than once.
    """
    counts = collections.Counter(name for spans in genes.values() for _, _, name in spans)
    repeated = {name for name, count in counts.items() if count > 1}
    return set(counts) - repeated, len(repeated)


def containing_genes(index: dict[str, tuple[list[tuple[int, int, str]], list[int], int]],
                     chrom: str, start: int, end: int) -> list[str]:
    """List the protein coding genes whose span holds a record whole.

    Args:
        index: The gene span index.
        chrom: The contig of the record, with the chr prefix.
        start: The one based start of the record.
        end: The one based end of the record.

    Returns:
        The gene names, shortest gene first.
    """
    entry = index.get(chrom)
    if entry is None:
        return []
    spans, starts, longest = entry
    first = bisect.bisect_left(starts, start - longest)
    found = []
    for gene_start, gene_end, name in spans[first:]:
        if gene_start > start:
            break
        if end <= gene_end:
            found.append((gene_end - gene_start, name))
    found.sort()
    return [name for _, name in found]


def select_pathogenic(vcf: Path, index: dict, unambiguous: set[str]) -> tuple[list[dict], dict]:
    """Apply the selection filters to the dbVar records and report what each filter kept.

    Args:
        vcf: The dbVar nstd102 GRCh38 variant call VCF.
        index: The gene span index of the protein coding genes.
        unambiguous: The gene names that the annotation gives exactly once.

    Returns:
        One entry for each record that passes every filter, and the count that each filter kept.
    """
    counts = collections.Counter()
    kept: list[dict] = []
    with gzip.open(vcf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            counts["records_in_file"] += 1
            fields = line.rstrip("\n").split("\t")
            info = scorer.parse_info(fields[7])
            terms = clinical_terms(info.get("CLNSIG", ""))
            if not any(term in PATHOGENIC_TERMS for term in terms):
                continue
            counts["pathogenic_or_likely_pathogenic"] += 1
            if fields[0] not in PRIMARY_CONTIGS:
                continue
            counts["on_a_primary_contig"] += 1
            pos = int(fields[1])
            svtype = info.get("SVTYPE", "NA").upper().split(":")[0]
            end, length = record_span(pos, info, svtype)
            if length is None or length < MIN_SPIKE_LENGTH:
                continue
            counts["at_least_50_base_pairs"] += 1
            chrom = f"chr{fields[0]}"
            names = containing_genes(index, chrom, pos, end)
            if not names:
                continue
            counts["inside_a_protein_coding_gene"] += 1
            names = [name for name in names if name in unambiguous]
            if not names:
                continue
            counts["gene_name_is_unambiguous"] += 1
            kept.append({
                "spike_id": fields[2],
                "chrom": chrom,
                "pos": pos,
                "end": end,
                "svlen": length,
                "svtype": svtype,
                "gene": names[0],
                "gene_count": len(names),
                "clnsig": ";".join(terms),
                "clinvar_accession": urllib.parse.unquote(info.get("CLNACC", "")).replace('"', ""),
                "imprecise": "IMPRECISE" in info,
            })
    counts["distinct_loci"] = len({(item["chrom"], item["pos"], item["end"], item["svtype"])
                                   for item in kept})
    counts["distinct_genes"] = len({item["gene"] for item in kept})
    return kept, dict(counts)


def draw_spikes(records: list[dict], count: int = SPIKE_COUNT, seed: int = SEED) -> list[dict]:
    """Draw the benchmark spikes, at most one for each gene.

    One locus can carry several dbVar submissions and one gene can carry many loci. The draw keeps
    one record for each locus and one locus for each gene, so no gene dominates the benchmark.

    Args:
        records: The records that passed every filter.
        count: The number of spikes to draw.
        seed: The seed of the draw.

    Returns:
        The drawn spikes, ordered by contig and position.
    """
    by_locus: dict[tuple, dict] = {}
    for record in records:
        key = (record["chrom"], record["pos"], record["end"], record["svtype"])
        if key not in by_locus or record["spike_id"] < by_locus[key]["spike_id"]:
            by_locus[key] = record
    by_gene: dict[str, list[dict]] = {}
    for record in by_locus.values():
        by_gene.setdefault(record["gene"], []).append(record)
    generator = random.Random(seed)
    one_per_gene = []
    for gene in sorted(by_gene):
        choices = sorted(by_gene[gene], key=lambda item: item["spike_id"])
        one_per_gene.append(choices[generator.randrange(len(choices))])
    drawn = one_per_gene if len(one_per_gene) <= count else generator.sample(one_per_gene, count)
    drawn.sort(key=lambda item: (item["chrom"], item["pos"]))
    return drawn


# --- Step 3. Build the smallest honest inputs. --------------------------------------------------


def write_gene_subset(gencode: Path, target: Path, wanted: set[str]) -> int:
    """Write a GENCODE subset that holds the gene and exon lines of the benchmark genes.

    Args:
        gencode: The full GENCODE annotation, gzip compressed.
        target: The output path, gzip compressed.
        wanted: The gene names to keep.

    Returns:
        The number of lines written.
    """
    written = 0
    with gzip.open(gencode, "rt") as source, gzip.open(target, "wt") as sink:
        sink.write("##description: GENCODE 50 basic annotation, the gene and exon lines of the "
                   "benchmark genes only\n")
        for line in source:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            if len(fields) < 9 or fields[2] not in ("gene", "exon"):
                continue
            match = re.search(r'gene_name "([^"]+)"', fields[8])
            if not match or match.group(1) not in wanted:
                continue
            sink.write(line)
            written += 1
    return written


def gene_regions(gencode_subset: Path, wanted: set[str]) -> list[tuple[str, int, int]]:
    """Return the padded span of each benchmark gene.

    Args:
        gencode_subset: The GENCODE subset, gzip compressed.
        wanted: The gene names.

    Returns:
        One region for each gene as a contig, a padded start and a padded end.
    """
    spans: dict[tuple[str, str], tuple[int, int]] = {}
    with gzip.open(gencode_subset, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            if fields[2] != "gene":
                continue
            match = re.search(r'gene_name "([^"]+)"', fields[8])
            if not match or match.group(1) not in wanted:
                continue
            key = (fields[0], match.group(1))
            start, end = int(fields[3]), int(fields[4])
            if key in spans:
                spans[key] = (min(spans[key][0], start), max(spans[key][1], end))
            else:
                spans[key] = (start, end)
    regions = [(chrom, max(1, start - GNOMAD_PAD), end + GNOMAD_PAD)
               for (chrom, _), (start, end) in spans.items()]
    regions.sort()
    return regions


def trim_gnomad_info(line: str) -> str:
    """Keep only the gnomAD INFO keys that the scorer reads.

    Args:
        line: One gnomAD SV sites data line.

    Returns:
        The same record with a reduced INFO column and no sample column.
    """
    fields = line.split("\t")
    info = scorer.parse_info(fields[7])
    parts = [f"{key}={info[key]}" for key in GNOMAD_KEPT_KEYS if key in info and info[key] != ""]
    return "\t".join(fields[:7] + [";".join(parts) or "."])


def write_gnomad_subset(target: Path, lines: list[str]) -> int:
    """Write the gnomAD SV subset that the benchmark loci need.

    Args:
        target: The output path, gzip compressed.
        lines: The gnomAD data lines of the benchmark regions.

    Returns:
        The number of records written.
    """
    header = ["##fileformat=VCFv4.2",
              "##source=gnomAD SV v4.1 sites, the benchmark loci only",
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    with gzip.open(target, "wt") as handle:
        handle.write("\n".join(header) + "\n")
        for line in lines:
            handle.write(trim_gnomad_info(line) + "\n")
    return len(lines)


def lift_background(giab: Path, chain: Path, target: Path) -> dict:
    """Lift the GIAB HG002 v0.6 Tier 1 call set to GRCh38 and write it as the background.

    The rule is the rule of `scripts/track1_sv_benchmark.py`, which the receipt
    `t1_liftover_control_v1` already measured on the same file.

    Args:
        giab: The GIAB v0.6 Tier 1 VCF on GRCh37.
        chain: The UCSC hg19 to hg38 chain file.
        target: The output VCF path, plain text.

    Returns:
        The counts of the lift.
    """
    from pyliftover import LiftOver

    lifter = LiftOver(str(chain))
    _, records = common.read_vcf(giab)
    lifted, counts = bench.lift_truth_records(records, lifter)
    order = {name: position for position, name in enumerate(common.PRIMARY_CONTIGS)}
    lifted.sort(key=lambda record: (order[record.chrom], record.pos))
    header = list(COMBINED_HEADER)
    with target.open("w") as handle:
        handle.write("\n".join(header) + "\n")
        for record in lifted:
            handle.write("\t".join([
                record.chrom, str(record.pos), record.identifier, record.ref, record.alt,
                record.qual, record.filter_field, scorer_info(record.info),
                "GT", record.tail[1] if len(record.tail) > 1 else "./."]) + "\n")
    return counts


def scorer_info(info: dict[str, str]) -> str:
    """Join an INFO mapping into a VCF INFO column.

    Args:
        info: The key to value mapping.

    Returns:
        The INFO column.
    """
    return ";".join(f"{key}={value}" for key, value in info.items()) or "."


# --- Step 4. Place one spike and run the shipped scorer. ---------------------------------------


def spike_line(spike: dict) -> str:
    """Render one pathogenic record as a VCF data line.

    The record carries the symbolic alternate allele of its class, a PASS filter, the heterozygous
    genotype and a genotype quality above the filter floor of the scorer. The dbVar confidence
    interval keys are dropped, because the scorer reads none of them.

    Args:
        spike: The selected pathogenic record.

    Returns:
        The data line, with no trailing newline.
    """
    length = spike["svlen"]
    signed = -length if spike["svtype"] == "DEL" else length
    info = f"SVTYPE={spike['svtype']};SVLEN={signed};END={spike['end']}"
    return "\t".join([spike["chrom"], str(spike["pos"]), spike["spike_id"], "N",
                      f"<{spike['svtype']}>", ".", "PASS", info, "GT:GQ",
                      f"{SPIKE_GENOTYPE}:{SPIKE_GENOTYPE_QUALITY}"])


def combined_lines(background: list[str], spike: dict | None) -> list[str]:
    """Build the combined call set of the background and at most one spike.

    Args:
        background: The background data lines, in coordinate order.
        spike: The pathogenic record to add, or None for the negative control.

    Returns:
        The header lines and the data lines, with the spike in coordinate order.
    """
    lines = list(COMBINED_HEADER)
    if spike is None:
        return lines + background
    added = spike_line(spike)
    position = 0
    order = {name: index for index, name in enumerate(common.PRIMARY_CONTIGS)}
    key = (order.get(spike["chrom"], 99), spike["pos"])
    for position, line in enumerate(background):
        fields = line.split("\t", 2)
        if (order.get(fields[0], 99), int(fields[1])) > key:
            return lines + background[:position] + [added] + background[position:]
    return lines + background + [added]


def run_scorer(sv_vcf: Path, mei_vcf: Path, gtf: Path, gnomad: Path, intervals: Path,
               out_dir: Path, genes: str, run_id: str) -> dict:
    """Run the shipped scorer once and read its outputs.

    Args:
        sv_vcf: The combined structural variant call set.
        mei_vcf: The mobile element call set, which holds no record in this benchmark.
        gtf: The GENCODE subset.
        gnomad: The gnomAD SV subset.
        intervals: The known BUB1B interval file.
        out_dir: The output directory of this run.
        genes: The value of the `--genes` option.
        run_id: The run identifier.

    Returns:
        The gate file, the annotated rows, the candidate rows and the submission rows.

    Raises:
        RuntimeError: When the scorer returns a FAIL verdict.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    argv = ["annotate_sv_mei.py", "--sv-vcf", str(sv_vcf), "--mei-vcf", str(mei_vcf),
            "--gtf", str(gtf), "--gnomad-sv", str(gnomad), "--intervals", str(intervals),
            "--out-dir", str(out_dir), "--genes", genes, "--run-id", run_id,
            "--min-frequency-records", "1"]
    saved = sys.argv
    sys.argv = argv
    try:
        with redirect_stdout(io.StringIO()):
            code = scorer.main()
    finally:
        sys.argv = saved
    if code != 0:
        raise RuntimeError(f"The scorer returned a FAIL verdict for {run_id}")
    return {
        "gates": json.loads((out_dir / "class_a_gates.json").read_text()),
        "annotated": read_table(out_dir / "class_c_sv_mei_annotated.tsv"),
        "candidates": read_table(out_dir / "class_b2_sv_mei_ranked.tsv"),
        "submission": read_table(out_dir / "class_b2_submission_rows.csv", ","),
    }


def scorer_default_genes() -> list[str]:
    """Read the `--genes` default of the shipped scorer from its own argument parser.

    Returns:
        The default gene names, in the order of the option.
    """
    probe = ["annotate_sv_mei.py", "--sv-vcf", "x", "--mei-vcf", "x", "--gtf", "x",
             "--gnomad-sv", "x", "--intervals", "x", "--out-dir", "x"]
    saved = sys.argv
    sys.argv = probe
    try:
        parsed = scorer.build_arguments()
    finally:
        sys.argv = saved
    return list(scorer.gene_names(parsed.genes))


def read_table(path: Path, delimiter: str = "\t") -> list[dict]:
    """Read one delimited table.

    Args:
        path: The table file.
        delimiter: The field delimiter.

    Returns:
        One row for each record.
    """
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def rank_of(candidates: list[dict], call_id: str) -> int | None:
    """Return the rank of one call in the candidate table.

    Args:
        candidates: The rows of the ranked candidate table.
        call_id: The identifier of the call.

    Returns:
        The rank, or None when the call is not in the table.
    """
    for row in candidates:
        if row["call_id"] == call_id:
            return int(row["rank"])
    return None


def in_row_budget(submission: list[dict], chrom: str, pos: int) -> bool:
    """Say whether one call holds a row of the submission budget.

    Args:
        submission: The rows of the proposed submission file.
        chrom: The contig of the call.
        pos: The position of the call.

    Returns:
        True when a submission row carries this coordinate.
    """
    return any(row["chrom_1"] == chrom and int(row["pos_1"]) == pos for row in submission)


def drop_reason(row: dict) -> str:
    """Name the filter that removed one annotated call from the candidate table.

    The filters run in the order of `annotate_sv_mei.select_candidates`.

    Args:
        row: One row of the annotated table.

    Returns:
        The filter name, or an empty string when no filter removed the call.
    """
    call = {"genotype": row["genotype"], "svtype": row["svtype"],
            "svlen": int(row["svlen"]) if row["svlen"] else 0,
            "source": row["source"],
            "genotype_quality": int(row["genotype_quality"]) if row["genotype_quality"] else None,
            "retroseq_level": int(row["retroseq_level"]) if row["retroseq_level"] else None}
    if not scorer.passes_quality(call):
        return "caller_quality_filter"
    if not row["mva_gene"]:
        return "no_disease_gene_overlap"
    frequency = row["gnomad_popmax_af"]
    if frequency and float(frequency) > scorer.MAX_GNOMAD_AF:
        return "gnomad_frequency_above_0.01"
    if int(row["score"]) < scorer.CANDIDATE_SCORE_FLOOR:
        return "below_candidate_score_floor_25"
    return ""


def measure_spike(spike: dict, outputs: dict) -> dict:
    """Read the rank, the score and the budget outcome of one spike.

    Args:
        spike: The pathogenic record that was added.
        outputs: The parsed outputs of the scorer run.

    Returns:
        One result row.
    """
    annotated = {row["call_id"]: row for row in outputs["annotated"]}
    row = annotated.get(spike["spike_id"], {})
    candidates = outputs["candidates"]
    rank = rank_of(candidates, spike["spike_id"])
    score = int(row["score"]) if row.get("score") else None
    competing = [int(item["score"]) for item in candidates if item["call_id"] != spike["spike_id"]]
    return {
        "spike_id": spike["spike_id"],
        "gene": spike["gene"],
        "chrom": spike["chrom"],
        "pos": spike["pos"],
        "end": spike["end"],
        "svlen": spike["svlen"],
        "svtype": spike["svtype"],
        "clnsig": spike["clnsig"],
        "clinvar_accession": spike["clinvar_accession"],
        "imprecise": int(bool(spike["imprecise"])),
        "reached_candidate_table": int(rank is not None),
        "rank": rank if rank is not None else "",
        "candidate_rows": len(candidates),
        "score": score if score is not None else "",
        "score_components": row.get("score_components", ""),
        "gene_effect": row.get("gene_effect", ""),
        "gnomad_popmax_af": row.get("gnomad_popmax_af", ""),
        "clears_candidate_floor": int(score is not None and score >= scorer.CANDIDATE_SCORE_FLOOR),
        "clears_submission_floor": int(score is not None and score >= scorer.SUBMISSION_SCORE_FLOOR),
        "in_row_budget": int(in_row_budget(outputs["submission"], spike["chrom"], spike["pos"])),
        "submission_rows": len(outputs["submission"]),
        "drop_reason": "" if rank is not None else (drop_reason(row) if row else "call_not_read"),
        "known_status": row.get("known_status", ""),
        "background_calls_in_gene": len(competing),
        "top_background_score": max(competing) if competing else "",
    }


# --- Step 5. The summary and the record. --------------------------------------------------------


def summarise(results: list[dict]) -> dict:
    """Compute the reported fractions with their Wilson intervals.

    Args:
        results: The per-spike result rows.

    Returns:
        The three fractions and the supporting counts.
    """
    total = len(results)
    rank_one = sum(1 for row in results if row["rank"] == 1)
    budget = sum(1 for row in results if row["in_row_budget"] == 1)
    dropped = sum(1 for row in results if row["reached_candidate_table"] == 0)
    scored = [row for row in results if row["score"] != ""]
    ranks = collections.Counter(row["rank"] for row in results if row["rank"] != "")
    competing = [row for row in results if row["background_calls_in_gene"] > 0]
    margins = sorted(row["score"] - row["top_background_score"] for row in competing
                     if row["score"] != "" and row["top_background_score"] != "")
    return {
        "spikes": total,
        "rank_1_fraction": stats.wilson(rank_one, total),
        "in_row_budget_fraction": stats.wilson(budget, total),
        "dropped_before_scoring_fraction": stats.wilson(dropped, total),
        "clears_candidate_floor_fraction": stats.wilson(
            sum(1 for row in results if row["clears_candidate_floor"] == 1), total),
        "clears_submission_floor_fraction": stats.wilson(
            sum(1 for row in results if row["clears_submission_floor"] == 1), total),
        "rank_distribution": {str(key): value for key, value in sorted(ranks.items())},
        "score_distribution": {str(key): value for key, value in sorted(
            collections.Counter(row["score"] for row in scored).items())},
        "drop_reasons": dict(collections.Counter(
            row["drop_reason"] for row in results if row["drop_reason"])),
        "mean_candidate_rows": round(
            sum(row["candidate_rows"] for row in results) / total, 3) if total else None,
        "svtype_counts": dict(collections.Counter(row["svtype"] for row in results)),
        "gene_effect_counts": dict(collections.Counter(row["gene_effect"] for row in results)),
        "clnsig_counts": dict(collections.Counter(row["clnsig"] for row in results)),
        "known_status_counts": dict(collections.Counter(row["known_status"] for row in results)),
        "imprecise_spikes": sum(row["imprecise"] for row in results),
        "runs_with_a_background_candidate": len(competing),
        "background_candidate_rows_total": sum(row["background_calls_in_gene"] for row in results),
        "score_margin_over_top_background_call": {
            "runs": len(margins),
            "minimum": margins[0] if margins else None,
            "median": statistics.median(margins) if margins else None,
            "maximum": margins[-1] if margins else None,
            "at_or_below_zero": sum(1 for value in margins if value <= 0),
        },
    }


def gate_verdicts(summary: dict, negative: dict, default_genes: list[str]) -> dict:
    """Apply every gate condition to the result.

    Args:
        summary: The summary fractions.
        negative: The negative control result.
        default_genes: The `--genes` default of the scorer, as a list.

    Returns:
        One verdict for each gate, with the observed value.
    """
    verdicts = {}
    rank_one = summary["rank_1_fraction"]
    budget = summary["in_row_budget_fraction"]
    dropped = summary["dropped_before_scoring_fraction"]
    verdicts["g1_benchmark_has_power"] = {
        "observed": summary["spikes"],
        "verdict": "PASS" if summary["spikes"] >= GATES["g1_benchmark_has_power"]["floor"] else "FAIL",
    }
    verdicts["g2_default_gene_set_unchanged"] = {
        "observed": default_genes,
        "verdict": "PASS" if default_genes == list(scorer.MVA_GENES) and default_genes[0] == "BUB1B"
                   else "FAIL",
    }
    verdicts["g3_rank_one_fraction"] = {
        "observed": rank_one["estimate"], "interval": [rank_one["low"], rank_one["high"]],
        "verdict": "PASS" if (rank_one["estimate"] is not None
                              and rank_one["estimate"] >= GATES["g3_rank_one_fraction"]["floor"]
                              and rank_one["low"] >= GATES["g3_rank_one_fraction"]["interval_floor"])
                   else "FAIL",
    }
    verdicts["g4_budget_fraction"] = {
        "observed": budget["estimate"], "interval": [budget["low"], budget["high"]],
        "verdict": "PASS" if (budget["estimate"] is not None
                              and budget["estimate"] >= GATES["g4_budget_fraction"]["floor"]
                              and budget["low"] >= GATES["g4_budget_fraction"]["interval_floor"])
                   else "FAIL",
    }
    verdicts["g5_pre_scoring_drop"] = {
        "observed": dropped["estimate"], "interval": [dropped["low"], dropped["high"]],
        "verdict": "PASS" if (dropped["estimate"] is not None
                              and dropped["estimate"] <= GATES["g5_pre_scoring_drop"]["ceiling"])
                   else "FAIL",
    }
    verdicts["g6_negative_control_clean"] = {
        "observed": negative["default_gene_set"]["submission_rows"],
        "verdict": "PASS" if negative["default_gene_set"]["submission_rows"] == 0 else "FAIL",
    }
    return verdicts


def negative_control(background: list[str], paths: dict, genes: list[str], work: Path) -> dict:
    """Run the scorer on the background call set with no spike.

    Args:
        background: The background data lines.
        paths: The input paths of the benchmark.
        genes: The benchmark gene names.
        work: The working directory.

    Returns:
        The result of the default gene set arm and of the per-gene arm.
    """
    sv_vcf = work / "negative.vcf"
    sv_vcf.write_text("\n".join(combined_lines(background, None)) + "\n")
    out_dir = work / "negative_default"
    outputs = run_scorer(sv_vcf, paths["mei"], paths["gtf"], paths["gnomad"], paths["intervals"],
                         out_dir, ",".join(scorer.MVA_GENES), "negative_default")
    default_arm = {
        "genes": list(scorer.MVA_GENES),
        "candidate_rows": len(outputs["candidates"]),
        "submission_rows": len(outputs["submission"]),
        "top_score": max((int(row["score"]) for row in outputs["candidates"]), default=None),
        "scores": sorted((int(row["score"]) for row in outputs["candidates"]), reverse=True),
        "rows": [{"rank": int(row["rank"]), "score": int(row["score"]), "gene": row["mva_gene"],
                  "gene_effect": row["gene_effect"], "svtype": row["svtype"],
                  "svlen": row["svlen"], "genotype": row["genotype"],
                  "gnomad_popmax_af": row["gnomad_popmax_af"],
                  "score_components": row["score_components"]}
                 for row in outputs["candidates"]],
    }
    shutil.rmtree(out_dir, ignore_errors=True)
    per_gene = []
    for position, gene in enumerate(genes, start=1):
        out_dir = work / f"negative_{position}"
        outputs = run_scorer(sv_vcf, paths["mei"], paths["gtf"], paths["gnomad"],
                             paths["intervals"], out_dir, gene, f"negative_{gene}")
        per_gene.append({
            "gene": gene,
            "candidate_rows": len(outputs["candidates"]),
            "submission_rows": len(outputs["submission"]),
            "top_score": max((int(row["score"]) for row in outputs["candidates"]), default=None),
        })
        shutil.rmtree(out_dir, ignore_errors=True)
        if position % 25 == 0:
            print(f"  negative control {position} of {len(genes)}")
    sv_vcf.unlink(missing_ok=True)
    with_rows = [item for item in per_gene if item["candidate_rows"] > 0]
    return {
        "default_gene_set": default_arm,
        "per_gene": {
            "genes_tested": len(per_gene),
            "genes_with_at_least_one_candidate": len(with_rows),
            "genes_with_a_submission_row": sum(1 for item in per_gene if item["submission_rows"]),
            "candidate_rows_total": sum(item["candidate_rows"] for item in per_gene),
            "candidate_row_distribution": {str(key): value for key, value in sorted(
                collections.Counter(item["candidate_rows"] for item in per_gene).items())},
            "top_score_distribution": {str(key): value for key, value in sorted(
                collections.Counter(item["top_score"] for item in with_rows).items(),
                key=lambda item: str(item[0]))},
            "submission_rows_total": sum(item["submission_rows"] for item in per_gene),
        },
        "per_gene_rows": per_gene,
    }


def write_results_table(path: Path, results: list[dict]) -> None:
    """Write the per-spike result table.

    Args:
        path: The output path.
        results: The result rows.
    """
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in results:
            writer.writerow(row)


def step_fetch() -> dict:
    """Run the fetch step.

    Returns:
        The source manifest.
    """
    manifest = fetch_sources()
    print(json.dumps({key: value["sha256"] for key, value in manifest.items()
                      if "sha256" in value}, indent=2))
    return manifest


def step_select() -> dict:
    """Run the selection step and build every benchmark input.

    Returns:
        The selection counts, the input checksums and the transfer counts.
    """
    genes = read_protein_coding_genes(CACHE / SOURCES["gencode"]["name"])
    unambiguous, repeated = unique_gene_names(genes)
    records, counts = select_pathogenic(CACHE / SOURCES["dbvar_nstd102"]["name"],
                                        gene_span_index(genes), unambiguous)
    counts["protein_coding_genes_in_gencode_50"] = sum(len(v) for v in genes.values())
    counts["gene_names_given_more_than_once"] = repeated
    spikes = draw_spikes(records)
    counts["spikes_drawn"] = len(spikes)
    (CACHE / "spikes.json").write_text(json.dumps(spikes, indent=2) + "\n")
    (CACHE / "selection_counts.json").write_text(json.dumps(counts, indent=2) + "\n")
    print(json.dumps(counts, indent=2))

    # The seven MVA genes join the subset so that the negative control can run with the shipped
    # default gene set. They carry no spike.
    wanted = {spike["gene"] for spike in spikes} | set(scorer.MVA_GENES)
    gtf = CACHE / "gencode_v50_benchmark_genes.gtf.gz"
    gtf_lines = write_gene_subset(CACHE / SOURCES["gencode"]["name"], gtf, wanted)
    regions = gene_regions(gtf, wanted)
    print(f"gene subset lines {gtf_lines}, regions {len(regions)}")

    index = read_tabix_index(CACHE / SOURCES["gnomad_sv_index"]["name"])
    lines, transfer = tabix_query(GNOMAD_SV_URL, index, regions)
    gnomad = CACHE / "gnomad_v4_1_sv_benchmark_loci.vcf.gz"
    transfer["records_written"] = write_gnomad_subset(gnomad, lines)
    print(json.dumps(transfer, indent=2))

    background = CACHE / "background_hg002_grch38.vcf"
    lift = lift_background(CACHE / SOURCES["giab_v06"]["name"],
                           CACHE / SOURCES["chain"]["name"], background)
    print(json.dumps(lift, indent=2))

    report = {
        "selection_counts": counts,
        "gene_subset": {"name": gtf.name, "lines": gtf_lines, "genes": len(wanted),
                        "sha256": sha256_of(gtf),
                        "rule": ("The subset holds the gene and exon lines of the spike genes and of "
                                 "the seven MVA genes. The seven MVA genes carry no spike and they "
                                 "are present so that the negative control can run with the shipped "
                                 "default gene set. The subset holds no other gene, so the "
                                 "informational `genes_overlapped` column of the scorer names fewer "
                                 "genes than a genome wide annotation would name. That column takes "
                                 "no part in the score and no part in the rank.")},
        "gnomad_subset": {"name": gnomad.name, "sha256": sha256_of(gnomad), **transfer,
                          "kept_info_keys": list(GNOMAD_KEPT_KEYS),
                          "rule": (f"A tabix range query over https reads every gnomAD SV v4.1 site "
                                   f"that overlaps a benchmark gene span padded by {GNOMAD_PAD} base "
                                   f"pairs, which is the regulatory window of the scorer. The INFO "
                                   f"column keeps the keys that annotate_sv_mei.read_gnomad_sv "
                                   f"reads and drops every other key.")},
        "background": {"name": background.name, "sha256": sha256_of(background),
                       "lift_counts": lift},
    }
    (CACHE / "input_report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def step_run() -> dict:
    """Run the scorer on every spike and on the negative control, then write the record.

    Returns:
        The summary, the negative control and the gate verdicts.
    """
    spikes = json.loads((CACHE / "spikes.json").read_text())
    input_report = json.loads((CACHE / "input_report.json").read_text())
    paths = {
        "gtf": CACHE / "gencode_v50_benchmark_genes.gtf.gz",
        "gnomad": CACHE / "gnomad_v4_1_sv_benchmark_loci.vcf.gz",
        "intervals": KNOWN_INTERVALS,
        "mei": CACHE / "empty_mei.vcf",
    }
    paths["mei"].write_text("\n".join(COMBINED_HEADER) + "\n")
    background = [line for line in (CACHE / "background_hg002_grch38.vcf").read_text().splitlines()
                  if line and not line.startswith("#")]
    work = CACHE / "runs"
    work.mkdir(parents=True, exist_ok=True)

    results = []
    for position, spike in enumerate(spikes, start=1):
        sv_vcf = work / "combined.vcf"
        sv_vcf.write_text("\n".join(combined_lines(background, spike)) + "\n")
        out_dir = work / "out"
        outputs = run_scorer(sv_vcf, paths["mei"], paths["gtf"], paths["gnomad"],
                             paths["intervals"], out_dir, spike["gene"],
                             f"spike_{spike['spike_id']}")
        results.append(measure_spike(spike, outputs))
        shutil.rmtree(out_dir, ignore_errors=True)
        sv_vcf.unlink(missing_ok=True)
        if position % 25 == 0:
            print(f"  spike {position} of {len(spikes)}")

    summary = summarise(results)
    negative = negative_control(background, paths, sorted({s["gene"] for s in spikes}), work)
    shutil.rmtree(work, ignore_errors=True)
    paths["mei"].unlink(missing_ok=True)

    default_genes = scorer_default_genes()
    verdicts = gate_verdicts(summary, negative, default_genes)

    write_results_table(RESULTS_TABLE, results)
    manifest = json.loads((CACHE / "source_manifest.json").read_text())
    receipt = {
        "receipt_id": RECEIPT_ID,
        "created_utc": utc_now(),
        "track": "t1_call",
        "data_class": "A",
        "defect": "R2-2 of reports/jury/defect_ledger_track1_round2_2026-09-17.md",
        "question": ("Does a known pathogenic structural variant that falls inside a protein coding "
                     "gene reach rank 1 of the candidate table and a row of the seven row "
                     "submission budget when the shipped scorer reads it inside a public background "
                     "call set?"),
        "object_under_test": {
            "path": "pipelines/sv_mei_annotate/bin/annotate_sv_mei.py",
            "point_terms": [50, 25, 20, 15, 12, 10, 8, 7],
            "candidate_score_floor": scorer.CANDIDATE_SCORE_FLOOR,
            "submission_score_floor": scorer.SUBMISSION_SCORE_FLOOR,
            "max_submission_rows": scorer.MAX_SUBMISSION_ROWS,
            "max_candidate_rows": scorer.MAX_CANDIDATE_ROWS,
            "default_genes": default_genes,
            "change_in_this_pass": ("A `--genes` option replaces the module constant MVA_GENES. Its "
                                    "default is that tuple, so a run that omits the option behaves "
                                    "as the shipped run behaves."),
        },
        "tools": {
            "python": platform.python_version(),
            "httpx": httpx.__version__,
            "liftover": "pyliftover 0.4.1",
            "tabix": ("No htslib. The script reads the published tabix index and the BGZF blocks "
                      "with the Python standard library and HTTP range requests."),
        },
        "parameters": {
            "spike_count": SPIKE_COUNT,
            "seed": SEED,
            "spike_genotype": SPIKE_GENOTYPE,
            "spike_genotype_quality": SPIKE_GENOTYPE_QUALITY,
            "min_spike_length_bp": MIN_SPIKE_LENGTH,
            "pathogenic_terms": list(PATHOGENIC_TERMS),
            "gnomad_region_pad_bp": GNOMAD_PAD,
            "spike_arm": ("Every spike enters the structural variant arm. No spike enters the "
                          "mobile element arm, so no spike earns the intronic mobile element term "
                          "of 20 points and every spike meets the 50 base pair floor of that arm."),
            "min_frequency_records": 1,
            "scorer_command": ("annotate_sv_mei.py --sv-vcf <background plus one spike> --mei-vcf "
                               "<header only> --gtf gencode_v50_benchmark_genes.gtf.gz --gnomad-sv "
                               "gnomad_v4_1_sv_benchmark_loci.vcf.gz --intervals intervals.bed "
                               "--out-dir <run> --genes <gene of the spike> "
                               "--min-frequency-records 1"),
        },
        "inputs": manifest,
        "derived_inputs": input_report,
        "gate_conditions": GATES,
        "gates": verdicts,
        "results": summary,
        "negative_control": negative,
        "outputs": {"results_table": str(RESULTS_TABLE.relative_to(ROOT)),
                    "results_table_sha256": sha256_of(RESULTS_TABLE),
                    "results_table_rows": len(results),
                    "record": "docs/benchmarks/t1_sv_causal_benchmark_v1.md",
                    "test": "harness/test_track1_sv_causal_benchmark.py"},
        "limits": [
            "The benchmark measures scoring and ranking. It does not measure calling from reads, so "
            "a causal event that a caller misses is outside it.",
            "Every spike is placed at its database coordinate with a heterozygous genotype and a "
            "genotype quality of 99. The genotype and the quality are choices of the authors.",
            "Every spike falls whole inside one protein coding gene, so the benchmark measures no "
            "event that spans a gene boundary and no break end.",
            "The gene of the spike is given to the scorer. The benchmark therefore measures the "
            "ranking and not gene selection.",
            "The background is one sample, GIAB HG002. A second sample was not tested.",
            "The scorer reads the nine AF_<population> keys of its GNOMAD_POPULATIONS tuple and it "
            "does not read AF_rmi, so its population maximum frequency ignores the remaining "
            "individuals group of gnomAD v4.1.",
        ],
        "verdict": "PASS" if all(item["verdict"] == "PASS" for item in verdicts.values()) else "FAIL",
    }
    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"gates": verdicts, "results": summary,
                      "negative_control": {"default_gene_set": negative["default_gene_set"],
                                           "per_gene": negative["per_gene"]}}, indent=2))
    return receipt


def main() -> int:
    """Run the requested step.

    Returns:
        Zero when every gate passed, one otherwise.
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("step", choices=("fetch", "select", "run", "all"), default="all", nargs="?")
    args = parser.parse_args()
    if args.step in ("fetch", "all"):
        step_fetch()
    if args.step in ("select", "all"):
        step_select()
    if args.step in ("run", "all"):
        receipt = step_run()
        return 0 if receipt["verdict"] == "PASS" else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
