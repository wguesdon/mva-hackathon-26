"""Resolve known MVA alleles from public literature to GRCh38 coordinates with deterministic gates (PRD S3).

Input: literature/known_alleles_seed.yaml. All content is public, data class A.

Gates for each allele:
    quote      The quote must occur in the literature pack abstract, the checked full-text quote file,
               or the PubMed abstract, after whitespace normalization.
    agreement  HGVS alleles: Ensembl VEP REST and VariantValidator must describe the same GRCh38 haplotype.
               VEP shifts indels right and VariantValidator shifts them left, so the gate applies both records
               to a reference window and compares the edited sequences. The row keeps the left-aligned record.
    sequence   A duplication with expected_inserted_sequence must insert a rotation of that sequence.
               rsIDs: NCBI dbSNP (SPDI) and Ensembl VEP must return the same GRCh38 VCF record.
    reference  A resolver warning about a reference mismatch fails the allele.
    insertion  A target site duplication must occur once in the last 300 bp of the named intron (MANE exon order).

Status values: resolved, interval_resolved, intron_level, linked, sequence_mismatch, unresolved, conflict.

Outputs:
    reports/known_alleles/known_alleles_grch38_<date>.tsv
    reports/known_alleles/clinvar_plp_mva_genes_<date>.tsv
    tracking/receipts/s3_known_alleles_<date>.json

Usage:
    uv run python scripts/resolve_known_alleles.py
"""

from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import re
import sys
import time
import urllib.parse
from pathlib import Path

import httpx
import yaml

ROOT = Path(__file__).resolve().parents[1]
SEED = ROOT / "literature/known_alleles_seed.yaml"
MANE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.summary.txt.gz"
VEP = "https://rest.ensembl.org"
VV = "https://rest.variantvalidator.org/VariantValidator/variantvalidator/GRCh38"
REFSNP = "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CLINVAR_GENES = ("BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD2L1BP", "MAD1L1", "BUB1", "SMC5")
FIELDS = ("allele_id", "gene", "kind", "variant", "transcript", "context", "pmid", "quote_gate", "status",
          "chrom", "pos", "ref", "alt", "interval_start", "interval_end", "vep_vcf", "vv_vcf", "dbsnp_vcf", "note")


def norm(text: str) -> str:
    """Collapse whitespace for quote matching.

    Args:
        text: Any text.

    Returns:
        The text with single spaces.
    """
    return " ".join(text.split())


def get_json(client: httpx.Client, url: str, pause: float = 0.4) -> dict | list | None:
    """GET a JSON resource with a short pause and one retry.

    Args:
        client: The HTTP client.
        url: The URL.
        pause: Seconds to wait after the call.

    Returns:
        The decoded JSON, or None on failure.
    """
    for _ in range(2):
        try:
            response = client.get(url, headers={"Accept": "application/json", "Content-Type": "application/json"})
            time.sleep(pause)
            if response.status_code == 200:
                return response.json()
            if response.status_code in (400, 404):
                try:
                    return response.json()
                except ValueError:
                    return None
        except httpx.HTTPError:
            time.sleep(2)
    return None


def mane_table(client: httpx.Client) -> dict[str, dict]:
    """Read the MANE Select transcripts by gene symbol.

    Args:
        client: The HTTP client.

    Returns:
        Gene symbol to RefSeq and Ensembl transcript identifiers.
    """
    raw = gzip.decompress(client.get(MANE_URL).content).decode()
    table = {}
    for row in csv.DictReader(io.StringIO(raw), delimiter="\t"):
        if row["MANE_status"] == "MANE Select":
            table[row["symbol"]] = {"refseq": row["RefSeq_nuc"], "ensembl": row["Ensembl_nuc"].split(".")[0]}
    return table


def pack_text(pack: Path) -> str:
    """Join all abstract files of the literature pack.

    Args:
        pack: The pack directory.

    Returns:
        The normalized text.
    """
    return norm(" ".join(path.read_text() for path in sorted((pack / "abstracts").glob("*.md"))))


def quote_gate(client: httpx.Client, allele: dict, abstracts: str, fulltext: dict[str, bool]) -> str:
    """Check the quote of one allele against its declared source.

    Args:
        client: The HTTP client.
        allele: The seed entry.
        abstracts: The normalized pack abstracts.
        fulltext: Quote to found flag from the full-text quote file.

    Returns:
        PASS or FAIL with the source.
    """
    quote = norm(allele.get("quote") or "")
    source = allele.get("quote_source")
    if not quote:
        return "FAIL:no_quote"
    if source == "abstract":
        return "PASS:abstract" if quote in abstracts else "FAIL:abstract"
    if source == "fulltext":
        return "PASS:fulltext" if fulltext.get(allele["quote"]) else "FAIL:fulltext"
    if source == "pubmed":
        text = client.get(f"{EUTILS}/efetch.fcgi", params={"db": "pubmed", "id": allele["pmid"],
                                                           "rettype": "abstract", "retmode": "text"}).text
        time.sleep(0.4)
        return "PASS:pubmed" if quote in norm(text) else "FAIL:pubmed"
    return "FAIL:unknown_source"


def vep_hgvs(client: httpx.Client, hgvs: str) -> tuple[str | None, str]:
    """Resolve an HGVS string with Ensembl VEP REST.

    Args:
        client: The HTTP client.
        hgvs: The transcript HGVS.

    Returns:
        The VCF string chrom-pos-ref-alt, or None, and a note.
    """
    data = get_json(client, f"{VEP}/vep/human/hgvs/{urllib.parse.quote(hgvs, safe=':')}?refseq=1&vcf_string=1", 0.2)
    if isinstance(data, list) and data and data[0].get("vcf_string"):
        return data[0]["vcf_string"], ""
    message = data.get("error", "no result") if isinstance(data, dict) else "no result"
    return None, f"VEP: {str(message)[:160]}"


def vv_hgvs(client: httpx.Client, hgvs: str) -> tuple[str | None, str]:
    """Resolve an HGVS string with VariantValidator.

    Args:
        client: The HTTP client.
        hgvs: The transcript HGVS.

    Returns:
        The VCF string chrom-pos-ref-alt, or None, and a note with validation warnings.
    """
    data = get_json(client, f"{VV}/{urllib.parse.quote(hgvs, safe='')}/mane_select", 1.0)
    if not isinstance(data, dict):
        return None, "VariantValidator: no result"
    notes, vcf = [], None
    for key, value in data.items():
        if not isinstance(value, dict):
            continue
        notes.extend(value.get("validation_warnings") or [])
        loci = (value.get("primary_assembly_loci") or {}).get("grch38") or {}
        record = loci.get("vcf")
        if record and vcf is None:
            vcf = f"{record['chr'].replace('chr', '')}-{record['pos']}-{record['ref']}-{record['alt']}"
    return vcf, ("VariantValidator: " + " | ".join(notes)[:200]) if notes else ""


def dbsnp_rs(client: httpx.Client, rsid: str) -> tuple[str | None, str]:
    """Resolve an rsID to a GRCh38 VCF record through NCBI dbSNP SPDI placements.

    Args:
        client: The HTTP client.
        rsid: The rsID.

    Returns:
        The VCF string of the first alternate allele, or None, and a note.
    """
    data = get_json(client, f"{REFSNP}/{rsid.removeprefix('rs')}", 0.4)
    if not isinstance(data, dict):
        return None, "dbSNP: no result"
    for placement in data.get("primary_snapshot_data", {}).get("placements_with_allele", []):
        traits = placement.get("placement_annot", {}).get("seq_id_traits_by_assembly", [])
        if not placement.get("is_ptlp") or not traits or not traits[0].get("assembly_name", "").startswith("GRCh38"):
            continue
        chrom = str(int(placement["seq_id"].split(".")[0].removeprefix("NC_")))
        for allele in placement["alleles"]:
            spdi = allele["allele"]["spdi"]
            if spdi["deleted_sequence"] != spdi["inserted_sequence"] and len(spdi["deleted_sequence"]) == 1 \
                    and len(spdi["inserted_sequence"]) == 1:
                return f"{chrom}-{spdi['position'] + 1}-{spdi['deleted_sequence']}-{spdi['inserted_sequence']}", ""
    return None, "dbSNP: no GRCh38 SNV placement"


def vep_rs(client: httpx.Client, rsid: str) -> str | None:
    """Resolve an rsID with Ensembl VEP REST.

    Args:
        client: The HTTP client.
        rsid: The rsID.

    Returns:
        The VCF string, or None.
    """
    data = get_json(client, f"{VEP}/vep/human/id/{rsid}?vcf_string=1", 0.2)
    return data[0].get("vcf_string") if isinstance(data, list) and data else None


def exon_start(exon: dict) -> int:
    """Sort key for Ensembl exon records.

    Args:
        exon: An exon record.

    Returns:
        The start coordinate.
    """
    return exon["start"]


def transcript_exons(client: httpx.Client, ensembl_id: str) -> tuple[str, int, list[dict]]:
    """Read the exons of a transcript in transcript order.

    Args:
        client: The HTTP client.
        ensembl_id: The Ensembl transcript identifier without version.

    Returns:
        The chromosome, the strand and the ordered exon records with start and end.
    """
    data = get_json(client, f"{VEP}/lookup/id/{ensembl_id}?expand=1", 0.2)
    exons = sorted(data["Exon"], key=exon_start, reverse=data["strand"] == -1)
    return data["seq_region_name"], data["strand"], exons


def intron_bounds(exons: list[dict], strand: int, number: int) -> tuple[int, int]:
    """Return the genomic bounds of an intron between exon N and exon N+1.

    Args:
        exons: The exons in transcript order.
        strand: 1 or -1.
        number: The intron number, 1-based.

    Returns:
        The lowest and highest intron coordinates.
    """
    left, right = exons[number - 1], exons[number]
    if strand == 1:
        return left["end"] + 1, right["start"] - 1
    return right["end"] + 1, left["start"] - 1


def region_sequence(client: httpx.Client, chrom: str, start: int, end: int, strand: int) -> str:
    """Fetch a GRCh38 sequence in the transcript orientation.

    Args:
        client: The HTTP client.
        chrom: The chromosome.
        start: The first coordinate.
        end: The last coordinate.
        strand: 1 or -1.

    Returns:
        The upper-case sequence.
    """
    data = get_json(client, f"{VEP}/sequence/region/human/{chrom}:{start}..{end}:{strand}?coord_system_version=GRCh38", 0.2)
    return data["seq"].upper()


def apply_vcf(window: str, window_start: int, vcf: str) -> str | None:
    """Apply one VCF record to a reference window.

    Args:
        window: The plus-strand reference sequence.
        window_start: The coordinate of the first window base.
        vcf: The chrom-pos-ref-alt string.

    Returns:
        The edited sequence, or None when the record does not fit the window or its REF differs.
    """
    _, pos, ref, alt = vcf.split("-", 3)
    offset = int(pos) - window_start
    if offset < 0 or offset + len(ref) > len(window) or window[offset:offset + len(ref)] != ref:
        return None
    return window[:offset] + alt + window[offset + len(ref):]


def same_haplotype(client: httpx.Client, first: str, second: str) -> bool:
    """Check that two VCF records give the same sequence on a shared reference window.

    Args:
        client: The HTTP client.
        first: A chrom-pos-ref-alt string.
        second: A chrom-pos-ref-alt string.

    Returns:
        True when both records fit the reference and produce the same edited sequence.
    """
    chrom_a, pos_a = first.split("-")[:2]
    chrom_b, pos_b = second.split("-")[:2]
    if chrom_a != chrom_b:
        return False
    start = min(int(pos_a), int(pos_b)) - 30
    end = max(int(pos_a), int(pos_b)) + 60
    window = region_sequence(client, chrom_a, start, end, 1)
    edited_a, edited_b = apply_vcf(window, start, first), apply_vcf(window, start, second)
    return edited_a is not None and edited_a == edited_b


def is_rotation(inserted: str, expected: str) -> bool:
    """Check that an inserted sequence is a cyclic rotation of the expected sequence.

    Args:
        inserted: The inserted bases.
        expected: The expected duplicated bases.

    Returns:
        True for a rotation of the same length.
    """
    return len(inserted) == len(expected) and inserted.upper() in (expected.upper() * 2)


def vcf_parts(vcf: str | None) -> tuple[str, str, str, str]:
    """Split a chrom-pos-ref-alt string.

    Args:
        vcf: The VCF string or None.

    Returns:
        Four fields, empty when absent.
    """
    return tuple(vcf.split("-", 3)) if vcf else ("", "", "", "")


def resolve(client: httpx.Client, allele: dict, mane: dict, resolved_rs: dict) -> dict:
    """Resolve one seed allele.

    Args:
        client: The HTTP client.
        allele: The seed entry.
        mane: The MANE Select table.
        resolved_rs: rsID rows already resolved in this run.

    Returns:
        The output row without the quote gate.
    """
    row = {field: "" for field in FIELDS}
    row.update(allele_id=allele["id"], gene=allele["gene"], kind=allele["kind"], variant=allele["variant"],
               context=allele["context"], pmid=allele["pmid"], note=allele.get("note", ""))
    tx = allele.get("transcript") or mane[allele["gene"]]["refseq"]
    kind = allele["kind"]
    if kind == "hgvs_c":
        hgvs = f"{tx}:{allele['variant']}"
        row["transcript"] = tx
        vep, vep_note = vep_hgvs(client, hgvs)
        mane_tx = mane[allele["gene"]]["refseq"]
        if vep is None and "Could not get a Transcript object" in vep_note and tx != mane_tx:
            vep, fallback_note = vep_hgvs(client, f"{mane_tx}:{allele['variant']}")
            vep_note = f"VEP lacks {tx}; VEP used {mane_tx}. {fallback_note}".strip()
        vv, vv_note = vv_hgvs(client, hgvs)
        row.update(vep_vcf=vep or "", vv_vcf=vv or "")
        notes = [note for note in (vep_note, vv_note, row["note"]) if note]
        mismatch = re.search(r"does not (match|agree)|Variant reference", " ".join((vep_note, vv_note)), re.IGNORECASE)
        if vep and vv and not mismatch and (vep == vv or same_haplotype(client, vep, vv)):
            row["status"] = "resolved"
            row["chrom"], row["pos"], row["ref"], row["alt"] = vcf_parts(vv)
            expected = allele.get("expected_inserted_sequence")
            inserted = row["alt"][len(row["ref"]):] if len(row["alt"]) > len(row["ref"]) else ""
            if expected and not is_rotation(inserted, expected):
                row["status"] = "sequence_mismatch"
                notes.insert(0, f"The inserted sequence {inserted} is not a rotation of the published {expected}; "
                                f"the published numbering does not fit {tx}")
        elif vep and vv:
            row["status"] = "conflict"
        else:
            row["status"] = "unresolved"
        row["note"] = " ; ".join(notes)
    elif kind == "rsid":
        dbsnp, db_note = dbsnp_rs(client, allele["variant"])
        vep = vep_rs(client, allele["variant"])
        row.update(dbsnp_vcf=dbsnp or "", vep_vcf=vep or "")
        if dbsnp and vep and dbsnp == vep:
            row["status"] = "resolved"
            row["chrom"], row["pos"], row["ref"], row["alt"] = vcf_parts(vep)
        else:
            row["status"] = "conflict" if dbsnp and vep else "unresolved"
        row["note"] = " ; ".join(note for note in (db_note, row["note"]) if note)
        resolved_rs[allele["variant"]] = row
    elif kind == "ssid":
        row["status"] = "linked"
        row["note"] = ("No public API resolves ss identifiers. PMID 28611924 full text links the 44 kb upstream "
                       "G>A allele to rs576524605; use that row. " + row["note"]).strip()
    elif kind in ("insertion_tsd", "insertion_interval"):
        ensembl = mane[allele["gene"]]["ensembl"]
        row["transcript"] = mane[allele["gene"]]["refseq"]
        chrom, strand, exons = transcript_exons(client, ensembl)
        low, high = intron_bounds(exons, strand, int(allele["intron"]))
        row["chrom"] = chrom
        if kind == "insertion_interval":
            row.update(status="intron_level", interval_start=str(low), interval_end=str(high))
            row["note"] = f"MANE intron {allele['intron']} spans {high - low + 1} bp. The paper gives no position. " + row["note"]
        else:
            if strand == 1:
                window_start, window_end = max(low, high - 299), high
            else:
                window_start, window_end = low, min(high, low + 299)
            sequence = region_sequence(client, chrom, window_start, window_end, strand)
            hits = [match.start() for match in re.finditer(f"(?={allele['tsd']})", sequence)]
            if len(hits) == 1:
                offset = hits[0]
                length = len(allele["tsd"])
                if strand == 1:
                    first, last = window_start + offset, window_start + offset + length - 1
                else:
                    last = window_end - offset
                    first = last - length + 1
                distance = (high - last) if strand == 1 else (first - low)
                row.update(status="interval_resolved", interval_start=str(first), interval_end=str(last))
                row["note"] = (f"The target site duplication occurs once in the last 300 bp of MANE intron {allele['intron']}, "
                               f"{distance} bp from the acceptor end. " + row["note"])
            else:
                row.update(status="unresolved", interval_start=str(low), interval_end=str(high))
                row["note"] = f"The target site duplication occurs {len(hits)} times in the last 300 bp of the intron. " + row["note"]
    return row


def clinvar_plp(client: httpx.Client, gene: str) -> list[dict]:
    """Read ClinVar pathogenic and likely pathogenic records of one gene with GRCh38 locations.

    Args:
        client: The HTTP client.
        gene: The gene symbol.

    Returns:
        One row per record.
    """
    term = f"{gene}[gene] AND (clinsig_pathogenic[prop] OR clinsig_likely_pathogenic[prop])"
    ids = client.get(f"{EUTILS}/esearch.fcgi", params={"db": "clinvar", "term": term, "retmode": "json",
                                                       "retmax": 1000}).json()["esearchresult"]["idlist"]
    time.sleep(0.4)
    rows = []
    for start in range(0, len(ids), 100):
        batch = ",".join(ids[start:start + 100])
        data = client.get(f"{EUTILS}/esummary.fcgi", params={"db": "clinvar", "id": batch, "retmode": "json"}).json()
        time.sleep(0.4)
        for uid in data["result"].get("uids", []):
            record = data["result"][uid]
            genes = {item.get("symbol") for item in record.get("genes", [])}
            if gene not in genes:
                continue
            variation = (record.get("variation_set") or [{}])[0]
            loc = next((item for item in variation.get("variation_loc", []) if item.get("assembly_name") == "GRCh38"), {})
            classification = (record.get("germline_classification") or {}).get("description", "")
            traits = "; ".join(trait.get("trait_name", "") for trait in (record.get("germline_classification") or {}).get("trait_set", []))
            rows.append({"clinvar_uid": uid, "gene": gene, "title": record.get("title", ""),
                         "classification": classification, "review_status": (record.get("germline_classification") or {}).get("review_status", ""),
                         "variant_type": variation.get("variant_type", ""), "chrom": loc.get("chr", ""),
                         "start": loc.get("start", ""), "stop": loc.get("stop", ""),
                         "spdi": variation.get("canonical_spdi", ""), "traits": traits[:300]})
    return rows


def write_tsv(path: Path, rows: list[dict], fields: tuple[str, ...] | list[str]) -> str:
    """Write rows to a TSV and return its SHA-256.

    Args:
        path: The output path.
        rows: The rows.
        fields: The column order.

    Returns:
        The hexadecimal digest.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    """Resolve all seed alleles, pull ClinVar records and write the receipt.

    Returns:
        Zero when every quote gate passes, one otherwise.
    """
    seed = yaml.safe_load(SEED.read_text())
    stamp = dt.datetime.now(dt.UTC).strftime("%Y_%m_%d")
    abstracts = pack_text(ROOT / seed["pack"])
    with (ROOT / seed["fulltext_quotes"]).open(newline="") as handle:
        fulltext = {row["quote"]: row["found_in_europepmc_fulltext"] == "True" for row in csv.DictReader(handle, delimiter="\t")}
    rows, resolved_rs = [], {}
    with httpx.Client(timeout=httpx.Timeout(60, connect=20), headers={"User-Agent": "mva-hackathon-2026"}) as client:
        mane = mane_table(client)
        for allele in seed["alleles"]:
            row = resolve(client, allele, mane, resolved_rs)
            row["quote_gate"] = quote_gate(client, allele, abstracts, fulltext)
            if row["quote_gate"].startswith("FAIL"):
                row["status"] = "unresolved"
            rows.append(row)
            print(f"{row['allele_id']:<34} {row['quote_gate']:<14} {row['status']:<18} "
                  f"{row['chrom']}:{row['pos'] or row['interval_start']} {row['ref']}>{row['alt']}", flush=True)
        clinvar = []
        for gene in CLINVAR_GENES:
            gene_rows = clinvar_plp(client, gene)
            clinvar.extend(gene_rows)
            print(f"clinvar {gene}: {len(gene_rows)} P/LP records", flush=True)
    out = ROOT / "reports/known_alleles"
    allele_path = out / f"known_alleles_grch38_{stamp}.tsv"
    clinvar_path = out / f"clinvar_plp_mva_genes_{stamp}.tsv"
    allele_sha = write_tsv(allele_path, rows, FIELDS)
    clinvar_sha = write_tsv(clinvar_path, clinvar, ["clinvar_uid", "gene", "title", "classification", "review_status",
                                                    "variant_type", "chrom", "start", "stop", "spdi", "traits"])
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    receipt = {"receipt_id": f"s3_known_alleles_{stamp}", "created_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
               "data_class": "A", "seed": str(SEED.relative_to(ROOT)),
               "seed_sha256": hashlib.sha256(SEED.read_bytes()).hexdigest(),
               "mane": MANE_URL, "resolvers": ["Ensembl VEP REST", "VariantValidator REST", "NCBI dbSNP refsnp API", "Ensembl sequence REST"],
               "status_counts": counts, "quote_gate_failures": [row["allele_id"] for row in rows if row["quote_gate"].startswith("FAIL")],
               "outputs": {str(allele_path.relative_to(ROOT)): allele_sha, str(clinvar_path.relative_to(ROOT)): clinvar_sha},
               "clinvar_records": len(clinvar)}
    receipt_path = ROOT / "tracking/receipts" / f"s3_known_alleles_{stamp}.json"
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt["status_counts"]), "quote failures:", receipt["quote_gate_failures"])
    return 1 if receipt["quote_gate_failures"] else 0


if __name__ == "__main__":
    sys.exit(main())
