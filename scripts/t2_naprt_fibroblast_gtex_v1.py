#!/usr/bin/env python3
"""Read the NAPRT expression of cultured fibroblasts from the GTEx Portal. Data class A.

Amendment 5 of ``reports/track2/preclinical_protocol_gate_v1.md`` adds a nicotinic acid arm to
stage 1 of the preclinical protocol. Nicotinic acid reaches the NAD+ pool through nicotinate
phosphoribosyltransferase, NAPRT, on the Preiss-Handler route, and the stage 1 cell type is the
primary human dermal fibroblast. Conditions C18 and C19 of that amendment ask whether the transcript
of *NAPRT* is present in that cell type, and they fix the threshold before this script runs.

The script reads the per sample expression values of *NAPRT* in the GTEx tissue
"Cells - Cultured fibroblasts", computes the median and the quartiles of transcripts per million,
and returns EXPRESSED when the median is 1.0 or more and NOT SHOWN otherwise. A transcript is not an
enzyme activity, so the verdict does not replace the NO ROUTE rule of section 12.3 of the gate.

Usage:
    uv run python scripts/t2_naprt_fibroblast_gtex_v1.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import t2_http_retry  # noqa: E402

EXPRESSION = "https://gtexportal.org/api/v2/expression/geneExpression"
GENE = "https://gtexportal.org/api/v2/reference/gene"
GENE_SYMBOL = "NAPRT"
TISSUE = "Cells_Cultured_fibroblasts"
# GTEx V10 quantifies expression on GENCODE v39. The reference endpoint returns v26 unless it is
# asked for v39, and a v26 identifier returns no V10 expression record.
GENCODE_VERSION = "v39"
GENOME_BUILD = "GRCh38/hg38"
GATE = ROOT / "reports" / "track2" / "preclinical_protocol_gate_v1.md"
OUT = ROOT / "tracking" / "receipts" / "t2_naprt_fibroblast_gtex_v1.json"

# Condition C19 of amendment 5 of the gate, fixed before this script ran.
EXPRESSED_MEDIAN_TPM = 1.0


def gene_record() -> dict:
    """Read the GTEx reference record of the gene.

    Returns:
        The reference fields that the portal returns for NAPRT.

    Raises:
        SystemExit: When the portal returns no record or more than one record for the symbol.
    """
    params = {"geneId": GENE_SYMBOL, "gencodeVersion": GENCODE_VERSION, "genomeBuild": GENOME_BUILD}
    payload = t2_http_retry.get(GENE, params, timeout=60.0).json()
    records = [row for row in payload["data"] if row.get("geneSymbol") == GENE_SYMBOL]
    if len(records) != 1:
        raise SystemExit(f"the portal returned {len(records)} records for {GENE_SYMBOL}")
    return records[0]


def expression_record(gencode_id: str) -> tuple[dict, str]:
    """Read the per sample expression values of the gene in cultured fibroblasts.

    Args:
        gencode_id: The versioned GENCODE identifier that the reference record gives.

    Returns:
        The single expression record of the tissue and the SHA-256 of the raw response body.

    Raises:
        SystemExit: When the endpoint returns a number of records other than one.
    """
    params = {"gencodeId": gencode_id, "tissueSiteDetailId": TISSUE}
    response = t2_http_retry.get(EXPRESSION, params, timeout=90.0)
    digest = hashlib.sha256(response.content).hexdigest()
    records = response.json()["data"]
    if len(records) != 1:
        raise SystemExit(f"the endpoint returned {len(records)} records and the run expects one")
    return records[0], digest


def tpm_summary(values: list[float]) -> dict:
    """Summarise a list of expression values in transcripts per million.

    Args:
        values: The per sample expression values.

    Returns:
        The sample count, the minimum, the quartiles, the maximum and the share of samples at or
        above the threshold of condition C19.
    """
    ordered = sorted(values)
    quartiles = statistics.quantiles(ordered, n=4, method="inclusive")
    at_or_above = sum(1 for value in ordered if value >= EXPRESSED_MEDIAN_TPM)
    return {
        "samples": len(ordered),
        "minimum_tpm": round(ordered[0], 4),
        "first_quartile_tpm": round(quartiles[0], 4),
        "median_tpm": round(statistics.median(ordered), 4),
        "third_quartile_tpm": round(quartiles[2], 4),
        "maximum_tpm": round(ordered[-1], 4),
        "samples_at_or_above_threshold": at_or_above,
        "share_at_or_above_threshold": round(at_or_above / len(ordered), 4),
    }


def verdict_of(median_tpm: float) -> str:
    """Apply the threshold of condition C19.

    Args:
        median_tpm: The median transcripts per million across the samples.

    Returns:
        EXPRESSED when the median reaches the threshold and NOT SHOWN otherwise.
    """
    return "EXPRESSED" if median_tpm >= EXPRESSED_MEDIAN_TPM else "NOT SHOWN"


def repository_commit() -> str:
    """Give the current commit of the repository.

    Returns:
        The short commit hash.
    """
    out = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                         capture_output=True, text=True, check=True)
    return out.stdout.strip()


def main() -> int:
    """Read the portal, apply the threshold and write the receipt.

    Returns:
        Zero when the receipt is written.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()

    reference = gene_record()
    record, digest = expression_record(reference["gencodeId"])
    summary = tpm_summary([float(value) for value in record["data"]])
    verdict = verdict_of(summary["median_tpm"])

    receipt = {
        "receipt_id": "t2_naprt_fibroblast_gtex_v1",
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "data_class": "A",
        "scope": "Public GTEx expression values of one gene in one tissue. No patient data.",
        "gate": str(GATE.relative_to(ROOT)),
        "gate_conditions": ["C18", "C19"],
        "repository_commit": repository_commit(),
        "script": "scripts/t2_naprt_fibroblast_gtex_v1.py",
        "source": {
            "name": "GTEx Portal",
            "endpoint": EXPRESSION,
            "dataset_id": record["datasetId"],
            "release": record["datasetId"].replace("gtex_", "").upper(),
            "tissue_site_detail_id": TISSUE,
            "unit": record.get("unit"),
            "response_sha256": digest,
        },
        "gene": {
            "symbol": reference["geneSymbol"],
            "gencode_id": reference["gencodeId"],
            "gencode_version": reference.get("gencodeVersion"),
            "description": reference.get("description"),
        },
        "threshold_median_tpm": EXPRESSED_MEDIAN_TPM,
        "summary": summary,
        "verdict": verdict,
        "what_the_verdict_covers": (
            "The presence of the NAPRT transcript in cultured fibroblasts of GTEx donors. It is not "
            "an enzyme activity and not a measurement of the stage 1 lines, so the NO ROUTE rule of "
            "section 12.3 of the gate, which reads the NAD+ metabolite endpoint under nicotinic "
            "acid, still decides the arm."),
    }
    args.out.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"{GENE_SYMBOL} in {TISSUE}: median {summary['median_tpm']} TPM over "
          f"{summary['samples']} samples, verdict {verdict}")
    print(f"receipt: {args.out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
