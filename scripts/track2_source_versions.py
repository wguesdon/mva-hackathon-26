"""Record the release identifier of every source that the Track 2 runner reads.

The Track 2 report states that the source endpoints return no release identifier, so a later access
cannot be pinned to the access of a recorded run. That statement holds for one of the four sources
and not for the other three. This module asks each source for its own version and writes the answer
into a receipt, so a run carries the release of every table that it read.

| Source | Version endpoint | What it returns |
|---|---|---|
| ChEMBL | `https://www.ebi.ac.uk/chembl/api/data/status.json` | The database version and its release date. |
| Open Targets | The `meta` field of the GraphQL API | The data version as a year and a month, and the API version. |
| PubMed through the NCBI E-utilities | `einfo.fcgi` on the `pubmed` database | The build identifier, the record count and the last update time. |
| PubTator3 | None | The API publishes no release identifier. The receipt records the access time and the PubMed build of the same day. |

Usage:
    uv run python scripts/track2_source_versions.py
    uv run python scripts/track2_source_versions.py --receipt-id t2_source_versions_roberts_v1
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
from pathlib import Path

import httpx

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402

RECEIPT_DIR = common.ROOT / "tracking" / "receipts"
CHEMBL_STATUS = "https://www.ebi.ac.uk/chembl/api/data/status.json"
OPEN_TARGETS_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
EINFO = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi"
PUBTATOR_ROOT = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/"

META_QUERY = '{"query":"{meta{name apiVersion{x y z} dataVersion{year month iteration}}}"}'


def chembl_version(client: httpx.Client) -> dict:
    """Read the ChEMBL release.

    Args:
        client: The HTTP client.

    Returns:
        The release identifier, the release date and the counts that the endpoint reports.
    """
    payload = client.get(CHEMBL_STATUS).json()
    return {
        "source": "chembl",
        "endpoint": CHEMBL_STATUS,
        "release": payload.get("chembl_db_version"),
        "release_date": payload.get("chembl_release_date"),
        "activities": payload.get("activities"),
        "distinct_compounds": payload.get("disinct_compounds"),
        "targets": payload.get("targets"),
        "licence": "CC BY-SA 3.0 Unported",
    }


def open_targets_version(client: httpx.Client) -> dict:
    """Read the Open Targets data release and API version.

    Args:
        client: The HTTP client.

    Returns:
        The data version, the API version and the endpoint.
    """
    payload = client.post(OPEN_TARGETS_GRAPHQL, content=META_QUERY,
                          headers={"Content-Type": "application/json"}).json()
    meta = (payload.get("data") or {}).get("meta") or {}
    api = meta.get("apiVersion") or {}
    data = meta.get("dataVersion") or {}
    return {
        "source": "open_targets",
        "endpoint": OPEN_TARGETS_GRAPHQL,
        "release": f"{data.get('year')}.{data.get('month')}",
        "api_version": f"{api.get('x')}.{api.get('y')}.{api.get('z')}",
        "name": meta.get("name"),
        "licence": "CC0 1.0 for the data, Apache License 2.0 for the code",
    }


def pubmed_version(client: httpx.Client) -> dict:
    """Read the PubMed build that the NCBI E-utilities serve.

    Args:
        client: The HTTP client.

    Returns:
        The build identifier, the record count and the last update time.
    """
    payload = client.get(EINFO, params={"db": "pubmed", "retmode": "json"}).json()
    info = (payload.get("einforesult") or {}).get("dbinfo") or [{}]
    record = info[0] if info else {}
    return {
        "source": "pubmed_eutils",
        "endpoint": EINFO,
        "release": record.get("dbbuild"),
        "record_count": record.get("count"),
        "last_update": record.get("lastupdate"),
        "licence": "Public domain for the NCBI generated content",
    }


def pubtator_version(client: httpx.Client, pubmed: dict) -> dict:
    """Record what PubTator3 publishes about its own release.

    Args:
        client: The HTTP client.
        pubmed: The PubMed build record, which the receipt keeps as the closest proxy.

    Returns:
        The access record.
    """
    response = client.get(PUBTATOR_ROOT)
    return {
        "source": "pubtator3",
        "endpoint": PUBTATOR_ROOT,
        "release": None,
        "release_note": ("The PubTator3 API publishes no release identifier at any endpoint that "
                         "this pipeline calls, so a later access returns the current annotation of "
                         "PubMed and not the annotation of this run. The PubMed build of the same "
                         "access is the closest proxy that a source supplies."),
        "pubmed_build_at_access": pubmed.get("release"),
        "http_status": response.status_code,
        "licence": "Public domain for the NCBI generated content",
    }


def collect() -> list[dict]:
    """Read the release of every source.

    Returns:
        One record per source, in the order that the receipt lists them.
    """
    with httpx.Client(timeout=60.0, follow_redirects=True,
                      headers={"User-Agent": common.USER_AGENT}) as client:
        pubmed = pubmed_version(client)
        return [pubtator_version(client, pubmed), pubmed, open_targets_version(client),
                chembl_version(client)]


def main() -> int:
    """Write the source version receipt.

    Returns:
        0 on success, 1 when a source did not answer.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--receipt-id", default="t2_source_versions_v1",
                        help="The receipt identifier.")
    parser.add_argument("--note", default="", help="A note for the receipt.")
    args = parser.parse_args()

    try:
        sources = collect()
    except (httpx.HTTPError, ValueError, KeyError) as error:
        print(f"FAIL: a source did not answer: {type(error).__name__}: {error}")
        return 1

    receipt = {
        "receipt_id": args.receipt_id,
        "note": args.note,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "track": "t2_rank",
        "control": "release_identifier_of_every_source_of_a_run",
        "queue_item": "R17",
        "decision_record": ("docs/decisions/"
                            "0009_track2_one_entry_point_and_a_prospective_disease.md"),
        "data_class_of_run": "A",
        "data_class_of_receipt": "A",
        "receipt_scope": ("Public release identifiers and public record counts only. No patient "
                          "data enters this run."),
        "script": "scripts/track2_source_versions.py",
        "sources": sources,
        "sources_with_a_release_identifier": sum(1 for row in sources if row.get("release")),
        "sources_total": len(sources),
    }
    path = common.write_json(RECEIPT_DIR / f"{args.receipt_id}.json", receipt)
    for row in sources:
        print(f"  {row['source']}: release {row.get('release') or 'none published'}")
    print(f"{receipt['sources_with_a_release_identifier']} of {receipt['sources_total']} sources "
          "publish a release identifier.")
    print(f"receipt: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
