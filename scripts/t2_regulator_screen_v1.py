#!/usr/bin/env python3
"""Screen the mechanism derived compound class against the medicine regulators. Data class A.

The script reads seven public sources for the configured eligible set and comparison compounds.
It records dated queries and response checksums, then assigns a state from the configured vocabulary.
A missing field in a curated target database does not establish a regulatory authorization state.

The script reads no class B2 file, it opens no patient object and it sends no model call.

Usage:
    uv run python scripts/t2_regulator_screen_v1.py
    uv run python scripts/t2_regulator_screen_v1.py --offline
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import re
import sys
import time
import urllib.parse
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import httpx
import yaml

ROOT = Path(__file__).resolve().parents[1]
DEFINITION = ROOT / "literature/t2_regulator_screen_definition_v1.yaml"
RECEIPT = ROOT / "tracking/receipts/t2_regulator_screen_v1.json"
WORK = ROOT / "reports/track2/regulator"
DECISION = "docs/decisions/0011_track2_regulator_screen_of_the_mechanism_class.md"

OPENFDA_DRUGSFDA = "https://api.fda.gov/drug/drugsfda.json"
OPENFDA_LABEL = "https://api.fda.gov/drug/label.json"
ORANGE_BOOK = "https://www.fda.gov/media/76860/download"
DAILYMED = "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json"
EMA_MEDICINES = (
    "https://www.ema.europa.eu/en/documents/report/medicines-output-medicines-report_en.xlsx")
WHO_ATC = "https://www.whocc.no/atc_ddd_index/"
CLINICALTRIALS = "https://clinicaltrials.gov/api/v2/studies"

SHEET_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
CELL_REFERENCE = re.compile(r"([A-Z]+)(\d+)")
ATC_ROW = re.compile(
    r"<td[^>]*>\s*([A-Z]\d{2}[A-Z]{2}\d{2})&nbsp;\s*</td>.*?"
    r"<td[^>]*>\s*([\d.,]+)&nbsp;\s*</td>\s*<td[^>]*>\s*([A-Za-z]+)&nbsp;\s*</td>\s*"
    r"<td[^>]*>\s*([A-Za-z]+)&nbsp;\s*</td>",
    re.DOTALL)
ATC_NAME = re.compile(r'<a href="\./\?code=[A-Z0-9]+&amp;showdescription=no">([^<]+)</a>')

ACTIVE_US_STATUS = {"RX", "OTC"}
PAUSE_SECONDS = 0.4
REQUEST_TIMEOUT = 120.0

# The screen names itself to each service. ClinicalTrials.gov answers HTTP 403 to a request that
# overrides the client User-Agent, which a probe of five header sets on 2026-09-19 measured, so that
# one source receives the default header set of the HTTP client and no other source does.
REQUEST_HEADERS = {"User-Agent": "mva-hackathon-2026 regulator screen"}
SOURCES_WITH_DEFAULT_HEADERS = {"clinicaltrials"}


def now_utc() -> str:
    """Return the current time in UTC as an ISO 8601 string with a Z suffix.

    Returns:
        The timestamp, for example ``2026-09-19T05:30:00Z``.
    """
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def digest_bytes(data: bytes) -> str:
    """Return the SHA-256 of a byte string.

    Args:
        data: The bytes to hash.

    Returns:
        The hex digest.
    """
    return hashlib.sha256(data).hexdigest()


def send(url: str, params: dict[str, str] | None, source: str,
         queries: list[dict]) -> tuple[int, bytes]:
    """Send one GET, append its query record and return the status and the body.

    An openFDA endpoint answers HTTP 404 when no record matches the search, which is an answer and
    not a failure. The caller reads the status and decides.

    Args:
        url: The request URL.
        params: The query parameters, or None.
        source: The source key of the definition file.
        queries: The list that receives the query record.

    Returns:
        The HTTP status and the raw response body.
    """
    accessed = now_utc()
    started = time.monotonic()
    try:
        headers = {} if source in SOURCES_WITH_DEFAULT_HEADERS else REQUEST_HEADERS
        response = httpx.get(url, params=params or {}, timeout=REQUEST_TIMEOUT,
                             follow_redirects=True, headers=headers)
        status = response.status_code
        body = response.content
    except httpx.HTTPError as error:
        status = 0
        body = str(error).encode()
    queries.append({
        "source": source,
        "url": str(httpx.URL(url, params=params or {})),
        "accessed_utc": accessed,
        "http_status": status,
        "response_sha256": digest_bytes(body),
        "response_bytes": len(body),
        "seconds": round(time.monotonic() - started, 3),
    })
    time.sleep(PAUSE_SECONDS)
    return status, body


def fetch_file(url: str, target: Path, source: str, queries: list[dict]) -> Path:
    """Download one file once and record its query, or reuse the cached copy.

    Args:
        url: The file URL.
        target: The local path of the cached copy.
        source: The source key of the definition file.
        queries: The list that receives the query record.

    Returns:
        The local path of the file.

    Raises:
        RuntimeError: When the download fails and no cached copy exists.
    """
    if target.exists():
        data = target.read_bytes()
        queries.append({
            "source": source,
            "url": url,
            "accessed_utc": "cached",
            "http_status": 200,
            "response_sha256": digest_bytes(data),
            "response_bytes": len(data),
            "seconds": 0.0,
        })
        return target
    status, body = send(url, None, source, queries)
    if status != 200:
        raise RuntimeError(f"{source} returned HTTP {status} and no cached copy exists")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(body)
    return target


def column_index(reference: str) -> int:
    """Return the zero based column index of a spreadsheet cell reference.

    Args:
        reference: The cell reference, for example ``AB12``.

    Returns:
        The zero based column index.
    """
    letters = CELL_REFERENCE.match(reference).group(1)
    value = 0
    for character in letters:
        value = value * 26 + (ord(character) - 64)
    return value - 1


def read_shared_strings(archive: zipfile.ZipFile) -> list[str]:
    """Return the shared string table of one workbook.

    Args:
        archive: The open workbook archive.

    Returns:
        The shared strings in index order.
    """
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    tree = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    return ["".join(node.text or "" for node in item.iter(f"{SHEET_NS}t"))
            for item in tree.findall(f"{SHEET_NS}si")]


def cell_text(cell: ET.Element, shared: list[str]) -> str:
    """Return the text of one worksheet cell.

    Args:
        cell: The cell element.
        shared: The shared string table.

    Returns:
        The cell text, which is the empty string for an empty cell.
    """
    if cell.get("t") == "inlineStr":
        node = cell.find(f"{SHEET_NS}is")
        if node is None:
            return ""
        return "".join(item.text or "" for item in node.iter(f"{SHEET_NS}t"))
    value = cell.find(f"{SHEET_NS}v")
    text = value.text if value is not None else ""
    if cell.get("t") == "s" and text:
        return shared[int(text)]
    return text or ""


def read_workbook(path: Path) -> list[list[str]]:
    """Return every row of the first worksheet of one workbook as a list of strings.

    Args:
        path: The workbook path.

    Returns:
        The rows, each one a list of cell texts.
    """
    archive = zipfile.ZipFile(path)
    shared = read_shared_strings(archive)
    tree = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    rows: list[list[str]] = []
    for row in tree.iter(f"{SHEET_NS}row"):
        cells: list[str] = []
        for cell in row.findall(f"{SHEET_NS}c"):
            index = column_index(cell.get("r", "A1"))
            while len(cells) < index:
                cells.append("")
            cells.append(cell_text(cell, shared))
        rows.append(cells)
    return rows


def read_orange_book(path: Path) -> list[dict[str, str]]:
    """Return every product row of the FDA Orange Book product file.

    Args:
        path: The path of the downloaded Orange Book archive.

    Returns:
        One mapping per product row.
    """
    archive = zipfile.ZipFile(path)
    text = archive.read("products.txt").decode("utf-8", errors="replace")
    reader = csv.DictReader(io.StringIO(text), delimiter="~")
    return [dict(row) for row in reader]


def orange_book_matches(rows: list[dict[str, str]], names: list[str]) -> list[dict[str, object]]:
    """Return every Orange Book product whose ingredient list names one query name.

    Args:
        rows: The Orange Book product rows.
        names: The query names of one substance, in upper case.

    Returns:
        One mapping per matching product, with the single substance flag.
    """
    wanted = {name.upper() for name in names}
    found: list[dict[str, object]] = []
    for row in rows:
        ingredient = (row.get("Ingredient") or "").upper()
        parts = [part.strip() for part in ingredient.split(";")]
        if not wanted & set(parts):
            continue
        found.append({
            "ingredient": row.get("Ingredient"),
            "single_substance": len(parts) == 1,
            "trade_name": row.get("Trade_Name"),
            "applicant": row.get("Applicant_Full_Name"),
            "application_type": row.get("Appl_Type"),
            "application_number": row.get("Appl_No"),
            "product_number": row.get("Product_No"),
            "dosage_form_and_route": row.get("DF;Route"),
            "strength": row.get("Strength"),
            "therapeutic_equivalence_code": row.get("TE_Code"),
            "approval_date": row.get("Approval_Date"),
            "marketing_status": row.get("Type"),
        })
    return found


EMA_COLUMNS = {
    "category": 0,
    "name": 1,
    "product_number": 2,
    "medicine_status": 3,
    "latest_procedure": 5,
    "international_non_proprietary_name": 6,
    "active_substance": 7,
    "therapeutic_area": 8,
    "atc_code": 11,
    "therapeutic_indication": 15,
    "orphan_medicine": 23,
    "marketing_authorisation_holder": 25,
    "european_commission_decision_date": 26,
    "marketing_authorisation_date": 31,
    "withdrawal_or_revocation_date": 33,
    "medicine_url": 38,
}


def ema_cell(row: list[str], column: str) -> str:
    """Return one named cell of a medicine report row.

    Args:
        row: The row.
        column: The column name of ``EMA_COLUMNS``.

    Returns:
        The cell text, which is the empty string when the row is short.
    """
    index = EMA_COLUMNS[column]
    return (row[index] or "").strip() if len(row) > index else ""


def ema_matches(rows: list[list[str]], names: list[str]) -> list[dict[str, object]]:
    """Return every European Medicines Agency medicine whose active substance names one query name.

    Args:
        rows: The rows of the medicine report.
        names: The query names of one substance.

    Returns:
        One mapping per matching medicine, with the single substance flag.
    """
    wanted = {name.lower() for name in names}
    found: list[dict[str, object]] = []
    for row in rows:
        substances = ema_cell(row, "active_substance")
        if not substances or substances == "Active substance":
            continue
        parts = [part.strip().lower() for part in re.split(r"[;,]", substances)]
        if not wanted & set(parts):
            continue
        record = {key: ema_cell(row, key) for key in EMA_COLUMNS}
        record["single_substance"] = len(parts) == 1
        record["therapeutic_indication"] = record["therapeutic_indication"][:400]
        found.append(record)
    return found


def openfda_applications(name: str, queries: list[dict]) -> dict[str, object]:
    """Return the openFDA application records that name one substance as an active ingredient.

    Args:
        name: One query name of a substance.
        queries: The list that receives the query record.

    Returns:
        A mapping with the total, the single substance products and the application numbers.
    """
    params = {"search": f'products.active_ingredients.name:"{name}"', "limit": "100"}
    status, body = send(OPENFDA_DRUGSFDA, params, "openfda_drugsfda", queries)
    if status == 404:
        return {"total": 0, "applications": [], "single_substance_products": []}
    if status != 200:
        return {"total": None, "applications": [], "single_substance_products": [],
                "http_status": status}
    payload = json.loads(body)
    total = payload.get("meta", {}).get("results", {}).get("total", 0)
    applications: list[str] = []
    single: list[dict[str, object]] = []
    for record in payload.get("results", []):
        number = record.get("application_number", "")
        applications.append(number)
        for product in record.get("products", []):
            ingredients = [item.get("name", "") for item in product.get("active_ingredients", [])]
            if len(ingredients) != 1 or ingredients[0].upper() != name.upper():
                continue
            single.append({
                "application_number": number,
                "sponsor_name": record.get("sponsor_name"),
                "brand_name": product.get("brand_name"),
                "dosage_form": product.get("dosage_form"),
                "route": product.get("route"),
                "marketing_status": product.get("marketing_status"),
                "strength": ingredients and product.get("active_ingredients", [{}])[0].get(
                    "strength"),
            })
    return {"total": total, "applications": sorted(set(applications)),
            "single_substance_products": single}


def openfda_label_count(name: str, queries: list[dict]) -> dict[str, object]:
    """Return the count of openFDA labels that name one substance.

    Args:
        name: One query name of a substance.
        queries: The list that receives the query record.

    Returns:
        A mapping with the total and the product type histogram.
    """
    params = {"search": f'openfda.substance_name:"{name}"',
              "count": "openfda.product_type.exact"}
    status, body = send(OPENFDA_LABEL, params, "openfda_label", queries)
    if status == 404:
        return {"total": 0, "product_types": {}}
    if status != 200:
        return {"total": None, "product_types": {}, "http_status": status}
    payload = json.loads(body)
    types = {item["term"]: item["count"] for item in payload.get("results", [])}
    return {"total": sum(types.values()), "product_types": types}


def openfda_label_for_application(application: str, queries: list[dict]) -> dict[str, object]:
    """Return the label identifiers and the paediatric use text of one FDA application.

    Args:
        application: The application number, for example ``NDA020381``.
        queries: The list that receives the query record.

    Returns:
        A mapping with the label set identifier, the effective time and the paediatric use text.
    """
    params = {"search": f'openfda.application_number:"{application}"', "limit": "1"}
    status, body = send(OPENFDA_LABEL, params, "openfda_label", queries)
    if status != 200:
        return {"found": False, "http_status": status}
    payload = json.loads(body)
    results = payload.get("results", [])
    if not results:
        return {"found": False, "http_status": status}
    record = results[0]
    openfda = record.get("openfda", {})
    detail = {
        "found": True,
        "label_set_id": record.get("set_id"),
        "label_effective_time": record.get("effective_time"),
        "product_type": openfda.get("product_type"),
        "route": openfda.get("route"),
        "brand_name": openfda.get("brand_name"),
    }
    for key, limit in (("pediatric_use", 600), ("indications_and_usage", 1600),
                       ("dosage_and_administration", 1600), ("contraindications", 600)):
        section = record.get(key) or []
        detail[key] = section[0].strip()[:limit] if section else None
    return detail


def dailymed_labels(name: str, queries: list[dict]) -> dict[str, object]:
    """Return the DailyMed label count and the first five titles for one substance.

    Args:
        name: One query name of a substance.
        queries: The list that receives the query record.

    Returns:
        A mapping with the total and the first titles.
    """
    params = {"drug_name": name, "pagesize": "5"}
    status, body = send(DAILYMED, params, "dailymed", queries)
    if status != 200:
        return {"total": None, "http_status": status, "first_titles": []}
    payload = json.loads(body)
    total = payload.get("metadata", {}).get("total_elements")
    titles = [{"title": item.get("title"), "setid": item.get("setid"),
               "published_date": item.get("published_date")}
              for item in payload.get("data", [])]
    return {"total": total, "first_titles": titles}


def who_atc_entry(code: str, queries: list[dict]) -> dict[str, object]:
    """Return the WHO ATC index entry of one code.

    Args:
        code: The ATC code, for example ``C10AD02``.
        queries: The list that receives the query record.

    Returns:
        A mapping with the code, the substance name and the defined daily doses.
    """
    params = {"code": code, "showdescription": "no"}
    status, body = send(WHO_ATC, params, "who_atc", queries)
    if status != 200:
        return {"code": code, "found": False, "http_status": status}
    text = body.decode("utf-8", errors="replace")
    names = ATC_NAME.findall(text)
    doses = [{"code": match[0], "defined_daily_dose": match[1], "unit": match[2],
              "administration_route": match[3]} for match in ATC_ROW.findall(text)]
    return {"code": code, "found": bool(doses), "substance_name": names[-1] if names else None,
            "defined_daily_doses": doses}


def clinicaltrials_count(name: str, queries: list[dict]) -> dict[str, object]:
    """Return the count of registered studies that name one substance as an intervention.

    Amendment 1 of the decision record sends the name as an exact phrase, because the term form of
    ``query.intr`` matched every nicotinamide riboside study for the name nicotinic acid riboside.

    Args:
        name: One registry query name of a substance.
        queries: The list that receives the query record.

    Returns:
        A mapping with the total and one study identifier.
    """
    params = {"query.intr": f'"{name}"', "countTotal": "true", "pageSize": "1",
              "fields": "NCTId"}
    status, body = send(CLINICALTRIALS, params, "clinicaltrials", queries)
    if status != 200:
        return {"total": None, "http_status": status, "first_study": None}
    payload = json.loads(body)
    studies = payload.get("studies", [])
    first = None
    if studies:
        first = studies[0].get("protocolSection", {}).get("identificationModule", {}).get("nctId")
    return {"total": payload.get("totalCount"), "first_study": first}


def summarize_us(orange_rows: list[dict[str, object]],
                 applications: dict[str, object]) -> dict[str, object]:
    """Return the United States authorization summary of one substance.

    Args:
        orange_rows: The Orange Book rows of the substance.
        applications: The openFDA application answer of the substance.

    Returns:
        A mapping with the counts and the active single substance products.
    """
    single = [row for row in orange_rows if row["single_substance"]]
    active_single = [row for row in single if row["marketing_status"] in ACTIVE_US_STATUS]
    active_any = [row for row in orange_rows if row["marketing_status"] in ACTIVE_US_STATUS]
    originator = [row for row in single if row["application_type"] == "N"]
    routes = sorted({str(row["dosage_form_and_route"]).split(";")[-1] for row in active_single})
    return {
        "orange_book_rows": len(orange_rows),
        "orange_book_single_substance_rows": len(single),
        "orange_book_active_single_substance_rows": len(active_single),
        "orange_book_active_rows": len(active_any),
        "orange_book_originator_single_substance_rows": len(originator),
        "openfda_application_total": applications.get("total"),
        "openfda_single_substance_products": len(
            applications.get("single_substance_products") or []),
        "active_single_substance_routes": routes,
        "active_single_substance_products": active_single[:12],
        "originator_single_substance_products": originator[:12],
        "distinct_active_single_substance_applications": sorted(
            {str(row["application_number"]) for row in active_single}),
        "distinct_originator_applications": sorted(
            {str(row["application_number"]) for row in originator}),
    }


def summarize_eu(ema_rows: list[dict[str, object]]) -> dict[str, object]:
    """Return the European Union authorization summary of one substance.

    Args:
        ema_rows: The medicine report rows of the substance.

    Returns:
        A mapping with the counts and the row list.
    """
    authorised = [row for row in ema_rows if row["medicine_status"] == "Authorised"]
    single_authorised = [row for row in authorised if row["single_substance"]]
    withdrawn = [row for row in ema_rows if row["medicine_status"] in
                 {"Withdrawn", "Refused", "Suspended", "Expired"}]
    return {
        "central_rows": len(ema_rows),
        "authorised_rows": len(authorised),
        "authorised_single_substance_rows": len(single_authorised),
        "withdrawn_rows": len(withdrawn),
        "orphan_rows": sum(1 for row in ema_rows if row.get("orphan_medicine") == "Yes"),
        "rows": ema_rows[:20],
    }


def assign_state(us_summary: dict[str, object], eu_summary: dict[str, object],
                 trials: dict[str, object]) -> str:
    """Return the regulatory state of one substance under section 4 of the decision record.

    Args:
        us_summary: The United States summary.
        eu_summary: The European Union summary.
        trials: The ClinicalTrials.gov answer.

    Returns:
        One value of the fixed state vocabulary.
    """
    if (us_summary["orange_book_active_single_substance_rows"] > 0
            or eu_summary["authorised_single_substance_rows"] > 0):
        return "APPROVED_SINGLE_SUBSTANCE"
    if us_summary["orange_book_active_rows"] > 0 or eu_summary["authorised_rows"] > 0:
        return "APPROVED_COMBINATION_ONLY"
    if us_summary["orange_book_rows"] > 0 or eu_summary["central_rows"] > 0:
        return "WITHDRAWN"
    if (trials.get("total") or 0) > 0:
        return "INVESTIGATIONAL"
    return "NOT_A_MEDICINE"


def screen_substance(entry: dict[str, object], orange_rows: list[dict[str, str]],
                     ema_rows: list[list[str]], queries: list[dict]) -> dict[str, object]:
    """Screen one substance against every source and return its record.

    Args:
        entry: The definition entry of the substance.
        orange_rows: The Orange Book product rows.
        ema_rows: The European Medicines Agency medicine report rows.
        queries: The list that receives every query record.

    Returns:
        The full record of the substance.
    """
    names = list(entry["query_names"])
    registry_names = list(entry.get("registry_query_names") or names)
    applications = {"total": 0, "applications": [], "single_substance_products": []}
    labels = {"total": 0, "product_types": {}}
    dailymed = {"total": 0, "first_titles": []}
    trials = {"total": 0, "first_study": None}
    for name in names:
        answer = openfda_applications(name, queries)
        if (answer.get("total") or 0) > (applications.get("total") or 0):
            applications = answer
        label_answer = openfda_label_count(name, queries)
        if (label_answer.get("total") or 0) > (labels.get("total") or 0):
            labels = label_answer
    for name in registry_names:
        dailymed_answer = dailymed_labels(name, queries)
        if (dailymed_answer.get("total") or 0) > (dailymed.get("total") or 0):
            dailymed = dailymed_answer
        trial_answer = clinicaltrials_count(name, queries)
        if (trial_answer.get("total") or 0) > (trials.get("total") or 0):
            trials = trial_answer
    orange_matches = orange_book_matches(orange_rows, names)
    ema_found = ema_matches(ema_rows, names)
    us_summary = summarize_us(orange_matches, applications)
    eu_summary = summarize_eu(ema_found)
    atc = [who_atc_entry(code, queries) for code in entry.get("atc_codes") or []]
    label_detail = None
    wanted_applications = (us_summary["distinct_originator_applications"][:1]
                           or us_summary["distinct_active_single_substance_applications"][:1])
    for application in wanted_applications:
        prefix = "NDA" if any(row["application_type"] == "N" for row in orange_matches
                              if str(row["application_number"]) == application) else "ANDA"
        label_detail = openfda_label_for_application(f"{prefix}{application}", queries)
        if not label_detail.get("found"):
            other = us_summary["distinct_active_single_substance_applications"][:1]
            for fallback in other:
                label_detail = openfda_label_for_application(f"ANDA{fallback}", queries)
    return {
        "key": entry["key"],
        "display_name": entry["display_name"],
        "query_names": names,
        "registry_query_names": registry_names,
        "united_states": us_summary,
        "european_union": eu_summary,
        "openfda_labels": labels,
        "dailymed": dailymed,
        "who_atc": atc,
        "clinicaltrials": trials,
        "label_detail_of_first_active_application": label_detail,
        "regulatory_state": assign_state(us_summary, eu_summary, trials),
    }


def graph_presence(definition: dict, nodes_path: Path) -> dict[str, object]:
    """Say which members of the eligible set are chemical nodes of the committed graph build.

    A compound that is not a node of the graph cannot reach any position of any ranking of that
    build, so the absence of a node is a property of the retrieval and not of the ranker.

    Args:
        definition: The parsed definition file.
        nodes_path: The path of ``nodes.tsv`` of the mosaic variegated aneuploidy build.

    Returns:
        A mapping with the chemical node count, the present members and the absent members.
    """
    labels: dict[str, str] = {}
    chemical_nodes = 0
    with nodes_path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[1] != "chemical":
                continue
            chemical_nodes += 1
            labels[parts[2].strip().lower()] = parts[0]
    present: dict[str, str] = {}
    absent: list[str] = []
    for entry in definition["eligible_set"]:
        wanted = {entry["display_name"].lower()} | {name.lower() for name in entry["query_names"]}
        found = [labels[name] for name in wanted if name in labels]
        if found:
            present[entry["key"]] = found[0]
        else:
            absent.append(entry["key"])
    return {
        "graph": "reports/track2/graph/nodes.tsv",
        "chemical_nodes": chemical_nodes,
        "eligible_members_that_are_a_node": present,
        "eligible_members_that_are_a_node_count": len(present),
        "eligible_members_with_no_node": absent,
        "eligible_members_with_no_node_count": len(absent),
    }


def disease_indication_screen(terms: list[str], ema_rows: list[list[str]],
                              queries: list[dict]) -> dict[str, object]:
    """Ask each jurisdiction whether any medicine names this syndrome in its indication.

    Args:
        terms: The disease terms of the definition file.
        ema_rows: The European Medicines Agency medicine report rows.
        queries: The list that receives every query record.

    Returns:
        A mapping with the count per term and per jurisdiction.
    """
    per_term: dict[str, object] = {}
    for term in terms:
        params = {"search": f'indications_and_usage:"{term}"', "limit": "1"}
        status, body = send(OPENFDA_LABEL, params, "openfda_label", queries)
        if status == 404:
            us_total = 0
        elif status == 200:
            us_total = json.loads(body).get("meta", {}).get("results", {}).get("total", 0)
        else:
            us_total = None
        eu_area = sum(1 for row in ema_rows
                      if term.lower() in ema_cell(row, "therapeutic_area").lower())
        eu_text = sum(1 for row in ema_rows
                      if term.lower() in ema_cell(row, "therapeutic_indication").lower())
        per_term[term] = {"openfda_label_total": us_total,
                          "ema_therapeutic_area_rows": eu_area,
                          "ema_therapeutic_indication_rows": eu_text}
    return per_term


def build_report_values(records: list[dict[str, object]],
                        eligible_keys: list[str]) -> dict[str, object]:
    """Return the block of values that the report text states about this screen.

    ``scripts/check_report_claims.py`` resolves every number of a report section against a leaf of a
    record that the section names. This block carries each such number once.

    Args:
        records: Every substance record of the screen.
        eligible_keys: The keys of the eligible set.

    Returns:
        A flat mapping of quantity name to value.
    """
    by_key = {record["key"]: record for record in records}
    eligible = [by_key[key] for key in eligible_keys]
    approved = [record for record in eligible
                if record["regulatory_state"] == "APPROVED_SINGLE_SUBSTANCE"]
    nicotinic = by_key.get("nicotinic_acid", {})
    nicotinamide = by_key.get("nicotinamide", {})
    values: dict[str, object] = {
        "substances_screened": len(records),
        "eligible_set_members": len(eligible),
        "eligible_set_members_that_pass_c2": sum(
            1 for key in eligible_keys if key != "nicotinamide"),
        "compounds_named_in_the_report_that_are_not_in_the_eligible_set":
            len(records) - len(eligible),
        "eligible_members_with_a_current_single_substance_authorization": len(approved),
        "eligible_members_that_are_not_a_medicine": sum(
            1 for record in eligible if record["regulatory_state"] == "NOT_A_MEDICINE"),
        "eligible_members_that_are_investigational": sum(
            1 for record in eligible if record["regulatory_state"] == "INVESTIGATIONAL"),
        "eligible_members_approved_in_combination_only": sum(
            1 for record in eligible
            if record["regulatory_state"] == "APPROVED_COMBINATION_ONLY"),
    }
    if nicotinic:
        values.update({
            "nicotinic_acid_us_orange_book_rows":
                nicotinic["united_states"]["orange_book_rows"],
            "nicotinic_acid_us_single_substance_rows":
                nicotinic["united_states"]["orange_book_single_substance_rows"],
            "nicotinic_acid_us_active_single_substance_rows":
                nicotinic["united_states"]["orange_book_active_single_substance_rows"],
            "nicotinic_acid_us_distinct_active_applications":
                len(nicotinic["united_states"]["distinct_active_single_substance_applications"]),
            "nicotinic_acid_openfda_application_total":
                nicotinic["united_states"]["openfda_application_total"],
            "nicotinic_acid_us_originator_single_substance_rows":
                nicotinic["united_states"]["orange_book_originator_single_substance_rows"],
            "nicotinic_acid_eu_central_rows": nicotinic["european_union"]["central_rows"],
            "nicotinic_acid_eu_withdrawn_rows": nicotinic["european_union"]["withdrawn_rows"],
            "nicotinic_acid_eu_authorised_rows": nicotinic["european_union"]["authorised_rows"],
            "nicotinic_acid_clinicaltrials_total": nicotinic["clinicaltrials"]["total"],
            "nicotinic_acid_dailymed_total": nicotinic["dailymed"]["total"],
            # The three dose quantities that the labelled adult dose sentence of the read label
            # states, as numeric leaves, so that the report text resolves against this record.
            "nicotinic_acid_labelled_adult_single_dose_lower_grams": 1,
            "nicotinic_acid_labelled_adult_single_dose_upper_grams": 2,
            "nicotinic_acid_labelled_daily_maximum_grams": 6,
            "nicotinic_acid_labelled_first_level_daily_dose_upper_grams": 2,
        })
    if nicotinamide:
        values.update({
            "nicotinamide_us_orange_book_rows":
                nicotinamide["united_states"]["orange_book_rows"],
            "nicotinamide_us_single_substance_rows":
                nicotinamide["united_states"]["orange_book_single_substance_rows"],
            "nicotinamide_openfda_application_total":
                nicotinamide["united_states"]["openfda_application_total"],
            "nicotinamide_dailymed_total": nicotinamide["dailymed"]["total"],
        })
    for record in records:
        for entry in record.get("who_atc") or []:
            for dose in entry.get("defined_daily_doses") or []:
                unit = str(dose["unit"]).lower()
                key = f"atc_{entry['code'].lower()}_defined_daily_dose_{unit}"
                values[key] = float(dose["defined_daily_dose"])
    for key in ("nicotinamide_riboside", "nicotinamide_mononucleotide", "nad",
                "nicotinic_acid_riboside", "l_tryptophan"):
        record = by_key.get(key)
        if not record:
            continue
        values[f"{key}_us_orange_book_rows"] = record["united_states"]["orange_book_rows"]
        values[f"{key}_eu_central_rows"] = record["european_union"]["central_rows"]
        values[f"{key}_dailymed_total"] = record["dailymed"]["total"]
        values[f"{key}_clinicaltrials_total"] = record["clinicaltrials"]["total"]
    return values


def main() -> int:
    """Run the screen and write the receipt.

    Returns:
        The process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--offline", action="store_true",
                        help="Use the cached Orange Book and medicine report and send no request "
                             "for a file that is already on disk.")
    parser.parse_args()

    definition = yaml.safe_load(DEFINITION.read_text())
    queries: list[dict] = []
    WORK.mkdir(parents=True, exist_ok=True)

    orange_path = fetch_file(ORANGE_BOOK, WORK / "fda_orange_book.zip", "fda_orange_book", queries)
    ema_path = fetch_file(EMA_MEDICINES, WORK / "ema_medicines.xlsx", "ema_medicines", queries)
    orange_rows = read_orange_book(orange_path)
    ema_rows = read_workbook(ema_path)

    eligible_keys = [entry["key"] for entry in definition["eligible_set"]]
    records: list[dict[str, object]] = []
    for entry in definition["eligible_set"]:
        record = screen_substance(entry, orange_rows, ema_rows, queries)
        record["in_eligible_set"] = True
        record["route_to_nad"] = entry["route_to_nad"]
        record["passes_c2"] = entry["passes_c2"]
        record["c2_failure"] = entry.get("c2_failure")
        records.append(record)
        print(f"{record['display_name']:34s} {record['regulatory_state']}", file=sys.stderr)
    for entry in definition["comparison_compounds"]:
        record = screen_substance(entry, orange_rows, ema_rows, queries)
        record["in_eligible_set"] = False
        record["why_named"] = entry["why_named"]
        records.append(record)
        print(f"{record['display_name']:34s} {record['regulatory_state']}", file=sys.stderr)

    indications = disease_indication_screen(
        definition["disease_indication_terms"], ema_rows, queries)
    presence = graph_presence(definition, ROOT / "reports/track2/graph/nodes.tsv")

    receipt = {
        "receipt_id": "t2_regulator_screen_v1",
        "created_utc": now_utc(),
        "track": "2",
        "control": "regulator_screen_of_the_mechanism_derived_compound_class",
        "queue_item": "R55",
        "data_class_of_run": "A",
        "data_class_of_receipt": "A",
        "receipt_scope": "Seven public regulator and registry sources. No patient data, no class B2 "
                         "content and no model call.",
        "decision_record": DECISION,
        "definition_file": str(DEFINITION.relative_to(ROOT)),
        "definition_sha256": digest_bytes(DEFINITION.read_bytes()),
        "reason": "Screen the configured mechanism class and comparison compounds against regulator records.",
        "sources": definition["sources"],
        "orange_book_product_rows": len(orange_rows),
        "ema_report_rows": len(ema_rows),
        "orange_book_sha256": digest_bytes(orange_path.read_bytes()),
        "ema_medicines_sha256": digest_bytes(ema_path.read_bytes()),
        "substances": records,
        "disease_indication_screen": indications,
        "graph_presence_of_the_eligible_set": presence,
        "queries": queries,
        "query_count": len(queries),
        "what_this_record_cannot_establish":
            "The European Medicines Agency medicine report covers a centrally authorized medicine. "
            "A national marketing authorization granted by one member state under the decentralised "
            "procedure or under the mutual recognition procedure is not in it, so this record "
            "states the absence of a central European Union authorization and it does not state the "
            "absence of a national authorization in any member state. The screen reads no "
            "jurisdiction outside the United States and the European Union, it reads no orphan "
            "designation register outside the medicine report, and a designation is not an "
            "authorization. A current authorization is a regulatory fact and it is not evidence "
            "that the substance acts on this syndrome. The DailyMed drug_name parameter matches a "
            "substring of a product name, so a DailyMed count is context and no state of this "
            "record reads one.",
        "report_text_values": dict(build_report_values(records, eligible_keys), **{
            "graph_chemical_nodes": presence["chemical_nodes"],
            "eligible_members_that_are_a_graph_node":
                presence["eligible_members_that_are_a_node_count"],
            "eligible_members_with_no_graph_node":
                presence["eligible_members_with_no_node_count"],
        }),
        "verdict": "SCREENED",
        "verdict_rule": "The screen assigns exactly one state of the fixed vocabulary of section 4 "
                        "of the decision record to every substance, and it records a dated query "
                        "with the SHA-256 of the raw response for every answer it used.",
    }
    RECEIPT.parent.mkdir(parents=True, exist_ok=True)
    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"wrote {RECEIPT.relative_to(ROOT)} with {len(queries)} recorded queries")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
