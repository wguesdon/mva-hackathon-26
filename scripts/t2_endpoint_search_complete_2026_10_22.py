#!/usr/bin/env python3
"""Retrieve every page of the registered public endpoint queries.

This script preserves the earlier capped search. It checks retrieval completeness
separately from scientific interpretation and reads no individual patient data.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import time

import httpx
import yaml

from literature_pack import clean

ROOT = Path(__file__).resolve().parents[1]
URL = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
PAGE_SIZE = 100
CUTOFF = '2026-10-22'
COMPOUND = re.compile(r'\bnicotinic acid\b|\bniacin\b', re.I)
BUBR1 = re.compile(r'\bbubr1\b|\bbub1b\b', re.I)
TARGET = re.compile(r'spindle|checkpoint|mitotic|chromosom|missegregat|aneuploid|sirt2|fibroblast|naprt|keratinocyt|endothel|blood mononuclear', re.I)


def digest(path: Path) -> str:
    """Return the SHA-256 digest of a public file.

    Args:
        path: The public input or output file.

    Returns:
        The hexadecimal digest.
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_registration(path: Path, commit: str) -> None:
    """Require the registered decision text to match its committed version.

    Args:
        path: The registration file.
        commit: The commit that records the registration.
    """
    relative = path.resolve().relative_to(ROOT)
    result = subprocess.run(['git', 'show', f'{commit}:{relative}'], cwd=ROOT,
                            capture_output=True, check=True)
    if hashlib.sha256(result.stdout).hexdigest() != digest(path):
        raise ValueError('The registration differs from the named commit.')


def retrieve_topic(client: httpx.Client, topic: dict, out: Path) -> tuple[list[dict], dict]:
    """Retrieve all cursor pages and validate the reported source-record total.

    Args:
        client: The HTTP client.
        topic: The original registered query topic.
        out: The new output directory.

    Returns:
        All records and the completeness record for the topic.
    """
    query = f"({topic['query']}) AND (SRC:MED OR SRC:PMC OR SRC:PPR) AND FIRST_PDATE:[* TO {CUTOFF}]"
    cursor = '*'
    cursors = set()
    records = []
    pages = []
    expected = None
    source_ids = set()
    while True:
        if cursor in cursors:
            raise ValueError('A query repeated its request cursor.')
        cursors.add(cursor)
        params = {'query': query, 'format': 'json', 'resultType': 'core',
                  'pageSize': PAGE_SIZE, 'cursorMark': cursor}
        response = client.get(URL, params=params)
        response.raise_for_status()
        payload = response.json()
        raw_path = out / 'raw' / f"{topic['id']}_page_{len(pages) + 1:04d}.json"
        raw_path.write_bytes(response.content)
        reported = int(payload['hitCount'])
        if expected is None:
            expected = reported
        if reported != expected:
            raise ValueError('The reported hit count changed during pagination.')
        batch = payload.get('resultList', {}).get('result', [])
        next_cursor = payload.get('nextCursorMark')
        pages.append({'page': len(pages) + 1, 'request_url': str(response.url),
                      'request_cursor': cursor, 'next_cursor': next_cursor,
                      'reported_hit_count': reported, 'returned': len(batch),
                      'retrieved_utc': datetime.now(timezone.utc).isoformat(),
                      'path': str(raw_path.relative_to(ROOT)), 'sha256': digest(raw_path)})
        for item in batch:
            key = f"{item['source']}:{item['id']}"
            if key in source_ids:
                raise ValueError('A query returned a duplicate source identifier.')
            source_ids.add(key)
            records.append(item)
        if len(records) > expected:
            raise ValueError('The retrieved total exceeds the reported hit count.')
        if len(records) == expected:
            termination = 'The retrieved unique source-record count equals the stable hit count.'
            break
        if not batch or not next_cursor or next_cursor in cursors:
            raise ValueError('Pagination ended before the reported hit count was retrieved.')
        cursor = next_cursor
        time.sleep(0.4)
    return records, {'topic': topic['id'], 'query': query, 'hit_count': expected,
                     'retrieved': len(records), 'unique_source_ids': len(source_ids),
                     'pages': pages, 'termination': termination, 'verdict': 'PASS'}


def identity(record: dict) -> str:
    """Return the deterministic identity used for cross-topic deduplication.

    Args:
        record: A public Europe PMC core record.

    Returns:
        A PMID when available, otherwise a source-specific identifier.
    """
    return f"PMID:{record['pmid']}" if record.get('pmid') else f"{record['source']}:{record['id']}"


def targeted_reasons(text: str) -> list[str]:
    """Select records for targeted assessment using registered lexical terms.

    Args:
        text: The combined title and abstract.

    Returns:
        The lexical selection reasons, without a scientific inclusion verdict.
    """
    reasons = []
    if BUBR1.search(text):
        reasons.append('The title or abstract names BubR1 or BUB1B.')
    if COMPOUND.search(text) and TARGET.search(text):
        reasons.append('The title or abstract names nicotinic acid or niacin and a target term.')
    return reasons


def merge_records(topics: list[tuple[str, list[dict]]]) -> list[dict]:
    """Deduplicate records while retaining every topic and source identity.

    Args:
        topics: Topic names and their retrieved records.

    Returns:
        The merged public record list.
    """
    merged = {}
    for topic, records in topics:
        for record in records:
            key = identity(record)
            source_id = f"{record['source']}:{record['id']}"
            if key not in merged:
                merged[key] = {'paper_id': key, 'topics': [], 'source_ids': [],
                               'title': clean(record.get('title')),
                               'abstract': clean(record.get('abstractText')),
                               'pmid': record.get('pmid'), 'pmcid': record.get('pmcid'),
                               'doi': record.get('doi'),
                               'first_publication_date': record.get('firstPublicationDate'),
                               'publication_types': record.get('pubTypeList', {}).get('pubType', []),
                               'assessment_status': 'UNASSESSED'}
            row = merged[key]
            if topic not in row['topics']:
                row['topics'].append(topic)
            if source_id not in row['source_ids']:
                row['source_ids'].append(source_id)
            if not row['abstract'] and record.get('abstractText'):
                row['abstract'] = clean(record['abstractText'])
    for row in merged.values():
        row['targeted_assessment_reasons'] = targeted_reasons(row['title'] + ' ' + row['abstract'])
    return [merged[key] for key in sorted(merged)]


def write_targeted_packet(rows: list[dict], path: Path) -> None:
    """Write a packet without assigning scientific inclusion decisions.

    Args:
        rows: The deduplicated public records.
        path: The new Markdown output file.
    """
    lines = ['# Targeted assessment packet', '',
             'The lexical filter selects these records for assessment. It does not establish inclusion.', '']
    for row in rows:
        if not row['targeted_assessment_reasons']:
            continue
        lines.extend([f"## {row['paper_id']}", '', row['title'], '',
                      row['abstract'] or 'The retrieved record has no abstract.', ''])
    path.write_text('\n'.join(lines) + '\n')


def main() -> None:
    """Run the committed retrieval specification and record every page."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration', type=Path, required=True)
    parser.add_argument('--registration-commit', required=True)
    args = parser.parse_args()
    verify_registration(args.registration, args.registration_commit)
    queries = ROOT / 'literature/task2_bubr1_endpoint_queries.yaml'
    config = yaml.safe_load(queries.read_text())
    stamp = datetime.now(timezone.utc).strftime('%Y_%m_%d_%H%M%S')
    out = ROOT / 'reports/literature' / f'track2_endpoint_complete_{stamp}'
    (out / 'raw').mkdir(parents=True, exist_ok=False)
    topics = []
    checks = []
    status = 'FAIL'
    error = None
    try:
        with httpx.Client(timeout=90.0, follow_redirects=True) as client:
            for topic in config['topics']:
                records, check = retrieve_topic(client, topic, out)
                topics.append((topic['id'], records))
                checks.append(check)
                print(f"{topic['id']}: {check['retrieved']} of {check['hit_count']}", flush=True)
        rows = merge_records(topics)
        (out / 'records.json').write_text(json.dumps(rows, indent=2) + '\n')
        write_targeted_packet(rows, out / 'targeted_packet.md')
        status = 'PASS'
    except (httpx.HTTPError, ValueError, KeyError) as failure:
        rows = []
        error = str(failure)
    receipt = {'created_utc': datetime.now(timezone.utc).isoformat(), 'data_class': 'A',
               'script': str(Path(__file__).relative_to(ROOT)), 'script_sha256': digest(Path(__file__)),
               'registration': str(args.registration), 'registration_commit': args.registration_commit,
               'registration_sha256': digest(args.registration), 'query_file': str(queries.relative_to(ROOT)),
               'query_sha256': digest(queries), 'publication_date_cutoff': CUTOFF,
               'date_field': 'FIRST_PDATE', 'page_size': PAGE_SIZE, 'retrieval_gate': status,
               'error': error, 'topics': checks, 'deduplicated_records': len(rows),
               'records_with_abstract': sum(bool(row['abstract']) for row in rows),
               'targeted_lexical_selection': sum(bool(row['targeted_assessment_reasons']) for row in rows),
               'scientifically_assessed_records': 0, 'deduplication': 'PMID when present; otherwise source:id.',
               'tools': {'python': platform.python_version(), 'httpx': httpx.__version__, 'PyYAML': yaml.__version__,
                         'uv': subprocess.run(['uv', '--version'], capture_output=True, text=True, check=True).stdout.strip()},
               'environment': 'Local uv environment; no container executes this retrieval.',
               'lock_sha256': digest(ROOT / 'uv.lock'),
               'files': [{'path': str(path.relative_to(ROOT)), 'sha256': digest(path)} for path in sorted(out.rglob('*')) if path.is_file()],
               'limitations': ['Complete retrieval of specified database queries does not establish exhaustive scientific evidence.',
                              'The lexical selection is not evidence adjudication. All records remain scientifically unassessed in this receipt.',
                              'The publication-date cutoff does not freeze the database indexing date.']}
    receipt_path = ROOT / 'tracking/receipts' / f't2_endpoint_search_complete_{stamp}.json'
    receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
    print(f'Receipt: {receipt_path}')
    print(f'Output: {out}')
    if status != 'PASS':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
