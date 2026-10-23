"""Test public literature pagination without network access or patient fixtures."""
from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import httpx

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import t2_endpoint_search_complete_2026_10_22 as search  # noqa: E402


class FakeClient:
    """Return recorded synthetic pages and retain the requested cursors."""

    def __init__(self, payloads: list[dict]) -> None:
        """Store synthetic API payloads.

        Args:
            payloads: Pages returned in request order.
        """
        self.payloads = payloads
        self.cursors: list[str] = []

    def get(self, url: str, params: dict) -> httpx.Response:
        """Return the next synthetic response.

        Args:
            url: The requested endpoint.
            params: The query and cursor parameters.

        Returns:
            A successful JSON response with the original request attached.
        """
        self.cursors.append(params['cursorMark'])
        payload = self.payloads[len(self.cursors) - 1]
        return httpx.Response(200, json=payload, request=httpx.Request('GET', url, params=params))


def record(identifier: str, abstract: str | None = 'A public synthetic abstract.') -> dict:
    """Return a synthetic source record.

    Args:
        identifier: The synthetic accession.
        abstract: The abstract, or None to omit it.

    Returns:
        The public metadata fixture.
    """
    result = {'source': 'MED', 'id': identifier, 'pmid': identifier, 'title': 'BubR1 research'}
    if abstract is not None:
        result['abstractText'] = abstract
    return result


def page(total: int, records: list[dict], cursor: str | None = None) -> dict:
    """Return a synthetic cursor page.

    Args:
        total: The source-reported total.
        records: The returned records.
        cursor: The next cursor, when present.

    Returns:
        A page in the Europe PMC response structure.
    """
    result = {'hitCount': total, 'resultList': {'result': records}}
    if cursor is not None:
        result['nextCursorMark'] = cursor
    return result


class CompleteSearchTests(unittest.TestCase):
    """Reject incomplete retrieval while preserving assessment uncertainty."""

    def setUp(self) -> None:
        """Prepare a temporary public output directory."""
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        (self.root / 'raw').mkdir()
        self.topic = {'id': 'synthetic_topic', 'query': 'BubR1'}
        self.root_patch = patch.object(search, 'ROOT', self.root)
        self.root_patch.start()
        self.sleep_patch = patch.object(search.time, 'sleep')
        self.sleep_patch.start()

    def tearDown(self) -> None:
        """Remove synthetic output and restore the module environment."""
        self.sleep_patch.stop()
        self.root_patch.stop()
        self.temporary.cleanup()

    def test_all_cursor_pages_are_required(self) -> None:
        """Complete traversal must retain records beyond the first page."""
        client = FakeClient([page(3, [record('1'), record('2')], 'second'),
                             page(3, [record('3')])])
        records, result = search.retrieve_topic(client, self.topic, self.root)
        self.assertEqual(client.cursors, ['*', 'second'])
        self.assertEqual([r['id'] for r in records], ['1', '2', '3'])
        self.assertEqual(result['verdict'], 'PASS')
        self.assertEqual(len(result['pages']), 2)

    def test_changing_hit_count_is_rejected(self) -> None:
        """A moving total must not produce a completeness verdict."""
        client = FakeClient([page(2, [record('1')], 'second'), page(3, [record('2')])])
        with self.assertRaisesRegex(ValueError, 'hit count changed'):
            search.retrieve_topic(client, self.topic, self.root)

    def test_duplicate_source_records_are_rejected(self) -> None:
        """A duplicate cannot count toward the reported source total."""
        client = FakeClient([page(2, [record('1')], 'second'), page(2, [record('1')])])
        with self.assertRaisesRegex(ValueError, 'duplicate source identifier'):
            search.retrieve_topic(client, self.topic, self.root)

    def test_early_cursor_termination_is_rejected(self) -> None:
        """Missing continuation cannot turn partial results into a complete set."""
        client = FakeClient([page(3, [record('1')])])
        with self.assertRaisesRegex(ValueError, 'before the reported hit count'):
            search.retrieve_topic(client, self.topic, self.root)

    def test_missing_abstract_is_not_an_assessment(self) -> None:
        """Metadata retrieval can finish while the scientific record remains unassessed."""
        client = FakeClient([page(1, [record('1', None)])])
        records, result = search.retrieve_topic(client, self.topic, self.root)
        merged = search.merge_records([('synthetic_topic', records)])
        self.assertEqual(result['verdict'], 'PASS')
        self.assertEqual(merged[0]['abstract'], '')
        self.assertEqual(merged[0]['assessment_status'], 'UNASSESSED')
        self.assertTrue(merged[0]['targeted_assessment_reasons'])

    def test_cross_topic_deduplication_preserves_provenance(self) -> None:
        """Duplicate PMID retrievals must retain both topics without doubling the record."""
        merged = search.merge_records([('first', [record('1', None)]),
                                       ('second', [record('1', 'Available abstract.')])])
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0]['topics'], ['first', 'second'])
        self.assertEqual(merged[0]['abstract'], 'Available abstract.')
        self.assertEqual(merged[0]['assessment_status'], 'UNASSESSED')


if __name__ == '__main__':
    unittest.main()
