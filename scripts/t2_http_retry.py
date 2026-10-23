#!/usr/bin/env python3
"""Send one HTTP GET and repeat it when the service asks the caller to wait.

The NCBI E-utilities return HTTP 429 when a host sends more than three requests per second without
an API key. A retrieval of this project can run at the same time as a literature query, so a single
GET has to survive that answer instead of ending the run. This module holds the one function that
both literature scripts call.
"""

from __future__ import annotations

import time

import httpx

RETRY_STATUS = (429, 500, 502, 503, 504)
DEFAULT_ATTEMPTS = 6
DEFAULT_PAUSE = 4.0


def get(url: str, params: dict[str, str], timeout: float = 90.0,
        attempts: int = DEFAULT_ATTEMPTS, pause: float = DEFAULT_PAUSE) -> httpx.Response:
    """Return the response of one GET, after up to `attempts` tries.

    Args:
        url: The request URL.
        params: The query parameters.
        timeout: The read timeout in seconds.
        attempts: The number of tries before the function raises.
        pause: The first wait in seconds. Each further wait is twice the one before it.

    Returns:
        The response with a status below 400.

    Raises:
        httpx.HTTPStatusError: When every try returns a status that this module retries.
    """
    wait = pause
    response = None
    for index in range(attempts):
        response = httpx.get(url, params=params, timeout=timeout, follow_redirects=True)
        if response.status_code not in RETRY_STATUS:
            response.raise_for_status()
            return response
        if index + 1 < attempts:
            time.sleep(wait)
            wait *= 2
    assert response is not None
    response.raise_for_status()
    return response
