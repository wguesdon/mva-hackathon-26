"""Copy stdin to stdout with private identifiers from PRIVACY_DENY_TOKENS replaced by [REDACTED].

Use it on every command output about a patient run that reaches an LLM session, for example:
    aws s3 cp s3://.../head/nextflow-head_failed.log - | uv run python scripts/redact_stream.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "harness"))
import privacy_guard  # noqa: E402


def main() -> int:
    """Redact the stream line by line.

    Returns:
        Zero.
    """
    tokens = privacy_guard.deny_tokens()
    if not tokens:
        print("[redact] PRIVACY_DENY_TOKENS is empty; refusing to pass the stream.", file=sys.stderr)
        return 2
    pattern = re.compile("|".join(re.escape(token) for token in tokens), re.IGNORECASE)
    for line in sys.stdin:
        sys.stdout.write(pattern.sub("[REDACTED]", line))
        sys.stdout.flush()
    return 0


if __name__ == "__main__":
    sys.exit(main())
