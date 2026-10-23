"""Reject tracked runtime data, binary artifacts and common credential patterns."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_ROOTS = {'data', 'raw', 'derived', 'inputs', 'outputs', 'out', 'reports', 'tracking',
                   'submission', 'research_notes', '.venv', '.claude', '.codex', 'aws'}
FORBIDDEN_SUFFIXES = {'.pdf', '.xlsx', '.html', '.jsonl', '.log', '.bam', '.cram', '.vcf', '.bcf',
                      '.fastq', '.fq', '.gz', '.zip', '.pem', '.key'}
PATTERNS = {
    'private_key': re.compile(r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----'),
    'aws_access_key': re.compile(r'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b'),
    'github_token': re.compile(r'\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,})\b'),
    'literal_secret': re.compile(r'''(?i)(?:api_key|secret_access_key|access_token)\s*[:=]\s*["'][A-Za-z0-9_/-]{24,}["']'''),
    'aws_account': re.compile(r'arn:aws:[^\s"\']*:\d{12}:'),
    'home_path': re.compile(r'/(?:home|Users)/[A-Za-z0-9_.-]+/'),
}


def check_path(name: str) -> list[str]:
    """Check one tracked path without printing its contents.

    Args:
        name: A repository-relative path.

    Returns:
        The failed check names.
    """
    path = Path(name)
    failures = []
    if path.parts[0] in FORBIDDEN_ROOTS or path.suffix in FORBIDDEN_SUFFIXES:
        failures.append('excluded_file_type_or_directory')
    if path.name == 'SESSION_SUMMARY.md' or (path.name.startswith('.env') and path.name != '.env.example'):
        failures.append('local_configuration')
    full = ROOT / path
    if full.is_symlink():
        return failures + ['symlink']
    try:
        text = full.read_text()
    except (UnicodeError, OSError):
        return failures + ['not_readable_text']
    if '\0' in text:
        failures.append('binary_content')
    for label, pattern in PATTERNS.items():
        if pattern.search(text):
            failures.append(label)
    return failures


def main() -> None:
    """Check all tracked files and exit nonzero when a file fails."""
    result = subprocess.run(['git', 'ls-files', '-z'], cwd=ROOT, capture_output=True, check=True)
    names = result.stdout.decode().strip('\0').split('\0')
    failures = []
    for name in names:
        for label in check_path(name):
            failures.append(f'{name}: {label}')
    if failures:
        print('\n'.join(failures))
        raise SystemExit(1)
    print(f'PASS: {len(names)} tracked source files meet the repository scope checks.')


if __name__ == '__main__':
    main()
