#!/usr/bin/env python3
"""Refuse selected shell operations on restricted paths. Data class A.

This optional command guard implements the public dispatcher's conservative path rules.
It is separate from the competition permission for qualifying interactive subscriptions.
The guard checks explicit paths and expands glob patterns before it returns a verdict.
It examines directory entries without opening restricted file contents.

A PreToolUse hook can pass a command in a JSON payload. A script can call decide
before executing a supplied command. The result contains a verdict, path, and reason.

Run:
    uv run python harness/class_b2_command_guard.py --command "grep -c x writeups/a/_class_b2_r.qmd"
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: A file whose base name starts with one of these strings carries its data class in its name.
CLASS_B2_BASENAME_PREFIXES = ("_class_b2_", "class_b2_")

#: A file with one of these base names is class B2 wherever it sits.
CLASS_B2_BASENAMES = (
    "allele_characterization_v2.md",
    "class_b2_phase_rows.md",
    "slot1_b0.csv",
)

#: A path that holds one of these directory sequences is class B2 or class C.
CLASS_B2_DIRECTORY_PARTS = (
    ("derived", "llm_safe"),
    ("data",),
    ("raw",),
)

#: A file with one of these suffixes carries the genome and is class C.
CLASS_C_SUFFIXES = (".vcf", ".vcf.gz", ".bam", ".cram", ".fastq", ".fastq.gz", ".fq", ".fq.gz")

#: The readers that this control approves. Each one prints a verdict, a path or a public stub, and
#: no content of a class B2 file.
APPROVED_READERS = (
    "harness/privacy_guard.py",
    "harness/class_b2_command_guard.py",
    "scripts/list_class_b2_prefix.sh",
    "scripts/read_class_a_gate_file.py",
    "scripts/check_report_claims.py",
    "scripts/stage_jury_review.sh",
    "writeups/build_public.sh",
)

#: A reader that walks a directory tree. A directory that holds a class B2 file is refused for one
#: of these, because the reader reaches the file without naming it.
RECURSIVE_READERS = ("grep", "rg", "ag", "ack", "find", "fgrep", "egrep")

#: A character that makes a token a shell glob pattern. The shell expands the token before the
#: reader runs, so the command reaches a class B2 file without naming it.
GLOB_METACHARACTERS = ("*", "?", "[")

#: A separator that ends one command segment and starts the next one.
SEGMENT_SEPARATORS = (";", "&&", "||", "|", "&")

ALLOW = "ALLOW"
REFUSE = "REFUSE"


def split_segments(tokens: list[str]) -> list[list[str]]:
    """Split a token list into the command segments of a shell line.

    Args:
        tokens: The tokens of one shell line.

    Returns:
        One token list for each segment, with the separators removed.
    """
    segments: list[list[str]] = [[]]
    for token in tokens:
        if token in SEGMENT_SEPARATORS:
            segments.append([])
            continue
        segments[-1].append(token)
    return [segment for segment in segments if segment]


def executable_of(segment: list[str]) -> str:
    """Return the executable of one command segment.

    Args:
        segment: The tokens of one segment.

    Returns:
        The first token that is not an environment assignment, or an empty string.
    """
    for token in segment:
        if "=" in token and not token.startswith("-") and "/" not in token.split("=")[0]:
            continue
        return token
    return ""


def is_approved_reader(segment: list[str]) -> bool:
    """State whether a segment runs a reader that this control approves.

    Args:
        segment: The tokens of one segment.

    Returns:
        True when the segment runs an approved reader.
    """
    text = " ".join(segment)
    return any(reader in text for reader in APPROVED_READERS)


def is_recursive_reader(segment: list[str]) -> bool:
    """State whether a segment runs a reader that walks a directory tree.

    Args:
        segment: The tokens of one segment.

    Returns:
        True when the segment names a recursive reader with a recursive option.
    """
    name = Path(executable_of(segment)).name
    if name in ("find", "rg", "ag", "ack"):
        return True
    if name not in RECURSIVE_READERS:
        return False
    return any(token in ("-r", "-R", "--recursive", "-rn", "-rl", "-ri") or
               (token.startswith("-") and not token.startswith("--") and
                ("r" in token[1:] or "R" in token[1:]))
               for token in segment[1:])


def is_class_b2_path(candidate: str) -> bool:
    """State whether a path names a class B2 or a class C file by its own name.

    Args:
        candidate: A path as the command wrote it.

    Returns:
        True when the path carries a class B2 or class C name.
    """
    cleaned = candidate.split("://")[-1]
    pure = Path(cleaned)
    name = pure.name
    if any(name.startswith(prefix) for prefix in CLASS_B2_BASENAME_PREFIXES):
        return True
    if name in CLASS_B2_BASENAMES:
        return True
    if any(name.endswith(suffix) for suffix in CLASS_C_SUFFIXES):
        return True
    parts = tuple(part for part in pure.parts if part not in ("/", "."))
    for sequence in CLASS_B2_DIRECTORY_PARTS:
        width = len(sequence)
        for start in range(len(parts) - width + 1):
            if parts[start:start + width] == sequence:
                return True
    return False


def is_glob_token(token: str) -> bool:
    """State whether a token is a shell glob pattern.

    Args:
        token: One token of a command segment.

    Returns:
        True when the token holds a glob metacharacter.
    """
    return any(character in token for character in GLOB_METACHARACTERS)


def glob_match_class_b2(token: str, root: Path, recursive: bool = False) -> str:
    """Return the first class B2 path that a glob token expands to.

    The shell expands a glob before the reader runs, so a pattern such as
    ``writeups/track1_methods/*.qmd`` names no class B2 file and still reaches one. This function
    expands the pattern against directory entries without opening files.
    A directory of the expansion is walked for a recursive reader only, because such a
    reader reaches every file under the directory.

    Args:
        token: One token of a command segment.
        root: The directory that a relative pattern expands against.
        recursive: True when the segment runs a reader that walks a directory tree.

    Returns:
        The first class B2 path of the expansion, or an empty string.
    """
    if not is_glob_token(token):
        return ""
    cleaned = token.split("://")[-1]
    if cleaned.startswith("/"):
        base = Path("/")
        pattern = cleaned.lstrip("/")
    else:
        base = root
        pattern = cleaned
    try:
        matches = sorted(base.glob(pattern))
    except (IndexError, NotImplementedError, OSError, ValueError):
        return ""
    for match in matches:
        if is_class_b2_path(str(match)):
            return str(match)
        if not recursive:
            continue
        found = directory_holds_class_b2(match)
        if found:
            return found
    return ""


def is_path_candidate(token: str, root: Path) -> bool:
    """State whether a token is a path and not a pattern.

    Args:
        token: One token of a command segment.
        root: The directory that a relative path resolves against.

    Returns:
        True when the token names a path.
    """
    if token.startswith("-"):
        return False
    if os.sep in token or token.startswith("s3://"):
        return True
    return (root / token).exists()


def directory_holds_class_b2(directory: Path) -> str:
    """Return the first class B2 path under a directory, or an empty string.

    The walk reads directory entries and it opens no file.

    Args:
        directory: The directory to walk.

    Returns:
        The path of the first class B2 file that the walk finds, or an empty string.
    """
    if not directory.is_dir():
        return ""
    for current, _directories, files in os.walk(directory):
        for name in files:
            if is_class_b2_path(name):
                return str(Path(current) / name)
    return ""


def decide(command: str, root: Path | None = None) -> tuple[str, str]:
    """Return the verdict of control 9 on one shell command.

    Args:
        command: The command line as the caller would run it.
        root: The directory that a relative path resolves against. The repository root by default.

    Returns:
        The verdict ALLOW or REFUSE, and the reason.
    """
    base = root if root is not None else ROOT
    try:
        tokens = shlex.split(command)
    except ValueError:
        tokens = command.split()
    for segment in split_segments(tokens):
        if is_approved_reader(segment):
            continue
        recursive = is_recursive_reader(segment)
        for token in segment[1:]:
            if token.startswith("-"):
                continue
            expanded = glob_match_class_b2(token, base, recursive)
            if expanded:
                return REFUSE, (
                    f"{token} is a glob that expands to the class B2 path {expanded}, so the "
                    "command reaches that file without naming it. Control 9 of "
                    "docs/llm_data_privacy.md refuses it."
                )
            if not is_path_candidate(token, base):
                continue
            if is_class_b2_path(token):
                return REFUSE, (
                    f"{token} is a class B2 or class C path and the command is not an approved "
                    "reader. Control 9 of docs/llm_data_privacy.md refuses it."
                )
            if not recursive:
                continue
            found = directory_holds_class_b2(base / token)
            if found:
                return REFUSE, (
                    f"{token} is a directory that holds the class B2 path {found}, and a recursive "
                    "reader reaches it without naming it. Control 9 of docs/llm_data_privacy.md "
                    "refuses it."
                )
    return ALLOW, "No token of the command names a class B2 or a class C path."


def read_hook_command(stream: object) -> str:
    """Read the command of a PreToolUse hook payload.

    Args:
        stream: The stream that carries the JSON payload.

    Returns:
        The command string, or an empty string when the payload carries none.
    """
    try:
        payload = json.load(stream)
    except (json.JSONDecodeError, ValueError):
        return ""
    tool_input = payload.get("tool_input") or {}
    return str(tool_input.get("command") or "")


def main(argv: list[str] | None = None) -> int:
    """Run the guard on one command and print the verdict.

    Args:
        argv: The argument list. ``sys.argv[1:]`` by default.

    Returns:
        0 when the command passes, 2 when the control refuses it.
    """
    parser = argparse.ArgumentParser(description="Control 9 of docs/llm_data_privacy.md.")
    parser.add_argument("--command", help="The command line to check.")
    parser.add_argument("--hook", action="store_true",
                        help="Read a PreToolUse hook payload from standard input.")
    parser.add_argument("--root", default=str(ROOT), help="The directory that a path resolves against.")
    arguments = parser.parse_args(argv)
    command = arguments.command or ""
    if arguments.hook:
        command = read_hook_command(sys.stdin)
    if not command:
        print("ALLOW: the payload carries no command.")
        return 0
    verdict, reason = decide(command, Path(arguments.root))
    if verdict == REFUSE:
        print(f"REFUSE: {reason}", file=sys.stderr)
        return 2
    print(f"ALLOW: {reason}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
