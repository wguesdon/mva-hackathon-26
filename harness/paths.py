"""Locate project files and the local environment file of the MVA harness.

The harness runs on a machine that holds no patient data. Patient data lives only in the
dedicated S3 bucket. This module never resolves a path to that data.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Credential and private identifier files that the guards read only to keep their values out of logs.
ENV_FILES = (
    ROOT / ".env",
)

# Directories that must never be sent to a model, whatever the arm.
PRIVATE_DIRECTORIES = ("data", "raw", "reference", "reports/llm_safe")


def read_env_file(path: Path) -> dict[str, str]:
    """Read simple KEY=VALUE lines from an environment file.

    Args:
        path: The environment file.

    Returns:
        The parsed values, or an empty mapping when the file is absent.
    """
    values: dict[str, str] = {}
    if not path.is_file():
        return values
    for line in path.read_text().splitlines():
        if "=" in line and not line.lstrip().startswith("#"):
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip().strip("\"'")
    return values


def local_env() -> dict[str, str]:
    """Return the values of the project `.env` file.

    Returns:
        The parsed values of `ROOT/.env`.
    """
    return read_env_file(ROOT / ".env")
