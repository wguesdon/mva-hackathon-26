#!/usr/bin/env bash
# Download the pinned release files that scripts/track2_run_disease.py reads from disk.
#
# The runner refuses to start when a file is absent or when its SHA-256 differs from the value that
# literature/track2_runner_resources.yaml pins. This script downloads each file and checks its
# digest. A digest that does not match is a failure and the script stops.
#
# Usage:
#     scripts/track2_fetch_resources.sh

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RELEASE="v2026-09-01"
BASE="https://github.com/obophenotype/human-phenotype-ontology/releases/download/${RELEASE}"
TARGET="${ROOT}/tracking/cache/hpo_v2026_09_01"

declare -A EXPECTED=(
  ["genes_to_disease.txt"]="7937a1b92779dc8614e39ace623b443c7c71dd15e4de11d5e52b14e331ab5978"
  ["phenotype.hpoa"]="e89aa39c8f97bf5a52c5d160f9250a603681d630ade8ec2679d84ac5aece1f72"
)

mkdir -p "${TARGET}"

for name in "${!EXPECTED[@]}"; do
  path="${TARGET}/${name}"
  if [[ -f "${path}" ]]; then
    have="$(sha256sum "${path}" | cut -d' ' -f1)"
    if [[ "${have}" == "${EXPECTED[$name]}" ]]; then
      echo "ok   ${name}: release ${RELEASE} is already on disk."
      continue
    fi
    echo "warn ${name}: the file on disk has the digest ${have:0:16}. Downloading again."
  fi
  echo "get  ${name} from ${BASE}/${name}"
  curl --fail --location --silent --show-error --output "${path}" "${BASE}/${name}"
  have="$(sha256sum "${path}" | cut -d' ' -f1)"
  if [[ "${have}" != "${EXPECTED[$name]}" ]]; then
    echo "FAIL ${name}: the download has the SHA-256 ${have} and release ${RELEASE} has ${EXPECTED[$name]}."
    echo "     The Human Phenotype Ontology replaced the file of this release, or the download is incomplete."
    exit 1
  fi
  echo "ok   ${name}: ${have:0:16}"
done

echo "Every pinned release file of release ${RELEASE} is on disk with the recorded digest."
