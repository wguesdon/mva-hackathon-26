"""Resolve a disease identifier into the row that the Track 2 pipeline needs.

`scripts/track2_transfer_pipeline.py` reads a disease row that holds an OMIM identifier, a causal
gene with its NCBI Gene identifier, a MeSH identifier, a PubTator3 accession and a graph node
identifier. A person wrote each row of `literature/track2_transfer_diseases.yaml` by hand, and the
second external jury review named that hand step as the first reason that the method does not yet
transfer without a person. This module removes the step.

The resolution chain is fixed by ADR 0009 and no step of it writes an identifier by hand.

| Field | Source |
|---|---|
| The causal gene and its NCBI Gene identifier | Human Phenotype Ontology `genes_to_disease` |
| The disease name | Human Phenotype Ontology `phenotype.hpoa` |
| The PubTator3 accession and the MeSH identifier | PubTator3 entity autocomplete |

The module refuses to guess. When the autocomplete endpoint returns no disease hit, or returns two
hits that the query cannot separate, `resolve` raises `ResolutionError` with the identifier that a
person must supply instead.

Usage:
    uv run python scripts/track2_disease_resolver.py --disease OMIM:268300
    uv run python scripts/track2_disease_resolver.py --disease MESH:C535687
    uv run python scripts/track2_disease_resolver.py --disease "Roberts syndrome"
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import re
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track2_common as common  # noqa: E402

AUTOCOMPLETE_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/entity/autocomplete/"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
DEFAULT_HPO_DIR = common.ROOT / "tracking" / "cache" / "hpo_v2026_09_01"
DEFAULT_ROW_DIR = common.ROOT / "literature" / "track2_run_disease"
HPO_RELEASE = "2026-09-01"
HPO_SOURCE = "https://github.com/obophenotype/human-phenotype-ontology/releases/tag/v2026-09-01"

# The association types of `genes_to_disease` that name a causal gene of a monogenic disease.
CAUSAL_ASSOCIATIONS = ("MENDELIAN",)

# The number of autocomplete hits that the resolver asks for.
AUTOCOMPLETE_LIMIT = 10

_OMIM = re.compile(r"^(?:OMIM:)?(\d{6})$", re.IGNORECASE)
_MESH = re.compile(r"^(?:MESH:)?([CD]\d{6})$", re.IGNORECASE)
_SLUG_STRIP = re.compile(r"[^a-z0-9]+")


class ResolutionError(Exception):
    """A disease identifier did not resolve to exactly one row."""


def slugify(label: str) -> str:
    """Build the directory slug of a disease label.

    Args:
        label: The disease name, for example `Roberts-SC phocomelia syndrome`.

    Returns:
        The slug, for example `roberts_sc_phocomelia_syndrome`.
    """
    return _SLUG_STRIP.sub("_", label.strip().lower()).strip("_")


def classify(identifier: str) -> tuple[str, str]:
    """Return the kind of an identifier and its normalized body.

    Args:
        identifier: The value that the command line gave.

    Returns:
        The kind, which is `omim`, `mesh` or `name`, and the normalized body.
    """
    text = identifier.strip()
    omim = _OMIM.match(text)
    if omim:
        return "omim", f"OMIM:{omim.group(1)}"
    mesh = _MESH.match(text)
    if mesh:
        return "mesh", mesh.group(1).upper()
    return "name", text


def read_gene_rows(path: Path, omim: str) -> list[dict]:
    """Read the causal gene rows of one OMIM identifier.

    Args:
        path: The `genes_to_disease.txt` file of the pinned release.
        omim: The OMIM identifier, for example `OMIM:268300`.

    Returns:
        One record per causal gene, with the symbol, the NCBI Gene identifier, the association type
        and the source column of the row.

    Raises:
        ResolutionError: When the file is missing.
    """
    if not path.is_file():
        raise ResolutionError(
            f"the Human Phenotype Ontology file {path} is missing. "
            f"Run `make track2-resources` to download release {HPO_RELEASE}.")
    rows: list[dict] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("disease_id") != omim:
                continue
            if row.get("association_type") not in CAUSAL_ASSOCIATIONS:
                continue
            rows.append({"symbol": row["gene_symbol"], "ncbi_gene": row["ncbi_gene_id"],
                         "association_type": row["association_type"], "source": row["source"]})
    return rows


def read_disease_name(path: Path, omim: str) -> str:
    """Read the disease name of one OMIM identifier.

    Args:
        path: The `phenotype.hpoa` file of the pinned release.
        omim: The OMIM identifier.

    Returns:
        The disease name, or an empty string when the file holds no row for the identifier.

    Raises:
        ResolutionError: When the file is missing.
    """
    if not path.is_file():
        raise ResolutionError(
            f"the Human Phenotype Ontology file {path} is missing. "
            f"Run `make track2-resources` to download release {HPO_RELEASE}.")
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("database_id"):
                continue
            fields = line.split("\t")
            if len(fields) > 1 and fields[0] == omim:
                return fields[1].strip()
    return ""


def autocomplete_disease(client, query: str) -> list[dict]:
    """Return the disease hits of the PubTator3 entity autocomplete endpoint.

    Args:
        client: The HTTP client.
        query: The disease name or the MeSH identifier.

    Returns:
        The hits whose biotype is `disease`, in endpoint order.
    """
    payload = common.get_json(client, AUTOCOMPLETE_URL,
                              {"query": query, "limit": AUTOCOMPLETE_LIMIT})
    hits = payload if isinstance(payload, list) else []
    return [hit for hit in hits if hit.get("biotype") == "disease"]


def mesh_entry_terms(client, mesh_id: str) -> list[str]:
    """Return the entry terms of a MeSH concept.

    Args:
        client: The HTTP client.
        mesh_id: The MeSH identifier, for example `C535687`.

    Returns:
        The entry terms in the order that the record gives them.

    Raises:
        ResolutionError: When the record holds no entry term.
    """
    uid = mesh_uid_of(mesh_id)
    if not uid:
        raise ResolutionError(
            f"{mesh_id} is neither a MeSH descriptor nor a MeSH supplementary concept record.")
    payload = common.get_json(client, ESUMMARY_URL,
                              {"db": "mesh", "id": uid, "retmode": "json"})
    record = (payload.get("result") or {}).get(uid) or {}
    terms = [str(term).strip() for term in record.get("ds_meshterms", []) if str(term).strip()]
    if not terms:
        raise ResolutionError(
            f"the MeSH record {mesh_id} returned no entry term, so the resolver has no name to "
            "query PubTator3 with. Give the disease name instead.")
    return terms


def mesh_concept_name(client, mesh_id: str) -> str:
    """Return the preferred term of a MeSH concept.

    The PubTator3 autocomplete endpoint searches names and not identifiers, so a MeSH identifier
    needs its name before the endpoint can answer.

    Args:
        client: The HTTP client.
        mesh_id: The MeSH identifier, for example `C535687`.

    Returns:
        The first entry term of the record.

    Raises:
        ResolutionError: When the record holds no entry term.
    """
    return mesh_entry_terms(client, mesh_id)[0]


def pick_disease_hit(hits: list[dict], query: str, mesh_id: str = "") -> dict:
    """Choose the one disease hit that the query names.

    Args:
        hits: The disease hits of the autocomplete endpoint.
        query: The query that produced them.
        mesh_id: The MeSH identifier that the caller asked for, when it has one.

    Returns:
        The chosen hit.

    Raises:
        ResolutionError: When no hit matches, or when two hits match equally well.
    """
    if mesh_id:
        exact = [hit for hit in hits if str(hit.get("db_id", "")).upper() == mesh_id]
        if len(exact) == 1:
            return exact[0]
        raise ResolutionError(
            f"the PubTator3 autocomplete endpoint returned {len(exact)} entities with the MeSH "
            f"identifier {mesh_id} for the query {query!r}. Give the disease name instead.")
    if not hits:
        raise ResolutionError(
            f"the PubTator3 autocomplete endpoint returned no disease for the query {query!r}. "
            "PubTator3 recognises no entity for this disease, so the method cannot build a "
            "literature layer for it. Give the name of a MeSH concept that PubTator3 recognises.")
    wanted = query.strip().lower()
    exact = [hit for hit in hits if str(hit.get("name", "")).strip().lower() == wanted]
    if len(exact) == 1:
        return exact[0]
    if len(hits) == 1:
        return hits[0]
    names = ", ".join(f"{hit.get('name')} ({hit.get('db_id')})" for hit in hits[:5])
    raise ResolutionError(
        f"the PubTator3 autocomplete endpoint returned {len(hits)} diseases for the query "
        f"{query!r} and none of them carries that exact name: {names}. "
        "Give the MeSH identifier of the one you want, for example MESH:C535687.")


def resolve(identifier: str, hpo_dir: Path = DEFAULT_HPO_DIR) -> dict:
    """Resolve a disease identifier into the pipeline row.

    Args:
        identifier: An OMIM identifier, a MeSH identifier or a disease name.
        hpo_dir: The directory of the pinned Human Phenotype Ontology release.

    Returns:
        The disease row, with a `resolution` block that names the source of every field.

    Raises:
        ResolutionError: When the identifier does not resolve to exactly one disease with at least
            one causal gene.
    """
    kind, body = classify(identifier)
    genes_path = hpo_dir / "genes_to_disease.txt"
    hpoa_path = hpo_dir / "phenotype.hpoa"
    omim = body if kind == "omim" else ""
    label = body if kind == "name" else ""
    mesh_id = body if kind == "mesh" else ""

    with common.make_client() as client:
        if kind == "mesh":
            concept = mesh_concept_name(client, mesh_id)
            hit = pick_disease_hit(autocomplete_disease(client, concept), concept, mesh_id)
            label = str(hit.get("name", "")).strip()
        else:
            if kind == "omim":
                label = read_disease_name(hpoa_path, omim)
                if not label:
                    raise ResolutionError(
                        f"{omim} has no row in {hpoa_path}. The Human Phenotype Ontology release "
                        f"{HPO_RELEASE} annotates no such disease, so the resolver has no name to "
                        "query PubTator3 with. Give the disease name instead.")
            hits = autocomplete_disease(client, label)
            if not hits and "-" in label:
                hits = autocomplete_disease(client, label.split("-")[0].strip())
            hit = pick_disease_hit(hits, label)

        accession = str(hit.get("_id", ""))
        database = str(hit.get("db", ""))
        mesh_id = str(hit.get("db_id", "")).upper()
        if database != "ncbi_mesh" or not mesh_id:
            raise ResolutionError(
                f"the PubTator3 entity {accession or identifier!r} resolves to the database "
                f"{database or 'none'} and not to MeSH. The mention rule reads MeSH entry terms, "
                "so the method cannot filter this disease node. Give a disease with a MeSH "
                "concept.")
        pubtator_name = str(hit.get("name", "")).strip() or label
        entry_terms = mesh_entry_terms(client, mesh_id)

    node = f"DISEASE:MESH:{mesh_id}"
    if not omim:
        omim = match_omim_by_name(hpoa_path, pubtator_name, label, *entry_terms)
        if not omim:
            raise ResolutionError(
                f"no disease of the Human Phenotype Ontology release {HPO_RELEASE} carries the "
                f"name {pubtator_name!r} or one of the {len(entry_terms)} entry terms of the MeSH "
                f"record {mesh_id}, so the resolver has no causal gene. Give the OMIM identifier "
                "instead, for example OMIM:268300.")

    genes = read_gene_rows(genes_path, omim)
    if not genes:
        raise ResolutionError(
            f"{omim} has no MENDELIAN gene row in {genes_path}. The method seeds the graph from a "
            "causal gene, so it cannot run on a disease with no such row.")

    hpo_label = read_disease_name(hpoa_path, omim) or pubtator_name
    return {
        "slug": slugify(hpo_label),
        "label": hpo_label,
        "pubtator_label": pubtator_name,
        "omim": omim,
        "genes": genes,
        "source": genes[0]["source"],
        "mesh_id": mesh_id,
        "mesh_uid": mesh_uid_of(mesh_id),
        "node": node,
        "accession": accession,
        "regime_expectation": "unknown",
        "reason": f"Resolved from the identifier {identifier} by scripts/track2_disease_resolver.py.",
        "resolution": {
            "identifier_given": identifier,
            "identifier_kind": kind,
            "resolved_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "hpo_release": HPO_RELEASE,
            "hpo_source": HPO_SOURCE,
            "genes_to_disease": str(genes_path.relative_to(common.ROOT)),
            "genes_to_disease_sha256": common.sha256_of(genes_path),
            "phenotype_hpoa": str(hpoa_path.relative_to(common.ROOT)),
            "phenotype_hpoa_sha256": common.sha256_of(hpoa_path),
            "autocomplete_endpoint": AUTOCOMPLETE_URL,
            "autocomplete_hit": {"_id": accession, "db": database, "db_id": mesh_id,
                                 "name": pubtator_name},
            "hand_written_fields": [],
        },
    }


def match_omim_by_name(hpoa_path: Path, *names: str) -> str:
    """Return the OMIM identifier whose disease name equals one of the given names.

    The comparison reads case and punctuation as noise, because the Human Phenotype Ontology writes
    "Roberts-SC phocomelia syndrome" where MeSH writes "Roberts-Sc Phocomelia Syndrome".

    Args:
        hpoa_path: The `phenotype.hpoa` file of the pinned release.
        *names: The names to try, in order of preference.

    Returns:
        The OMIM identifier of the first name that matches, or an empty string when none does.
    """
    wanted = [slugify(name) for name in names if name and name.strip()]
    if not wanted:
        return ""
    by_key: dict[str, str] = {}
    with hpoa_path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("database_id"):
                continue
            fields = line.split("\t")
            if len(fields) < 2 or not fields[0].startswith("OMIM:"):
                continue
            by_key.setdefault(slugify(fields[1]), fields[0])
    for key in wanted:
        if key in by_key:
            return by_key[key]
    return ""


def mesh_uid_of(mesh_id: str) -> str:
    """Return the NCBI unique identifier of a MeSH identifier.

    Args:
        mesh_id: The MeSH identifier, for example `C535687`.

    Returns:
        The unique identifier, or an empty string when the prefix is neither `C` nor `D`.
    """
    body = mesh_id.split(":")[-1].strip().upper()
    if body.startswith("D"):
        return f"68{body[1:]}"
    if body.startswith("C"):
        return f"67{body[1:]}"
    return ""


def write_row(row: dict, directory: Path = DEFAULT_ROW_DIR) -> Path:
    """Write the resolved row to a file that a reader checks against the sources.

    Args:
        row: The resolved disease row.
        directory: The output directory.

    Returns:
        The written path.
    """
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{row['slug']}.yaml"
    header = ("# The disease row that scripts/track2_disease_resolver.py resolved. Every field "
              "carries its\n# source in the `resolution` block. No field of this file is written "
              "by hand.\n")
    path.write_text(header + yaml.safe_dump({"diseases": [row]}, sort_keys=False), encoding="utf-8")
    return path


def main() -> int:
    """Resolve one disease identifier from the command line.

    Returns:
        0 on success, 3 when the identifier does not resolve.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--disease", required=True,
                        help="An OMIM identifier, a MeSH identifier or a disease name.")
    parser.add_argument("--hpo-dir", type=Path, default=DEFAULT_HPO_DIR,
                        help="The directory of the pinned Human Phenotype Ontology release.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_ROW_DIR,
                        help="The directory that receives the resolved row.")
    parser.add_argument("--no-write", action="store_true", help="Print the row and write no file.")
    args = parser.parse_args()

    try:
        row = resolve(args.disease, args.hpo_dir)
    except ResolutionError as error:
        print(f"FAIL: {error}")
        return 3
    print(f"disease: {row['label']} ({row['omim']})")
    print(f"  genes: {', '.join(gene['symbol'] for gene in row['genes'])}")
    print(f"  mesh: {row['mesh_id']}  node: {row['node']}  accession: {row['accession']}")
    if not args.no_write:
        print(f"  row: {write_row(row, args.out_dir)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
