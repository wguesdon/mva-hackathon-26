"""Shared helpers for the Track 2 knowledge graph fetch and build scripts.

Every function here works on public class A data only: gene symbols, MeSH identifiers,
ChEMBL identifiers, Reactome identifiers and PMIDs. No patient data reaches this module.

The PubTator3 API names an entity with an accession such as `@GENE_BUB1B` or
`@DISEASE_Mosaic_variegated_aneuploidy_syndrome`. The helpers below split that accession into
its biotype and its body, and map it to the stable node identifier used by the graph contract.
"""

from __future__ import annotations

import atexit
import hashlib
import json
import math
import os
import re
import time
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit

import httpx

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "reports" / "track2" / "raw"
GRAPH_DIR = ROOT / "reports" / "track2" / "graph"

USER_AGENT = "mva-hackathon-2026-track2/1.0 (public data only)"

# The environment variable that names the file where `dump_counters` writes the API call counter.
CALL_LOG_VARIABLE = "TRACK2_CALL_LOG"

# The API call counter of this process, keyed by host.
_CALL_COUNTS: dict[str, dict[str, int]] = {}
DEFAULT_TIMEOUT = 60.0
DEFAULT_PAUSE_S = 0.34
MAX_ATTEMPTS = 6

# NCBI answers a burst with HTTP 429. The first wait after that status, in seconds.
THROTTLE_BASE_S = 8.0

# NCBI documents a limit of three requests in one second for a client with no API key, and
# `DEFAULT_PAUSE_S` paces this pipeline at that ceiling with no margin. A second process on the same
# machine shares the limit, so one process alone cannot stay under it. A throttle status therefore
# raises a pause that every later call of the process pays, and a run of answered calls lowers it
# again. The status codes below are the two that NCBI returns when it refuses a burst.
THROTTLE_STATUS = (429, 503)
THROTTLE_PAUSE_STEP_S = 0.5
THROTTLE_PAUSE_MAX_S = 8.0
THROTTLE_PAUSE_DECAY = 0.98

# A throttle status is not a failure of the request, so it carries a larger attempt budget than a
# malformed response or a dropped connection.
MAX_THROTTLE_ATTEMPTS = 14

# The pause that every call of this process pays on top of the pace of its caller.
_EXTRA_PAUSE_S = 0.0

# PubTator3 biotype -> graph node kind.
KIND_BY_BIOTYPE = {
    "GENE": "gene",
    "DISEASE": "disease",
    "CHEMICAL": "chemical",
    "SPECIES": "species",
    "VARIANT": "variant",
}

# Graph node kind -> node identifier prefix.
PREFIX_BY_KIND = {
    "gene": "GENE",
    "disease": "DISEASE",
    "chemical": "CHEMICAL",
    "pathway": "PATHWAY",
}

# PubTator3 database name -> identifier namespace used inside a node identifier.
NAMESPACE_BY_DB = {
    "ncbi_mesh": "MESH",
    "ncbi_gene": "NCBIGene",
    "ncbi_taxonomy": "NCBITaxon",
    "litvar": "LITVAR",
}

_ACCESSION = re.compile(r"^@([A-Z]+)_(.+)$")
_HIGHLIGHT_MARKED = re.compile(r"@<m>([A-Z]+_[^<]*)</m>\s*")
_HIGHLIGHT_TOKEN = re.compile(r"@[A-Z]+_\S+\s*")
_HIGHLIGHT_SURFACE = re.compile(r"@@@(.*?)@@@", re.DOTALL)
_HIGHLIGHT_PAIR = re.compile(r"@<m>([A-Z]+_[^<]*)</m>\s*@[A-Z]+_(\S+)")
_TAG = re.compile(r"</?m>")


def split_accession(accession: str) -> tuple[str, str]:
    """Split a PubTator3 accession into its biotype and its body.

    Args:
        accession: An accession such as `@GENE_BUB1B`.

    Returns:
        The upper case biotype and the body, for example `("GENE", "BUB1B")`.
        Both values are empty when the accession does not parse.
    """
    match = _ACCESSION.match(accession.strip())
    if not match:
        return "", ""
    return match.group(1), match.group(2)


def accession_kind(accession: str) -> str:
    """Return the graph node kind of a PubTator3 accession.

    Args:
        accession: An accession such as `@DISEASE_Aneuploidy`.

    Returns:
        The node kind, or an empty string when the biotype is unknown.
    """
    biotype, _ = split_accession(accession)
    return KIND_BY_BIOTYPE.get(biotype, "")


def accession_label(accession: str) -> str:
    """Return a human readable label for a PubTator3 accession.

    Args:
        accession: An accession such as `@DISEASE_Mosaic_variegated_aneuploidy_syndrome`.

    Returns:
        The body of the accession with underscores replaced by spaces.
    """
    _, body = split_accession(accession)
    return body.replace("_", " ")


def node_id_for_accession(accession: str, identifier: str = "", database: str = "") -> str:
    """Build the stable graph node identifier of a PubTator3 accession.

    A gene keeps its symbol because the symbol is stable and readable. A disease or a chemical
    uses its source database identifier when the entity endpoint resolved it. An unresolved
    entity falls back to the PubTator3 accession body, so the identifier is never invented.

    Args:
        accession: An accession such as `@CHEMICAL_Nocodazole`.
        identifier: The database identifier returned by the entity endpoint, for example `D015739`.
        database: The database name returned by the entity endpoint, for example `ncbi_mesh`.

    Returns:
        The node identifier, or an empty string when the accession does not parse.
    """
    kind = accession_kind(accession)
    prefix = PREFIX_BY_KIND.get(kind)
    if not prefix:
        return ""
    _, body = split_accession(accession)
    if kind == "gene":
        return f"{prefix}:{body}"
    if identifier:
        namespace = NAMESPACE_BY_DB.get(database, "PUBTATOR")
        if ":" in identifier:
            return f"{prefix}:{identifier}"
        return f"{prefix}:{namespace}:{identifier}"
    return f"{prefix}:PUBTATOR:{body}"


def plain_sentence(text_hl: str) -> str:
    """Strip the PubTator3 annotation markup from a highlighted passage.

    A highlighted passage marks each annotated entity as an accession, an identifier and the
    surface form, for example `@<m>GENE_BUB1B</m> @GENE_701 @@@BUB1B@@@`. This function keeps
    only the surface form, so the result reads as the sentence that the authors wrote.

    Args:
        text_hl: The `text_hl` field of a PubTator3 search result.

    Returns:
        The sentence without markup and without repeated whitespace.
    """
    if not text_hl:
        return ""
    text = _HIGHLIGHT_MARKED.sub("", text_hl)
    text = _HIGHLIGHT_TOKEN.sub("", text)
    text = _HIGHLIGHT_SURFACE.sub(r"\1", text)
    text = _TAG.sub("", text)
    return " ".join(text.split())


def highlight_identifiers(text_hl: str) -> dict[str, str]:
    """Read the database identifier of each marked entity in a highlighted passage.

    Args:
        text_hl: The `text_hl` field of a PubTator3 search result.

    Returns:
        A mapping of accession to database identifier, for example
        `{"@DISEASE_Mosaic_variegated_aneuploidy_syndrome": "MESH:C536987"}`.
    """
    found: dict[str, str] = {}
    for body, identifier in _HIGHLIGHT_PAIR.findall(text_hl or ""):
        found.setdefault(f"@{body}", identifier)
    return found


def edge_key(source_id: str, target_id: str, relation: str) -> tuple[str, str, str]:
    """Build the deduplication key of an edge.

    Args:
        source_id: The source node identifier.
        target_id: The target node identifier.
        relation: The relation type.

    Returns:
        The three values as a tuple.
    """
    return source_id, target_id, relation


def pair_query(accession_1: str, accession_2: str, relation: str = "ANY") -> str:
    """Build the PubTator3 search text that retrieves the publications of one entity pair.

    Args:
        accession_1: The first accession, for example `@GENE_BUB1B`.
        accession_2: The second accession.
        relation: The relation type, or `ANY` for every relation type.

    Returns:
        The value of the `text` query parameter.
    """
    return f"relations:{relation}|{accession_1}|{accession_2}"


def tier2_extra_chemicals(records_by_gene: dict[str, list[dict]], neighbours: list[str],
                          min_publications: int) -> set[str]:
    """Collect the chemical accessions that a derived neighbour gene links to strongly.

    A compound that a spindle assembly checkpoint neighbour is linked to is the repurposing
    signal that Track 2 looks for. The threshold removes the long tail of one-off mentions.

    Args:
        records_by_gene: The relation records of each gene, keyed by symbol.
        neighbours: The derived neighbour gene symbols.
        min_publications: The smallest publication count of a kept chemical relation.

    Returns:
        The chemical accessions to admit into the graph.
    """
    chemicals: set[str] = set()
    for symbol in neighbours:
        accession = f"@GENE_{symbol}"
        for record in records_by_gene.get(symbol, []):
            if int(record.get("publications", 0) or 0) < min_publications:
                continue
            source, target = record.get("source", ""), record.get("target", "")
            other = target if source == accession else source
            if other.startswith("@CHEMICAL_"):
                chemicals.add(other)
    return chemicals


def select_graph(records_by_gene: dict[str, list[dict]], seeds: list[str], neighbours: list[str],
                 tier2_chemical_min: int,
                 extra_chemicals: frozenset[str] = frozenset()) -> tuple[set[str], list[dict]]:
    """Choose the nodes and the relation records of the literature layer.

    The first tier keeps every relation of the seed genes, so the MVA neighbourhood is complete.
    The second tier keeps a relation of a derived neighbour gene only when both of its endpoints
    already belong to the node set. That rule closes the subgraph without letting a hub gene pull
    in the whole of PubMed.

    `extra_chemicals` admits a chemical accession at the second tier whatever its publication count
    is. ADR 0016 uses it for the members of a named chemical class, which the public chemical
    vocabulary enumerates. The default is empty, so a build with no option reproduces every earlier
    receipt.

    Args:
        records_by_gene: The relation records of each gene, keyed by symbol.
        seeds: The seed gene symbols.
        neighbours: The derived neighbour gene symbols.
        tier2_chemical_min: The smallest publication count of a chemical admitted at the second tier.
        extra_chemicals: The chemical accessions that enter whatever their publication count is.

    Returns:
        The node accessions and the selected relation records.
    """
    nodes: set[str] = set()
    edges: list[dict] = []
    seen: set[tuple[str, str, str]] = set()
    for symbol in seeds:
        for record in records_by_gene.get(symbol, []):
            source, target = record.get("source", ""), record.get("target", "")
            nodes.add(source)
            nodes.add(target)
            key = edge_key(source, target, record.get("type", ""))
            if key not in seen:
                seen.add(key)
                edges.append(record)
    nodes.update(f"@GENE_{symbol}" for symbol in neighbours)
    nodes.update(tier2_extra_chemicals(records_by_gene, neighbours, tier2_chemical_min))
    nodes.update(extra_chemicals)
    for symbol in neighbours:
        for record in records_by_gene.get(symbol, []):
            source, target = record.get("source", ""), record.get("target", "")
            if source not in nodes or target not in nodes:
                continue
            key = edge_key(source, target, record.get("type", ""))
            if key not in seen:
                seen.add(key)
                edges.append(record)
    return nodes, edges


def pair_total(records: list[dict]) -> dict[tuple[str, str], int]:
    """Sum the publication count of every relation type between the same two entities.

    The relations endpoint reports a count for each relation type. A specificity weight describes
    the entity pair, not the relation type, so the counts of a pair are added together first.

    Args:
        records: The selected relation records.

    Returns:
        The total co-mention count of each unordered entity pair.
    """
    totals: dict[tuple[str, str], int] = {}
    for record in records:
        source, target = record.get("source", ""), record.get("target", "")
        key = tuple(sorted((source, target)))
        totals[key] = totals.get(key, 0) + int(record.get("publications", 0) or 0)
    return totals


def npmi(pair_count: int, count_a: int, count_b: int, corpus_size: int) -> float | None:
    """Compute the normalized pointwise mutual information of one entity pair.

    The measure compares the observed co-mention count against the count that the two entity
    totals predict. A pair that co-occurs only because both entities are ubiquitous scores near
    zero. A pair that occurs together far more often than chance scores near one.

    The value is `log(p(a,b) / (p(a) p(b))) / -log(p(a,b))`, with each probability taken over the
    corpus size.

    Args:
        pair_count: The number of publications that mention both entities.
        count_a: The number of publications that mention the first entity.
        count_b: The number of publications that mention the second entity.
        corpus_size: The number of records in the corpus.

    Returns:
        The value in the range -1 to 1, or None when an input makes the measure undefined.
    """
    if pair_count <= 0 or count_a <= 0 or count_b <= 0 or corpus_size <= 0:
        return None
    if pair_count > count_a or pair_count > count_b or pair_count >= corpus_size:
        return None
    joint = pair_count / corpus_size
    expected = (count_a / corpus_size) * (count_b / corpus_size)
    if expected <= 0:
        return None
    return (math.log(joint / expected)) / (-math.log(joint))


def consequence_accessions(raw_dir: Path) -> set[str]:
    """Read the consequence entities that the consequence fetch selected.

    Args:
        raw_dir: The raw cache directory.

    Returns:
        The consequence accessions, or an empty set when the fetch has not run.
    """
    cache = read_json(raw_dir / "pubtator_consequences.json") or {}
    return {entry["accession"] for entry in cache.get("consequences", [])}


def select_consequence_layer(consequences: list[dict], records_by_accession: dict[str, list[dict]],
                             nodes_in_graph: set[str],
                             chemical_min: int) -> tuple[set[str], list[dict]]:
    """Choose the nodes and the relation records of the consequence layer.

    A consequence relation is kept when the other endpoint is already a node of the graph, or
    when it is a compound that reaches the publication threshold. The compound rule is the point
    of this layer. It is what lets a ranker travel from the disease state that BUBR1 loss
    produces to a compound that acts on that state.

    Args:
        consequences: The selected consequence records, each with an `accession`.
        records_by_accession: The relation records of each consequence, keyed by accession.
        nodes_in_graph: The accessions that the literature layer already holds.
        chemical_min: The smallest publication count of an admitted compound.

    Returns:
        The accessions that the layer adds and the selected relation records.
    """
    added: set[str] = set()
    selected: list[dict] = []
    seen: set[tuple[str, str, str]] = set()
    for entry in consequences:
        accession = entry["accession"]
        added.add(accession)
        for record in records_by_accession.get(accession, []):
            source, target = record.get("source", ""), record.get("target", "")
            other = target if source == accession else source
            known = other in nodes_in_graph or other in added
            strong = (other.startswith("@CHEMICAL_")
                      and int(record.get("publications", 0) or 0) >= chemical_min)
            if not known and not strong:
                continue
            if not accession_kind(other):
                continue
            added.add(other)
            key = edge_key(source, target, record.get("type", ""))
            if key in seen:
                continue
            seen.add(key)
            selected.append(record)
    return added, selected


def base_selection(raw_dir: Path, chemical_min: int) -> tuple[set[str], list[dict]]:
    """Run the two tier selection of the literature layer from the cache.

    Args:
        raw_dir: The raw cache directory.
        chemical_min: The smallest publication count of a second tier chemical.

    Returns:
        The node accessions and the selected relation records. Both are empty when the relation
        cache is absent.
    """
    index = read_json(raw_dir / "pubtator_relations_index.json") or {}
    seeds = index.get("seed_genes", [])
    neighbours = index.get("derived_neighbours", [])
    if not seeds:
        return set(), []
    records = load_relation_cache(raw_dir, seeds + neighbours, "pubtator_relations")
    return select_graph(records, seeds, neighbours, chemical_min)


def consequence_selection(raw_dir: Path, chemical_min: int) -> tuple[set[str], list[dict]]:
    """Run the consequence layer selection from the cache.

    Args:
        raw_dir: The raw cache directory.
        chemical_min: The smallest publication count of an admitted compound.

    Returns:
        The accessions that the layer adds and the selected relation records. Both are empty when
        the consequence fetch has not run.
    """
    cache = read_json(raw_dir / "pubtator_consequences.json") or {}
    consequences = cache.get("consequences", [])
    if not consequences:
        return set(), []
    base_nodes, _ = base_selection(raw_dir, chemical_min)
    by_accession = {}
    for entry in consequences:
        body = split_accession(entry["accession"])[1]
        by_accession[entry["accession"]] = read_json(
            raw_dir / "pubtator_consequences" / f"{body}.json") or []
    return select_consequence_layer(consequences, by_accession, base_nodes, chemical_min)


def consequence_nodes(raw_dir: Path, chemical_min: int) -> set[str]:
    """Return the accessions that the consequence layer adds.

    Args:
        raw_dir: The raw cache directory.
        chemical_min: The smallest publication count of an admitted compound.

    Returns:
        The accessions, or an empty set when the consequence fetch has not run.
    """
    return consequence_selection(raw_dir, chemical_min)[0]


def consequence_records(raw_dir: Path, chemical_min: int) -> list[dict]:
    """Return the relation records that the consequence layer adds.

    Args:
        raw_dir: The raw cache directory.
        chemical_min: The smallest publication count of an admitted compound.

    Returns:
        The records, or an empty list when the consequence fetch has not run.
    """
    return consequence_selection(raw_dir, chemical_min)[1]


def load_relation_cache(raw_dir: Path, symbols: list[str], subdir: str) -> dict[str, list[dict]]:
    """Read the cached relation file of each gene.

    Args:
        raw_dir: The raw cache directory.
        symbols: The gene symbols to load.
        subdir: The cache subdirectory that holds one file per gene.

    Returns:
        The relation records of each gene that has a cache file.
    """
    records: dict[str, list[dict]] = {}
    for symbol in symbols:
        cached = read_json(raw_dir / subdir / f"{symbol}.json")
        if isinstance(cached, list):
            records[symbol] = cached
    return records


def make_client() -> httpx.Client:
    """Create the HTTP client used by every Track 2 fetch script.

    Returns:
        A client with the project user agent and a generous timeout.
    """
    return httpx.Client(
        timeout=DEFAULT_TIMEOUT,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
        follow_redirects=True,
    )


def retry_delay(error: Exception, attempt: int) -> float:
    """Return the number of seconds to wait before the next attempt.

    NCBI answers a burst with HTTP 429. That status needs a long wait, not the short wait that
    suits a dropped connection. The server sometimes names the wait in a `Retry-After` header,
    and that value is used when it is present.

    Args:
        error: The error that the attempt raised.
        attempt: The zero based attempt number.

    Returns:
        The wait in seconds.
    """
    response = getattr(error, "response", None)
    status = getattr(response, "status_code", 0)
    if status in (429, 503):
        header = (response.headers or {}).get("Retry-After", "") if response is not None else ""
        if header.strip().isdigit():
            return min(float(header.strip()), 120.0)
        return min(THROTTLE_BASE_S * (2 ** attempt), 120.0)
    return 1.5 * (attempt + 1)


def is_throttle(error: Exception) -> bool:
    """Return whether an error is a refusal of a burst rather than a failed request.

    Args:
        error: The error that an attempt raised.

    Returns:
        True when the server answered with a throttle status.
    """
    response = getattr(error, "response", None)
    return getattr(response, "status_code", 0) in THROTTLE_STATUS


def note_throttle() -> None:
    """Raise the process wide pause after a server refused a burst."""
    global _EXTRA_PAUSE_S  # noqa: PLW0603
    _EXTRA_PAUSE_S = min(_EXTRA_PAUSE_S + THROTTLE_PAUSE_STEP_S, THROTTLE_PAUSE_MAX_S)


def note_success() -> None:
    """Lower the process wide pause after a server answered a call."""
    global _EXTRA_PAUSE_S  # noqa: PLW0603
    _EXTRA_PAUSE_S = _EXTRA_PAUSE_S * THROTTLE_PAUSE_DECAY
    if _EXTRA_PAUSE_S < 0.01:
        _EXTRA_PAUSE_S = 0.0


def extra_pause() -> float:
    """Return the process wide pause in seconds.

    Returns:
        The pause that every call of this process pays on top of the pace of its caller.
    """
    return _EXTRA_PAUSE_S


def pay_extra_pause() -> None:
    """Wait for the process wide pause before a request."""
    if _EXTRA_PAUSE_S > 0.0:
        time.sleep(_EXTRA_PAUSE_S)


def attempt_budget(attempt: int, throttled: int) -> bool:
    """Return whether another attempt is allowed.

    Args:
        attempt: The number of attempts already made.
        throttled: The number of those attempts that a throttle status ended.

    Returns:
        True when the call may be sent again.
    """
    if throttled and throttled == attempt:
        return attempt < MAX_THROTTLE_ATTEMPTS
    return attempt < MAX_ATTEMPTS


def record_call(url: str, byte_count: int, failed: bool = False) -> None:
    """Add one HTTP call to the counter of its host.

    The counter answers the cost question of the Track 2 report: how many calls a run sends to
    each source and how many bytes it retrieves. It holds no request body and no response body.

    Args:
        url: The request URL.
        byte_count: The size of the response body in bytes.
        failed: True when the attempt raised an error.
    """
    host = urlsplit(url).netloc or "unknown"
    entry = _CALL_COUNTS.setdefault(host, {"calls": 0, "bytes": 0, "failed": 0})
    entry["calls"] += 1
    entry["bytes"] += int(byte_count)
    entry["failed"] += 1 if failed else 0


def counters() -> dict[str, dict[str, int]]:
    """Return the call counter of this process.

    Returns:
        The number of calls, the number of response bytes and the number of failed attempts of
        each host.
    """
    return {host: dict(values) for host, values in sorted(_CALL_COUNTS.items())}


def dump_counters() -> None:
    """Write the call counter to the file that `TRACK2_CALL_LOG` names.

    The function does nothing when the variable is unset. It is registered with `atexit`, so a
    fetch script needs no change to record its own API usage.
    """
    target = os.environ.get(CALL_LOG_VARIABLE, "")
    if not target or not _CALL_COUNTS:
        return
    path = Path(target)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(counters(), indent=1, sort_keys=True))


atexit.register(dump_counters)


def response_size(response: httpx.Response) -> int:
    """Return the size of a response body in bytes.

    Args:
        response: The HTTP response.

    Returns:
        The number of bytes, or 0 when the body cannot be read.
    """
    try:
        return len(response.content)
    except (httpx.HTTPError, RuntimeError):
        return 0


def get_json(client: httpx.Client, url: str, params: dict[str, Any] | None = None) -> Any:
    """Send a GET request and decode the JSON body, with a retry on a transient failure.

    Args:
        client: The HTTP client.
        url: The request URL.
        params: The query parameters.

    Returns:
        The decoded JSON body.

    Raises:
        httpx.HTTPError: When every attempt failed.
    """
    last: Exception | None = None
    attempt = 0
    throttled = 0
    while attempt_budget(attempt, throttled):
        try:
            pay_extra_pause()
            response = client.get(url, params=params)
            record_call(url, response_size(response))
            response.raise_for_status()
            note_success()
            return response.json()
        except (httpx.HTTPError, json.JSONDecodeError) as error:
            last = error
            record_call(url, 0, failed=True)
            if is_throttle(error):
                throttled += 1
                note_throttle()
            time.sleep(retry_delay(error, attempt))
            attempt += 1
    raise httpx.HTTPError(f"GET failed after {attempt} attempts, of which {throttled} met a "
                          f"throttle status: {url} {params} ({last})")


def post_json(client: httpx.Client, url: str, payload: dict[str, Any]) -> Any:
    """Send a POST request with a JSON body and decode the JSON response.

    Args:
        client: The HTTP client.
        url: The request URL.
        payload: The JSON body.

    Returns:
        The decoded JSON body.

    Raises:
        httpx.HTTPError: When every attempt failed.
    """
    last: Exception | None = None
    for attempt in range(MAX_ATTEMPTS):
        try:
            response = client.post(url, json=payload)
            record_call(url, response_size(response))
            response.raise_for_status()
            return response.json()
        except (httpx.HTTPError, json.JSONDecodeError) as error:
            last = error
            record_call(url, 0, failed=True)
            time.sleep(1.5 * (attempt + 1))
    raise httpx.HTTPError(f"POST failed after {MAX_ATTEMPTS} attempts: {url} ({last})")


def write_json(path: Path, payload: Any) -> Path:
    """Write a JSON file and create the parent directory.

    Args:
        path: The output file.
        payload: The value to serialize.

    Returns:
        The output path.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=1, sort_keys=False))
    return path


def read_json(path: Path) -> Any:
    """Read a JSON file.

    Args:
        path: The input file.

    Returns:
        The decoded value, or None when the file is absent.
    """
    if not path.is_file():
        return None
    return json.loads(path.read_text())


def sha256_of(path: Path) -> str:
    """Compute the SHA-256 digest of a file.

    Args:
        path: The file to hash.

    Returns:
        The hexadecimal digest, or an empty string when the file is absent.
    """
    if not path.is_file():
        return ""
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def write_tsv(path: Path, header: list[str], rows: list[dict[str, Any]]) -> Path:
    """Write a tab separated file with a header row.

    A tab, a carriage return and a newline inside a value become a space, so every record
    stays on one line and the column count never changes.

    Args:
        path: The output file.
        header: The column names in order.
        rows: The records, keyed by column name.

    Returns:
        The output path.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)]
    for row in rows:
        values = []
        for column in header:
            value = str(row.get(column, ""))
            values.append(value.replace("\t", " ").replace("\r", " ").replace("\n", " "))
        lines.append("\t".join(values))
    path.write_text("\n".join(lines) + "\n")
    return path
