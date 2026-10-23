"""Build the R12 target interval list for the targeted local reassembly. Data class A.

The list holds the seven MVA genes on GRCh38. Each gene contributes its full gene body, 50 kb
upstream of the transcription start and 10 kb downstream of the transcription end. The direction
follows the strand of the gene. The list also holds the two known BUB1B intervals of the A5 work,
and each of those two carries the label `known`.

The gene model is the GENCODE 50 basic annotation that queue item R1 pinned. The script records the
SHA-256 of that file, so a later run can prove that it read the same model.

Usage:
    uv run python scripts/build_r12_targets.py \
        --gtf inputs/gencode.v50.basic.annotation.gtf.gz \
        --known-intervals inputs/intervals.bed \
        --out-bed out/r12_targets.bed \
        --out-manifest out/r12_targets_manifest.json
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

# The structural annotation set of literature/track1_gene_panel_v1.yaml and of
# literature/track1_gene_panel_v2.yaml. harness/test_track1_gene_panel.py fails when a tuple no
# longer equals the operating set that the matching panel derives. The shipped default stays at v1,
# because the recorded R12 target list used the v1 set.
PANEL_SETS = {
    "v1": ("BUB1B", "CEP57", "TRIP13", "BUB1", "BUB3", "CDC20", "PTPA"),
    "v2": ("BUB1B", "CEP57", "TRIP13", "CENATAC", "MAD1L1", "MAD2L1BP", "BUB1", "BUB3", "CDC20",
           "PTPA"),
}
SHIPPED_PANEL_VERSION = "v1"
MVA_GENES = PANEL_SETS[SHIPPED_PANEL_VERSION]
UPSTREAM_BP = 50_000
DOWNSTREAM_BP = 10_000


def file_sha256(path: Path) -> str:
    """Return the SHA-256 of one file.

    Args:
        path: The file path.

    Returns:
        The digest in lower case hexadecimal.
    """
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def attribute_value(attributes: str, key: str) -> str:
    """Read one attribute of a GTF attribute column.

    Args:
        attributes: The ninth column of a GTF line.
        key: The attribute name.

    Returns:
        The value without quotes, or an empty string when the key is absent.
    """
    marker = key + ' "'
    start = attributes.find(marker)
    if start < 0:
        return ""
    start += len(marker)
    end = attributes.find('"', start)
    return attributes[start:end]


def read_gene_records(gtf_path: Path, genes: tuple[str, ...]) -> dict[str, dict]:
    """Read the gene body of each requested gene from a GENCODE GTF.

    A gene that appears more than once keeps the record on the primary chromosome with the widest
    span. A record on a scaffold or on a patch is skipped.

    Args:
        gtf_path: The gzip-compressed GENCODE GTF.
        genes: The gene symbols to read.

    Returns:
        A dictionary from gene symbol to a record with the contig, the start, the end, the strand
        and the GENCODE gene identifier. The start is 1-based and inclusive.
    """
    wanted = set(genes)
    found: dict[str, dict] = {}
    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            name = attribute_value(fields[8], "gene_name")
            if name not in wanted:
                continue
            contig = fields[0]
            if "_" in contig:
                continue
            record = {"gene": name, "contig": contig, "start": int(fields[3]), "end": int(fields[4]),
                      "strand": fields[6], "gene_id": attribute_value(fields[8], "gene_id"),
                      "gene_type": attribute_value(fields[8], "gene_type")}
            previous = found.get(name)
            if previous is None or (record["end"] - record["start"]) > (previous["end"] - previous["start"]):
                found[name] = record
    return found


def padded_interval(record: dict, upstream: int, downstream: int) -> tuple[int, int]:
    """Add a strand-aware pad to one gene body.

    Args:
        record: A gene record with a start, an end and a strand.
        upstream: The pad upstream of the transcription start, in base pairs.
        downstream: The pad downstream of the transcription end, in base pairs.

    Returns:
        The padded start and the padded end. The start is 1-based and inclusive.
    """
    if record["strand"] == "-":
        start = record["start"] - downstream
        end = record["end"] + upstream
    else:
        start = record["start"] - upstream
        end = record["end"] + downstream
    return max(1, start), end


def read_known_intervals(path: Path) -> list[dict]:
    """Read the known BUB1B intervals from a BED file.

    Args:
        path: The BED file of the A5 work.

    Returns:
        One record per interval, with a 0-based start and an exclusive end.
    """
    records = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        fields = line.split("\t")
        records.append({"contig": fields[0], "start": int(fields[1]), "end": int(fields[2]),
                        "name": fields[3] if len(fields) > 3 else "known_interval"})
    return records


def contig_rank(contig: str) -> tuple[int, str]:
    """Sort key for a GRCh38 primary contig.

    Args:
        contig: The contig name with a chr prefix.

    Returns:
        A number for an autosome and a large number for any other contig, with the name as a tie
        break.
    """
    name = contig.removeprefix("chr")
    return (int(name) if name.isdigit() else 99), name


def merge_intervals(intervals: list[dict]) -> list[dict]:
    """Merge overlapping or touching intervals of the same contig.

    Args:
        intervals: Records with a contig, a 0-based start, an exclusive end and a name.

    Returns:
        The merged records, sorted by contig and start. A merged record joins the names with a
        plus sign and keeps the `known` flag when any input carried it.
    """
    ordered = sorted(intervals, key=lambda item: (contig_rank(item["contig"]), item["start"], item["end"]))
    merged: list[dict] = []
    for item in ordered:
        last = merged[-1] if merged else None
        if last is not None and last["contig"] == item["contig"] and item["start"] <= last["end"]:
            last["end"] = max(last["end"], item["end"])
            if item["name"] not in last["name"].split("+"):
                last["name"] = last["name"] + "+" + item["name"]
            last["known"] = last["known"] or item["known"]
        else:
            merged.append(dict(item))
    return merged


def build_targets(gtf_path: Path, known_path: Path) -> tuple[list[dict], dict]:
    """Build the merged target list and the gene table.

    Args:
        gtf_path: The GENCODE 50 basic annotation.
        known_path: The BED file with the two known BUB1B intervals.

    Returns:
        The merged interval records and a table with one entry per gene.
    """
    genes = read_gene_records(gtf_path, MVA_GENES)
    missing = [name for name in MVA_GENES if name not in genes]
    if missing:
        raise SystemExit(f"The gene model holds no record for {', '.join(missing)}.")
    gene_table = {}
    intervals = []
    for name in MVA_GENES:
        record = genes[name]
        start, end = padded_interval(record, UPSTREAM_BP, DOWNSTREAM_BP)
        gene_table[name] = {
            "gene_id": record["gene_id"], "gene_type": record["gene_type"], "strand": record["strand"],
            "gene_body": f"{record['contig']}:{record['start']}-{record['end']}",
            "gene_body_bp": record["end"] - record["start"] + 1,
            "target": f"{record['contig']}:{start}-{end}", "target_bp": end - start + 1,
        }
        intervals.append({"contig": record["contig"], "start": start - 1, "end": end,
                          "name": name, "known": False})
    for item in read_known_intervals(known_path):
        intervals.append({"contig": item["contig"], "start": item["start"], "end": item["end"],
                          "name": item["name"], "known": True})
    return merge_intervals(intervals), gene_table


def write_outputs(merged: list[dict], gene_table: dict, gtf_path: Path, known_path: Path,
                  out_bed: Path, out_manifest: Path) -> dict:
    """Write the target BED file and the manifest.

    Args:
        merged: The merged interval records.
        gene_table: The per-gene table.
        gtf_path: The gene model file, for the checksum.
        known_path: The known interval file, for the checksum.
        out_bed: The BED output path.
        out_manifest: The manifest output path.

    Returns:
        The manifest content.
    """
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{item['contig']}\t{item['start']}\t{item['end']}\t{item['name']}\t"
             f"{'known' if item['known'] else 'gene'}" for item in merged]
    out_bed.write_text("\n".join(lines) + "\n")
    manifest = {
        "target_list_id": "r12_targets_v1",
        "data_class": "A",
        "gene_model": {
            "name": "GENCODE 50 basic annotation, GRCh38",
            "source": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/"
                      "gencode.v50.basic.annotation.gtf.gz",
            "path": str(gtf_path),
            "sha256": file_sha256(gtf_path),
            "pinned_by": "r1_sv_mei_v1",
        },
        "known_intervals": {
            "name": "The two BUB1B literature intervals of the A5 and S6 work",
            "path": str(known_path),
            "sha256": file_sha256(known_path),
            "count": len(read_known_intervals(known_path)),
        },
        "pad_bp": {"upstream": UPSTREAM_BP, "downstream": DOWNSTREAM_BP,
                   "rule": "The pad follows the strand of the gene."},
        "genes": gene_table,
        "intervals": merged,
        "interval_count": len(merged),
        "total_bp": sum(item["end"] - item["start"] for item in merged),
        "bed_sha256": hashlib.sha256(out_bed.read_bytes()).hexdigest(),
    }
    out_manifest.parent.mkdir(parents=True, exist_ok=True)
    out_manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> int:
    """Build the target list and print the manifest.

    Returns:
        Zero on success.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--known-intervals", type=Path, required=True)
    parser.add_argument("--out-bed", type=Path, required=True)
    parser.add_argument("--out-manifest", type=Path, required=True)
    args = parser.parse_args()
    merged, gene_table = build_targets(args.gtf, args.known_intervals)
    manifest = write_outputs(merged, gene_table, args.gtf, args.known_intervals,
                             args.out_bed, args.out_manifest)
    print(json.dumps({"interval_count": manifest["interval_count"], "total_bp": manifest["total_bp"],
                      "bed": str(args.out_bed)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
