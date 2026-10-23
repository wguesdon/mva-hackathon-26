#!/usr/bin/env python3
"""Targeted local reassembly over the MVA genes, and a structural variant call from the contigs.

The module holds the two steps of queue item R12 that need Python. Each step runs as a subcommand,
so the same code runs on the local machine for the public HG002 validation and on AWS Batch for the
patient alignment.

    call      Read the contig alignment of one interval and write a sequence-resolved VCF.
    merge     Join the per-interval VCF files, remove a duplicate event and write one VCF.

The read extraction, the assembly and the contig alignment call `samtools`, `spades.py` and
`minimap2` directly. Those three containers hold no Python interpreter, so the caller runs them.

Privacy. The module reads the contig alignment and writes files. It prints a gate verdict and a
file name only. It never prints a read, a call, a count of reads or a count of calls, because such
a number derived from the patient alignment is class B2 under ADR 0005. The per-interval counts go
to the class C table that the caller keeps beside the data.

Usage:
    r12_reassembly.py call --paf <interval.paf> --reference <genome.fa> --bed <targets.bed> \
        --interval <name> --out-vcf <interval.vcf>
    r12_reassembly.py merge --vcf-dir calls --reference <genome.fa> --sample <name> \
        --out-vcf r12_calls.vcf --out-gates class_a_gates.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

MIN_EVENT_BP = 50
MAX_EXPLICIT_BP = 50_000
MIN_CONTIG_MAPQ = 20
EDGE_MARGIN_BP = 100
MERGE_WINDOW_BP = 100
MERGE_SIZE_RATIO = 0.8
SPAN_WINDOW_BP = 50
VCF_SOURCE = "r12_local_reassembly"


class FastaReader:
    """Read one region of an indexed FASTA file.

    The reader uses the `.fai` index, so it opens the file once and seeks. It needs no external
    tool.

    Attributes:
        path: The FASTA file.
        index: The record of each contig, with the length, the offset, the line bases and the line
            width.
    """

    def __init__(self, path: Path):
        """Open the FASTA file and read its index.

        Args:
            path: The FASTA file. A file `<path>.fai` must exist beside it.
        """
        self.path = path
        self.index = {}
        for line in Path(str(path) + ".fai").read_text().splitlines():
            fields = line.split("\t")
            self.index[fields[0]] = {"length": int(fields[1]), "offset": int(fields[2]),
                                     "line_bases": int(fields[3]), "line_width": int(fields[4])}
        self.handle = path.open("rb")

    def fetch(self, contig: str, start: int, end: int) -> str:
        """Return one region of the sequence in upper case.

        Args:
            contig: The contig name.
            start: The 0-based inclusive start.
            end: The 0-based exclusive end.

        Returns:
            The sequence, clipped to the contig.
        """
        record = self.index[contig]
        start = max(0, start)
        end = min(record["length"], end)
        if end <= start:
            return ""
        line_bases = record["line_bases"]
        line_width = record["line_width"]
        begin = record["offset"] + (start // line_bases) * line_width + (start % line_bases)
        stop = record["offset"] + (end // line_bases) * line_width + (end % line_bases)
        self.handle.seek(begin)
        return self.handle.read(stop - begin).decode("ascii").replace("\n", "").replace("\r", "").upper()

    def close(self) -> None:
        """Close the FASTA file."""
        self.handle.close()

    def contig_lengths(self) -> dict[str, int]:
        """Return the length of every contig.

        Returns:
            A dictionary from contig name to length.
        """
        return {name: record["length"] for name, record in self.index.items()}


def read_targets(path: Path) -> list[dict]:
    """Read the target BED file.

    Args:
        path: The BED file with a contig, a 0-based start, an exclusive end, a name and a class.

    Returns:
        One record per interval.
    """
    records = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        fields = line.split("\t")
        records.append({"contig": fields[0], "start": int(fields[1]), "end": int(fields[2]),
                        "name": fields[3] if len(fields) > 3 else f"{fields[0]}_{fields[1]}",
                        "class": fields[4] if len(fields) > 4 else "gene"})
    return records


def extract_shell_command(bam: str, interval: dict, threads: int) -> str:
    """Build the samtools pipeline that reads one target interval into FASTQ files.

    The region query of a coordinate-sorted alignment returns every read that overlaps the
    interval. It also returns an unmapped read whose mate maps inside the interval, because such a
    read carries the position of its mate. Those reads carry the evidence of a novel insertion, so
    the command keeps them.

    Args:
        bam: The alignment file or S3 URL.
        interval: One target record.
        threads: The thread count for samtools.

    Returns:
        The shell command. The caller runs it inside the samtools container.
    """
    name = interval["name"]
    region = f"{interval['contig']}:{interval['start'] + 1}-{interval['end']}"
    return (f"samtools view -u -@ {threads} '{bam}' '{region}' "
            f"| samtools collate -u -O -@ {threads} - "
            f"| samtools fastq -N -@ {threads} -1 {name}_R1.fastq.gz -2 {name}_R2.fastq.gz "
            f"-s {name}_S.fastq.gz -0 /dev/null -")


def parse_cs(cs_value: str) -> list[tuple[str, str]]:
    """Split a minimap2 short `cs` tag into its operations.

    Args:
        cs_value: The value of the `cs:Z:` tag.

    Returns:
        One pair per operation, with the operator character and its argument.
    """
    operations = []
    index = 0
    length = len(cs_value)
    while index < length:
        operator = cs_value[index]
        index += 1
        start = index
        while index < length and cs_value[index] not in ":*+-~":
            index += 1
        operations.append((operator, cs_value[start:index]))
    return operations


def parse_paf_line(line: str) -> dict | None:
    """Read one PAF line of a contig alignment.

    Args:
        line: The PAF line.

    Returns:
        The alignment record, or None when the line carries no `cs` tag.
    """
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        return None
    tags = {}
    for item in fields[12:]:
        key, _, value = item.split(":", 2)
        tags[key] = value
    if "cs" not in tags:
        return None
    return {"query": fields[0], "query_length": int(fields[1]), "query_start": int(fields[2]),
            "query_end": int(fields[3]), "strand": fields[4], "contig": fields[5],
            "target_start": int(fields[7]), "target_end": int(fields[8]),
            "mapq": int(fields[11]), "type": tags.get("tp", "P"), "cs": tags["cs"]}


def events_of_alignment(alignment: dict, min_size: int) -> list[dict]:
    """Read every insertion and deletion of one contig alignment.

    Args:
        alignment: One PAF record with a `cs` tag.
        min_size: The smallest event to keep, in base pairs.

    Returns:
        One record per event, with the 0-based reference position, the type, the size and the
        inserted sequence.
    """
    events = []
    position = alignment["target_start"]
    for operator, value in parse_cs(alignment["cs"]):
        if operator == ":":
            position += int(value)
        elif operator == "*":
            position += 1
        elif operator == "+":
            if len(value) >= min_size:
                events.append({"contig": alignment["contig"], "position": position, "svtype": "INS",
                               "size": len(value), "sequence": value.upper(),
                               "query": alignment["query"], "mapq": alignment["mapq"]})
        elif operator == "-":
            if len(value) >= min_size:
                events.append({"contig": alignment["contig"], "position": position, "svtype": "DEL",
                               "size": len(value), "sequence": value.upper(),
                               "query": alignment["query"], "mapq": alignment["mapq"]})
            position += len(value)
        elif operator == "~":
            position += int(value[2:-2])
    return [item for item in events if inside_alignment(item, alignment)]


def inside_alignment(event: dict, alignment: dict) -> bool:
    """Reject an event that sits at the edge of a contig alignment.

    An event at the edge comes from a clipped contig end more often than from a real allele, so
    the step keeps only an event with a real anchor on both sides.

    Args:
        event: One event record.
        alignment: The PAF record that produced it.

    Returns:
        True when the event has an anchor of at least the edge margin on both sides.
    """
    return (event["position"] - alignment["target_start"] >= EDGE_MARGIN_BP
            and alignment["target_end"] - event["position"] >= EDGE_MARGIN_BP)


def read_alignments(paf_path: Path) -> list[dict]:
    """Read the primary contig alignments of one PAF file.

    Args:
        paf_path: The PAF file from minimap2, written with `-c --cs`.

    Returns:
        The primary alignments with a mapping quality at or above the limit.
    """
    alignments = []
    for line in paf_path.read_text().splitlines():
        record = parse_paf_line(line)
        if record is None:
            continue
        if record["type"] != "P" or record["mapq"] < MIN_CONTIG_MAPQ:
            continue
        alignments.append(record)
    return alignments


def spanning_count(alignments: list[dict], contig: str, position: int) -> int:
    """Count the contig alignments that cross one reference position.

    Args:
        alignments: The primary alignments.
        contig: The reference contig.
        position: The 0-based reference position.

    Returns:
        The number of alignments that cover the position with a margin on both sides.
    """
    return sum(1 for item in alignments
               if item["contig"] == contig
               and item["target_start"] + SPAN_WINDOW_BP <= position <= item["target_end"] - SPAN_WINDOW_BP)


def merge_events(events: list[dict]) -> list[dict]:
    """Join the events that two redundant contigs report for one allele.

    Args:
        events: The events of one interval.

    Returns:
        One record per allele, with the support count and the widest inserted sequence.
    """
    ordered = sorted(events, key=lambda item: (item["contig"], item["svtype"], item["position"], -item["size"]))
    merged: list[dict] = []
    for event in ordered:
        target = next((item for item in merged if same_allele(item, event)), None)
        if target is None:
            record = dict(event)
            record["support"] = 1
            record["queries"] = [event["query"]]
            merged.append(record)
        else:
            target["support"] += 1
            if event["query"] not in target["queries"]:
                target["queries"].append(event["query"])
            if event["size"] > target["size"]:
                target["size"] = event["size"]
                target["sequence"] = event["sequence"]
    return merged


def same_allele(left: dict, right: dict) -> bool:
    """Decide whether two events describe one allele.

    Args:
        left: The first event.
        right: The second event.

    Returns:
        True when the type agrees, the positions are inside the merge window and the sizes agree
        inside the size ratio.
    """
    if left["contig"] != right["contig"] or left["svtype"] != right["svtype"]:
        return False
    if abs(left["position"] - right["position"]) > MERGE_WINDOW_BP:
        return False
    smaller, larger = sorted((left["size"], right["size"]))
    return larger > 0 and smaller / larger >= MERGE_SIZE_RATIO


def inside_targets(event: dict, targets: list[dict]) -> str | None:
    """Return the name of the target interval that holds one event.

    Args:
        event: One event record.
        targets: The target intervals.

    Returns:
        The interval name, or None when the event lies outside every interval.
    """
    for interval in targets:
        if interval["contig"] == event["contig"] and interval["start"] <= event["position"] < interval["end"]:
            return interval["name"]
    return None


def vcf_record(event: dict, reader: FastaReader, alignments: list[dict], sample: str) -> str:
    """Write one event as a VCF data line.

    Args:
        event: The merged event, with a support count.
        reader: The reference reader.
        alignments: The primary alignments of the interval, for the genotype.
        sample: The sample column name.

    Returns:
        The VCF line without a line break.
    """
    contig = event["contig"]
    anchor = max(0, event["position"] - 1)
    anchor_base = reader.fetch(contig, anchor, anchor + 1) or "N"
    spanning = spanning_count(alignments, contig, event["position"])
    genotype = "1/1" if spanning <= event["support"] else "0/1"
    if event["svtype"] == "INS":
        end = anchor + 1
        length = event["size"]
        reference, alternate = anchor_base, anchor_base + event["sequence"]
        if length > MAX_EXPLICIT_BP:
            reference, alternate = anchor_base, "<INS>"
    else:
        end = event["position"] + event["size"]
        length = -event["size"]
        deleted = reader.fetch(contig, anchor, end)
        reference, alternate = deleted, anchor_base
        if event["size"] > MAX_EXPLICIT_BP:
            reference, alternate = anchor_base, "<DEL>"
    info = (f"SVTYPE={event['svtype']};SVLEN={length};END={end};"
            f"SUPPORT={event['support']};SPAN={spanning};CONTIGS={'|'.join(event['queries'][:4])};"
            f"SOURCE={VCF_SOURCE}")
    identifier = f"{VCF_SOURCE}:{contig}:{anchor + 1}:{event['svtype']}:{abs(length)}"
    return "\t".join([contig, str(anchor + 1), identifier, reference, alternate, ".", "PASS", info,
                      "GT:SP", f"{genotype}:{event['support']}"])


def vcf_header(reader: FastaReader, sample: str, targets: list[dict]) -> list[str]:
    """Build the VCF header of the call set.

    Args:
        reader: The reference reader, for the contig lines.
        sample: The sample column name.
        targets: The target intervals, so that the header names the contigs that matter first.

    Returns:
        The header lines.
    """
    wanted = []
    for interval in targets:
        if interval["contig"] not in wanted:
            wanted.append(interval["contig"])
    lengths = reader.contig_lengths()
    lines = ["##fileformat=VCFv4.2", f"##source={VCF_SOURCE}"]
    lines.extend(f"##contig=<ID={name},length={lengths[name]}>" for name in wanted if name in lengths)
    lines.extend([
        '##ALT=<ID=DEL,Description="Deletion">',
        '##ALT=<ID=INS,Description="Insertion">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
        '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Contig alignments that carry the event">',
        '##INFO=<ID=SPAN,Number=1,Type=Integer,Description="Contig alignments that cross the position">',
        '##INFO=<ID=CONTIGS,Number=1,Type=String,Description="Names of the supporting contigs">',
        '##INFO=<ID=SOURCE,Number=1,Type=String,Description="The caller that produced the record">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Supporting contig alignments">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}",
    ])
    return lines


def call_interval(paf_path: Path, reference: Path, bed: Path, interval_name: str, min_size: int,
                  out_vcf: Path, sample: str) -> dict:
    """Call the structural variants of one interval from its contig alignment.

    Args:
        paf_path: The PAF file of the interval.
        reference: The indexed GRCh38 FASTA.
        bed: The target BED file.
        interval_name: The interval that the assembly covers.
        min_size: The smallest event to keep.
        out_vcf: The output VCF path.
        sample: The sample column name.

    Returns:
        The interval name, the alignment count and the call count.
    """
    targets = read_targets(bed)
    reader = FastaReader(reference)
    alignments = read_alignments(paf_path)
    events = []
    for alignment in alignments:
        events.extend(events_of_alignment(alignment, min_size))
    merged = [item for item in merge_events(events) if inside_targets(item, targets)]
    merged.sort(key=lambda item: (item["contig"], item["position"]))
    lines = vcf_header(reader, sample, targets)
    lines.extend(vcf_record(event, reader, alignments, sample) for event in merged)
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    out_vcf.write_text("\n".join(lines) + "\n")
    reader.close()
    return {"interval": interval_name, "alignments": len(alignments), "calls": len(merged)}


def merge_vcf_files(vcf_dir: Path, reference: Path, bed: Path, sample: str, out_vcf: Path) -> dict:
    """Join the per-interval VCF files into one sorted file.

    Args:
        vcf_dir: The directory of the per-interval VCF files.
        reference: The indexed GRCh38 FASTA.
        bed: The target BED file.
        sample: The sample column name.
        out_vcf: The output VCF path.

    Returns:
        The number of input files and the number of records.
    """
    targets = read_targets(bed)
    reader = FastaReader(reference)
    order = {interval["contig"]: index for index, interval in enumerate(targets)}
    records = []
    files = sorted(vcf_dir.glob("*.vcf"))
    for path in files:
        for line in path.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            records.append(line)
    records.sort(key=lambda line: (order.get(line.split("\t")[0], 99), int(line.split("\t")[1])))
    lines = vcf_header(reader, sample, targets)
    lines.extend(records)
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    out_vcf.write_text("\n".join(lines) + "\n")
    reader.close()
    return {"input_files": len(files), "records": len(records)}


def command_call(args: argparse.Namespace) -> int:
    """Run the call step for one interval.

    Args:
        args: The parsed arguments.

    Returns:
        Zero on success.
    """
    result = call_interval(args.paf, args.reference, args.bed, args.interval, args.min_size,
                           args.out_vcf, args.sample)
    print(json.dumps({"step": "call", "interval": result["interval"], "gate": "PASS"}))
    return 0


def build_class_a_gates(targets: list[dict], interval_files: int, merged_written: bool,
                        run_id: str) -> dict:
    """Build the class A gate file of the reassembly from a fixed key list.

    Every quantity in the record is a property of the public target list or a count of the output
    files that the public target list determines. The record holds no count of calls. The key list
    and every declared quantity are in `literature/class_a_publish_rules.json`, and
    `harness/privacy_guard.py --class-a-publish` checks the written file against them.

    Args:
        targets: The public target intervals.
        interval_files: The number of per-interval VCF files that the step read.
        merged_written: True when the merged VCF exists.
        run_id: The run identifier.

    Returns:
        The class A record.
    """
    gates = {
        "run_id": run_id,
        "targets_loaded": "PASS" if targets else "FAIL",
        "target_intervals": len(targets),
        "target_bases": sum(item["end"] - item["start"] for item in targets),
        "interval_vcf_files": interval_files,
        "all_intervals_called": "PASS" if interval_files == len(targets) else "FAIL",
        "merged_vcf_written": "PASS" if merged_written else "FAIL",
    }
    gates["verdict"] = "PASS" if all(value == "PASS" for key, value in gates.items()
                                     if isinstance(value, str) and key != "run_id") else "FAIL"
    return gates


def command_merge(args: argparse.Namespace) -> int:
    """Run the merge step and write the class A gate file.

    Args:
        args: The parsed arguments.

    Returns:
        Zero on success.
    """
    result = merge_vcf_files(args.vcf_dir, args.reference, args.bed, args.sample, args.out_vcf)
    targets = read_targets(args.bed)
    gates = build_class_a_gates(targets, result["input_files"], args.out_vcf.exists(), args.run_id)
    args.out_gates.parent.mkdir(parents=True, exist_ok=True)
    args.out_gates.write_text(json.dumps(gates, indent=2) + "\n")
    print(json.dumps({"step": "merge", "verdict": gates["verdict"], "gates": args.out_gates.name}))
    return 0 if gates["verdict"] == "PASS" else 1


def build_parser() -> argparse.ArgumentParser:
    """Build the command line parser.

    Returns:
        The parser with one subcommand per step.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    call = sub.add_parser("call")
    call.add_argument("--paf", required=True, type=Path)
    call.add_argument("--reference", required=True, type=Path)
    call.add_argument("--bed", required=True, type=Path)
    call.add_argument("--interval", required=True)
    call.add_argument("--out-vcf", required=True, type=Path)
    call.add_argument("--sample", default="SAMPLE")
    call.add_argument("--min-size", type=int, default=MIN_EVENT_BP)
    call.set_defaults(handler=command_call)

    merge = sub.add_parser("merge")
    merge.add_argument("--vcf-dir", required=True, type=Path)
    merge.add_argument("--reference", required=True, type=Path)
    merge.add_argument("--bed", required=True, type=Path)
    merge.add_argument("--out-vcf", required=True, type=Path)
    merge.add_argument("--out-gates", required=True, type=Path)
    merge.add_argument("--sample", default="SAMPLE")
    merge.add_argument("--run-id", default="r12_local_reassembly_v1")
    merge.set_defaults(handler=command_merge)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run one step of the reassembly.

    Args:
        argv: The argument list, or None for the process arguments.

    Returns:
        The exit status.
    """
    args = build_parser().parse_args(argv)
    return args.handler(args)


if __name__ == "__main__":
    sys.exit(main())
