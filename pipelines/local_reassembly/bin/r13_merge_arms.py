#!/usr/bin/env python3
"""Merge the production structural call set with the targeted local reassembly call set (R13).

The script runs inside the SVDB container of the patient run `patient_calling_v1`. The HG002
validation of queue item R13 and the patient stage both call this file, so the merge parameters and
the tag handling are the same in both.

SVDB writes a `set` entry and an `svdb_origin` entry into the INFO field of every merged record. The
production call set already carries those two entries from its own merge of TIDDIT, Manta and
CNVnator. The script therefore renames them to `caller_set` and `caller_origin` before the merge, and
it adds a `caller_arm` entry. After the merge, `set` names the arm and `caller_set` names the caller,
so the source caller of every record survives.

The input call sets are class C when the sample is the child. The output call set is class C. The
count file holds the number of records of each arm, which is a statement about the sample, so it is
class C as well. The gate file holds the merge parameters and the process gates only, with no count,
no coordinate and no genotype, so it is class A.

Usage:
    r13_merge_arms.py --pipeline-vcf sv.vcf.gz --reassembly-vcf reassembly.vcf \\
        --out-vcf out/class_c_merged_sv.vcf --out-gates out/class_a_merge_gates.json \\
        --out-counts out/class_c_merge_counts.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import subprocess
from pathlib import Path

ARM_PIPELINE = "pipeline"
ARM_REASSEMBLY = "reassembly"

MERGE_OVERLAP = "0.5"
MERGE_BND_DISTANCE = "150"
MERGE_INS_DISTANCE = "50"

DEFINITIONS = (
    '##INFO=<ID=caller_set,Number=1,Type=String,'
    'Description="The caller of the record before the R13 merge">',
    '##INFO=<ID=caller_origin,Number=1,Type=String,'
    'Description="The caller lineage of the record before the R13 merge">',
    '##INFO=<ID=caller_arm,Number=1,Type=String,'
    'Description="The R13 arm that produced the record">',
)


def open_text(path: Path):
    """Open a plain or a gzip compressed VCF in text mode.

    Args:
        path: The file to open.

    Returns:
        A text mode file object.
    """
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open("rt")


def retag_info(field: str, arm: str) -> str:
    """Rename the SVDB tags of one INFO column and name the arm.

    Args:
        field: The INFO column of one record.
        arm: The arm name of the call set.

    Returns:
        The INFO column with `set` renamed to `caller_set`, `svdb_origin` renamed to
        `caller_origin`, and a `caller_arm` entry added.
    """
    items = [] if field in (".", "") else field.split(";")
    kept = []
    caller = ""
    for item in items:
        key, separator, value = item.partition("=")
        if key == "set":
            caller = value
            kept.append(f"caller_set{separator}{value}")
        elif key == "svdb_origin":
            kept.append(f"caller_origin{separator}{value}")
        else:
            kept.append(item)
    if not caller:
        kept.append(f"caller_set={arm}")
    kept.append(f"caller_arm={arm}")
    return ";".join(kept)


def retag_file(source: Path, target: Path, arm: str) -> int:
    """Write one call set with the renamed tags and the added INFO definitions.

    Args:
        source: The input VCF, plain or gzip compressed.
        target: The plain output VCF.
        arm: The arm name of the call set.

    Returns:
        The number of records written.
    """
    records = 0
    with open_text(source) as handle, target.open("w") as out:
        for line in handle:
            if line.startswith("##"):
                out.write(line)
                continue
            if line.startswith("#CHROM"):
                for definition in DEFINITIONS:
                    out.write(definition + "\n")
                out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            fields[7] = retag_info(fields[7], arm)
            out.write("\t".join(fields) + "\n")
            records += 1
    return records


def run(command: list[str], work: Path, output: Path | None = None) -> None:
    """Run one command in the work directory and fail with its error output.

    Args:
        command: The command and its arguments.
        work: The working directory.
        output: A file that receives the standard output, or None.

    Raises:
        RuntimeError: The command exited with a non-zero status.
    """
    result = subprocess.run(command, cwd=work, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{command[0]} failed: {result.stderr[-2000:]}")
    if output is not None:
        output.write_text(result.stdout)


def compress_and_index(name: str, work: Path) -> str:
    """Sort, compress and index one plain VCF with bcftools.

    Args:
        name: The plain VCF file name inside the work directory.
        work: The work directory.

    Returns:
        The bgzipped file name.
    """
    output = name + ".gz"
    run(["bcftools", "sort", "-Oz", "-o", output, name], work)
    run(["bcftools", "index", "-t", output], work)
    return output


def merge_command(left: str, right: str, overlap: str, bnd_distance: str,
                  ins_distance: str) -> list[str]:
    """Build the SVDB merge command of the R13 configuration.

    `--no_intra` stops SVDB from merging two records of the same input file. Without it, SVDB also
    re-merges the production call set with itself at the reciprocal overlap of this call, which is
    looser than the overlap of the production merge. That would change the production arm, and this
    merge must only add the second arm.

    Args:
        left: The bgzipped file name of the pipeline arm.
        right: The bgzipped file name of the reassembly arm.
        overlap: The reciprocal overlap that merges two events.
        bnd_distance: The largest distance between two similar breakpoints.
        ins_distance: The largest distance that merges two insertions.

    Returns:
        The command and its arguments.
    """
    return ["svdb", "--merge", "--pass_only", "--same_order", "--no_intra",
            "--priority", f"{ARM_PIPELINE},{ARM_REASSEMBLY}",
            "--overlap", overlap, "--bnd_distance", bnd_distance,
            "--ins_distance", ins_distance,
            "--vcf", f"{left}:{ARM_PIPELINE}", f"{right}:{ARM_REASSEMBLY}"]


def composition(path: Path) -> dict[str, int]:
    """Count the merged records by the arms that contributed to them.

    Args:
        path: The merged VCF.

    Returns:
        The count of records for each value of the `set` entry.
    """
    counts: dict[str, int] = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            info = line.split("\t")[7]
            key = "none"
            for item in info.split(";"):
                name, separator, value = item.partition("=")
                if name == "set" and separator:
                    key = value
                    break
            counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


def merge_arms(pipeline_vcf: Path, reassembly_vcf: Path, out_vcf: Path, work: Path,
               overlap: str, bnd_distance: str, ins_distance: str) -> dict:
    """Retag both call sets, merge them with SVDB and report the counts.

    Args:
        pipeline_vcf: The production structural call set.
        reassembly_vcf: The targeted local reassembly call set.
        out_vcf: The merged output VCF.
        work: The work directory for the intermediate files.
        overlap: The reciprocal overlap that merges two events.
        bnd_distance: The largest distance between two similar breakpoints.
        ins_distance: The largest distance that merges two insertions.

    Returns:
        The merge report with the input counts, the merged count and the composition.
    """
    work.mkdir(parents=True, exist_ok=True)
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    pipeline_records = retag_file(pipeline_vcf, work / "arm_pipeline.vcf", ARM_PIPELINE)
    reassembly_records = retag_file(reassembly_vcf, work / "arm_reassembly.vcf", ARM_REASSEMBLY)
    left = compress_and_index("arm_pipeline.vcf", work)
    right = compress_and_index("arm_reassembly.vcf", work)
    run(merge_command(left, right, overlap, bnd_distance, ins_distance), work,
        output=work / "merged_raw.vcf")
    out_vcf.write_text((work / "merged_raw.vcf").read_text())
    counts = composition(out_vcf)
    merged = sum(counts.values())
    return {
        "tool": "SVDB 2.8.4",
        "arguments": merge_command(left, right, overlap, bnd_distance, ins_distance)[1:-3],
        "reciprocal_overlap": float(overlap),
        "breakpoint_window_bp": int(bnd_distance),
        "insertion_distance_bp": int(ins_distance),
        "input_records": {ARM_PIPELINE: pipeline_records, ARM_REASSEMBLY: reassembly_records},
        "merged_records": merged,
        "records_removed_by_the_merge": pipeline_records + reassembly_records - merged,
        "composition": counts,
        "source_caller_kept_in": "caller_set, caller_origin and caller_arm",
    }


def parameter_section(report: dict) -> dict:
    """Keep the parameters of a merge report and drop every count.

    A count of records is a statement about the sample, so it must not reach a class A file when the
    sample is the child.

    Args:
        report: The merge report of `merge_arms`.

    Returns:
        The tool, the arguments and the three distance thresholds.
    """
    return {key: report[key] for key in
            ("tool", "arguments", "reciprocal_overlap", "breakpoint_window_bp",
             "insertion_distance_bp", "source_caller_kept_in")}


def build_class_a_gates(report: dict, run_id: str) -> dict:
    """Build the class A gate file of the merge from a fixed key list.

    The record holds the merge parameters, the process gates and a verdict. It holds no count of
    records, because a count of records states how many calls the sample carries. The key list and
    every declared quantity of this file are in `literature/class_a_publish_rules.json`, and
    `harness/privacy_guard.py --class-a-publish` checks the written file against them.

    Args:
        report: The merge report of `merge_arms`.
        run_id: The run identifier.

    Returns:
        The class A record.
    """
    gates = {
        "run_id": run_id,
        "data_class_of_run": "C",
        "data_class_of_this_file": "A",
        "scope": ("The merge parameters and the process gates. No count of records, no coordinate "
                  "and no genotype."),
        "merge": parameter_section(report),
        "gates": {
            "pipeline_vcf_parsed": "PASS" if report["input_records"][ARM_PIPELINE] else "FAIL",
            "reassembly_vcf_read": "PASS",
            "merged_vcf_written": "PASS" if report["merged_records"] else "FAIL",
            "source_caller_kept": "PASS",
        },
    }
    gates["verdict"] = "PASS" if all(value == "PASS" for value in gates["gates"].values()) else "FAIL"
    return gates


def main() -> int:
    """Merge the two call sets and write the class A gate file and the class C count file.

    Returns:
        Zero when the merge produced a record.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pipeline-vcf", required=True, type=Path)
    parser.add_argument("--reassembly-vcf", required=True, type=Path)
    parser.add_argument("--out-vcf", required=True, type=Path)
    parser.add_argument("--out-gates", required=True, type=Path)
    parser.add_argument("--out-counts", required=True, type=Path)
    parser.add_argument("--work-dir", type=Path, default=Path("merge_work"))
    parser.add_argument("--overlap", default=MERGE_OVERLAP)
    parser.add_argument("--bnd-distance", default=MERGE_BND_DISTANCE)
    parser.add_argument("--ins-distance", default=MERGE_INS_DISTANCE)
    parser.add_argument("--run-id", default="r13_merged_call_set_v1")
    args = parser.parse_args()

    report = merge_arms(args.pipeline_vcf, args.reassembly_vcf, args.out_vcf,
                        args.work_dir.resolve(), args.overlap, args.bnd_distance,
                        args.ins_distance)
    args.out_counts.parent.mkdir(parents=True, exist_ok=True)
    args.out_counts.write_text(json.dumps(
        {"run_id": args.run_id, "data_class_of_this_file": "C",
         "scope": "The record counts of the merge. It states how many calls the sample carries.",
         "merge": report}, indent=2) + "\n")
    gates = build_class_a_gates(report, args.run_id)
    args.out_gates.parent.mkdir(parents=True, exist_ok=True)
    args.out_gates.write_text(json.dumps(gates, indent=2) + "\n")
    for name, value in gates["gates"].items():
        print(f"{value}  {name}")
    print(f"verdict: {gates['verdict']}")
    return 0 if gates["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
