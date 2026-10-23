#!/usr/bin/env python3
"""Freeze the MVA gene panel from a versioned ontology release. Data class A.

Panel version v1 answers defect R2-9 of `reports/jury/defect_ledger_track1_round2_2026-09-17.md`
and defect 6 of the third jury round. Both state that the three gene sets of this work are a manual
choice with no versioned source and no written inclusion rule.

Panel version v2 answers the gene coverage gap that the round 4 panel check found. That check
reported CENATAC and MAD1L1 as panel genes that no operating set holds, so a causal variant in
either gene earns no gene points anywhere in the submission. The cause is not the panel rule. It is
the set derivation. Two mechanisms produced the omission and v2 removes both.

    Mechanism 1, an enumerated disease list. `derive_sets` of v1 admits a gene to the established
    disease gene set only when one of its links names OMIM:257300, OMIM:614114 or OMIM:617598. That
    enumeration was written when the OMIM phenotypic series PS257300 held three subtypes. OMIM has
    since added mosaic variegated aneuploidy syndrome 4 as OMIM:620153 and syndrome 7 as
    OMIM:620189, and the frozen release carries both links. A gene whose only link is one of the
    two new identifiers therefore enters the panel and cannot enter any operating set. Every
    operating set is the established set plus a declared addition, so the enumeration is a hard
    gate in front of all three sets. v2 replaces the enumeration by the property that the
    enumeration approximated: the release links the gene to an OMIM phenotype of the mosaic
    variegated aneuploidy surface form group with the association type MENDELIAN.

    Mechanism 2, a surface form inside the code. The v1 inclusion rule matches the disease name
    against one regular expression, ``mosaic variegated aneuploidy``. A curator who names one
    disease of the group differently defeats it. The frozen release names two such diseases. It
    links SLF2 to OMIM:620184 and SMC5 to OMIM:620185 and names both diseases Atelis syndrome, and
    it links TRIM37 to OMIM:253250 and names that disease Mulibrey nanism. The v1 expression
    matched none of the three, so SLF2 and TRIM37 reached no set at all and SMC5 entered the forced
    genotyping list through a manual addition whose written reason states that the release links
    SMC5 to no matching disease. v2 moves every accepted surface form into
    `literature/track1_mva_disease_surface_forms_v1.yaml`, where each form carries its accession,
    its access date and its licence, so a new subtype name is a data change with a citation.

Five closure gates make the omission unable to recur. The builder refuses to write a panel that
fails one, and `harness/test_track1_gene_panel.py` runs the same gates.

    G1  No panel gene is outside every operating set.
    G2  Every operating set member is a panel gene or a declared addition with a reason and an
        evidence field.
    G3  Every declared consumer constant equals the operating set that the panel derives for it.
    G4  Every gene of the evidence file whose disease identifier is a disease of the matching set
        is a panel gene.
    G5  Every gene of the evidence file that carries a published case is in at least one operating
        set.

The panel is the source of record. Each operating set is the panel plus a declared restriction and
a declared addition, and the panel file states every difference. A Nextflow task runs the scorer
inside a container that holds no copy of `literature/`, so each consumer keeps its gene tuple in
the code and the test holds the tuple to the panel. Each consumer carries one tuple per panel
version and its shipped default stays at v1, because the run that produced the submission used the
v1 sets and queue row R30 blocks every patient rerun.

Usage:
    uv run python scripts/build_track1_gene_panel.py write
    uv run python scripts/build_track1_gene_panel.py check
    uv run python scripts/build_track1_gene_panel.py write --panel-version v2
    uv run python scripts/build_track1_gene_panel.py check --panel-version v2
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
RELEASE_DIR = ROOT / "tracking" / "cache" / "hpo_v2026_09_01"
SURFACE_FORM_FILE = ROOT / "literature" / "track1_mva_disease_surface_forms_v1.yaml"
EVIDENCE_FILE = ROOT / "literature" / "track1_mva_gene_evidence_v1.yaml"
RELEASE = "v2026-09-01"
DISEASE_PATTERN = re.compile(r"mosaic variegated aneuploidy", re.IGNORECASE)

ESTABLISHED_DISEASES = ("OMIM:257300", "OMIM:614114", "OMIM:617598")

MANUAL_ADDITIONS = {
    "SMC5": {
        "set": "forced_genotyping_target_list",
        "evidence": "PMID 41374403",
        "reason": ("The literature seed resolved two biallelic *SMC5* alleles in a family with "
                   "near-tetraploidy and mosaic variegated aneuploidy. The frozen release links "
                   "SMC5 to no disease whose name matches the panel rule."),
    },
    "BUB1": {
        "set": "structural_annotation_set",
        "evidence": "ORPHA:1052 of the frozen release",
        "reason": ("A spindle assembly checkpoint component. It is tier 2 of the panel, so the "
                   "frozen release supports it and the panel rule admits it."),
    },
    "BUB3": {
        "set": "structural_annotation_set",
        "evidence": "ORPHA:1052 of the frozen release",
        "reason": ("A spindle assembly checkpoint component. It is tier 2 of the panel, so the "
                   "frozen release supports it and the panel rule admits it."),
    },
    "CDC20": {
        "set": "structural_annotation_set",
        "evidence": "None in the frozen release",
        "reason": ("A spindle assembly checkpoint component and the target of the checkpoint "
                   "complex. The frozen release links CDC20 to no disease whose name matches the "
                   "panel rule, so this member rests on the manual addition rule alone."),
    },
    "PTPA": {
        "set": "structural_annotation_set",
        "evidence": "None in the frozen release",
        "reason": ("A regulator of protein phosphatase 2A, which acts on the checkpoint. The "
                   "frozen release links PTPA to no disease whose name matches the panel rule, so "
                   "this member rests on the manual addition rule alone."),
    },
}

SET_RULES = {
    "established_disease_gene_set": (
        "Every tier 1 gene of the panel whose disease is one of OMIM:257300, OMIM:614114 and "
        "OMIM:617598, which are mosaic variegated aneuploidy syndrome 1, 2 and 3. These are the "
        "three genes that a published case series reports with a biallelic genotype and a "
        "confirmed diagnosis."),
    "forced_genotyping_target_list": (
        "The established disease gene set plus every gene of the manual addition list whose set "
        "field names this list. The addition rule is a published biallelic genotype in a person "
        "with a mosaic variegated aneuploidy phenotype, with the PMID recorded beside the gene."),
    "structural_annotation_set": (
        "The established disease gene set plus every gene of the manual addition list whose set "
        "field names this list. The addition rule is membership of the spindle assembly "
        "checkpoint, because a structural event in a checkpoint component is the class of second "
        "hit that this arm looks for. *SMC5* is a structural maintenance of chromosomes protein "
        "and not a checkpoint component, so the rule excludes it from this set and admits it to "
        "the forced genotyping list."),
}

# The declared additions of panel version v2. A gene of this list reaches an operating set through
# the written addition rule of that set and not through the established clause. The `sets` field
# holds every set that the gene enters, because one gene can satisfy two addition rules.
ADDITIONS_V2 = {
    "SMC5": {
        "sets": ["forced_genotyping_target_list"],
        "evidence": "OMIM:620185 of the frozen release and PMID 41374403",
        "reason": ("A panel gene of the related chromosomal instability group. The frozen release "
                   "links SMC5 to OMIM:620185, Atelis syndrome 2, with the association type "
                   "MENDELIAN, and PMID 41374403 reports two biallelic SMC5 alleles in a family "
                   "with near-tetraploidy and mosaic variegated aneuploidy. SMC5 is a structural "
                   "maintenance of chromosomes protein and not a spindle assembly checkpoint "
                   "component, so the structural annotation rule excludes it."),
    },
    "SLF2": {
        "sets": ["forced_genotyping_target_list"],
        "evidence": "OMIM:620184 of the frozen release and PMID 36333305",
        "reason": ("A panel gene of the related chromosomal instability group. The frozen release "
                   "links SLF2 to OMIM:620184, Atelis syndrome 1, with the association type "
                   "MENDELIAN, and PMID 36333305 reports 7 patients from 5 families with a "
                   "biallelic SLF2 genotype and mosaic variegated hyperploidy. SLF2 is a "
                   "component of the RAD18-SLF1/2-SMC5/6 pathway and not a spindle assembly "
                   "checkpoint component, so the structural annotation rule excludes it. Panel "
                   "version v1 held SLF2 in no set, because its one regular expression did not "
                   "match the disease name Atelis syndrome 1."),
    },
    "TRIM37": {
        "sets": ["forced_genotyping_target_list"],
        "evidence": "OMIM:253250 of the frozen release and PMID 41702694",
        "reason": ("A panel gene of the related chromosomal instability group. The frozen release "
                   "links TRIM37 to OMIM:253250, Mulibrey nanism, with the association type "
                   "MENDELIAN, and PMID 41702694 reports mosaic variegated aneuploidy in 7 to 36 "
                   "percent of metaphases in 2 siblings with a biallelic TRIM37 genotype. TRIM37 "
                   "is a tripartite motif protein and not a spindle assembly checkpoint "
                   "component, so the structural annotation rule excludes it."),
    },
    "MAD2L1BP": {
        "sets": ["forced_genotyping_target_list", "structural_annotation_set"],
        "evidence": "PMID 37796616",
        "reason": ("A published mosaic aneuploidy gene that the frozen release links to no "
                   "disease. PMID 37796616 reports a homozygous nonsense allele in 2 siblings "
                   "with high grade mosaic variegated aneuploidy and a second homozygous "
                   "truncation in a third patient. The protein silences the spindle assembly "
                   "checkpoint by disassembling MAD2 containing complexes, so it is a checkpoint "
                   "component and the structural annotation rule admits it as well."),
    },
    "CEP192": {
        "sets": ["forced_genotyping_target_list"],
        "evidence": "PMID 37981762",
        "reason": ("A published mosaic aneuploidy gene that the frozen release links to no "
                   "disease. PMID 37981762 reports a biallelic CEP192 genotype in 2 siblings with "
                   "mosaic variegated aneuploidy and tetraploidy. CEP192 is a centrosomal "
                   "scaffold protein and not a spindle assembly checkpoint component, so the "
                   "structural annotation rule excludes it."),
    },
    "BUB1": {
        "sets": ["structural_annotation_set"],
        "evidence": "ORPHA:1052 of the frozen release",
        "reason": ("A spindle assembly checkpoint component and a tier 2 panel gene on the "
                   "Orphanet association. No published case carries a biallelic BUB1 genotype "
                   "with a mosaic variegated aneuploidy diagnosis, so the evidence file scores it "
                   "disputed and it enters this set on the checkpoint criterion alone."),
    },
    "BUB3": {
        "sets": ["structural_annotation_set"],
        "evidence": "ORPHA:1052 of the frozen release",
        "reason": ("A spindle assembly checkpoint component and a tier 2 panel gene on the "
                   "Orphanet association. No published case carries a biallelic BUB3 genotype "
                   "with a mosaic variegated aneuploidy diagnosis, so the evidence file scores it "
                   "disputed and it enters this set on the checkpoint criterion alone."),
    },
    "CDC20": {
        "sets": ["structural_annotation_set"],
        "evidence": "PMID 33094908",
        "reason": ("A spindle assembly checkpoint component and the target of the checkpoint "
                   "complex. The frozen release links CDC20 to no disease of the matching set. "
                   "PMID 33094908 reports one de novo monoallelic missense allele with random "
                   "chromosome number instability, which does not fit the recessive inheritance "
                   "of the panel diseases, so this member rests on the checkpoint criterion."),
    },
    "PTPA": {
        "sets": ["structural_annotation_set"],
        "evidence": "None in the frozen release and no published human case",
        "reason": ("A regulator of protein phosphatase 2A, which acts on the checkpoint. The "
                   "frozen release links PTPA to no disease of the matching set and no published "
                   "human case carries a PTPA genotype with a phenotype of this class, so this "
                   "member rests on the checkpoint criterion alone."),
    },
}

SET_RULES_V2 = {
    "established_disease_gene_set": (
        "Every gene that the frozen release links with the association type MENDELIAN to an OMIM "
        "phenotype whose name matches a surface form group of "
        "literature/track1_mva_disease_surface_forms_v1.yaml whose "
        "admits_to_the_established_disease_gene_set field is true. The order is the ascending OMIM "
        "entry number of the matching phenotype of each gene. Panel version v1 wrote this rule as "
        "the enumeration OMIM:257300, OMIM:614114 and OMIM:617598, and that enumeration is the "
        "mechanism that held CENATAC and MAD1L1 out of every operating set."),
    "forced_genotyping_target_list": (
        "The established disease gene set plus every gene of the declared addition list whose sets "
        "field names this list. The addition rule is a published biallelic genotype in a person "
        "with a mosaic variegated aneuploidy or a mosaic variegated hyperploidy phenotype, with "
        "the PMID recorded beside the gene in "
        "literature/track1_mva_gene_evidence_v1.yaml."),
    "structural_annotation_set": (
        "The established disease gene set plus every gene of the declared addition list whose sets "
        "field names this list. The addition rule is membership of the spindle assembly "
        "checkpoint, because a structural event in a checkpoint component is the class of second "
        "hit that this arm looks for. SMC5, SLF2, CEP192 and TRIM37 are not checkpoint components, "
        "so the rule excludes each one from this set and admits each one to the forced genotyping "
        "list. MAD2L1BP silences the checkpoint, so the rule admits it to both."),
    "benchmark_mva_gene_set": (
        "The established disease gene set. scripts/build_spikein_suite.py reads it to decide which "
        "gene of the public spike-in suite carries the mva_clinvar and the mva_literature group "
        "label and which gene is eligible for the recessive_clinvar group. Panel version v1 held "
        "this list as an undeclared third copy of the three established genes, so the panel test "
        "did not check it and the suite carried no CENATAC and no MAD1L1 scenario."),
}

CONSUMERS_V1 = {
    "pipelines/sv_mei_annotate/bin/annotate_sv_mei.py": {
        "constant": "MVA_GENES",
        "set": "structural_annotation_set",
    },
    "scripts/build_s6_targets.py": {
        "constant": "GENES",
        "set": "forced_genotyping_target_list",
    },
}

CONSUMERS_V2 = {
    "pipelines/sv_mei_annotate/bin/annotate_sv_mei.py": {
        "constant": "PANEL_SETS",
        "set": "structural_annotation_set",
    },
    "scripts/build_s6_targets.py": {
        "constant": "PANEL_SETS",
        "set": "forced_genotyping_target_list",
    },
    "scripts/build_r12_targets.py": {
        "constant": "PANEL_SETS",
        "set": "structural_annotation_set",
    },
    "scripts/track1_sv_benchmark.py": {
        "constant": "PANEL_SETS",
        "set": "structural_annotation_set",
    },
    "scripts/build_spikein_suite.py": {
        "constant": "PANEL_SETS",
        "set": "benchmark_mva_gene_set",
    },
}


def panel_file(version: str) -> Path:
    """Give the panel file of one panel version.

    Args:
        version: The panel version, ``v1`` or ``v2``.

    Returns:
        The panel file path.
    """
    return ROOT / "literature" / f"track1_gene_panel_{version}.yaml"


def receipt_file(version: str) -> Path:
    """Give the receipt of one panel version.

    Args:
        version: The panel version, ``v1`` or ``v2``.

    Returns:
        The receipt path.
    """
    return ROOT / "tracking" / "receipts" / f"t1_gene_panel_{version}.json"


def verify_release() -> dict[str, str]:
    """Check the two release files against the checksum file of the release.

    Returns:
        A map of file name to its SHA-256.

    Raises:
        SystemExit: When a file is missing or its checksum does not match the record.
    """
    checksums = RELEASE_DIR / "sha256sums.txt"
    if not checksums.exists():
        raise SystemExit(f"No checksum file at {checksums}.")
    recorded: dict[str, str] = {}
    for line in checksums.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2:
            recorded[parts[1]] = parts[0]
    verified: dict[str, str] = {}
    for name in ("genes_to_disease.txt", "phenotype.hpoa"):
        path = RELEASE_DIR / name
        if not path.exists():
            raise SystemExit(f"No release file at {path}.")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if recorded.get(name) != digest:
            raise SystemExit(f"{name} does not match the checksum of the release record.")
        verified[name] = digest
    return verified


def load_surface_forms() -> dict:
    """Read the accepted disease surface forms of panel version v2.

    Returns:
        The surface form document.

    Raises:
        SystemExit: When the file is absent.
    """
    if not SURFACE_FORM_FILE.exists():
        raise SystemExit(f"No surface form file at {SURFACE_FORM_FILE}.")
    with SURFACE_FORM_FILE.open() as handle:
        return yaml.safe_load(handle)


def load_evidence() -> dict:
    """Read the verified gene evidence of panel version v2.

    Returns:
        The evidence document.

    Raises:
        SystemExit: When the file is absent.
    """
    if not EVIDENCE_FILE.exists():
        raise SystemExit(f"No evidence file at {EVIDENCE_FILE}. Run "
                         f"scripts/fetch_mva_gene_evidence.py write first.")
    with EVIDENCE_FILE.open() as handle:
        return yaml.safe_load(handle)


def read_matching_diseases() -> dict[str, str]:
    """Find every disease of the release whose name matches the v1 panel rule.

    Returns:
        A map of database identifier to disease name.
    """
    diseases: dict[str, str] = {}
    with (RELEASE_DIR / "phenotype.hpoa").open() as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("database_id"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            if DISEASE_PATTERN.search(fields[1]):
                diseases[fields[0]] = fields[1]
    return diseases


def read_matching_diseases_v2(surface_forms: dict) -> tuple[dict[str, str], dict[str, str]]:
    """Find every disease of the release whose name matches a declared surface form.

    Args:
        surface_forms: The surface form document.

    Returns:
        A map of database identifier to disease name, and a map of database identifier to the name
        of the surface form group that matched it.
    """
    patterns = {name: re.compile(group["pattern"], re.IGNORECASE)
                for name, group in surface_forms["groups"].items()}
    diseases: dict[str, str] = {}
    groups: dict[str, str] = {}
    with (RELEASE_DIR / "phenotype.hpoa").open() as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("database_id"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            for name, pattern in patterns.items():
                if pattern.search(fields[1]):
                    diseases[fields[0]] = fields[1]
                    groups.setdefault(fields[0], name)
                    break
    return diseases, groups


def read_panel(diseases: dict[str, str]) -> dict[str, dict]:
    """Map every matching disease to its genes in the release.

    Args:
        diseases: The matching diseases as a map of identifier to name.

    Returns:
        A map of gene symbol to its tier, its NCBI gene identifier and its links.
    """
    panel: dict[str, dict] = {}
    with (RELEASE_DIR / "genes_to_disease.txt").open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or fields[0] == "ncbi_gene_id":
                continue
            gene_id, symbol, association, disease_id, source = fields[:5]
            if disease_id not in diseases:
                continue
            entry = panel.setdefault(symbol, {"ncbi_gene_id": gene_id, "tier": 2, "links": []})
            entry["links"].append({"disease_id": disease_id, "disease_name": diseases[disease_id],
                                   "association_type": association, "source": source})
            if association == "MENDELIAN":
                entry["tier"] = 1
    return panel


def derive_sets(panel: dict[str, dict]) -> dict[str, list[str]]:
    """Derive the three operating sets of panel version v1 from the panel.

    Args:
        panel: The frozen panel.

    Returns:
        A map of set name to its members in the order that the code uses.
    """
    established = [symbol for symbol, entry in panel.items()
                   if entry["tier"] == 1
                   and any(link["disease_id"] in ESTABLISHED_DISEASES for link in entry["links"])]
    order = {"BUB1B": 0, "CEP57": 1, "TRIP13": 2}
    established.sort(key=lambda symbol: order.get(symbol, 99))
    forced = established + [gene for gene, spec in MANUAL_ADDITIONS.items()
                            if spec["set"] == "forced_genotyping_target_list"]
    structural = established + [gene for gene, spec in MANUAL_ADDITIONS.items()
                                if spec["set"] == "structural_annotation_set"]
    return {"established_disease_gene_set": established,
            "forced_genotyping_target_list": forced,
            "structural_annotation_set": structural}


def established_order_key(entry: dict, admitting: set[str], groups: dict[str, str]) -> int:
    """Give the ordering key of one established gene.

    Args:
        entry: The panel entry of the gene.
        admitting: The surface form groups that admit a gene to the established set.
        groups: The surface form group of every matching disease.

    Returns:
        The smallest OMIM entry number among the admitting MENDELIAN links of the gene.
    """
    numbers = [int(link["disease_id"].split(":")[1]) for link in entry["links"]
               if link["association_type"] == "MENDELIAN"
               and link["disease_id"].startswith("OMIM:")
               and groups.get(link["disease_id"]) in admitting]
    return min(numbers) if numbers else 10 ** 9


def derive_sets_v2(panel: dict[str, dict], surface_forms: dict,
                   groups: dict[str, str]) -> dict[str, list[str]]:
    """Derive the four operating sets of panel version v2 from the panel.

    Args:
        panel: The frozen panel.
        surface_forms: The surface form document.
        groups: The surface form group of every matching disease.

    Returns:
        A map of set name to its members in the order that the code uses.
    """
    admitting = {name for name, group in surface_forms["groups"].items()
                 if group["admits_to_the_established_disease_gene_set"]}
    established = [symbol for symbol, entry in panel.items()
                   if established_order_key(entry, admitting, groups) < 10 ** 9]
    established.sort(key=lambda symbol: established_order_key(panel[symbol], admitting, groups))
    forced = established + [gene for gene, spec in ADDITIONS_V2.items()
                            if "forced_genotyping_target_list" in spec["sets"]]
    structural = established + [gene for gene, spec in ADDITIONS_V2.items()
                                if "structural_annotation_set" in spec["sets"]]
    return {"established_disease_gene_set": established,
            "forced_genotyping_target_list": forced,
            "structural_annotation_set": structural,
            "benchmark_mva_gene_set": list(established)}


def consumer_constants(consumers: dict[str, dict], version: str) -> dict[str, tuple[str, ...]]:
    """Read the gene tuple that each declared consumer carries for one panel version.

    Args:
        consumers: The declared consumers of the panel version.
        version: The panel version, ``v1`` or ``v2``.

    Returns:
        A map of consumer path to the tuple that its constant holds, for the consumers that this
        process can import.
    """
    sys.path.insert(0, str(ROOT / "pipelines" / "sv_mei_annotate" / "bin"))
    sys.path.insert(0, str(ROOT / "scripts"))
    modules = {
        "pipelines/sv_mei_annotate/bin/annotate_sv_mei.py": "annotate_sv_mei",
        "scripts/build_s6_targets.py": "build_s6_targets",
        "scripts/build_r12_targets.py": "build_r12_targets",
        "scripts/track1_sv_benchmark.py": "track1_sv_benchmark",
        "scripts/build_spikein_suite.py": "build_spikein_suite",
    }
    found: dict[str, tuple[str, ...]] = {}
    for path in consumers:
        name = modules.get(path)
        if name is None:
            continue
        module = __import__(name)
        sets = getattr(module, "PANEL_SETS", None)
        if version == "v1" and sets is None:
            legacy = getattr(module, consumers[path]["constant"], None)
            if legacy is not None:
                found[path] = tuple(legacy)
            continue
        if sets is None or version not in sets:
            continue
        found[path] = tuple(sets[version])
    return found


def closure_gates(panel: dict[str, dict], sets: dict[str, list[str]], diseases: dict[str, str],
                  evidence: dict, consumers: dict[str, dict],
                  version: str) -> list[str]:
    """Run the five closure gates of panel version v2.

    Args:
        panel: The frozen panel.
        sets: The derived operating sets.
        diseases: The matching diseases.
        evidence: The verified evidence document.
        consumers: The declared consumers.
        version: The panel version.

    Returns:
        The failure messages, empty when every gate holds.
    """
    failures: list[str] = []
    in_sets: set[str] = set()
    for members in sets.values():
        in_sets |= set(members)
    silent = sorted(set(panel) - in_sets)
    if silent:
        failures.append(f"G1 a panel gene is outside every operating set: {silent}. A causal "
                        f"variant in such a gene earns no gene points anywhere in the submission.")
    for name, members in sets.items():
        for gene in members:
            if gene in panel:
                continue
            spec = ADDITIONS_V2.get(gene)
            if spec is None:
                failures.append(f"G2 {gene} is in {name} and it is neither a panel gene nor a "
                                f"declared addition")
            elif not spec.get("reason") or not spec.get("evidence"):
                failures.append(f"G2 the declared addition {gene} carries no reason or no "
                                f"evidence")
    constants = consumer_constants(consumers, version)
    for path, spec in consumers.items():
        expected = tuple(sets[spec["set"]])
        found = constants.get(path)
        if found is None:
            failures.append(f"G3 {path} was not read, so its constant was not checked")
        elif found != expected:
            failures.append(f"G3 {path} {spec['constant']}[{version}] is {found} and the panel "
                            f"derives {expected} for {spec['set']}")
    for symbol, entry in evidence["genes"].items():
        identifier = str(entry.get("disease_id") or "")
        if identifier and identifier in diseases and symbol not in panel:
            failures.append(f"G4 the evidence file links {symbol} to {identifier}, which is a "
                            f"disease of the matching set, and the panel does not hold it")
        if entry.get("first_pmid") and entry["individuals_in_the_first_report"] > 0 \
                and symbol not in in_sets:
            failures.append(f"G5 the evidence file carries a published case for {symbol} at PMID "
                            f"{entry['first_pmid']} and no operating set holds it")
    return failures


def build_document(verified: dict[str, str], diseases: dict[str, str],
                   panel: dict[str, dict], sets: dict[str, list[str]]) -> dict:
    """Build the panel document of version v1.

    Args:
        verified: The checksums of the release files.
        diseases: The matching diseases.
        panel: The frozen panel.
        sets: The three derived operating sets.

    Returns:
        The document that the panel file holds.
    """
    panel_genes = sorted(panel)
    return {
        "panel_id": "track1_gene_panel_v1",
        "inclusion_rule": (
            "A gene enters this panel when the Human Phenotype Ontology release " + RELEASE +
            " links it to a disease whose name matches the regular expression "
            "'mosaic variegated aneuploidy', case insensitive. The disease names come from "
            "phenotype.hpoa of that release and the gene to disease links from "
            "genes_to_disease.txt of the same release. The association_type of the link gives the "
            "tier: MENDELIAN is tier 1 and any other value is tier 2."),
        "source": {
            "name": "Human Phenotype Ontology",
            "release": RELEASE,
            "directory": "tracking/cache/hpo_v2026_09_01",
            "files": verified,
            "citation": "PMID 37953324",
        },
        "generated_by": "scripts/build_track1_gene_panel.py",
        "matching_diseases": dict(sorted(diseases.items())),
        "panel": {symbol: panel[symbol] for symbol in panel_genes},
        "set_rules": SET_RULES,
        "manual_additions": MANUAL_ADDITIONS,
        "operating_sets": sets,
        "differences_against_the_panel": {
            "in_the_panel_and_in_no_operating_set": sorted(
                set(panel_genes) - set().union(*[set(members) for members in sets.values()])),
            "in_an_operating_set_and_not_in_the_panel": sorted(
                set().union(*[set(members) for members in sets.values()]) - set(panel_genes)),
        },
        "consumers": CONSUMERS_V1,
        "test": "harness/test_track1_gene_panel.py",
    }


def build_document_v2(verified: dict[str, str], diseases: dict[str, str],
                      groups: dict[str, str], panel: dict[str, dict],
                      sets: dict[str, list[str]], surface_forms: dict,
                      evidence: dict) -> dict:
    """Build the panel document of version v2.

    Args:
        verified: The checksums of the release files.
        diseases: The matching diseases.
        groups: The surface form group of every matching disease.
        panel: The frozen panel.
        sets: The four derived operating sets.
        surface_forms: The surface form document.
        evidence: The verified evidence document.

    Returns:
        The document that the panel file holds.
    """
    panel_genes = sorted(panel)
    in_sets: set[str] = set()
    for members in sets.values():
        in_sets |= set(members)
    return {
        "panel_id": "track1_gene_panel_v2",
        "supersedes": "track1_gene_panel_v1",
        "inclusion_rule": (
            "A gene enters this panel when the Human Phenotype Ontology release " + RELEASE +
            " links it to a disease whose name matches a declared surface form of "
            "literature/track1_mva_disease_surface_forms_v1.yaml, case insensitive. The disease "
            "names come from phenotype.hpoa of that release and the gene to disease links from "
            "genes_to_disease.txt of the same release. The association_type of the link gives the "
            "tier: MENDELIAN is tier 1 and any other value is tier 2. Panel version v1 held one "
            "regular expression inside the builder, 'mosaic variegated aneuploidy'. That form "
            "missed OMIM:620184 and OMIM:620185, which the same release names Atelis syndrome 1 "
            "and Atelis syndrome 2, and OMIM:253250, which it names Mulibrey nanism."),
        "correction": (
            "The round 4 panel check reported CENATAC and MAD1L1 as panel genes that no operating "
            "set holds. The cause is the established clause of version v1, which enumerated "
            "OMIM:257300, OMIM:614114 and OMIM:617598 and therefore excluded every mosaic "
            "variegated aneuploidy subtype that OMIM added after subtype 3. Version v2 replaces "
            "the enumeration by the property that the release links the gene with the association "
            "type MENDELIAN to an OMIM phenotype of a surface form group that admits to the "
            "established set. Version v2 also declares three further consumers, because the gene "
            "tuple of scripts/build_spikein_suite.py, scripts/build_r12_targets.py and "
            "scripts/track1_sv_benchmark.py was an undeclared copy that no test checked."),
        "source": {
            "name": "Human Phenotype Ontology",
            "release": RELEASE,
            "directory": "tracking/cache/hpo_v2026_09_01",
            "files": verified,
            "citation": "PMID 37953324",
        },
        "surface_forms": {
            "file": str(SURFACE_FORM_FILE.relative_to(ROOT)),
            "sha256": hashlib.sha256(SURFACE_FORM_FILE.read_bytes()).hexdigest(),
            "groups": {name: {"pattern": group["pattern"],
                              "admits_to_the_established_disease_gene_set":
                                  group["admits_to_the_established_disease_gene_set"]}
                       for name, group in surface_forms["groups"].items()},
        },
        "evidence": {
            "file": str(EVIDENCE_FILE.relative_to(ROOT)),
            "sha256": hashlib.sha256(EVIDENCE_FILE.read_bytes()).hexdigest(),
            "evidence_class_rules": evidence["evidence_class_rules"],
            "genes": {symbol: {"disease": entry["disease"],
                               "disease_id": entry["disease_id"],
                               "inheritance": entry["inheritance"],
                               "families_in_the_first_report":
                                   entry["families_in_the_first_report"],
                               "individuals_in_the_first_report":
                                   entry["individuals_in_the_first_report"],
                               "first_pmid": entry["first_pmid"],
                               "evidence_class": entry["evidence_class"],
                               "clingen": entry["clingen"]}
                      for symbol, entry in sorted(evidence["genes"].items())},
        },
        "generated_by": "scripts/build_track1_gene_panel.py",
        "matching_diseases": dict(sorted(diseases.items())),
        "matching_disease_groups": dict(sorted(groups.items())),
        "panel": {symbol: panel[symbol] for symbol in panel_genes},
        "set_rules": SET_RULES_V2,
        "declared_additions": ADDITIONS_V2,
        "operating_sets": sets,
        "differences_against_the_panel": {
            "in_the_panel_and_in_no_operating_set": sorted(set(panel_genes) - in_sets),
            "in_an_operating_set_and_not_in_the_panel": sorted(in_sets - set(panel_genes)),
        },
        "closure_gates": {
            "G1": "No panel gene is outside every operating set.",
            "G2": "Every operating set member is a panel gene or a declared addition with a "
                  "reason and an evidence field.",
            "G3": "Every declared consumer constant equals the operating set that the panel "
                  "derives for it.",
            "G4": "Every gene of the evidence file whose disease identifier is a disease of the "
                  "matching set is a panel gene.",
            "G5": "Every gene of the evidence file that carries a published case is in at least "
                  "one operating set.",
        },
        "consumers": CONSUMERS_V2,
        "shipped_panel_version": {
            "version": "v1",
            "reason": "Every consumer keeps its shipped default at the v1 tuple, because the runs "
                      "that produced the submission used the v1 sets and queue row R30 blocks "
                      "every patient rerun until the compliance question closes. The v2 tuple is "
                      "in the code, the test holds it to this panel, and the change of the "
                      "shipped default is an owner decision that queue row R33 records.",
        },
        "test": "harness/test_track1_gene_panel.py",
    }


def load_committed(version: str) -> dict | None:
    """Read the committed panel file of one panel version.

    Args:
        version: The panel version.

    Returns:
        The committed document, or None when the file is absent.
    """
    path = panel_file(version)
    if not path.exists():
        return None
    with path.open() as handle:
        return yaml.safe_load(handle)


def panel_header(version: str) -> str:
    """Give the comment header of one panel file.

    Args:
        version: The panel version.

    Returns:
        The header text.
    """
    if version == "v1":
        return (
            "# The MVA gene panel of Track 1, frozen from one ontology release.\n"
            "#\n"
            "# scripts/build_track1_gene_panel.py writes this file and `check` refuses a file "
            "that the\n"
            "# frozen release no longer produces. harness/test_track1_gene_panel.py holds the two "
            "code\n"
            "# constants to the operating sets below, so a gene set cannot move in the code "
            "alone.\n"
            "#\n"
            "# Do not edit this file by hand. Edit the rule or the manual addition list in the "
            "script.\n"
        )
    return (
        "# The MVA gene panel of Track 1, version 2, frozen from one ontology release.\n"
        "#\n"
        "# This file supersedes literature/track1_gene_panel_v1.yaml. Version v1 stays in the\n"
        "# repository and it still reproduces, because every receipt of the recorded runs cites "
        "it.\n"
        "#\n"
        "# scripts/build_track1_gene_panel.py --panel-version v2 writes this file and `check`\n"
        "# refuses a file that the frozen release, the surface form file and the evidence file no\n"
        "# longer produce. Five closure gates run on every write and on every check, and\n"
        "# harness/test_track1_gene_panel.py runs the same gates.\n"
        "#\n"
        "# Do not edit this file by hand. Edit the surface form file, the gene seed or the\n"
        "# declared addition list in the script.\n"
    )


def write_panel(document: dict, version: str, gate_report: list[str]) -> None:
    """Write one panel file and its receipt.

    Args:
        document: The panel document.
        version: The panel version.
        gate_report: The closure gate failures, which is empty on a write.
    """
    path = panel_file(version)
    with path.open("w") as handle:
        handle.write(panel_header(version))
        yaml.safe_dump(document, handle, sort_keys=False, default_flow_style=False, width=98)
    receipt = {
        "receipt_id": f"t1_gene_panel_{version}",
        "created_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "track": "t1_call",
        "data_class": "A",
        "queue_item": "R29" if version == "v1" else "R32",
        "control": "a_versioned_gene_panel_and_a_written_inclusion_rule",
        "answers": ("Defect R2-9 of round 2 and defect 6 of round 3, which state that the three "
                    "gene sets carry no versioned source and no written inclusion rule.")
        if version == "v1" else
        ("The gene coverage gap of the round 4 panel check, which reported CENATAC and MAD1L1 as "
         "panel genes that no operating set holds."),
        "panel_file": str(path.relative_to(ROOT)),
        "panel_file_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "builder": "scripts/build_track1_gene_panel.py",
        "builder_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source": document["source"],
        "inclusion_rule": document["inclusion_rule"],
        "matching_diseases": document["matching_diseases"],
        "panel_size": len(document["panel"]),
        "panel_tier_1": sorted(symbol for symbol, entry in document["panel"].items()
                               if entry["tier"] == 1),
        "panel_tier_2": sorted(symbol for symbol, entry in document["panel"].items()
                               if entry["tier"] == 2),
        "operating_sets": document["operating_sets"],
        "operating_set_sizes": {name: len(members)
                                for name, members in document["operating_sets"].items()},
        "matching_disease_count": len(document["matching_diseases"]),
        "differences_against_the_panel": document["differences_against_the_panel"],
        "shipped_configuration_changed": False,
    }
    if version == "v1":
        receipt["note"] = (
            "The panel was frozen after the patient run, so it changes no shipped value. The "
            "structural annotation set and the forced genotyping target list hold the same "
            "members as the runs that this report records, and the difference against the "
            "panel is stated rather than repaired by moving a shipped set after the result.")
    else:
        v1_sets = (load_committed("v1") or {}).get("operating_sets", {})
        receipt.update({
            "supersedes": "track1_gene_panel_v1",
            "correction": document["correction"],
            "surface_forms": document["surface_forms"],
            "evidence": {"file": document["evidence"]["file"],
                         "sha256": document["evidence"]["sha256"]},
            "closure_gates": document["closure_gates"],
            "closure_gate_count": len(document["closure_gates"]),
            "closure_gate_failures": gate_report,
            "closure_gate_verdict": "PASS" if not gate_report else "FAIL",
            "surface_form_group_count": len(document["surface_forms"]["groups"]),
            "matching_disease_count_per_group": {
                group: sum(1 for name in document["matching_disease_groups"].values()
                           if name == group)
                for group in sorted(set(document["matching_disease_groups"].values()))},
            "consumer_count": len(document["consumers"]),
            "panel_gene_count_per_evidence_class": {
                name: sum(1 for symbol in document["panel"]
                          if document["evidence"]["genes"].get(symbol, {})
                          .get("evidence_class") == name)
                for name in sorted({entry["evidence_class"]
                                    for entry in document["evidence"]["genes"].values()})},
            "genes_added_to_each_set_against_v1": {
                name: sorted(set(members) - set(v1_sets.get(name, [])))
                for name, members in document["operating_sets"].items()},
            "genes_removed_from_each_set_against_v1": {
                name: sorted(set(v1_sets.get(name, [])) - set(members))
                for name, members in document["operating_sets"].items()},
            "evidence_class_of_each_panel_gene": {
                symbol: document["evidence"]["genes"][symbol]["evidence_class"]
                for symbol in sorted(document["panel"])
                if symbol in document["evidence"]["genes"]},
            "shipped_panel_version": document["shipped_panel_version"],
            "note": (
                "The corrected panel changes no shipped value. Every consumer keeps its shipped "
                "default at the v1 tuple, and the v2 tuple reaches a run through the panel "
                "version argument of that consumer. Queue row R33 records what the owner must "
                "authorize for a patient rerun under this panel."),
        })
    with receipt_file(version).open("w") as handle:
        json.dump(receipt, handle, indent=2)
        handle.write("\n")


def build(version: str) -> tuple[dict, dict[str, dict], dict[str, list[str]], list[str]]:
    """Build one panel version from the frozen release.

    Args:
        version: The panel version, ``v1`` or ``v2``.

    Returns:
        The panel document, the panel, the operating sets and the closure gate failures.

    Raises:
        SystemExit: When no disease of the release matches the rule.
    """
    verified = verify_release()
    if version == "v1":
        diseases = read_matching_diseases()
        if not diseases:
            raise SystemExit("No disease of the release matches the panel rule.")
        panel = read_panel(diseases)
        sets = derive_sets(panel)
        return build_document(verified, diseases, panel, sets), panel, sets, []
    surface_forms = load_surface_forms()
    evidence = load_evidence()
    diseases, groups = read_matching_diseases_v2(surface_forms)
    if not diseases:
        raise SystemExit("No disease of the release matches a declared surface form.")
    panel = read_panel(diseases)
    sets = derive_sets_v2(panel, surface_forms, groups)
    document = build_document_v2(verified, diseases, groups, panel, sets, surface_forms, evidence)
    failures = closure_gates(panel, sets, diseases, evidence, CONSUMERS_V2, "v2")
    return document, panel, sets, failures


def main() -> int:
    """Write or check one frozen gene panel version.

    Returns:
        Zero on success, one when the committed panel differs from the frozen release or a closure
        gate fails.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=("write", "check"), default="check", nargs="?")
    parser.add_argument("--panel-version", choices=("v1", "v2"), default="v1",
                        help="The panel version to write or to check. The default is v1.")
    args = parser.parse_args()
    version = args.panel_version
    document, panel, sets, failures = build(version)
    if failures:
        for message in failures:
            print(f"FAIL: {message}", file=sys.stderr)
        print(f"{len(failures)} closure gates of panel {version} failed.", file=sys.stderr)
        return 1
    if args.step == "write":
        write_panel(document, version, failures)
        sizes = ", ".join(f"{len(members)} {name}" for name, members in sets.items())
        print(f"{panel_file(version).relative_to(ROOT)}: {len(panel)} genes in the panel, {sizes}")
        return 0
    committed = load_committed(version)
    if committed is None:
        print(f"No panel file at {panel_file(version)}.", file=sys.stderr)
        return 1
    keys = ["inclusion_rule", "panel", "operating_sets", "set_rules", "matching_diseases",
            "differences_against_the_panel", "consumers"]
    keys.append("manual_additions" if version == "v1" else "declared_additions")
    if version == "v2":
        keys += ["surface_forms", "evidence", "matching_disease_groups", "closure_gates",
                 "correction"]
    for key in keys:
        if committed.get(key) != document.get(key):
            print(f"The committed panel {version} differs from the frozen release on {key}.",
                  file=sys.stderr)
            return 1
    if committed.get("source", {}).get("files") != document["source"]["files"]:
        print("The release files no longer match the checksums of the panel file.", file=sys.stderr)
        return 1
    print(f"OK: the committed panel {version} reproduces from {RELEASE} with {len(panel)} genes "
          f"and every closure gate passes.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
