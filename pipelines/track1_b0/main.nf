#!/usr/bin/env nextflow
// Track 1 B0 pipeline (ADR 0003).
// mode suite:   VEP on the HG002 background and the spike-in VCF, then B0 on 50 scenarios with the official scorer.
// mode patient: VEP on the S5 patient VCF and the S6 forced genotypes, then B0 and a slot 1 file with a format check.
// mode phenotype_patient: map the clinical phenotype document to HPO terms, then B0 with the phenotype layer on the
//                         already annotated patient calls of run b0_patient_v1. The term list and the rows are class B2.

nextflow.enable.dsl = 2

params.mode          = null
params.vep_cache     = null
params.literature    = null
params.outdir        = null
// suite
params.background    = null
params.spikes        = null
params.manifest      = null
params.min_mean_rank_points = 25
params.min_pair_top10       = 0.5
// patient
params.patient_vcf   = null
params.forced_vcf    = null
// phenotype_patient
params.annotated_patient_vcf = null
params.annotated_forced_vcf  = null
params.phenotype_document    = null
params.label                 = null
params.hpo_obo               = null
params.genes_to_phenotype    = null
// The mismatched phenotype control of docs/benchmarks/t1_phenotype_mismatch_v1.md returned
// GATE NOT MET on 2026-09-17, so the layer ships at weight 0.0. Receipt
// tracking/receipts/t1_phenotype_mismatch_v1.json.
params.phenotype_weight      = 0.0
params.max_rows              = 3
params.epcr_top              = 0.99
params.epcr_step             = 0.01
params.epcr_digits           = 2
params.fasta                 = null
params.fai                   = null

process VEP_ANNOTATE {
    tag "${label}"
    container 'docker.io/ensemblorg/ensembl-vep:release_116.2'
    cpus 32
    memory '64 GB'
    time '12h'

    input:
    tuple val(label), path(vcf)
    path cache_tar

    output:
    tuple val(label), path("${label}.vep.vcf.gz")

    script:
    """
    mkdir -p cache
    tar -xzf ${cache_tar} -C cache
    vep --offline --cache --dir_cache cache --cache_version 116 --assembly GRCh38 \\
        --input_file ${vcf} --output_file ${label}.vep.vcf.gz --vcf --compress_output bgzip \\
        --fork ${task.cpus} --buffer_size 20000 --no_stats \\
        --symbol --canonical --mane_select --af_gnomadg --af_gnomade --check_existing --flag_pick_allele_gene \\
        --fields "Allele,Consequence,IMPACT,SYMBOL,Gene,PICK,MANE_SELECT,CANONICAL,gnomADg_AF,gnomADe_AF,CLIN_SIG,Existing_variation"
    rm -rf cache
    """
}

process PUBLISH_ANNOTATION {
    tag "${label}"
    container 'docker.io/library/python:3.12-bookworm'
    cpus 1
    memory '2 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(label), path(annotated)

    output:
    path "annotation/*"

    script:
    """
    mkdir -p annotation
    cp ${annotated} annotation/
    python3 -c "import hashlib, json, sys; d = hashlib.sha256(open(sys.argv[1], 'rb').read()).hexdigest(); print(json.dumps({'file': sys.argv[1], 'sha256': d, 'data_class': 'A'}, indent=2))" ${annotated} > annotation/${label}_annotation_gate.json
    """
}

process B0_SUITE {
    container 'docker.io/library/python:3.12-bookworm'
    cpus 4
    memory '32 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    path background
    path spikes
    path manifest
    path literature

    output:
    path "suite/**"

    script:
    """
    mkdir -p suite
    b0_suite.py --background ${background} --spikes ${spikes} --manifest ${manifest} --literature ${literature} \\
        --min-mean-rank-points ${params.min_mean_rank_points} --min-pair-top10 ${params.min_pair_top10} --out-dir suite
    """
}

process B0_PATIENT {
    container 'docker.io/library/python:3.12-bookworm'
    cpus 4
    memory '32 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    path patient
    path forced
    path literature
    path guard
    path publish_rules

    output:
    path "slot1/**"

    script:
    """
    mkdir -p slot1
    b0_rank.py --vcf ${patient} --extra-vcf ${forced} --literature ${literature} \\
        --out slot1/slot1_b0.csv --summary slot1/b0_summary.json
    BIN_DIR=\$(dirname "\$(command -v b0_rank.py)")
    BIN_DIR="\$BIN_DIR" python3 - <<'PY'
import hashlib, json, os, sys
sys.path.insert(0, os.environ["BIN_DIR"])
import hackathon_evaluation as evaluation
path = "slot1/slot1_b0.csv"
budget = int("${params.max_rows}")
check = {"file": path, "sha256": hashlib.sha256(open(path, "rb").read()).hexdigest(),
         "data_class": "A", "scope": "The shape of the file, not its content.",
         "row_budget": budget}
try:
    submission = evaluation.load_submission(path)
    rows = submission.get("PROBAND01", [])
    # Control 7 of docs/llm_data_privacy.md. The number of rows is a quantity of the child,
    # because a value below the budget states how many candidates the ranker found. This file
    # carries the verdict against the budget and not the count.
    check.update(format_check="PASS", proband_ids=sorted(submission),
                 rows_present="PASS" if rows else "FAIL",
                 rows_within_budget="PASS" if len(rows) <= budget else "FAIL")
except Exception as error:
    check.update(format_check="FAIL", error=type(error).__name__ + ": " + str(error)[:200])
json.dump(check, open("slot1/format_check.json", "w"), indent=2)
PY
    python3 ${guard} --class-a-publish --rules ${publish_rules} slot1/format_check.json
    """
}

process EXTRACT_PHENOTYPE {
    tag "case_phenotype"
    container 'docker.io/library/python:3.12-bookworm'
    cpus 2
    memory '8 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    path document
    path obo
    path guard
    path publish_rules

    output:
    path "phenotype/case_hpo_terms.tsv", emit: terms
    path "phenotype/phenotype_extraction_gate.json", emit: gate

    script:
    """
    mkdir -p phenotype
    extract_case_phenotype.py --document ${document} --obo ${obo} \\
        --out-terms phenotype/case_hpo_terms.tsv \\
        --out-gate phenotype/phenotype_extraction_gate.json
    python3 ${guard} --class-a-publish --rules ${publish_rules} \\
        phenotype/phenotype_extraction_gate.json
    """
}

process B0_PATIENT_PHENOTYPE {
    tag "b0_phenotype"
    // The SVDB container of the patient run carries Python 3.13 and bcftools, so the ranker and the
    // R14 normalization check run in one task with one recorded image digest.
    container 'community.wave.seqera.io/library/bcftools_svdb:12db401acbacc624'
    cpus 4
    memory '32 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    path patient
    path forced
    path literature
    path obo
    path genes_to_phenotype
    path terms
    path fasta
    path fai
    path guard
    path publish_rules

    output:
    path "slot1/**"

    script:
    """
    mkdir -p slot1
    b0_rank.py --vcf ${patient} --extra-vcf ${forced} --literature ${literature} \\
        --hpo-obo ${obo} --genes-to-phenotype ${genes_to_phenotype} --case-phenotype ${terms} \\
        --phenotype-weight ${params.phenotype_weight} \\
        --max-rows ${params.max_rows} --epcr-top ${params.epcr_top} \\
        --epcr-step ${params.epcr_step} --epcr-digits ${params.epcr_digits} \\
        --reference ${fasta} --norm-backend bcftools --bcftools bcftools \\
        --norm-report slot1/normalization_report.json \\
        --phenotype-report slot1/phenotype_layer.json \\
        --out slot1/slot1_b0_phenotype.csv --summary slot1/b0_summary.json
    BIN_DIR=\$(dirname "\$(command -v b0_rank.py)")
    BIN_DIR="\$BIN_DIR" python3 - <<'PY'
import csv, hashlib, json, os, sys
sys.path.insert(0, os.environ["BIN_DIR"])
import hackathon_evaluation as evaluation
path = "slot1/slot1_b0_phenotype.csv"
budget = int("${params.max_rows}")
check = {"file": path, "sha256": hashlib.sha256(open(path, "rb").read()).hexdigest(),
         "data_class": "A", "scope": "The shape of the file, not its content.",
         "row_budget": budget}
try:
    submission = evaluation.load_submission(path)
    rows = submission.get("PROBAND01", [])
    values = [row["epcr"] for row in csv.DictReader(open(path))]
    # Control 7 of docs/llm_data_privacy.md. The row count and the list of EPCR values both
    # state how many rows the file of the child holds, so this file carries the verdicts and
    # neither of them.
    check.update(format_check="PASS", proband_ids=sorted(submission),
                 rows_present="PASS" if rows else "FAIL",
                 rows_within_budget="PASS" if len(rows) <= budget else "FAIL",
                 epcr_values_are_distinct=len(set(values)) == len(values),
                 epcr_minimum_above_the_structural_ceiling=min(float(v) for v in values) > 0.95
                 if values else None)
except Exception as error:
    check.update(format_check="FAIL", error=type(error).__name__ + ": " + str(error)[:200])
json.dump(check, open("slot1/format_check.json", "w"), indent=2)
PY
    python3 ${guard} --class-a-publish --rules ${publish_rules} slot1/format_check.json
    """
}

workflow {
    if (params.mode == 'phenotype_patient') {
        ['annotated_patient_vcf', 'annotated_forced_vcf', 'phenotype_document', 'hpo_obo',
         'genes_to_phenotype', 'literature', 'fasta', 'fai', 'outdir'].each { name ->
            if (!params[name]) { error "Set ${name}." }
        }
        guard = file("${projectDir}/../../harness/privacy_guard.py")
        publish_rules = file("${projectDir}/../../literature/class_a_publish_rules.json")
        extracted = EXTRACT_PHENOTYPE(file(params.phenotype_document), file(params.hpo_obo),
                                      guard, publish_rules)
        B0_PATIENT_PHENOTYPE(file(params.annotated_patient_vcf), file(params.annotated_forced_vcf),
                             file(params.literature), file(params.hpo_obo),
                             file(params.genes_to_phenotype), extracted.terms,
                             file(params.fasta), file(params.fai), guard, publish_rules)
        return
    }
    if (!params.mode || !params.vep_cache || !params.literature || !params.outdir) {
        error "Set mode, vep_cache, literature and outdir."
    }
    cache = file(params.vep_cache)
    literature = file(params.literature)
    if (params.mode == 'suite') {
        inputs = channel.of(['hg002_background', file(params.background)], ['spikes', file(params.spikes)])
        annotated = VEP_ANNOTATE(inputs, cache)
        background = annotated.filter { it[0] == 'hg002_background' }.map { it[1] }
        spikes = annotated.filter { it[0] == 'spikes' }.map { it[1] }
        B0_SUITE(background, spikes, file(params.manifest), literature)
    } else if (params.mode == 'annotate') {
        // The annotate mode runs VEP alone and publishes the annotated file. The held-out time
        // split suite of docs/benchmarks/t1_time_split_v1.md uses it, so that the frozen cache
        // annotates the new spikes with the command line of the suite run s4_b0_suite_v1 and the
        // 25.7 GB cache never leaves the region.
        if (!params.spikes) { error "Set spikes." }
        annotated = VEP_ANNOTATE(channel.of([params.label ?: 'spikes', file(params.spikes)]), cache)
        PUBLISH_ANNOTATION(annotated)
    } else if (params.mode == 'patient') {
        inputs = channel.of(['patient_calls', file(params.patient_vcf)], ['patient_forced', file(params.forced_vcf)])
        annotated = VEP_ANNOTATE(inputs, cache)
        patient = annotated.filter { it[0] == 'patient_calls' }.map { it[1] }
        forced = annotated.filter { it[0] == 'patient_forced' }.map { it[1] }
        B0_PATIENT(patient, forced, literature,
                   file("${projectDir}/../../harness/privacy_guard.py"),
                   file("${projectDir}/../../literature/class_a_publish_rules.json"))
    } else {
        error "mode must be suite, annotate or patient"
    }
}
