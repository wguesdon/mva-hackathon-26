#!/usr/bin/env nextflow
// R1: annotate, filter and rank the patient structural, copy number and mobile element calls.
//
// The input calls are class C, so the process runs on AWS Batch near the data (ADR 0005). Three
// publishDir patterns separate the outputs. The full annotated table follows outdir and stays class C.
// The ranked candidate table, the Markdown summary and the proposed submission rows go to the
// llm_safe prefix as class B2. The gate file follows outdir and is class A.

nextflow.enable.dsl = 2

params.sv_vcf      = null
params.mei_vcf     = null
params.gtf         = null
params.gnomad_sv   = null
params.intervals   = null
params.deny_tokens = ''
params.run_id      = 'r1_sv_mei_v1'
params.outdir      = null
params.llm_safe    = null

process SV_MEI_ANNOTATE {
    tag "${params.run_id}"
    container 'docker.io/library/python:3.12-bookworm'
    cpus 4
    memory '32 GB'
    publishDir params.outdir, mode: 'copy', pattern: 'out/class_c_*'
    publishDir params.outdir, mode: 'copy', pattern: 'out/class_a_*'
    publishDir params.llm_safe, mode: 'copy', pattern: 'out/class_b2_*'

    input:
    path sv_vcf
    path mei_vcf
    path gtf
    path gnomad_sv
    path intervals
    path guard
    path publish_rules

    output:
    path "out/class_c_sv_mei_annotated.tsv"
    path "out/class_b2_sv_mei_ranked.tsv"
    path "out/class_b2_sv_mei_ranked.md"
    path "out/class_b2_submission_rows.csv"
    path "out/class_b2_submission_rows.md"
    path "out/class_a_gates.json"

    script:
    """
    annotate_sv_mei.py --sv-vcf ${sv_vcf} --mei-vcf ${mei_vcf} --gtf ${gtf} \\
        --gnomad-sv ${gnomad_sv} --intervals ${intervals} \\
        --deny-tokens '${params.deny_tokens}' --run-id ${params.run_id} --out-dir out
    python3 ${guard} --class-a-publish --rules ${publish_rules} out/class_a_gates.json
    """
}

workflow {
    ['sv_vcf', 'mei_vcf', 'gtf', 'gnomad_sv', 'intervals', 'outdir', 'llm_safe'].each { name ->
        if (!params[name]) { error "Set ${name}." }
    }
    SV_MEI_ANNOTATE(file(params.sv_vcf), file(params.mei_vcf), file(params.gtf),
                    file(params.gnomad_sv), file(params.intervals),
                    file("${projectDir}/../../harness/privacy_guard.py"),
                    file("${projectDir}/../../literature/class_a_publish_rules.json"))
}
