#!/usr/bin/env nextflow
// Summarize the S6 patient outputs on AWS Batch, near the data (ADR 0005).
//
// The class C genotype table stays under the run results prefix. The class B2 summary goes to a
// separate llm_safe prefix, so that only the summary can reach a Bedrock model after owner review.
// No model reads the class C table.

nextflow.enable.dsl = 2

params.forced    = null
params.scan      = null
params.mei       = null
params.targets   = null
params.control   = null
params.outdir    = null
params.llm_safe  = null

process S6_SUMMARY {
    container 'docker.io/library/python:3.12-bookworm'
    cpus 2
    memory '4 GB'
    publishDir params.outdir, mode: 'copy', pattern: 'out/class_c_*'
    publishDir params.llm_safe, mode: 'copy', pattern: 'out/class_b_*'

    input:
    path forced
    path scan
    path mei
    path targets
    path control

    output:
    path "out/class_c_genotypes.tsv"
    path "out/class_b_summary.md"

    script:
    """
    summarize_s6_patient.py --forced ${forced} --scan ${scan} --mei ${mei} \\
        --targets ${targets} --control ${control} --out-dir out
    """
}

workflow {
    ['forced', 'scan', 'mei', 'targets', 'control', 'outdir', 'llm_safe'].each { name ->
        if (!params[name]) { error "Set ${name}." }
    }
    S6_SUMMARY(file(params.forced), file(params.scan), file(params.mei),
               file(params.targets), file(params.control))
}
