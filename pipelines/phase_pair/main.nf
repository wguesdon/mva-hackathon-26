#!/usr/bin/env nextflow
// Read backed phasing of the compound heterozygous pairs of a Track 1 submission.
//
// Defect R2-4 of the second jury round states that the B0 ranker reports a compound heterozygous
// pair with no phase. The disease of the case is recessive, so the two alleles must lie on opposite
// haplotypes. This workflow measures whether the reads support that state.
//
// The alignment and the call set are class C, so every process runs on AWS Batch near the data
// (ADR 0005). The method and the decision gates are in docs/benchmarks/t1_phase_v1.md. The method
// is validated on GIAB HG002 first, receipt tracking/receipts/t1_phase_control_v1.json.
//
// Two publishDir patterns separate the outputs. The phase result of the sample goes to the llm_safe
// prefix as class B2. The gate file follows outdir and is class A. The class A file holds process
// gates, tool versions and checksums, and no result of the sample.

nextflow.enable.dsl = 2

params.bam           = null
params.bai           = null
params.calls_vcf     = null
params.forced_vcf    = null
params.submission    = null
params.targets       = null
params.reference     = null
params.reference_fai = null
params.run_id        = 'phase_pair_patient_v1'
params.outdir        = null
params.llm_safe      = null

// Pinned images. The HG002 validation in scripts/run_track1_phase_control.sh uses the same three.
params.samtools_image = 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'
params.bcftools_image = 'quay.io/biocontainers/bcftools:1.24--h118bc1c_2'
params.whatshap_image = 'quay.io/biocontainers/whatshap:2.8--py312hf731ba3_0'
params.python_image   = 'docker.io/library/python:3.12-bookworm'

process SLICE_ALIGNMENT {
    tag "${params.run_id}"
    container params.samtools_image
    cpus 4
    memory '16 GB'

    input:
    path bam
    path bai
    path targets

    output:
    tuple path("targets.bam"), path("targets.bam.bai")

    script:
    """
    samtools view -b -M -L ${targets} -@ ${task.cpus} -o targets.bam ${bam}
    samtools index -@ ${task.cpus} targets.bam
    """
}

process HETEROZYGOUS_CALLS {
    tag "${params.run_id}"
    container params.bcftools_image
    cpus 2
    memory '8 GB'

    input:
    path calls
    path forced
    path targets

    output:
    tuple path("het.vcf.gz"), path("het.vcf.gz.tbi")

    script:
    // The DeepVariant output of the calling run writes a missing FILTER value, so the filter list
    // holds both PASS and the missing value. The B0 ranker applies the same rule.
    """
    for input in ${calls} ${forced}; do
        bcftools index -t -f "\${input}" || bcftools index -f "\${input}"
    done
    bcftools concat -a -d exact -Oz -o union.vcf.gz ${calls} ${forced}
    bcftools index -t -f union.vcf.gz
    bcftools view -R ${targets} -f PASS,. -i 'GT="het"' -Oz -o het.vcf.gz union.vcf.gz
    bcftools index -t -f het.vcf.gz
    """
}

process PHASE {
    tag "${params.run_id}"
    container params.whatshap_image
    cpus 4
    memory '16 GB'

    input:
    tuple path(bam), path(bai)
    tuple path(vcf), path(tbi)
    path reference
    path reference_fai

    output:
    path "phased.vcf"

    script:
    """
    whatshap phase --reference ${reference} --ignore-read-groups -o phased.vcf ${vcf} ${bam}
    """
}

process READ_PAIRS {
    tag "${params.run_id}"
    container params.python_image
    cpus 2
    memory '8 GB'
    publishDir params.outdir, mode: 'copy', pattern: 'out/class_a_*'
    publishDir params.llm_safe, mode: 'copy', pattern: 'out/class_b2_*'

    input:
    path phased
    path submission
    path reader
    path guard
    path publish_rules

    output:
    path "out/class_b2_phase_rows.tsv"
    path "out/class_b2_phase_rows.md"
    path "out/class_a_phase_gates.json"

    script:
    // The class A publish gate of control 7 runs here, before publishDir copies the gate file to
    // the class A prefix and before any agent reads it. A refusal stops the task, so the run fails.
    """
    python3 ${reader} --phased ${phased} --submission ${submission} --out-dir out \\
        --run-id ${params.run_id} \\
        --whatshap-image ${params.whatshap_image} \\
        --samtools-image ${params.samtools_image} \\
        --bcftools-image ${params.bcftools_image}
    python3 ${guard} --class-a-publish --rules ${publish_rules} out/class_a_phase_gates.json
    """
}

workflow {
    ['bam', 'bai', 'calls_vcf', 'forced_vcf', 'submission', 'targets', 'reference',
     'reference_fai', 'outdir', 'llm_safe'].each { name ->
        if (!params[name]) { error "Set ${name}." }
    }

    targets = file(params.targets)
    sliced = SLICE_ALIGNMENT(file(params.bam), file(params.bai), targets)
    het = HETEROZYGOUS_CALLS(file(params.calls_vcf), file(params.forced_vcf), targets)
    phased = PHASE(sliced, het, file(params.reference), file(params.reference_fai))
    READ_PAIRS(phased, file(params.submission),
               file("${projectDir}/bin/phase_pairs.py"),
               file("${projectDir}/../../harness/privacy_guard.py"),
               file("${projectDir}/../../literature/class_a_publish_rules.json"))
}
