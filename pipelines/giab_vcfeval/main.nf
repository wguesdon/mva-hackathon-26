#!/usr/bin/env nextflow
// Compare a query VCF with a GIAB truth set using RTG Tools vcfeval.
// It replaces the built-in rtg vcfeval of nf-core/raredisease 3.1.2, which fails (ADR 0002 notes).

nextflow.enable.dsl = 2

params.query_vcf  = null
params.truth_vcf  = null
params.truth_bed  = null
params.fasta      = null
params.sample     = null
params.outdir     = null

process RTG_FORMAT {
    container 'quay.io/biocontainers/rtg-tools:3.13--hdfd78af_0'
    cpus 4
    memory '16 GB'

    input:
    path fasta

    output:
    path 'reference.sdf'

    script:
    """
    rtg format -o reference.sdf ${fasta}
    """
}

process RTG_VCFEVAL {
    container 'quay.io/biocontainers/rtg-tools:3.13--hdfd78af_0'
    cpus 16
    memory '32 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    path sdf
    path query_vcf
    path truth_vcf
    path truth_bed

    output:
    path 'vcfeval/*'

    script:
    """
    rtg index -f vcf ${query_vcf} || true
    rtg index -f vcf ${truth_vcf} || true
    rtg vcfeval --baseline ${truth_vcf} --calls ${query_vcf} --bed-regions ${truth_bed} \\
        --template ${sdf} --output vcfeval --threads ${task.cpus} \\
        --sample ${params.sample},${params.sample} --ref-overlap --decompose
    """
}

workflow {
    if (!params.query_vcf || !params.truth_vcf || !params.truth_bed || !params.fasta || !params.sample || !params.outdir) {
        error "Set query_vcf, truth_vcf, truth_bed, fasta, sample and outdir."
    }
    sdf = RTG_FORMAT(file(params.fasta))
    RTG_VCFEVAL(sdf, file(params.query_vcf), file(params.truth_vcf), file(params.truth_bed))
}
