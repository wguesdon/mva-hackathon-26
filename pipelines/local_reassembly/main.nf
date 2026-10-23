#!/usr/bin/env nextflow
// R12: targeted local reassembly across the seven MVA genes, and a sequence-resolved structural
// variant call from the contigs.
//
// R13 adds one step. When `merge_sv_vcf` names the production structural call set of the same
// sample, the workflow merges that call set with the reassembly call set through SVDB 2.8.4, and it
// ranks the merged set together with the mobile element call set named by `mei_vcf`. When
// `merge_sv_vcf` is null the workflow is the R12 workflow and it ranks the reassembly call set
// alone. The amended gate of R13 is in docs/benchmarks/r13_merged_call_set_v1.md.
//
// The alignment is class C, so every process runs on AWS Batch near the data (ADR 0005). The
// workflow reads one target interval at a time, assembles the reads of that interval on their own,
// aligns the contigs to the whole GRCh38 reference and reads the insertions and the deletions from
// the contig alignment. The last process runs the R1 annotator on the new call set, so the ranking
// rule, the gnomAD frequency limit and the known interval labels are the ones of queue item R1.
//
// Three publishDir patterns separate the outputs. The call set and the annotated table follow
// outdir and stay class C. The ranked candidate table, the Markdown summary and the proposed
// submission rows go to the llm_safe prefix as class B2. The gate files follow outdir and are
// class A.

nextflow.enable.dsl = 2

params.bam         = null
params.bai         = null
params.targets     = null
params.reference   = null
params.reference_fai = null
params.gtf         = null
params.gnomad_sv   = null
params.intervals   = null
params.annotator   = null
params.sample      = 'SAMPLE'
params.deny_tokens = ''
params.min_size    = 50
params.run_id      = 'r12_local_reassembly_v1'
params.outdir      = null
params.llm_safe    = null

// Pinned images. The local HG002 validation in scripts/r12_run_local.py uses the same four.
params.samtools_image = 'quay.io/biocontainers/samtools:1.22--h96c455f_0'
params.spades_image   = 'quay.io/biocontainers/spades:3.15.5--h95f258a_1'
params.minimap2_image = 'quay.io/biocontainers/minimap2:2.28--he4a0461_0'
params.python_image   = 'docker.io/library/python:3.12-bookworm'
// The SVDB container of the patient run patient_calling_v1. Its header line reads SVDB_version=2.8.4.
params.svdb_image     = 'community.wave.seqera.io/library/bcftools_svdb:12db401acbacc624'

params.spades_kmers    = '21,33,55,77'
params.minimap2_preset = 'asm5'

// R13. When merge_sv_vcf is set, the workflow merges the production structural call set with the
// reassembly call set and it ranks the merged set together with the mobile element call set. When
// merge_sv_vcf is null, the workflow is the R12 workflow and it ranks the reassembly call set alone.
params.merge_sv_vcf        = null
params.mei_vcf             = null
params.merge_overlap       = '0.5'
params.merge_bnd_distance  = '150'
params.merge_ins_distance  = '50'

process EXTRACT_READS {
    tag "${params.run_id}"
    container params.samtools_image
    cpus 8
    memory '16 GB'

    input:
    path bam
    path bai
    path targets

    output:
    path "reads"

    script:
    // One task reads every interval, so the alignment is staged from S3 once and not once per gene.
    """
    mkdir -p reads
    while IFS=\$'\\t' read -r contig start end name rest; do
        [ -n "\${contig}" ] || continue
        region="\${contig}:\$((start + 1))-\${end}"
        samtools view -u -@ ${task.cpus} ${bam} "\${region}" \\
            | samtools collate -u -O -@ ${task.cpus} - \\
            | samtools fastq -N -@ ${task.cpus} \\
                -1 reads/\${name}_R1.fastq.gz -2 reads/\${name}_R2.fastq.gz \\
                -s reads/\${name}_S.fastq.gz -0 /dev/null -
    done < ${targets}
    """
}

process ASSEMBLE {
    tag "${name}"
    container params.spades_image
    cpus 8
    memory '32 GB'

    input:
    tuple val(name), path(reads)

    output:
    tuple val(name), path("${name}_contigs.fasta")

    script:
    """
    spades.py --pe1-1 ${reads}/${name}_R1.fastq.gz --pe1-2 ${reads}/${name}_R2.fastq.gz \\
        --pe1-s ${reads}/${name}_S.fastq.gz \\
        -k ${params.spades_kmers} -t ${task.cpus} -m ${task.memory.toGiga()} -o assembly_${name}
    cp assembly_${name}/contigs.fasta ${name}_contigs.fasta
    """
}

process ALIGN_CONTIGS {
    tag "${name}"
    container params.minimap2_image
    cpus 8
    memory '32 GB'

    input:
    tuple val(name), path(contigs)
    path reference

    output:
    tuple val(name), path("${name}.paf")

    script:
    """
    minimap2 -c --cs -x ${params.minimap2_preset} -t ${task.cpus} ${reference} ${contigs} > ${name}.paf
    """
}

process CALL_FROM_CONTIGS {
    tag "${name}"
    container params.python_image
    cpus 2
    memory '8 GB'

    input:
    tuple val(name), path(paf)
    path reference
    path reference_fai
    path targets

    output:
    path "${name}.vcf"

    script:
    """
    r12_reassembly.py call --paf ${paf} --reference ${reference} --bed ${targets} \\
        --interval ${name} --out-vcf ${name}.vcf --sample ${params.sample} \\
        --min-size ${params.min_size}
    """
}

process MERGE_CALLS {
    tag "${params.run_id}"
    container params.python_image
    cpus 2
    memory '8 GB'
    publishDir params.outdir, mode: 'copy', pattern: 'out/class_c_*'
    publishDir params.outdir, mode: 'copy', pattern: 'out/class_a_*'

    input:
    path vcf_files
    path reference
    path reference_fai
    path targets
    path guard
    path publish_rules

    output:
    path "out/class_c_r12_calls.vcf", emit: calls
    path "out/class_c_r12_empty_mei.vcf", emit: empty_mei
    path "out/class_a_reassembly_gates.json"

    script:
    """
    mkdir -p interval_calls out
    for item in ${vcf_files}; do cp \$item interval_calls/; done
    r12_reassembly.py merge --vcf-dir interval_calls --reference ${reference} --bed ${targets} \\
        --out-vcf out/class_c_r12_calls.vcf --out-gates out/class_a_reassembly_gates.json \\
        --sample ${params.sample} --run-id ${params.run_id}
    printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${params.sample}\\n' \\
        > out/class_c_r12_empty_mei.vcf
    python3 ${guard} --class-a-publish --rules ${publish_rules} out/class_a_reassembly_gates.json
    """
}

process MERGE_ARMS {
    tag "${params.run_id}"
    container params.svdb_image
    cpus 2
    memory '16 GB'
    publishDir params.outdir, mode: 'copy', pattern: 'out/class_c_*'
    publishDir params.outdir, mode: 'copy', pattern: 'out/class_a_*'

    input:
    path reassembly_vcf
    path pipeline_vcf
    path guard
    path publish_rules

    output:
    path "out/class_c_merged_sv.vcf", emit: calls
    path "out/class_c_merge_counts.json"
    path "out/class_a_merge_gates.json"

    script:
    """
    mkdir -p out
    r13_merge_arms.py --pipeline-vcf ${pipeline_vcf} --reassembly-vcf ${reassembly_vcf} \\
        --out-vcf out/class_c_merged_sv.vcf --out-gates out/class_a_merge_gates.json \\
        --out-counts out/class_c_merge_counts.json \\
        --work-dir merge_work --overlap ${params.merge_overlap} \\
        --bnd-distance ${params.merge_bnd_distance} --ins-distance ${params.merge_ins_distance} \\
        --run-id ${params.run_id}
    python3 ${guard} --class-a-publish --rules ${publish_rules} out/class_a_merge_gates.json
    """
}

process ANNOTATE_RANK {
    tag "${params.run_id}"
    container params.python_image
    cpus 4
    memory '32 GB'
    publishDir params.outdir, mode: 'copy', pattern: 'out/class_c_*'
    publishDir params.outdir, mode: 'copy', pattern: 'out/class_a_*'
    publishDir params.llm_safe, mode: 'copy', pattern: 'out/class_b2_*'

    input:
    path calls
    path empty_mei
    path annotator
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
    python3 ${annotator} --sv-vcf ${calls} --mei-vcf ${empty_mei} --gtf ${gtf} \\
        --gnomad-sv ${gnomad_sv} --intervals ${intervals} \\
        --deny-tokens '${params.deny_tokens}' --run-id ${params.run_id} --out-dir out
    python3 ${guard} --class-a-publish --rules ${publish_rules} out/class_a_gates.json
    """
}

workflow {
    ['bam', 'bai', 'targets', 'reference', 'reference_fai', 'gtf', 'gnomad_sv', 'intervals',
     'annotator', 'outdir', 'llm_safe'].each { name ->
        if (!params[name]) { error "Set ${name}." }
    }

    reference = file(params.reference)
    reference_fai = file(params.reference_fai)
    targets = file(params.targets)

    names = Channel
        .fromPath(params.targets)
        .splitCsv(sep: '\t')
        .map { row -> row[3] }

    reads = EXTRACT_READS(file(params.bam), file(params.bai), targets)
    contigs = ASSEMBLE(names.combine(reads))
    paf = ALIGN_CONTIGS(contigs, reference)
    calls = CALL_FROM_CONTIGS(paf, reference, reference_fai, targets)
    guard = file("${projectDir}/../../harness/privacy_guard.py")
    publish_rules = file("${projectDir}/../../literature/class_a_publish_rules.json")

    merged = MERGE_CALLS(calls.collect(), reference, reference_fai, targets, guard, publish_rules)

    if (params.merge_sv_vcf) {
        if (!params.mei_vcf) { error "Set mei_vcf beside merge_sv_vcf." }
        arms = MERGE_ARMS(merged.calls, file(params.merge_sv_vcf), guard, publish_rules)
        ANNOTATE_RANK(arms.calls, file(params.mei_vcf), file(params.annotator), file(params.gtf),
                      file(params.gnomad_sv), file(params.intervals), guard, publish_rules)
    }
    else {
        ANNOTATE_RANK(merged.calls, merged.empty_mei, file(params.annotator), file(params.gtf),
                      file(params.gnomad_sv), file(params.intervals), guard, publish_rules)
    }
}
