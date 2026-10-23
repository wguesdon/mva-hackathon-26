#!/usr/bin/env nextflow
// S6: forced genotyping at known MVA allele sites and an insertion evidence scan in BUB1B introns (ADR 0003).
// Run it first on GIAB HG002 as a class A control, then on the patient BAM (class C, outputs stay in S3).

nextflow.enable.dsl = 2

params.label     = null
params.bam       = null
params.bai       = null
params.mei_vcf   = null
params.fasta     = null
params.fai       = null
params.targets   = null
params.intervals = null
params.outdir    = null

process FORCED_GENOTYPE {
    tag "${params.label}"
    container 'quay.io/biocontainers/bcftools:1.24--h118bc1c_2'
    cpus 4
    memory '8 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    path bam
    path bai
    path fasta
    path fai
    path targets

    output:
    path "${params.label}_forced.vcf.gz"
    path "${params.label}_forced.vcf.gz.tbi"

    script:
    """
    bcftools sort -Oz -o targets.vcf.gz ${targets}
    bcftools index -t targets.vcf.gz
    bcftools query -f '%CHROM\\t%POS\\t%REF,%ALT\\n' targets.vcf.gz | bgzip -c > alleles.tsv.gz
    tabix -s1 -b2 -e2 alleles.tsv.gz
    # bcftools annotate -a needs an indexed input file, so each step writes and indexes a file
    # (run s6_hg002_control_v1 failed with "Could not retrieve index file for '-'").
    bcftools mpileup --threads ${task.cpus} -f ${fasta} -T targets.vcf.gz -a FORMAT/AD,FORMAT/DP -d 2000 -q 20 -Q 20 -Ou ${bam} \\
        | bcftools call -m -C alleles -T alleles.tsv.gz -Oz -o called.vcf.gz
    bcftools index -t called.vcf.gz
    bcftools annotate -a targets.vcf.gz -c ID -Oz -o annotated.vcf.gz called.vcf.gz
    echo "${params.label}" > sample_name.txt
    bcftools reheader -s sample_name.txt -o ${params.label}_forced.vcf.gz annotated.vcf.gz
    bcftools index -t ${params.label}_forced.vcf.gz
    """
}

process INSERTION_SCAN {
    tag "${params.label}"
    container 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'
    cpus 2
    memory '4 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    path bam
    path bai
    path intervals
    path mei_vcf

    output:
    path "${params.label}_insertion_scan.tsv"
    path "${params.label}_mei_in_intervals.txt"

    script:
    """
    printf 'interval\\tbin_start\\tbin_end\\treads\\tsoftclip_20bp\\tdiscordant_other_chrom\\tsupplementary_sa\\n' > ${params.label}_insertion_scan.tsv
    while read -r chrom start end name; do
        for ((bin = start; bin < end; bin += 50)); do
            region="\${chrom}:\$((bin + 1))-\$((bin + 50))"
            samtools view -F 0x904 -q 20 ${bam} "\${region}" \\
              | awk -v n="\${name}" -v s="\$bin" -v e="\$((bin + 50))" 'BEGIN{r=0;c=0;d=0;a=0}
                  {r++; lead=\$6; tail=\$6; clip=0
                   if (lead ~ /^[0-9]+S/) { sub(/S.*/, "", lead); if (lead + 0 >= 20) clip=1 }
                   if (tail ~ /[0-9]+S\$/) { sub(/S\$/, "", tail); sub(/^.*[MIDNHP=X]/, "", tail); if (tail + 0 >= 20) clip=1 }
                   c += clip
                   if (\$7 != "=" && \$7 != "*") d++
                   if (index(\$0, "SA:Z:") > 0) a++}
                  END{printf "%s\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\n", n, s, e, r, c, d, a}' >> ${params.label}_insertion_scan.tsv
        done
    done < ${intervals}
    : > ${params.label}_mei_in_intervals.txt
    zcat ${mei_vcf} | grep -v '^#' | while IFS=\$'\\t' read -r c p rest; do
        while read -r chrom start end name; do
            if [ "\$c" = "\$chrom" ] && [ "\$p" -gt "\$start" ] && [ "\$p" -le "\$end" ]; then
                printf '%s\\t%s\\t%s\\n' "\$name" "\$c" "\$p" >> ${params.label}_mei_in_intervals.txt
            fi
        done < ${intervals}
    done
    """
}

workflow {
    ['label', 'bam', 'bai', 'mei_vcf', 'fasta', 'fai', 'targets', 'intervals', 'outdir'].each { name ->
        if (!params[name]) { error "Set ${name}." }
    }
    FORCED_GENOTYPE(file(params.bam), file(params.bai), file(params.fasta), file(params.fai), file(params.targets))
    INSERTION_SCAN(file(params.bam), file(params.bai), file(params.intervals), file(params.mei_vcf))
}
