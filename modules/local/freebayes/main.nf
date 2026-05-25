process FREEBAYES {
    tag "${clone_id}"
    label 'process_high'
    publishDir "${params.outdir}/variant_calling/freebayes", mode: params.publish_dir_mode

    // FreeBayes: MIT license, fully open-source
    container 'quay.io/biocontainers/freebayes:1.3.7--hbfe0e7f_2'

    input:
        tuple val(clone_id), path(bam), path(bai)
        path fasta
        path fai
        path intervals  // may be []

    output:
        tuple val(clone_id), path("${clone_id}.freebayes.vcf.gz"), emit: vcf
        path "${clone_id}.freebayes.stats.tsv", emit: stats

    script:
    def regions_arg = (intervals && intervals.name != 'NO_FILE') ? "--targets ${intervals}" : ""
    """
    freebayes \\
        --fasta-reference ${fasta} \\
        --bam ${bam} \\
        ${regions_arg} \\
        ${params.freebayes_extra_args} \\
    | bcftools view -f 'PASS,.' \\
    | bgzip -c > ${clone_id}.freebayes.vcf.gz

    tabix -p vcf ${clone_id}.freebayes.vcf.gz

    bcftools stats ${clone_id}.freebayes.vcf.gz > ${clone_id}.freebayes.stats.tsv
    """
}
