process FILTER_MUTECT2 {
    tag "${clone_id}"
    label 'process_medium'
    publishDir "${params.outdir}/variant_calling/mutect2/filtered", mode: params.publish_dir_mode

    container 'broadinstitute/gatk:4.6.0.0'

    input:
        tuple val(clone_id), path(vcf), path(tbi)
        path stats
        path fasta
        path fai
        path dict

    output:
        tuple val(clone_id), path("${clone_id}.mutect2.filtered.vcf.gz"), path("${clone_id}.mutect2.filtered.vcf.gz.tbi"), emit: filtered_vcf
        path "${clone_id}.mutect2.filtering.stats", emit: filter_stats

    script:
    """
    gatk FilterMutectCalls \\
        -R ${fasta} \\
        -V ${vcf} \\
        --stats ${stats} \\
        --filtering-stats ${clone_id}.mutect2.filtering.stats \\
        -O ${clone_id}.mutect2.filtered.vcf.gz
    """
}
