process MERGE_MUTECT_STATS {
    tag "${clone_id}"
    label 'process_single'
    publishDir "${params.outdir}/variant_calling/mutect2/raw", mode: params.publish_dir_mode

    container 'broadinstitute/gatk:4.6.0.0'

    input:
        tuple val(clone_id), path(stats)

    output:
        tuple val(clone_id), path("${clone_id}.mutect2.vcf.gz.stats"), emit: stats

    script:
    def inputs = stats.collect { "-stats ${it}" }.join(' ')
    """
    gatk MergeMutectStats \\
        ${inputs} \\
        -O ${clone_id}.mutect2.vcf.gz.stats
    """
}
