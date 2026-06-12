process MERGE_VCFS {
    tag "${clone_id}"
    label 'process_single'
    publishDir "${params.outdir}/variant_calling/mutect2/raw", mode: params.publish_dir_mode

    container 'broadinstitute/gatk:4.6.0.0'

    input:
        tuple val(clone_id), path(vcfs), path(tbis)

    output:
        tuple val(clone_id), path("${clone_id}.mutect2.vcf.gz"), path("${clone_id}.mutect2.vcf.gz.tbi"), emit: vcf

    script:
    def inputs = vcfs.collect { "-I ${it}" }.join(' ')
    """
    gatk MergeVcfs \\
        ${inputs} \\
        --CREATE_INDEX true \\
        -O ${clone_id}.mutect2.vcf.gz
    """
}
