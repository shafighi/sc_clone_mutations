process MERGE_CALLER_VCFS {
    tag "${clone_id}"
    label 'process_low'
    publishDir "${params.outdir}/variant_analysis/merged_per_clone", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/bcftools:1.23.1--ha83d96e_0'

    input:
        tuple val(clone_id), val(callers), path(vcfs), path(tbis)

    output:
        tuple val(clone_id), path("${clone_id}.merged_callers.vcf.gz"), path("${clone_id}.merged_callers.vcf.gz.tbi"), emit: merged_vcf

    script:
    def vcf_list = vcfs instanceof List ? vcfs.join(' ') : vcfs
    """
    bcftools merge \\
        --force-samples \\
        --merge none \\
        --output-type z \\
        --output ${clone_id}.merged_callers.vcf.gz \\
        ${vcf_list}

    tabix -p vcf ${clone_id}.merged_callers.vcf.gz
    """
}
