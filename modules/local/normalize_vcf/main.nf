process NORMALIZE_VCF {
    tag "${clone_id} | ${caller}"
    label 'process_low'
    publishDir "${params.outdir}/variant_calling/${caller}/normalized", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/bcftools:1.23.1--ha83d96e_0'

    input:
        tuple val(clone_id), val(caller), path(vcf), path(tbi)
        path fasta

    output:
        tuple val(clone_id), val(caller), path("${clone_id}.${caller}.norm.vcf.gz"), path("${clone_id}.${caller}.norm.vcf.gz.tbi"), emit: normalized_vcf

    script:
    def tbi_arg = (tbi && tbi.name != 'NO_FILE') ? "" : ""  // tabix auto-detected
    """
    bcftools norm \\
        --fasta-ref ${fasta} \\
        --multiallelics -both \\
        --output-type z \\
        --output ${clone_id}.${caller}.norm.vcf.gz \\
        ${vcf}

    tabix -p vcf ${clone_id}.${caller}.norm.vcf.gz
    """
}
