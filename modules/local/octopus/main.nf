process OCTOPUS {
    tag "${clone_id}"
    label 'process_high'
    publishDir "${params.outdir}/variant_calling/octopus", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/octopus:0.7.4--h4ac6f70_2'

    input:
        tuple val(clone_id), path(tumor_cram), path(tumor_crai), path(normal_cram), path(normal_crai)
        path fasta
        path fai
        path intervals  // may be []

    output:
        tuple val(clone_id), path("${clone_id}.octopus.vcf.gz"), path("${clone_id}.octopus.vcf.gz.tbi"), emit: vcf
        path "${clone_id}.octopus.stats.tsv", emit: stats

    script:
    def has_normal    = normal_cram && normal_cram.name != 'NO_FILE'
    def normal_arg    = has_normal ? "--normal-sample ${normal_cram.baseName} --reads ${normal_cram}" : ""
    def intervals_arg = (intervals && intervals.name != 'NO_FILE') ? "--regions-file ${intervals}" : ""
    def caller        = has_normal ? "cancer" : "cancer"  // cancer caller works tumor-only too
    """
    octopus \\
        --reference ${fasta} \\
        --reads ${tumor_cram} \\
        ${normal_arg} \\
        --caller ${caller} \\
        ${intervals_arg} \\
        --threads ${task.cpus} \\
        --output ${clone_id}.octopus.vcf.gz \\
        ${params.octopus_extra_args}

    tabix -p vcf ${clone_id}.octopus.vcf.gz

    bcftools stats ${clone_id}.octopus.vcf.gz > ${clone_id}.octopus.stats.tsv
    """
}
