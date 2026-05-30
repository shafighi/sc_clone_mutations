process STRELKA2_SOMATIC {
    tag "${clone_id}"
    label 'process_high'
    publishDir "${params.outdir}/variant_calling/strelka2", mode: params.publish_dir_mode

    // Strelka2 Docker image — use the official Illumina release
    // Note: Strelka2 is free for research use. Source: https://github.com/Illumina/strelka
    container 'quay.io/biocontainers/strelka:2.9.10--hdfd78af_2'

    input:
        tuple val(clone_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
        path fasta
        path fai
        path intervals  // may be []

    output:
        tuple val(clone_id), path("${clone_id}.strelka2.vcf.gz"), path("${clone_id}.strelka2.vcf.gz.tbi"), emit: vcf
        path "${clone_id}.strelka2.stats.tsv", emit: stats

    script:
    def has_normal    = normal_bam && normal_bam.name != 'NO_FILE'
    def intervals_arg = (intervals && intervals.name != 'NO_FILE') ? "--callRegions ${intervals}.gz" : ""
    def normal_arg    = has_normal ? "--normalBam ${normal_bam}" : ""
    def config_script = has_normal ? "configureStrelkaSomaticWorkflow.py" : "configureStrelkaGermlineWorkflow.py"

    // NOTE: In tumor-only mode, Strelka2 uses its germline workflow.
    // True somatic tumor-only calling is not natively supported; this mode
    // calls variants relative to a germline prior. Consider using Mutect2
    // tumor-only mode as primary for somatic calls.
    """
    # Compress and index intervals BED if provided
    if [ -n "${intervals_arg}" ]; then
        bgzip -c ${intervals} > ${intervals}.gz
        tabix -p bed ${intervals}.gz
    fi

    # Configure Strelka2
    ${config_script} \\
        --tumorBam  ${tumor_bam} \\
        ${normal_arg} \\
        --referenceFasta ${fasta} \\
        ${intervals_arg} \\
        ${params.strelka2_extra_args} \\
        --runDir strelka2_${clone_id}

    # Run on available CPUs
    strelka2_${clone_id}/runWorkflow.py -m local -j ${task.cpus}

    # Collect and rename output VCF
    if [ "${has_normal}" = "true" ]; then
        bcftools concat \\
            strelka2_${clone_id}/results/variants/somatic.snvs.vcf.gz \\
            strelka2_${clone_id}/results/variants/somatic.indels.vcf.gz \\
            -a -O z -o ${clone_id}.strelka2.vcf.gz
    else
        cp strelka2_${clone_id}/results/variants/variants.vcf.gz ${clone_id}.strelka2.vcf.gz
    fi

    tabix -p vcf ${clone_id}.strelka2.vcf.gz

    # Minimal stats
    bcftools stats ${clone_id}.strelka2.vcf.gz > ${clone_id}.strelka2.stats.tsv
    """
}
