process SAMTOOLS_FLAGSTAT {
    tag "${clone_id}"
    label 'process_single'
    publishDir "${params.outdir}/pseudobulk/qc/flagstat", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/samtools:1.23.1--ha83d96e_0'

    input:
        tuple val(clone_id), path(bam), path(bai)

    output:
        path "${clone_id}.flagstat.txt", emit: flagstat

    script:
    """
    samtools flagstat -@ ${task.cpus} ${bam} > ${clone_id}.flagstat.txt
    """
}
