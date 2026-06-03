process SAMTOOLS_FLAGSTAT {
    tag "${clone_id}"
    label 'process_single'
    publishDir "${params.outdir}/pseudobulk/qc/flagstat", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/samtools:1.23.1--ha83d96e_0'

    input:
        tuple val(clone_id), path(cram), path(crai)
        path fasta

    output:
        path "${clone_id}.flagstat.txt", emit: flagstat

    script:
    """
    export REF_PATH=${fasta}
    samtools flagstat -@ ${task.cpus} ${cram} > ${clone_id}.flagstat.txt
    """
}
