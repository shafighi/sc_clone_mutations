process MARKDUPLICATES {
    tag "${clone_id}"
    label 'process_medium'
    publishDir "${params.outdir}/pseudobulk/bams", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/picard:3.2.0--hdfd78af_0'

    input:
        tuple val(clone_id), path(bam), path(bai)

    output:
        tuple val(clone_id), path("${clone_id}.markdup.bam"), path("${clone_id}.markdup.bam.bai"), emit: bam
        path "${clone_id}.markdup_metrics.txt", emit: metrics

    script:
    // NOTE: For pseudobulk somatic calling, consider whether to REMOVE duplicates
    // (better sensitivity) or keep them. Default here: mark but do not remove,
    // which lets callers decide. Set REMOVE_DUPLICATES=true to remove them.
    """
    picard MarkDuplicates \\
        -Xmx${task.memory.toGiga()}g \\
        INPUT=${bam} \\
        OUTPUT=${clone_id}.markdup.bam \\
        METRICS_FILE=${clone_id}.markdup_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        VALIDATION_STRINGENCY=LENIENT \\
        CREATE_INDEX=true

    mv ${clone_id}.markdup.bai ${clone_id}.markdup.bam.bai
    """
}
