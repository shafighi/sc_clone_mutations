process MARKDUPLICATES {
    tag "${clone_id}"
    label 'process_high'
    errorStrategy 'retry'
    maxRetries 2
    publishDir "${params.outdir}/pseudobulk/bams", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/picard:3.2.0--hdfd78af_0'

    input:
        tuple val(clone_id), path(cram), path(crai)
        path fasta
        path fai

    output:
        tuple val(clone_id), path("${clone_id}.markdup.cram"), path("${clone_id}.markdup.crai"), emit: bam
        path "${clone_id}.markdup_metrics.txt", emit: metrics

    script:
    """
    picard MarkDuplicates \\
        -Xmx${task.memory.toGiga()}g \\
        INPUT=${cram} \\
        OUTPUT=${clone_id}.markdup.cram \\
        METRICS_FILE=${clone_id}.markdup_metrics.txt \\
        REFERENCE_SEQUENCE=${fasta} \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        VALIDATION_STRINGENCY=LENIENT \\
        COMPRESSION_LEVEL=6 \\
        CREATE_INDEX=true
    """
}
