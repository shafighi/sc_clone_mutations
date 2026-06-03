process MERGE_BAMS {
    tag "${clone_id}"
    label 'process_medium'
    label 'process_long'
    errorStrategy 'retry'
    maxRetries 2
    publishDir "${params.outdir}/pseudobulk/bams", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/samtools:1.23.1--ha83d96e_0'

    input:
        tuple val(clone_id), path(bams), path(bais)

    output:
        tuple val(clone_id), path("${clone_id}.merged.bam"), path("${clone_id}.merged.bam.bai"), emit: merged_bam

    script:
    def bam_list = bams instanceof List ? bams.join(' ') : bams
    def n_bams   = bams instanceof List ? bams.size() : 1
    """
    # Input BAMs are already coordinate-sorted; merge preserves sort order
    if [ "${n_bams}" -eq 1 ]; then
        cp ${bam_list} ${clone_id}.merged.bam
    else
        samtools merge -@ ${task.cpus} -f ${clone_id}.merged.bam ${bam_list}
    fi

    samtools index ${clone_id}.merged.bam

    # Verify output
    samtools quickcheck ${clone_id}.merged.bam
    """
}
