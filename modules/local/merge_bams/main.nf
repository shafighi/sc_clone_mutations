process MERGE_BAMS {
    tag "${clone_id}"
    label 'process_medium'
    publishDir "${params.outdir}/pseudobulk/bams", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
        tuple val(clone_id), path(bams), path(bais)

    output:
        tuple val(clone_id), path("${clone_id}.merged.bam"), path("${clone_id}.merged.bam.bai"), emit: merged_bam

    script:
    def bam_list = bams instanceof List ? bams.join(' ') : bams
    def n_bams   = bams instanceof List ? bams.size() : 1
    """
    # Add clone-aware read-group tag to each input BAM, then merge
    # If only one BAM, skip merge and just sort+index
    if [ "${n_bams}" -eq 1 ]; then
        samtools sort -@ ${task.cpus} -o ${clone_id}.merged.bam ${bam_list}
    else
        samtools merge -@ ${task.cpus} -f - ${bam_list} \\
        | samtools sort -@ ${task.cpus} -o ${clone_id}.merged.bam
    fi

    samtools index ${clone_id}.merged.bam

    # Verify output
    samtools quickcheck ${clone_id}.merged.bam
    """
}
