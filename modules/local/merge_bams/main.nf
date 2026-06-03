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
        path fasta
        path fai

    output:
        tuple val(clone_id), path("${clone_id}.merged.cram"), path("${clone_id}.merged.cram.crai"), emit: merged_bam

    script:
    def bam_list = bams instanceof List ? bams.join(' ') : bams
    def n_bams   = bams instanceof List ? bams.size() : 1
    """
    # Input BAMs are already coordinate-sorted; merge preserves sort order.
    # Output as CRAM for ~75% size reduction vs BAM.
    if [ "${n_bams}" -eq 1 ]; then
        samtools view -@ ${task.cpus} -C -T ${fasta} -o ${clone_id}.merged.cram ${bam_list}
    else
        samtools merge -@ ${task.cpus} -l 6 --output-fmt CRAM --reference ${fasta} -f ${clone_id}.merged.cram ${bam_list}
    fi

    samtools index ${clone_id}.merged.cram

    samtools quickcheck ${clone_id}.merged.cram
    """
}
