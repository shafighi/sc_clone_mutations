process MOSDEPTH {
    tag "${clone_id}"
    label 'process_low'
    publishDir "${params.outdir}/pseudobulk/qc/mosdepth", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0'

    input:
        tuple val(clone_id), path(bam), path(bai)
        path fasta
        path fai

    output:
        path "${clone_id}.mosdepth.summary.txt",  emit: summary
        path "${clone_id}.mosdepth.global.dist.txt"
        path "${clone_id}.mosdepth.region.dist.txt", optional: true

    script:
    def intervals_arg = params.intervals ? "--by ${params.intervals}" : "--by 500"
    """
    mosdepth \\
        --threads ${task.cpus} \\
        --no-abbrev \\
        ${intervals_arg} \\
        ${clone_id} \\
        ${bam}
    """
}
