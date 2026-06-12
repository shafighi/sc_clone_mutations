process SPLIT_INTERVALS {
    tag "scatter=${params.mutect2_scatter_count}"
    label 'process_single'

    container 'broadinstitute/gatk:4.6.0.0'

    input:
        path fasta
        path fai
        path dict
        path intervals      // may be [] (whole genome)

    output:
        path 'intervals_scatter/*.interval_list', emit: intervals

    script:
    def l_arg = (intervals && intervals.name != 'NO_FILE') ? "-L ${intervals}" : ""
    """
    mkdir -p intervals_scatter
    gatk SplitIntervals \\
        -R ${fasta} \\
        ${l_arg} \\
        --scatter-count ${params.mutect2_scatter_count} \\
        --subdivision-mode INTERVAL_SUBDIVISION \\
        -O intervals_scatter
    """
}
