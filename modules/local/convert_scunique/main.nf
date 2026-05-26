process CONVERT_SCUNIQUE {
    tag "convert_rds"
    label 'process_low'
    container 'rocker/tidyverse:4.3.2'
    publishDir "${params.outdir}/converted_inputs", mode: params.publish_dir_mode

    input:
    path scunique_dir

    output:
    path "medicc2_tree.new",     emit: medicc2_tree
    path "scunique_events.tsv",  emit: scunique_events
    path "cell_metadata.csv",    emit: cell_metadata
    path "bam_manifest.csv",     emit: bam_manifest, optional: true

    script:
    def sample_arg = params.sample_id ?: ''
    def bam_arg    = params.bam_dir ?: ''
    """
    Rscript ${projectDir}/bin/convert_scunique_rds.R \\
        ${scunique_dir} \\
        . \\
        "${sample_arg}" \\
        "${bam_arg}"
    """
}
