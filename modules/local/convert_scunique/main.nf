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

    script:
    def sample_id = params.sample_id ?: ''
    """
    Rscript ${projectDir}/bin/convert_scunique_rds.R \\
        ${scunique_dir} \\
        . \\
        ${sample_id}
    """
}
