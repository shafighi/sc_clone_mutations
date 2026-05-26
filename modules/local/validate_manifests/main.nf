process VALIDATE_MANIFESTS {
    tag "input_validation"
    label 'process_single'
    publishDir "${params.outdir}/validation", mode: params.publish_dir_mode

    container 'ghcr.io/shafighi/scclone-python:1.0.0'

    input:
        path bam_manifest
        path cell_metadata
        path medicc2_tree
        path scunique_events
        path normal_manifest  // may be empty file

    output:
        path 'bam_manifest_validated.csv',  emit: bam_manifest
        path 'cell_metadata_validated.csv', emit: cell_metadata
        path 'validation_report.json',      emit: report

    script:
    def normal_arg = normal_manifest ? "--normal_manifest ${normal_manifest}" : ""
    """
    validate_inputs.py \\
        --bam_manifest     ${bam_manifest} \\
        --cell_metadata    ${cell_metadata} \\
        --medicc2_tree     ${medicc2_tree} \\
        --scunique_events  ${scunique_events} \\
        ${normal_arg} \\
        --out_bam_manifest bam_manifest_validated.csv \\
        --out_cell_metadata cell_metadata_validated.csv \\
        --out_report       validation_report.json \\
        --min_mapped_reads ${params.min_mapped_reads} \\
        --skip_bam_check
    """
}
