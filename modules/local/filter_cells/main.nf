process FILTER_CELLS {
    tag "filter_cells"
    label 'process_single'
    publishDir "${params.outdir}/pseudobulk/cell_qc", mode: params.publish_dir_mode

    container 'ghcr.io/TODO/scclone-python:1.0.0'

    input:
        path cell_clone_assignments
        path bam_manifest

    output:
        path 'filtered_manifest.csv', emit: filtered_manifest
        path 'cell_qc_summary.csv',   emit: qc_summary

    script:
    """
    filter_cells.py \\
        --assignments       ${cell_clone_assignments} \\
        --bam_manifest      ${bam_manifest} \\
        --min_mapped_reads  ${params.min_mapped_reads} \\
        --max_dup_rate      ${params.max_duplication_rate} \\
        --min_confidence    ${params.min_clone_confidence} \\
        --out_manifest      filtered_manifest.csv \\
        --out_qc_summary    cell_qc_summary.csv
    """
}
