process FILTER_CELLS {
    tag "filter_cells"
    label 'process_single'
    publishDir "${params.outdir}/pseudobulk/cell_qc", mode: params.publish_dir_mode

    container 'ghcr.io/shafighi/scclone-python:1.0.0'

    input:
        path cell_clone_assignments
        path bam_manifest

    output:
        path 'filtered_manifest.csv', emit: filtered_manifest
        path 'cell_qc_summary.csv',   emit: qc_summary
        path 'clone_reliability.csv', emit: reliability
        path 'clone_reliability.md',  emit: reliability_md

    script:
    """
    python3 ${projectDir}/bin/filter_cells.py \\
        --assignments           ${cell_clone_assignments} \\
        --bam_manifest          ${bam_manifest} \\
        --min_mapped_reads      ${params.min_mapped_reads} \\
        --max_dup_rate          ${params.max_duplication_rate} \\
        --min_confidence        ${params.min_clone_confidence} \\
        --min_cells_for_calling ${params.min_cells_for_calling} \\
        --min_cells_reliable    ${params.min_cells_reliable} \\
        --out_manifest          filtered_manifest.csv \\
        --out_qc_summary        cell_qc_summary.csv \\
        --out_reliability       clone_reliability.csv \\
        --out_reliability_md    clone_reliability.md
    """
}
