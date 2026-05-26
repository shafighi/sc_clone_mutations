process CUSTOM_REPORT {
    tag "custom_report"
    label 'process_single'
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode

    container 'ghcr.io/shafighi/scclone-python:1.0.0'

    input:
        path validation_report
        path clone_summary
        path consensus_table
        path cross_clone_matrix

    output:
        path 'sc_clone_mutations_report.html', emit: report
        path 'sc_clone_mutations_report.md',   emit: report_md

    script:
    """
    generate_report.py \\
        --validation_report  ${validation_report} \\
        --clone_summary      ${clone_summary} \\
        --consensus_table    ${consensus_table} \\
        --cross_clone_matrix ${cross_clone_matrix} \\
        --pipeline_version   ${workflow.manifest.version} \\
        --out_html           sc_clone_mutations_report.html \\
        --out_md             sc_clone_mutations_report.md
    """
}
