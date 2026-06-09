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
        path tree_png
        path clone_reliability
        path signatures,       stageAs: 'sig_exposures/*'
        path signature_plots,  stageAs: 'sig_plots/*'

    output:
        path 'sc_clone_mutations_report.html', emit: report
        path 'sc_clone_mutations_report.md',   emit: report_md

    script:
    def sig_arg  = signatures      ? "--signatures ${signatures.join(' ')}"           : ""
    def plot_arg = signature_plots ? "--signature_plots ${signature_plots.join(' ')}" : ""
    """
    python3 ${projectDir}/bin/generate_report.py \\
        --validation_report  ${validation_report} \\
        --clone_summary      ${clone_summary} \\
        --consensus_table    ${consensus_table} \\
        --cross_clone_matrix ${cross_clone_matrix} \\
        --clone_reliability  ${clone_reliability} \\
        --clone_tree_png     ${tree_png} \\
        ${sig_arg} \\
        ${plot_arg} \\
        --pipeline_version   ${workflow.manifest.version} \\
        --out_html           sc_clone_mutations_report.html \\
        --out_md             sc_clone_mutations_report.md
    """
}
