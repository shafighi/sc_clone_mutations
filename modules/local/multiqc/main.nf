process MULTIQC {
    tag "multiqc"
    label 'process_single'
    publishDir "${params.outdir}/reports", mode: params.publish_dir_mode

    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'

    input:
        path qc_files
        path multiqc_config

    output:
        path 'multiqc_report.html', emit: report
        path 'multiqc_data',        emit: data

    script:
    def config_arg = multiqc_config ? "--config ${multiqc_config}" : ""
    """
    multiqc \\
        ${config_arg} \\
        --title "sc_clone_mutations QC Report" \\
        --filename multiqc_report.html \\
        .
    """
}
