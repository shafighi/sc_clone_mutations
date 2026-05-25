/*
================================================================================
  subworkflows/reporting.nf
  Collects all QC outputs and generates MultiQC + custom HTML report
================================================================================
*/

include { MULTIQC        } from '../modules/local/multiqc/main'
include { CUSTOM_REPORT  } from '../modules/local/custom_report/main'

workflow REPORTING {
    take:
        ch_validation_report
        ch_clone_summary
        ch_pseudobulk_qc
        ch_caller_stats
        ch_consensus_table
        ch_cross_clone_matrix

    main:
        // Gather all QC inputs for MultiQC
        ch_multiqc_inputs = ch_pseudobulk_qc
            .mix(ch_caller_stats)
            .collect()

        MULTIQC(
            ch_multiqc_inputs,
            Channel.fromPath(params.multiqc_config, checkIfExists: true)
        )

        // Custom HTML / Markdown report
        CUSTOM_REPORT(
            ch_validation_report,
            ch_clone_summary,
            ch_consensus_table,
            ch_cross_clone_matrix
        )

    emit:
        multiqc_report = MULTIQC.out.report
        custom_report  = CUSTOM_REPORT.out.report
}
