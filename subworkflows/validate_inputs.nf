/*
================================================================================
  subworkflows/validate_inputs.nf
  Validates all input manifests and reference files before the main pipeline
================================================================================
*/

include { VALIDATE_MANIFESTS } from '../modules/local/validate_manifests/main'

workflow VALIDATE_INPUTS {
    take:
        ch_bam_manifest
        ch_cell_metadata
        ch_medicc2_tree
        ch_scunique_events
        ch_normal_manifest
        ch_fasta
        ch_fai

    main:
        VALIDATE_MANIFESTS(
            ch_bam_manifest,
            ch_cell_metadata,
            ch_medicc2_tree,
            ch_scunique_events,
            ch_normal_manifest
        )

    emit:
        bam_manifest_validated  = VALIDATE_MANIFESTS.out.bam_manifest
        cell_metadata_validated = VALIDATE_MANIFESTS.out.cell_metadata
        validation_report       = VALIDATE_MANIFESTS.out.report
}
