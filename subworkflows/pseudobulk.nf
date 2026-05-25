/*
================================================================================
  subworkflows/pseudobulk.nf

  Builds clone-level pseudobulk BAMs by:
    1. Filtering cells by QC metrics
    2. Merging per-cell BAMs for each clone
    3. Sorting, indexing, marking duplicates
    4. Computing coverage QC (mosdepth + samtools flagstat)

  Emits one (clone_id, bam, bai) tuple per clone.
================================================================================
*/

include { FILTER_CELLS       } from '../modules/local/filter_cells/main'
include { MERGE_BAMS         } from '../modules/local/merge_bams/main'
include { MARKDUPLICATES     } from '../modules/local/markduplicates/main'
include { MOSDEPTH           } from '../modules/local/mosdepth/main'
include { SAMTOOLS_FLAGSTAT  } from '../modules/local/samtools_flagstat/main'

workflow PSEUDOBULK {
    take:
        ch_cell_clone_assignments  // CSV: cell_id,clone_id,confidence,...
        ch_bam_manifest            // CSV: cell_id,bam_path,bai_path,sample_id,patient_id
        ch_fasta                   // reference FASTA (needed for coverage)
        ch_fai

    main:
        // Step 1: Filter cells by QC thresholds
        FILTER_CELLS(
            ch_cell_clone_assignments,
            ch_bam_manifest
        )
        // Output: filtered_manifest.csv (same schema as bam_manifest, subset of rows)

        // Step 2: Group cells by clone → one channel item per clone
        // Channel: [ clone_id, [ list of bam paths ], [ list of bai paths ] ]
        ch_clone_bam_groups = FILTER_CELLS.out.filtered_manifest
            .splitCsv(header: true)
            .map { row ->
                def bam = file(row.bam_path, checkIfExists: true)
                def bai = row.bai_path ? file(row.bai_path, checkIfExists: true)
                                       : file("${row.bam_path}.bai")
                tuple(row.clone_id, bam, bai)
            }
            .groupTuple(by: 0)

        // Step 3: Merge BAMs per clone
        MERGE_BAMS(ch_clone_bam_groups)

        // Step 4: Mark duplicates (picard MarkDuplicates)
        //         Skip if params.mark_duplicates == false
        if (params.mark_duplicates) {
            MARKDUPLICATES(MERGE_BAMS.out.merged_bam)
            ch_processed_bams = MARKDUPLICATES.out.bam
        } else {
            ch_processed_bams = MERGE_BAMS.out.merged_bam
        }

        // Step 5: Coverage QC
        MOSDEPTH(ch_processed_bams, ch_fasta, ch_fai)
        SAMTOOLS_FLAGSTAT(ch_processed_bams)

        // Collect all QC for MultiQC
        ch_qc_reports = MOSDEPTH.out.summary
            .mix(SAMTOOLS_FLAGSTAT.out.flagstat)

    emit:
        clone_bams  = ch_processed_bams   // [ clone_id, bam, bai ]
        qc_reports  = ch_qc_reports
}
