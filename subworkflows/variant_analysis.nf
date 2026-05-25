/*
================================================================================
  subworkflows/variant_analysis.nf

  Cross-clone and cross-caller variant comparison:
    1. Merge VCFs per clone across callers → caller-concordance table
    2. Compare presence/absence of variants across clones
    3. Build consensus callset (≥N callers)
    4. Optional VEP annotation
    5. Generate summary tables and plots
================================================================================
*/

include { COMPARE_VARIANTS   } from '../modules/local/compare_variants/main'
include { BUILD_CONSENSUS    } from '../modules/local/build_consensus/main'
include { ANNOTATE_VEP       } from '../modules/local/annotate_vep/main'

workflow VARIANT_ANALYSIS {
    take:
        ch_vcfs_per_clone      // [ clone_id, caller, vcf, tbi ]
        ch_clone_summary       // clone summary CSV
        ch_fasta
        ch_fai

    main:
        // Group all VCFs per clone: [ clone_id, [ [caller,vcf,tbi], ... ] ]
        ch_clone_vcf_groups = ch_vcfs_per_clone
            .map { clone_id, caller, vcf, tbi -> tuple(clone_id, tuple(caller, vcf, tbi)) }
            .groupTuple(by: 0)

        // Step 1: Per-clone caller concordance and cross-clone comparison
        COMPARE_VARIANTS(
            ch_vcfs_per_clone.collect(),
            ch_clone_summary
        )

        // Step 2: Build consensus callset
        BUILD_CONSENSUS(
            COMPARE_VARIANTS.out.variant_matrix,
            COMPARE_VARIANTS.out.per_clone_vcfs
        )

        // Step 3: Optional VEP annotation of consensus VCF
        if (params.annotate && params.vep_cache) {
            ANNOTATE_VEP(
                BUILD_CONSENSUS.out.consensus_vcf,
                params.vep_cache
            )
            ch_final_vcf = ANNOTATE_VEP.out.annotated_vcf
        } else {
            ch_final_vcf = BUILD_CONSENSUS.out.consensus_vcf
        }

    emit:
        consensus_table    = BUILD_CONSENSUS.out.consensus_table
        cross_clone_matrix = COMPARE_VARIANTS.out.variant_matrix
        final_vcf          = ch_final_vcf
        comparison_plots   = COMPARE_VARIANTS.out.plots
}
