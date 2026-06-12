/*
================================================================================
  subworkflows/mutation_calling.nf

  Runs somatic mutation callers on each clone pseudobulk CRAM:
    - Mutect2  (GATK4)
    - Octopus  (haplotype-aware, replaces Strelka2)
    - FreeBayes

  Supports:
    - Paired tumor-normal mode (when normal_manifest is provided)
    - Tumor-only mode          (when --tumor_only or no normal)

  Per-clone outputs: raw VCF, filtered/PASS VCF, normalized VCF
================================================================================
*/

include { SPLIT_INTERVALS     } from '../modules/local/split_intervals/main'
include { MUTECT2             } from '../modules/local/mutect2/main'
include { MERGE_VCFS          } from '../modules/local/merge_vcfs/main'
include { MERGE_MUTECT_STATS  } from '../modules/local/merge_mutect_stats/main'
include { FILTER_MUTECT2      } from '../modules/local/filter_mutect2/main'
include { OCTOPUS             } from '../modules/local/octopus/main'
include { FREEBAYES           } from '../modules/local/freebayes/main'
include { NORMALIZE_VCF as NORMALIZE_MUTECT2_VCF   } from '../modules/local/normalize_vcf/main'
include { NORMALIZE_VCF as NORMALIZE_OCTOPUS_VCF   } from '../modules/local/normalize_vcf/main'
include { NORMALIZE_VCF as NORMALIZE_FREEBAYES_VCF } from '../modules/local/normalize_vcf/main'
include { MERGE_CALLER_VCFS   } from '../modules/local/merge_caller_vcfs/main'

workflow MUTATION_CALLING {
    take:
        ch_clone_bams            // [ clone_id, bam, bai ]
        ch_normal_manifest       // CSV or empty channel
        ch_fasta
        ch_fai
        ch_dict
        ch_germline_resource
        ch_germline_resource_tbi
        ch_pon
        ch_pon_tbi
        ch_intervals

    main:
        def callers = params.callers instanceof List
            ? params.callers
            : params.callers.tokenize(',')*.trim()

        // Parse normal manifest if provided
        ch_normals = params.normal_manifest && !params.tumor_only
            ? Channel.fromPath(params.normal_manifest)
                .splitCsv(header: true)
                .map { row ->
                    def bam = file(row.bam_path, checkIfExists: true)
                    def bai = row.bai_path ? file(row.bai_path, checkIfExists: true)
                                           : file("${row.bam_path}.bai")
                    tuple(row.sample_id, bam, bai)
                }
            : Channel.empty()

        // Join clone BAMs with normals on sample_id
        // clone_bams schema: [ clone_id, bam, bai ] — need to add sample_id
        // For now, treat each clone independently in tumor-only mode
        // TODO: if a clone belongs to a patient with a matched normal, join here
        ch_tumor_normal_pairs = params.tumor_only
            ? ch_clone_bams.map { clone_id, bam, bai -> tuple(clone_id, bam, bai, [], []) }
            : ch_clone_bams.map { clone_id, bam, bai -> tuple(clone_id, bam, bai, [], []) }
            // TODO: implement proper sample-level normal joining when normal_manifest provided

        ch_vcfs = Channel.empty()

        // ── Mutect2 (scatter by interval chunk → gather per clone) ───────────
        if ('mutect2' in callers) {
            // Split the calling region into N chunks so each clone is called in
            // parallel jobs that each stay well under the per-job time cap.
            SPLIT_INTERVALS(ch_fasta, ch_fai, ch_dict, ch_intervals)

            // One Mutect2 job per (clone × chunk): combine each clone with every chunk.
            ch_mutect2_in = ch_tumor_normal_pairs
                .combine(SPLIT_INTERVALS.out.intervals.flatten())

            MUTECT2(
                ch_mutect2_in,
                ch_fasta,
                ch_fai,
                ch_dict,
                ch_germline_resource,
                ch_germline_resource_tbi,
                ch_pon,
                ch_pon_tbi
            )

            // Gather all chunks back per clone.
            MERGE_VCFS( MUTECT2.out.vcf.groupTuple(by: 0) )
            MERGE_MUTECT_STATS( MUTECT2.out.stats.groupTuple(by: 0) )

            // Join merged VCF + merged stats by clone_id (robust against process
            // ordering) before filtering.
            ch_filter_in = MERGE_VCFS.out.vcf.join(MERGE_MUTECT_STATS.out.stats)

            FILTER_MUTECT2(
                ch_filter_in,
                ch_fasta,
                ch_fai,
                ch_dict
            )
            NORMALIZE_MUTECT2_VCF(
                FILTER_MUTECT2.out.filtered_vcf.map { id, vcf, tbi -> tuple(id, 'mutect2', vcf, tbi) },
                ch_fasta
            )
            ch_vcfs = ch_vcfs.mix(NORMALIZE_MUTECT2_VCF.out.normalized_vcf)
        }

        // ── Octopus ───────────────────────────────────────────────────────────
        if ('octopus' in callers) {
            OCTOPUS(
                ch_tumor_normal_pairs,
                ch_fasta,
                ch_fai,
                ch_intervals
            )
            NORMALIZE_OCTOPUS_VCF(
                OCTOPUS.out.vcf.map { id, vcf, tbi -> tuple(id, 'octopus', vcf, tbi) },
                ch_fasta
            )
            ch_vcfs = ch_vcfs.mix(NORMALIZE_OCTOPUS_VCF.out.normalized_vcf)
        }

        // ── FreeBayes ────────────────────────────────────────────────────────
        if ('freebayes' in callers) {
            FREEBAYES(
                ch_clone_bams,
                ch_fasta,
                ch_fai,
                ch_intervals
            )
            NORMALIZE_FREEBAYES_VCF(
                FREEBAYES.out.vcf.map { id, vcf -> tuple(id, 'freebayes', vcf, []) },
                ch_fasta
            )
            ch_vcfs = ch_vcfs.mix(NORMALIZE_FREEBAYES_VCF.out.normalized_vcf)
        }

        // Collect caller-level statistics
        ch_caller_stats = Channel.empty()
        if ('mutect2'   in callers) ch_caller_stats = ch_caller_stats.mix(MERGE_MUTECT_STATS.out.stats.map { id, s -> s })
        if ('octopus'   in callers) ch_caller_stats = ch_caller_stats.mix(OCTOPUS.out.stats)
        if ('freebayes' in callers) ch_caller_stats = ch_caller_stats.mix(FREEBAYES.out.stats)

    emit:
        vcfs_per_clone = ch_vcfs            // [ clone_id, caller, vcf, tbi ]
        caller_stats   = ch_caller_stats
}
