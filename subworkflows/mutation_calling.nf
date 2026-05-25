/*
================================================================================
  subworkflows/mutation_calling.nf

  Runs somatic mutation callers on each clone pseudobulk BAM:
    - Mutect2  (GATK4)
    - Strelka2
    - FreeBayes

  Supports:
    - Paired tumor-normal mode (when normal_manifest is provided)
    - Tumor-only mode          (when --tumor_only or no normal)

  Per-clone outputs: raw VCF, filtered/PASS VCF, normalized VCF
================================================================================
*/

include { MUTECT2             } from '../modules/local/mutect2/main'
include { FILTER_MUTECT2      } from '../modules/local/filter_mutect2/main'
include { STRELKA2_SOMATIC    } from '../modules/local/strelka2/main'
include { FREEBAYES           } from '../modules/local/freebayes/main'
include { NORMALIZE_VCF as NORMALIZE_MUTECT2_VCF   } from '../modules/local/normalize_vcf/main'
include { NORMALIZE_VCF as NORMALIZE_STRELKA2_VCF  } from '../modules/local/normalize_vcf/main'
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

        // ── Mutect2 ──────────────────────────────────────────────────────────
        if ('mutect2' in callers) {
            MUTECT2(
                ch_tumor_normal_pairs,
                ch_fasta,
                ch_fai,
                ch_dict,
                ch_germline_resource,
                ch_germline_resource_tbi,
                ch_pon,
                ch_pon_tbi,
                ch_intervals
            )
            FILTER_MUTECT2(
                MUTECT2.out.vcf,
                MUTECT2.out.stats,
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

        // ── Strelka2 ─────────────────────────────────────────────────────────
        if ('strelka2' in callers) {
            STRELKA2_SOMATIC(
                ch_tumor_normal_pairs,
                ch_fasta,
                ch_fai,
                ch_intervals
            )
            NORMALIZE_STRELKA2_VCF(
                STRELKA2_SOMATIC.out.vcf.map { id, vcf, tbi -> tuple(id, 'strelka2', vcf, tbi) },
                ch_fasta
            )
            ch_vcfs = ch_vcfs.mix(NORMALIZE_STRELKA2_VCF.out.normalized_vcf)
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
        if ('mutect2'   in callers) ch_caller_stats = ch_caller_stats.mix(MUTECT2.out.stats)
        if ('strelka2'  in callers) ch_caller_stats = ch_caller_stats.mix(STRELKA2_SOMATIC.out.stats)
        if ('freebayes' in callers) ch_caller_stats = ch_caller_stats.mix(FREEBAYES.out.stats)

    emit:
        vcfs_per_clone = ch_vcfs            // [ clone_id, caller, vcf, tbi ]
        caller_stats   = ch_caller_stats
}
