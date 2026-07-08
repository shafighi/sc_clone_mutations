#!/usr/bin/env nextflow

/*
================================================================================
  sc_clone_mutations — main workflow
  Clone-aware somatic mutation analysis from single-cell DNA sequencing data

  Inputs:
    1. MEDICC2 phylogenetic tree + events
    2. scUnique unique/recent events per cell
    3. Per-cell BAM manifest
    4. (Optional) matched-normal BAM manifest
    5. Reference genome resources

  Stages:
    A. Input validation
    B. Clone assignment (tree-based / event-based / hybrid)
    C. Pseudobulk BAM construction per clone
    D. Somatic mutation calling (Mutect2 + Octopus + FreeBayes)
    E. VCF normalization, filtering, and cross-clone comparison
    F. Consensus callset and reporting
================================================================================
*/

nextflow.enable.dsl = 2

// ─── Help / version banners ───────────────────────────────────────────────────

def helpMessage() {
    log.info """
    ╔══════════════════════════════════════════════════════════════════════════╗
    ║            sc_clone_mutations  v${workflow.manifest.version}                   ║
    ║     Clone-aware somatic mutation analysis from scDNA-seq                ║
    ╚══════════════════════════════════════════════════════════════════════════╝

    Usage:
      nextflow run main.nf [options]

    Required:
      --bam_manifest       CSV  cell_id,bam_path,bai_path,sample_id,patient_id
      --cell_metadata      CSV  cell_id,sample_id,patient_id,[...extra columns]
      --medicc2_tree       Newick tree file from MEDICC2 (*_final_tree.new)
      --scunique_events    TSV  per-cell unique/recent events from scUnique
      --fasta              Reference genome FASTA

    Optional:
      --normal_manifest    CSV  sample_id,patient_id,bam_path,bai_path
      --medicc2_events     TSV  MEDICC2 events-per-branch TSV
      --germline_resource  VCF.gz  germline variants for Mutect2 (e.g. gnomAD)
      --panel_of_normals   VCF.gz  PoN for Mutect2
      --intervals          BED  target intervals (whole genome if omitted)
      --vep_cache          DIR  VEP offline cache (required if --annotate)

    Clone assignment:
      --clone_strategy     [internal_node | distance | event_profile | hybrid]
                           default: internal_node
      --min_cells_per_clone INT   default: 5
      --min_branch_length  FLOAT  default: 0.0
      --distance_threshold FLOAT  required for distance strategy
      --small_clone_action [drop | merge | flag]  default: merge

    Callers:
      --callers            Comma-separated list: mutect2,octopus,freebayes
      --tumor_only         Run callers without matched normal

    Outputs:
      --outdir             Output directory (default: results)

    Profiles:
      -profile docker       Run locally with Docker
      -profile singularity  Run on HPC with Singularity/Apptainer
      -profile slurm        Submit jobs to SLURM (combine with singularity)
      -profile test         Use synthetic test inputs

    Examples:
      # Local Docker run
      nextflow run main.nf -profile docker \\
        --bam_manifest bam_manifest.csv \\
        --cell_metadata cell_metadata.csv \\
        --medicc2_tree tree.new \\
        --scunique_events events.tsv \\
        --fasta hg38.fa

      # HPC SLURM + Singularity
      nextflow run main.nf -profile singularity,slurm \\
        --bam_manifest bam_manifest.csv \\
        --cell_metadata cell_metadata.csv \\
        --medicc2_tree tree.new \\
        --scunique_events events.tsv \\
        --fasta /ref/hg38/hg38.fa \\
        --germline_resource /ref/gnomad.hg38.vcf.gz

      # Smoke test
      nextflow run main.nf -profile test,docker
    """.stripIndent()
}

if (params.help) { helpMessage(); exit 0 }
if (params.version) { log.info "sc_clone_mutations v${workflow.manifest.version}"; exit 0 }

// ─── Import subworkflows ──────────────────────────────────────────────────────

include { VALIDATE_INPUTS      } from './subworkflows/validate_inputs'
include { CLONE_DEFINITION     } from './subworkflows/clone_definition'
include { PSEUDOBULK           } from './subworkflows/pseudobulk'
include { MUTATION_CALLING     } from './subworkflows/mutation_calling'
include { VARIANT_ANALYSIS     } from './subworkflows/variant_analysis'
include { REPORTING            } from './subworkflows/reporting'
include { CONVERT_SCUNIQUE     } from './modules/local/convert_scunique/main'
include { MUTATIONAL_SIGNATURES } from './modules/local/mutational_signatures/main'
include { AGGREGATE_SIGNATURES  } from './modules/local/aggregate_signatures/main'

// ─── Parameter validation ─────────────────────────────────────────────────────

def validateRequiredParams() {
    def errors = []
    if (!params.bam_manifest && !params.bam_dir) errors << "  --bam_manifest or --bam_dir is required"
    if (!params.fasta)           errors << "  --fasta is required"

    // Either provide pre-processed files OR a scUnique results directory
    if (params.scunique_dir) {
        // RDS auto-conversion mode — tree, events, metadata derived from dir
        log.info "Using scUnique RDS auto-conversion mode (--scunique_dir)"
    } else {
        if (!params.cell_metadata)   errors << "  --cell_metadata is required (or use --scunique_dir)"
        if (!params.medicc2_tree)    errors << "  --medicc2_tree is required (or use --scunique_dir)"
        if (!params.scunique_events) errors << "  --scunique_events is required (or use --scunique_dir)"
    }

    if (params.clone_strategy == 'distance' && !params.distance_threshold)
        errors << "  --distance_threshold is required when --clone_strategy distance"
    if (params.annotate && !params.vep_cache)
        errors << "  --vep_cache is required when --annotate is set"
    if (errors) {
        log.error "Parameter validation failed:\n${errors.join('\n')}"
        exit 1
    }
}

if (params.validate_params) validateRequiredParams()

// ─── Main workflow ────────────────────────────────────────────────────────────

workflow {

    // ── Log run parameters ──────────────────────────────────────────────────
    log.info """
    ╔══════════════════════════════════════════════════════════════════════════╗
    ║         sc_clone_mutations  v${workflow.manifest.version}
    ╠══════════════════════════════════════════════════════════════════════════╣
    ║  BAM manifest       : ${params.bam_manifest}
    ║  Input mode         : ${params.scunique_dir ? 'scUnique RDS dir' : 'pre-processed files'}
    ║  scUnique dir       : ${params.scunique_dir ?: 'N/A'}
    ║  Cell metadata      : ${params.cell_metadata ?: '(from scUnique dir)'}
    ║  MEDICC2 tree       : ${params.medicc2_tree ?: '(from scUnique dir)'}
    ║  scUnique events    : ${params.scunique_events ?: '(from scUnique dir)'}
    ║  Reference FASTA    : ${params.fasta}
    ║  Normal manifest    : ${params.normal_manifest ?: 'not provided (tumor-only)'}
    ║  Clone strategy     : ${params.clone_strategy}
    ║  Callers            : ${params.callers.join(', ')}
    ║  Output dir         : ${params.outdir}
    ╚══════════════════════════════════════════════════════════════════════════╝
    """.stripIndent()

    // ── Build reference channel ─────────────────────────────────────────────
    // IMPORTANT: reference/resource files must be VALUE channels
    // (Channel.value), never queue channels (Channel.fromPath). A single-item
    // *queue* channel used as a secondary process input is consumed in lockstep
    // with the primary channel and limits the process to ONE execution. That
    // silently collapsed scatter-gather MUTECT2 (24 interval chunks) down to a
    // single chunk — only chrM + chr1:1-130.7Mb was ever called. Value channels
    // are re-emitted for every task, so the scatter fans out correctly.
    ch_fasta = Channel.value(file(params.fasta, checkIfExists: true))
    ch_fai   = params.fai  ? Channel.value(file(params.fai,  checkIfExists: true))
                           : ch_fasta.map { fa -> file("${fa}.fai") }
    ch_dict  = params.dict ? Channel.value(file(params.dict, checkIfExists: true))
                           : ch_fasta.map { fa -> file("${fa.parent}/${fa.baseName}.dict") }

    ch_germline_resource = params.germline_resource
        ? Channel.value(file(params.germline_resource, checkIfExists: true))
        : Channel.value([])
    ch_germline_resource_tbi = params.germline_resource_tbi
        ? Channel.value(file(params.germline_resource_tbi, checkIfExists: true))
        : Channel.value([])
    ch_pon = params.panel_of_normals
        ? Channel.value(file(params.panel_of_normals, checkIfExists: true))
        : Channel.value([])
    ch_pon_tbi = params.panel_of_normals_tbi
        ? Channel.value(file(params.panel_of_normals_tbi, checkIfExists: true))
        : Channel.value([])
    ch_intervals = params.intervals
        ? Channel.value(file(params.intervals, checkIfExists: true))
        : Channel.value([])

    // ── Parse input manifests ────────────────────────────────────────────────
    // Normal manifest
    ch_normal_manifest = params.normal_manifest
        ? Channel.fromPath(params.normal_manifest, checkIfExists: true)
        : Channel.value([])

    // ── Determine input mode: scUnique RDS dir vs pre-processed files ────────
    if (params.scunique_dir) {
        // Auto-convert RDS files to pipeline formats
        ch_scunique_dir = Channel.fromPath(params.scunique_dir, type: 'dir', checkIfExists: true)
        CONVERT_SCUNIQUE(ch_scunique_dir)

        ch_medicc2_tree    = CONVERT_SCUNIQUE.out.medicc2_tree
        ch_scunique_events = CONVERT_SCUNIQUE.out.scunique_events
        ch_cell_metadata   = params.cell_metadata
            ? Channel.fromPath(params.cell_metadata, checkIfExists: true)
            : CONVERT_SCUNIQUE.out.cell_metadata
        ch_bam_manifest    = params.bam_manifest
            ? Channel.fromPath(params.bam_manifest, checkIfExists: true)
            : CONVERT_SCUNIQUE.out.bam_manifest
    } else {
        // Standard mode: user provides pre-processed files
        ch_cell_metadata = Channel
            .fromPath(params.cell_metadata, checkIfExists: true)
        ch_medicc2_tree = Channel
            .fromPath(params.medicc2_tree, checkIfExists: true)
        ch_scunique_events = Channel
            .fromPath(params.scunique_events, checkIfExists: true)
        ch_bam_manifest = Channel
            .fromPath(params.bam_manifest, checkIfExists: true)
    }

    ch_medicc2_events = params.medicc2_events
        ? Channel.fromPath(params.medicc2_events, checkIfExists: true)
        : Channel.value([])

    // ── Stage A: Validate inputs ─────────────────────────────────────────────
    VALIDATE_INPUTS(
        ch_bam_manifest,
        ch_cell_metadata,
        ch_medicc2_tree,
        ch_scunique_events,
        ch_normal_manifest,
        ch_fasta,
        ch_fai
    )

    // ── Stage B: Clone assignment ────────────────────────────────────────────
    CLONE_DEFINITION(
        VALIDATE_INPUTS.out.bam_manifest_validated,
        VALIDATE_INPUTS.out.cell_metadata_validated,
        ch_medicc2_tree,
        ch_medicc2_events,
        ch_scunique_events
    )

    // ── Stage C: Pseudobulk BAM construction ─────────────────────────────────
    PSEUDOBULK(
        CLONE_DEFINITION.out.cell_clone_assignments,
        VALIDATE_INPUTS.out.bam_manifest_validated,
        ch_fasta,
        ch_fai
    )

    // ── Stage D: Somatic mutation calling ────────────────────────────────────
    MUTATION_CALLING(
        PSEUDOBULK.out.clone_bams,        // [ clone_id, bam, bai ]
        ch_normal_manifest,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_resource_tbi,
        ch_pon,
        ch_pon_tbi,
        ch_intervals
    )

    // ── Stage E: Variant analysis ────────────────────────────────────────────
    VARIANT_ANALYSIS(
        MUTATION_CALLING.out.vcfs_per_clone,   // [ clone_id, caller, vcf, tbi ]
        CLONE_DEFINITION.out.clone_summary,
        ch_fasta,
        ch_fai
    )

    // ── Stage E2: Mutational signatures (opt-in: --run_signatures) ───────────
    // Per-clone SBS96 refit to COSMIC, with an audit that can return
    // 'insufficient_evidence'. Consumes the normalized per-caller VCFs; only the
    // configured caller is used so there is one fit per clone. Terminal stage —
    // does not affect any existing process hash.
    ch_sig_exposures = Channel.value([])
    ch_sig_plots     = Channel.value([])
    if (params.run_signatures) {
        if (!params.cosmic_signatures)
            error "ERROR: --run_signatures requires --cosmic_signatures (COSMIC SBS reference). See assets/cosmic/README.md"
        ch_cosmic = Channel.value(file(params.cosmic_signatures, checkIfExists: true))
        ch_sig_vcfs = MUTATION_CALLING.out.vcfs_per_clone
            .filter { clone_id, caller, vcf, tbi -> caller == params.sig_caller }
        MUTATIONAL_SIGNATURES(ch_sig_vcfs, ch_fasta, ch_fai, ch_cosmic)
        ch_sig_exposures = MUTATIONAL_SIGNATURES.out.exposures.collect()

        // Cross-clone combined exposure plot (one figure across all clones)
        AGGREGATE_SIGNATURES(MUTATIONAL_SIGNATURES.out.exposures.collect())
        ch_sig_plots = MUTATIONAL_SIGNATURES.out.plot
            .mix(AGGREGATE_SIGNATURES.out.plot)
            .collect()
    }

    // ── Stage F: Reporting ────────────────────────────────────────────────────
    REPORTING(
        VALIDATE_INPUTS.out.validation_report,
        CLONE_DEFINITION.out.clone_summary,
        PSEUDOBULK.out.qc_reports,
        MUTATION_CALLING.out.caller_stats,
        VARIANT_ANALYSIS.out.consensus_table,
        VARIANT_ANALYSIS.out.cross_clone_matrix,
        CLONE_DEFINITION.out.tree_plot_png,
        PSEUDOBULK.out.reliability,
        ch_sig_exposures,
        ch_sig_plots
    )
}

// ─── Workflow completion handler ──────────────────────────────────────────────

workflow.onComplete {
    def status = workflow.success ? "SUCCESS" : "FAILED"
    log.info """
    ══════════════════════════════════════════════════════════
    Pipeline ${status}
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Results     : ${params.outdir}
    ══════════════════════════════════════════════════════════
    """
    if (params.email && workflow.success) {
        sendMail(
            to:      params.email,
            subject: "sc_clone_mutations completed: ${workflow.runName}",
            body:    "Pipeline finished successfully. Results: ${params.outdir}"
        )
    }
}

workflow.onError {
    log.error "Pipeline error: ${workflow.errorMessage}"
}
