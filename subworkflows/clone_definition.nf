/*
================================================================================
  subworkflows/clone_definition.nf

  Assigns cells to clones using the MEDICC2 tree and scUnique events.

  Strategies (controlled by params.clone_strategy):
    internal_node  — cut tree at meaningful internal nodes (default)
    distance       — hierarchical clustering on tree pairwise distances
    event_profile  — cluster on scUnique event similarity (Jaccard)
    hybrid         — tree topology constrained by event similarity

  Outputs per-cell clone assignments + summary statistics.
================================================================================
*/

include { PARSE_MEDICC2_TREE } from '../modules/local/parse_medicc2_tree/main'
include { ASSIGN_CLONES      } from '../modules/local/assign_clones/main'
include { PLOT_CLONE_TREE    } from '../modules/local/plot_clone_tree/main'

workflow CLONE_DEFINITION {
    take:
        ch_bam_manifest         // validated BAM manifest CSV
        ch_cell_metadata        // validated cell metadata CSV
        ch_medicc2_tree         // Newick tree file
        ch_medicc2_events       // MEDICC2 branch-events TSV (may be empty)
        ch_scunique_events      // scUnique unique-events TSV

    main:
        // Step 1: Parse and annotate the MEDICC2 tree
        PARSE_MEDICC2_TREE(
            ch_medicc2_tree,
            ch_medicc2_events
        )

        // Step 2: Assign cells to clones
        ASSIGN_CLONES(
            PARSE_MEDICC2_TREE.out.tree_data,
            PARSE_MEDICC2_TREE.out.node_table,
            ch_scunique_events,
            ch_cell_metadata,
            ch_bam_manifest
        )

        // Step 3: Plot annotated clone tree (for QC and figures)
        PLOT_CLONE_TREE(
            PARSE_MEDICC2_TREE.out.tree_data,
            ASSIGN_CLONES.out.cell_clone_assignments
        )

    emit:
        cell_clone_assignments  = ASSIGN_CLONES.out.cell_clone_assignments
        clone_summary           = ASSIGN_CLONES.out.clone_summary
        clone_events            = ASSIGN_CLONES.out.clone_events
        tree_plot               = PLOT_CLONE_TREE.out.plot
}
