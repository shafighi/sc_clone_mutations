process ASSIGN_CLONES {
    tag "clone_assignment"
    label 'process_low'
    publishDir "${params.outdir}/clone_definition", mode: params.publish_dir_mode

    container 'ghcr.io/shafighi/scclone-python:1.0.0'

    input:
        path tree_data          // pickled dendropy Tree
        path node_table         // node metadata CSV
        path scunique_events    // scUnique unique events TSV
        path cell_metadata      // validated cell metadata CSV
        path bam_manifest       // validated BAM manifest CSV

    output:
        path 'cell_clone_assignments.csv', emit: cell_clone_assignments
        path 'clone_summary.csv',          emit: clone_summary
        path 'clone_events.csv',           emit: clone_events
        path 'clone_tree_annotated.new',   emit: annotated_tree

    script:
    def dist_arg  = params.distance_threshold ? "--distance_threshold ${params.distance_threshold}" : ""
    def max_arg   = params.max_clones         ? "--max_clones ${params.max_clones}"                  : ""
    """
    assign_clones.py \\
        --tree_pkl             ${tree_data} \\
        --node_table           ${node_table} \\
        --scunique_events      ${scunique_events} \\
        --cell_metadata        ${cell_metadata} \\
        --bam_manifest         ${bam_manifest} \\
        --strategy             ${params.clone_strategy} \\
        --min_cells_per_clone  ${params.min_cells_per_clone} \\
        --min_branch_length    ${params.min_branch_length} \\
        ${dist_arg} \\
        --event_similarity_thr ${params.event_similarity_thr} \\
        ${max_arg} \\
        --small_clone_action   ${params.small_clone_action} \\
        --out_assignments      cell_clone_assignments.csv \\
        --out_summary          clone_summary.csv \\
        --out_events           clone_events.csv \\
        --out_tree             clone_tree_annotated.new
    """
}
