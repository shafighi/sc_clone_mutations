process PARSE_MEDICC2_TREE {
    tag "parse_tree"
    label 'process_single'
    publishDir "${params.outdir}/clone_definition/tree", mode: params.publish_dir_mode

    container 'ghcr.io/shafighi/scclone-python:1.0.0'

    input:
        path medicc2_tree
        path medicc2_events   // may be empty []

    output:
        path 'tree_data.pkl',    emit: tree_data    // pickled dendropy Tree object
        path 'node_table.csv',   emit: node_table   // node metadata table
        path 'edge_table.csv',   emit: edge_table
        path 'pairwise_dist.csv',emit: pairwise_dist

    script:
    def events_arg = medicc2_events ? "--events_tsv ${medicc2_events}" : ""
    """
    parse_medicc2_tree.py \\
        --tree_file      ${medicc2_tree} \\
        ${events_arg} \\
        --out_pkl        tree_data.pkl \\
        --out_nodes      node_table.csv \\
        --out_edges      edge_table.csv \\
        --out_distances  pairwise_dist.csv
    """
}
