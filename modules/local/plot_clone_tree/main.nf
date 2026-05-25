process PLOT_CLONE_TREE {
    tag "plot_tree"
    label 'process_single'
    publishDir "${params.outdir}/clone_definition/figures", mode: params.publish_dir_mode

    container 'ghcr.io/TODO/scclone-python:1.0.0'

    input:
        path tree_data
        path cell_clone_assignments

    output:
        path 'clone_tree.pdf',  emit: plot
        path 'clone_tree.png',  emit: plot_png

    script:
    """
    plot_clone_tree.py \\
        --tree_pkl       ${tree_data} \\
        --assignments    ${cell_clone_assignments} \\
        --out_pdf        clone_tree.pdf \\
        --out_png        clone_tree.png
    """
}
