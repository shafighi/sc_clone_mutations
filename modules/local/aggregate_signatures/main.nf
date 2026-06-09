process AGGREGATE_SIGNATURES {
    tag "aggregate_signatures"
    label 'process_single'
    publishDir "${params.outdir}/signatures", mode: params.publish_dir_mode

    container 'ghcr.io/shafighi/scclone-python:1.0.0'

    input:
        path exposures, stageAs: 'exp/*'

    output:
        path 'signatures_combined.tsv'
        path 'signatures_combined.png', emit: plot, optional: true

    script:
    """
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    mkdir -p "\$MPLCONFIGDIR"

    python3 ${projectDir}/bin/aggregate_signatures.py \\
        --exposures ${exposures.join(' ')} \\
        --out_tsv   signatures_combined.tsv \\
        --out_png   signatures_combined.png
    """
}
