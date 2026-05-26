process COMPARE_VARIANTS {
    tag "cross_clone_comparison"
    label 'process_medium'
    publishDir "${params.outdir}/variant_analysis", mode: params.publish_dir_mode

    container 'ghcr.io/shafighi/scclone-python:1.0.0'

    input:
        path vcf_list     // all normalized VCFs as a collected list
        path clone_summary

    output:
        path 'variant_matrix.csv',         emit: variant_matrix
        path 'caller_concordance.csv',     emit: caller_concordance
        path 'private_shared_table.csv',   emit: private_shared
        path 'per_clone_vcfs',             emit: per_clone_vcfs
        path 'comparison_plots',           emit: plots

    script:
    """
    mkdir -p per_clone_vcfs comparison_plots

    compare_variants.py \\
        --vcf_list         ${vcf_list.join(' ')} \\
        --clone_summary    ${clone_summary} \\
        --out_matrix       variant_matrix.csv \\
        --out_concordance  caller_concordance.csv \\
        --out_private      private_shared_table.csv \\
        --out_per_clone    per_clone_vcfs/ \\
        --out_plots        comparison_plots/ \\
        --pass_only        ${params.consensus_pass_only}
    """
}
