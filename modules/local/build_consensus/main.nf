process BUILD_CONSENSUS {
    tag "build_consensus"
    label 'process_medium'
    publishDir "${params.outdir}/variant_analysis/consensus", mode: params.publish_dir_mode

    container 'ghcr.io/shafighi/scclone-python:1.0.0'

    input:
        path variant_matrix
        path per_clone_vcf_dir

    output:
        path 'consensus_mutations.vcf.gz',     emit: consensus_vcf
        path 'consensus_mutations.vcf.gz.tbi'
        path 'consensus_table.csv',            emit: consensus_table
        path 'consensus_summary.json',         emit: summary

    script:
    """
    build_consensus.py \\
        --variant_matrix     ${variant_matrix} \\
        --per_clone_vcf_dir  ${per_clone_vcf_dir} \\
        --min_callers        ${params.consensus_min_callers} \\
        --pass_only          ${params.consensus_pass_only} \\
        --out_vcf            consensus_mutations.vcf.gz \\
        --out_table          consensus_table.csv \\
        --out_summary        consensus_summary.json

    tabix -p vcf consensus_mutations.vcf.gz
    """
}
