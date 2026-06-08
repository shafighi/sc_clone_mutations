process MUTATIONAL_SIGNATURES {
    tag "${clone_id}"
    label 'process_low'
    publishDir "${params.outdir}/signatures/${clone_id}", mode: params.publish_dir_mode

    container 'ghcr.io/shafighi/scclone-python:1.0.0'

    input:
        tuple val(clone_id), val(caller), path(vcf), path(tbi)
        path fasta
        path fai
        path cosmic

    output:
        path "${clone_id}.sbs96_counts.tsv"
        path "${clone_id}.exposures.tsv",            emit: exposures
        path "${clone_id}.exposures_bootstrap.tsv"
        path "${clone_id}.signatures_audit.json",    emit: audit

    script:
    """
    python3 ${projectDir}/bin/fit_signatures.py \\
        --clone_id   ${clone_id} \\
        --vcf        ${vcf} \\
        --fasta      ${fasta} \\
        --cosmic     ${cosmic} \\
        --pass_only  ${params.consensus_pass_only} \\
        --tumor_only ${params.normal_manifest ? 'false' : 'true'} \\
        --min_snv    ${params.sig_min_snv} \\
        --min_cosine ${params.sig_min_cosine} \\
        --n_boot     ${params.sig_n_boot} \\
        --seed       ${params.sig_seed} \\
        --out_prefix ${clone_id}
    """
}
