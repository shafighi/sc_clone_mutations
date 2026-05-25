process ANNOTATE_VEP {
    tag "vep_annotation"
    label 'process_high_memory'
    publishDir "${params.outdir}/variant_analysis/annotated", mode: params.publish_dir_mode

    // VEP requires a local cache — use ensemblorg/ensembl-vep Docker image
    // Note: VEP is dual-licensed (Apache 2.0 / Perl Artistic License).
    // Cache must be pre-downloaded: vep_install -a c -s homo_sapiens -y GRCh38
    container 'ensemblorg/ensembl-vep:release_112.0'

    input:
        tuple val(meta_id), path(vcf), path(tbi)
        val vep_cache_dir

    output:
        tuple val(meta_id), path("*.vep.vcf.gz"), path("*.vep.vcf.gz.tbi"), emit: annotated_vcf
        path "*.vep_summary.html", emit: summary

    script:
    // TODO: adjust VEP flags to match your annotation requirements
    """
    vep \\
        --input_file        ${vcf} \\
        --output_file       ${meta_id}.vep.vcf.gz \\
        --format            vcf \\
        --vcf \\
        --compress_output   bgzip \\
        --cache \\
        --offline \\
        --dir_cache         ${vep_cache_dir} \\
        --species           ${params.vep_species} \\
        --assembly          ${params.vep_genome_build} \\
        --fork              ${task.cpus} \\
        --everything \\
        --stats_file        ${meta_id}.vep_summary.html

    tabix -p vcf ${meta_id}.vep.vcf.gz
    """
}
