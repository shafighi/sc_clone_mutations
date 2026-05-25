process MUTECT2 {
    tag "${clone_id}"
    label 'process_high'
    publishDir "${params.outdir}/variant_calling/mutect2/raw", mode: params.publish_dir_mode

    container 'broadinstitute/gatk:4.6.0.0'

    input:
        tuple val(clone_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
        path fasta
        path fai
        path dict
        path germline_resource      // may be []
        path germline_resource_tbi
        path pon                    // may be []
        path pon_tbi
        path intervals              // may be []

    output:
        tuple val(clone_id), path("${clone_id}.mutect2.vcf.gz"), path("${clone_id}.mutect2.vcf.gz.tbi"), emit: vcf
        path "${clone_id}.mutect2.vcf.gz.stats", emit: stats
        path "${clone_id}.f1r2.tar.gz",          emit: f1r2     // for orientation bias model

    script:
    def normal_arg     = (normal_bam && normal_bam.name != 'NO_FILE') ? "-I ${normal_bam} --normal-sample ${normal_bam.baseName}" : ""
    def germline_arg   = (germline_resource && germline_resource.name != 'NO_FILE') ? "--germline-resource ${germline_resource}" : ""
    def pon_arg        = (pon && pon.name != 'NO_FILE') ? "--panel-of-normals ${pon}" : ""
    def intervals_arg  = (intervals && intervals.name != 'NO_FILE') ? "-L ${intervals}" : ""
    """
    gatk Mutect2 \\
        -R ${fasta} \\
        -I ${tumor_bam} \\
        ${normal_arg} \\
        ${germline_arg} \\
        ${pon_arg} \\
        ${intervals_arg} \\
        --f1r2-tar-gz ${clone_id}.f1r2.tar.gz \\
        ${params.mutect2_extra_args} \\
        -O ${clone_id}.mutect2.vcf.gz

    gatk IndexFeatureFile -I ${clone_id}.mutect2.vcf.gz
    """
}
