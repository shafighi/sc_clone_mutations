process MUTECT2 {
    tag "${clone_id}:${intervals.baseName}"
    label 'process_high'

    container 'broadinstitute/gatk:4.6.0.0'

    input:
        tuple val(clone_id), path(tumor_cram), path(tumor_crai), path(normal_cram), path(normal_crai), path(intervals)
        path fasta
        path fai
        path dict
        path germline_resource      // may be []
        path germline_resource_tbi
        path pon                    // may be []
        path pon_tbi

    output:
        tuple val(clone_id), path("${clone_id}.*.mutect2.vcf.gz"), path("${clone_id}.*.mutect2.vcf.gz.tbi"), emit: vcf
        tuple val(clone_id), path("${clone_id}.*.mutect2.vcf.gz.stats"), emit: stats
        tuple val(clone_id), path("${clone_id}.*.f1r2.tar.gz"),          emit: f1r2     // for orientation bias model

    script:
    // One Mutect2 job per (clone x interval chunk). The chunk tag keeps output
    // names unique so the per-clone gather (MergeVcfs/MergeMutectStats) can stage
    // all chunks together without filename collisions.
    def chunk          = (intervals && intervals.name != 'NO_FILE') ? intervals.baseName : 'all'
    def prefix         = "${clone_id}.${chunk}"
    def normal_arg     = (normal_cram && normal_cram.name != 'NO_FILE') ? "-I ${normal_cram} --normal-sample ${normal_cram.baseName}" : ""
    def germline_arg   = (germline_resource && germline_resource.name != 'NO_FILE') ? "--germline-resource ${germline_resource}" : ""
    def pon_arg        = (pon && pon.name != 'NO_FILE') ? "--panel-of-normals ${pon}" : ""
    def intervals_arg  = (intervals && intervals.name != 'NO_FILE') ? "-L ${intervals}" : ""
    """
    gatk Mutect2 \\
        -R ${fasta} \\
        -I ${tumor_cram} \\
        ${normal_arg} \\
        ${germline_arg} \\
        ${pon_arg} \\
        ${intervals_arg} \\
        --f1r2-tar-gz ${prefix}.f1r2.tar.gz \\
        ${params.mutect2_extra_args} \\
        -O ${prefix}.mutect2.vcf.gz

    gatk IndexFeatureFile -I ${prefix}.mutect2.vcf.gz
    """
}
