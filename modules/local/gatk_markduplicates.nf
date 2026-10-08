process GATK_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_high'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/align" }, mode: params.publish_dir_mode,
        pattern: '*.{bam,bai,txt}'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.dedup.bam"), path("${meta.id}.dedup.bai"), emit: bam
    path "${meta.id}.dedup_metrics.txt", emit: metrics
    path "versions.yml", emit: versions

    script:
    """
    gatk MarkDuplicates \\
        -I ${bam} \\
        -O ${meta.id}.dedup.bam \\
        -M ${meta.id}.dedup_metrics.txt \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY SILENT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oE 'v[0-9.]+' | head -1 || echo 'unknown')
    END_VERSIONS
    """
}
