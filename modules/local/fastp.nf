process FASTP {
    tag "$meta.id"
    label 'process_medium'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/qc" }, mode: params.publish_dir_mode, pattern: '*'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_R1.clean.fastq.gz"), path("${meta.id}_R2.clean.fastq.gz"), emit: reads
    path "${meta.id}_fastp.html", emit: html
    path "${meta.id}_fastp.json", emit: json
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: '--detect_adapter_for_pe -c'
    """
    fastp \\
        -i ${r1} -I ${r2} \\
        -o ${meta.id}_R1.clean.fastq.gz \\
        -O ${meta.id}_R2.clean.fastq.gz \\
        -h ${meta.id}_fastp.html \\
        -j ${meta.id}_fastp.json \\
        -w ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e 's/fastp //g')
    END_VERSIONS
    """
}
