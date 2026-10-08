process GATK_BASERECALIBRATOR {
    tag "$meta.id"
    label 'process_high'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/align" }, mode: params.publish_dir_mode,
        pattern: '*.table'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai
    path  dict
    path  dbsnp
    path  mills   // optional NO_FILE

    output:
    tuple val(meta), path("${meta.id}.recal.table"), emit: table
    path "versions.yml", emit: versions

    when:
    dbsnp.name != 'NO_FILE'

    script:
    def use_mills = mills.name != 'NO_FILE'
    """
    # Ensure known-sites VCF has a local index (Nextflow often stages only the .vcf.gz)
    stage_known() {
        local src="\$1" dst="\$2"
        local real
        real=\$(readlink -f "\$src")
        cp -L "\$src" "\$dst"
        if [[ -f "\${real}.tbi" ]]; then
            cp -L "\${real}.tbi" "\${dst}.tbi"
        elif [[ -f "\${real}.idx" ]]; then
            cp -L "\${real}.idx" "\${dst}.idx"
        else
            gatk IndexFeatureFile -I "\$dst"
        fi
    }

    stage_known ${dbsnp} dbsnp.known.vcf.gz
    KS_ARGS=(--known-sites dbsnp.known.vcf.gz)

    if [[ "${use_mills}" == "true" ]]; then
        stage_known ${mills} mills.known.vcf.gz
        KS_ARGS+=(--known-sites mills.known.vcf.gz)
    fi

    gatk BaseRecalibrator \\
        -R ${fasta} -I ${bam} \\
        "\${KS_ARGS[@]}" \\
        -O ${meta.id}.recal.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oE 'v[0-9.]+' | head -1 || echo 'unknown')
    END_VERSIONS
    """
}

process GATK_APPLYBQSR {
    tag "$meta.id"
    label 'process_high'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/align" }, mode: params.publish_dir_mode,
        pattern: '*.{bam,bai}'

    input:
    tuple val(meta), path(bam), path(bai), path(table)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("${meta.id}.recal.bam"), path("${meta.id}.recal.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    """
    gatk ApplyBQSR \\
        -R ${fasta} -I ${bam} \\
        --bqsr-recal-file ${table} \\
        -O ${meta.id}.recal.bam \\
        --create-output-bam-index true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oE 'v[0-9.]+' | head -1 || echo 'unknown')
    END_VERSIONS
    """
}
