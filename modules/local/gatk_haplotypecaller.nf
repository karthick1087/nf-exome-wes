process GATK_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_high_memory'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/variants" }, mode: params.publish_dir_mode,
        pattern: '*.{vcf.gz,tbi}'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai
    path  dict
    path  bed   // optional NO_FILE

    output:
    tuple val(meta), path("${meta.id}.gatk.raw.vcf.gz"), path("${meta.id}.gatk.raw.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def bed_arg = bed.name != 'NO_FILE' ? "-L ${bed}" : ''
    """
    gatk HaplotypeCaller \\
        -R ${fasta} -I ${bam} \\
        -O ${meta.id}.gatk.raw.vcf.gz \\
        ${bed_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oE 'v[0-9.]+' | head -1 || echo 'unknown')
    END_VERSIONS
    """
}

process GATK_VARIANTFILTRATION {
    tag "$meta.id"
    label 'process_medium'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/variants" }, mode: params.publish_dir_mode,
        pattern: '*.{vcf.gz,tbi,txt}'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), path("${meta.id}.filtered.vcf.gz.tbi"), emit: filtered
    tuple val(meta), path("${meta.id}.pass.vcf.gz"), path("${meta.id}.pass.vcf.gz.tbi"), emit: pass
    path "${meta.id}_clinical_summary.txt", emit: summary
    path "versions.yml", emit: versions

    script:
    """
    gatk VariantFiltration \\
        -R ${fasta} -V ${vcf} -O ${meta.id}.filtered.vcf.gz \\
        --filter-name "FS60" --filter-expression "FS > 60.0" \\
        --filter-name "QD2"  --filter-expression "QD < 2.0" \\
        --filter-name "MQ40" --filter-expression "MQ < 40.0" \\
        --filter-name "SOR3"  --filter-expression "SOR > 3.0" \\
        --filter-name "DP10"  --filter-expression "DP < 10"

    bcftools view -f PASS -Oz -o ${meta.id}.pass.vcf.gz ${meta.id}.filtered.vcf.gz
    bcftools index -t ${meta.id}.pass.vcf.gz
    # filtered.vcf already indexed by GATK usually; ensure tbi
    bcftools index -f -t ${meta.id}.filtered.vcf.gz 2>/dev/null || true

    {
        echo "Clinical WES Pipeline Summary (GATK / Nextflow)"
        echo "=============================================="
        echo "Sample:     ${meta.name} (${meta.id})"
        echo "Date:       \$(date)"
        echo "Reference:  ${fasta}"
        echo ""
        echo "Outputs:"
        echo "  Raw VCF:      ${meta.id}.gatk.raw.vcf.gz"
        echo "  Filtered VCF: ${meta.id}.filtered.vcf.gz"
        echo "  PASS VCF:     ${meta.id}.pass.vcf.gz  (\$(bcftools view -H ${meta.id}.pass.vcf.gz | wc -l) variants)"
    } > ${meta.id}_clinical_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oE 'v[0-9.]+' | head -1 || echo 'unknown')
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}

process GATK_DICT {
    tag "$fasta"
    label 'process_low'

    input:
    path fasta

    output:
    path "*.dict", emit: dict
    path "versions.yml", emit: versions

    script:
    def base = fasta.name.replaceAll(/\.fa(sta)?$/, '')
    """
    # GATK expects hg38.dict not hg38.fasta.dict
    if [[ -f "${base}.dict" ]]; then
        cp -f "${base}.dict" ${base}.dict 2>/dev/null || true
    elif [[ -f "${fasta}.dict" ]]; then
        cp -f "${fasta}.dict" ${base}.dict
    else
        gatk CreateSequenceDictionary -R ${fasta} -O ${base}.dict
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oE 'v[0-9.]+' | head -1 || echo 'unknown')
    END_VERSIONS
    """
}

process GATK_COLLECTHSMETRICS {
    tag "$meta.id"
    label 'process_medium'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/qc" }, mode: params.publish_dir_mode,
        pattern: '*.txt'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai
    path  dict
    path  bed

    output:
    path "${meta.id}.hs_metrics.txt", optional: true, emit: metrics
    path "versions.yml", emit: versions

    when:
    bed.name != 'NO_FILE'

    script:
    """
    # Picard rejects a bare BED. Convert with the reference dictionary first.
    gatk BedToIntervalList \\
        -I ${bed} \\
        -O targets.interval_list \\
        -SD ${dict}

    gatk CollectHsMetrics \\
        -I ${bam} -R ${fasta} \\
        -BAIT_INTERVALS targets.interval_list \\
        -TARGET_INTERVALS targets.interval_list \\
        -O ${meta.id}.hs_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep -oE 'v[0-9.]+' | head -1 || echo 'unknown')
    END_VERSIONS
    """
}
