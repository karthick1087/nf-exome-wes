process BCFTOOLS_MPILEUP_CALL {
    tag "$meta.id"
    label 'process_high'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/variants" }, mode: params.publish_dir_mode,
        pattern: '*.{vcf.gz,tbi,txt}'

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai
    path  bed   // optional; may be empty list

    output:
    tuple val(meta), path("${meta.id}.raw.vcf.gz"), path("${meta.id}.raw.vcf.gz.tbi"), emit: raw
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), path("${meta.id}.filtered.vcf.gz.tbi"), emit: filtered
    path "${meta.id}_summary.txt", emit: summary
    path "versions.yml", emit: versions

    script:
    def bed_arg = bed.name != 'NO_FILE' ? "-R ${bed}" : ''
    def min_qual = params.min_qual
    def min_depth = params.min_depth
    """
    bcftools mpileup -Ou -f ${fasta} \\
        -a FORMAT/DP,FORMAT/AD,FORMAT/SP \\
        --threads ${task.cpus} \\
        ${bed_arg} ${bam} \\
    | bcftools call -mv -Oz --ploidy GRCh38 \\
        --threads ${task.cpus} -o ${meta.id}.raw.vcf.gz

    bcftools index -t ${meta.id}.raw.vcf.gz

    bcftools filter \\
        -i "QUAL>=${min_qual} && FORMAT/DP>=${min_depth}" \\
        -Oz -o ${meta.id}.filtered.vcf.gz ${meta.id}.raw.vcf.gz
    bcftools index -t ${meta.id}.filtered.vcf.gz

    {
        echo "WES → VCF Pipeline Summary (bcftools / Nextflow)"
        echo "=============================================="
        echo "Sample:          ${meta.name} (${meta.id})"
        echo "Date:            \$(date)"
        echo "Reference:       ${fasta}"
        echo "Min QUAL:        ${min_qual}"
        echo "Min Depth:       ${min_depth}"
        echo ""
        echo "Final VCF:       ${meta.id}.filtered.vcf.gz"
        echo "Variants (filt): \$(bcftools view -H ${meta.id}.filtered.vcf.gz | wc -l)"
        echo ""
        samtools flagstat ${bam}
    } > ${meta.id}_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}

process BCFTOOLS_VIEW_EXOME {
    tag "$meta.id"
    label 'process_low'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/variants" }, mode: params.publish_dir_mode,
        pattern: '*.{vcf.gz,tbi}'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  bed

    output:
    tuple val(meta), path("${meta.id}.exome.vcf.gz"), path("${meta.id}.exome.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    bed.name != 'NO_FILE'

    script:
    """
    bcftools view -R ${bed} -Oz -o ${meta.id}.exome.vcf.gz ${vcf}
    bcftools index -t ${meta.id}.exome.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
