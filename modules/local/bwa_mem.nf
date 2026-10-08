process BWA_MEM {
    tag "$meta.id"
    label 'process_high'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/align" }, mode: params.publish_dir_mode,
        pattern: '*.{bam,bai,flagstat.txt,stats.txt}'

    input:
    tuple val(meta), path(r1), path(r2)
    path  genome_files  // fasta + bwa indices only (do not also stage fasta separately)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam
    path "${meta.id}.flagstat.txt", emit: flagstat
    path "${meta.id}.stats.txt", emit: stats
    path "versions.yml", emit: versions

    script:
    def rg = meta.lane == 'results-gatk' \
        ? "@RG\\tID:${meta.id}\\tSM:${meta.name}\\tPL:ILLUMINA\\tLB:WES\\tPU:unit1" \
        : "@RG\\tID:${meta.id}\\tSM:${meta.name}\\tPL:ILLUMINA"
    """
    # Locate fasta among staged genome files
    FA=\$(ls *.fasta *.fa 2>/dev/null | head -1)
    if [[ -z "\$FA" ]]; then echo "ERROR: no fasta in genome files"; ls -la; exit 1; fi

    bwa mem -t ${task.cpus} \\
        -R "${rg}" \\
        "\$FA" ${r1} ${r2} \\
    | samtools view -@ ${task.cpus} -bS - \\
    | samtools sort -@ ${task.cpus} -m 2G -o ${meta.id}.sorted.bam -

    samtools index -@ ${task.cpus} ${meta.id}.sorted.bam
    samtools flagstat ${meta.id}.sorted.bam > ${meta.id}.flagstat.txt
    samtools stats ${meta.id}.sorted.bam > ${meta.id}.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -oE 'Version: [^ ]+' | head -1 || echo present)
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
