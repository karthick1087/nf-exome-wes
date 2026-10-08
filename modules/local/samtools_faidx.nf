process SAMTOOLS_FAIDX {
    tag "$fasta"
    label 'process_low'

    input:
    path fasta

    output:
    path "${fasta}.fai", emit: fai
    path "versions.yml", emit: versions

    script:
    """
    if [[ ! -f "${fasta}.fai" ]]; then
        samtools faidx ${fasta}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
