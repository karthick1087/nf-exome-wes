process BWA_INDEX {
    tag "$fasta"
    label 'process_high'

    input:
    path fasta

    output:
    path "${fasta}*", emit: index
    path "versions.yml", emit: versions

    script:
    """
    # Reuse prebuilt index files staged next to the fasta when present
    if [[ ! -f "${fasta}.bwt" ]]; then
        bwa index -a bwtsw ${fasta}
    fi
    if [[ ! -f "${fasta}.fai" ]]; then
        samtools faidx ${fasta}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
