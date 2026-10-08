/*
 * Stage reference + existing indices via symlinks (no multi-hour reindex).
 */
process PREPARE_REF {
    tag "$ref_path"
    label 'process_low'

    input:
    val ref_path

    output:
    path "genome/*", emit: genome_files
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p genome
    REF="${ref_path}"
    BASE="\$(basename "\${REF}")"
    PARENT="\$(dirname "\${REF}")"
    SIMPLE="\${BASE%.fasta}"
    SIMPLE="\${SIMPLE%.fa}"

    ln -sf "\${REF}" "genome/\${BASE}"

    for ext in fai bwt pac amb ann sa dict; do
        if [[ -e "\${REF}.\${ext}" ]]; then
            ln -sf "\${REF}.\${ext}" "genome/\${BASE}.\${ext}"
        fi
    done

    if [[ -e "\${PARENT}/\${SIMPLE}.dict" ]]; then
        ln -sf "\${PARENT}/\${SIMPLE}.dict" "genome/\${SIMPLE}.dict"
    elif [[ -e "genome/\${BASE}.dict" && ! -e "genome/\${SIMPLE}.dict" ]]; then
        ln -sf "\${BASE}.dict" "genome/\${SIMPLE}.dict"
    fi

    if [[ ! -e "genome/\${BASE}.fai" ]]; then
        samtools faidx "genome/\${BASE}"
    fi
    if [[ ! -e "genome/\${BASE}.bwt" ]]; then
        echo "WARNING: no prebuilt BWA index — building (slow)..."
        bwa index -a bwtsw "genome/\${BASE}"
    fi
    if [[ ! -e "genome/\${SIMPLE}.dict" ]]; then
        if command -v gatk >/dev/null 2>&1; then
            gatk CreateSequenceDictionary -R "genome/\${BASE}" -O "genome/\${SIMPLE}.dict"
        else
            samtools dict "genome/\${BASE}" -o "genome/\${SIMPLE}.dict"
        fi
    fi

    # materialize symlinks to real files for Nextflow staging of outputs
    # (keep symlinks — NF follows them when publishing inputs downstream)

    ls -la genome/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
