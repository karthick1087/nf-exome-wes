process BCFTOOLS_ANNOTATE_DBSNP {
    tag "$meta.id"
    label 'process_medium'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/annotation/bcftools" }, mode: params.publish_dir_mode,
        pattern: '*.{vcf.gz,tbi}'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  dbsnp

    output:
    tuple val(meta), path("${meta.id}.dbsnp.vcf.gz"), path("${meta.id}.dbsnp.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    dbsnp.name != 'NO_FILE'

    script:
    """
    # Stage annotation VCF with index (Nextflow may only stage the .vcf.gz symlink)
    ensure_ann() {
        local src="\$1" dst="\$2"
        local real
        real=\$(readlink -f "\$src")
        cp -L "\$src" "\$dst"
        if [[ -f "\${real}.tbi" ]]; then
            cp -L "\${real}.tbi" "\${dst}.tbi"
        elif [[ -f "\${real}.csi" ]]; then
            cp -L "\${real}.csi" "\${dst}.csi"
        elif [[ -f "\${src}.tbi" ]]; then
            cp -L "\${src}.tbi" "\${dst}.tbi"
        else
            bcftools index -t "\$dst" || tabix -p vcf "\$dst"
        fi
    }
    ensure_ann ${dbsnp} dbsnp.ann.vcf.gz
    bcftools annotate -a dbsnp.ann.vcf.gz -c ID -Oz -o ${meta.id}.dbsnp.vcf.gz ${vcf} || cp ${vcf} ${meta.id}.dbsnp.vcf.gz
    bcftools index -t ${meta.id}.dbsnp.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}

process BCFTOOLS_ANNOTATE_CLINVAR {
    tag "$meta.id"
    label 'process_medium'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/annotation/bcftools" }, mode: params.publish_dir_mode,
        pattern: '*.{vcf.gz,tbi}'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  clinvar

    output:
    tuple val(meta), path("${meta.id}.clinvar.vcf.gz"), path("${meta.id}.clinvar.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    clinvar.name != 'NO_FILE'

    script:
    """
    ensure_ann() {
        local src="\$1" dst="\$2"
        local real
        real=\$(readlink -f "\$src")
        cp -L "\$src" "\$dst"
        if [[ -f "\${real}.tbi" ]]; then
            cp -L "\${real}.tbi" "\${dst}.tbi"
        elif [[ -f "\${real}.csi" ]]; then
            cp -L "\${real}.csi" "\${dst}.csi"
        elif [[ -f "\${src}.tbi" ]]; then
            cp -L "\${src}.tbi" "\${dst}.tbi"
        else
            bcftools index -t "\$dst" || tabix -p vcf "\$dst"
        fi
    }
    ensure_ann ${clinvar} clinvar.ann.vcf.gz
    bcftools annotate -a clinvar.ann.vcf.gz \\
        -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNHGVS,INFO/CLNDISDB,INFO/GENEINFO,INFO/ORIGIN \\
        -Oz -o ${meta.id}.clinvar.vcf.gz ${vcf}
    bcftools index -t ${meta.id}.clinvar.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}

process BCFTOOLS_ANNOTATE_GNOMAD {
    tag "$meta.id"
    label 'process_medium'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/annotation/bcftools" }, mode: params.publish_dir_mode,
        pattern: '*.{vcf.gz,tbi}'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  gnomad

    output:
    tuple val(meta), path("${meta.id}.gnomad.vcf.gz"), path("${meta.id}.gnomad.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    gnomad.name != 'NO_FILE'

    script:
    """
    ensure_ann() {
        local src="\$1" dst="\$2"
        local real
        real=\$(readlink -f "\$src")
        cp -L "\$src" "\$dst"
        if [[ -f "\${real}.tbi" ]]; then
            cp -L "\${real}.tbi" "\${dst}.tbi"
        elif [[ -f "\${real}.csi" ]]; then
            cp -L "\${real}.csi" "\${dst}.csi"
        elif [[ -f "\${src}.tbi" ]]; then
            cp -L "\${src}.tbi" "\${dst}.tbi"
        else
            bcftools index -t "\$dst" || tabix -p vcf "\$dst"
        fi
    }
    ensure_ann ${gnomad} gnomad.ann.vcf.gz
    bcftools annotate -a gnomad.ann.vcf.gz \\
        -c INFO/AF,INFO/AF_popmax,INFO/AF_eas,INFO/AF_nfe,INFO/AF_sas \\
        -Oz -o ${meta.id}.gnomad.vcf.gz ${vcf}
    bcftools index -t ${meta.id}.gnomad.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}

process BCFTOOLS_FINALIZE_ANNOT {
    tag "$meta.id"
    label 'process_low'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/annotation/bcftools" }, mode: params.publish_dir_mode,
        pattern: '*.{vcf.gz,tbi}'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}.fully_annotated.vcf.gz"), path("${meta.id}.fully_annotated.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${meta.id}.step0.vcf.gz"), path("${meta.id}.step0.vcf.gz.tbi"), emit: step0
    path "versions.yml", emit: versions

    script:
    """
    cp -f ${vcf} ${meta.id}.fully_annotated.vcf.gz
    bcftools index -t ${meta.id}.fully_annotated.vcf.gz
    # step0 mirror (same layout as legacy annotate_vcf.sh)
    cp -f ${vcf} ${meta.id}.step0.vcf.gz
    bcftools index -t ${meta.id}.step0.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}

process SNPEFF_ANNOTATE {
    tag "$meta.id"
    label 'process_high_memory'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/annotation/snpeff" }, mode: params.publish_dir_mode,
        pattern: '*.{csv,txt}'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path snpeff_jar, stageAs: 'snpEff.jar'
    path snpeff_data, stageAs: 'snpeff_data'

    output:
    tuple val(meta), path("${meta.id}.snpeff.vcf"), emit: vcf
    path "${meta.id}_snpeff.csv", optional: true, emit: csv
    path "${meta.id}_snpeff.genes.txt", optional: true, emit: genes
    path "versions.yml", emit: versions

    when:
    !params.skip_snpeff && !params.skip_annotation

    script:
    def db = params.snpeff_db
    def mem = params.snpeff_mem
    def java = params.java ?: 'java'
    // The SnpEff image is 5.2 (Java + snpEff + gzip). It cannot read a 4.3 database
    // and it has no bcftools, bgzip, or tabix. A staged jar runs on the image JRE.
    // BCFTOOLS_PACK_SNPEFF compresses the plain VCF.
    """
    gzip -dc ${vcf} > ${meta.id}.input.vcf

    if [[ ! -s snpEff.jar && ! -d snpeff_data ]]; then
        echo "SnpEff jar and data directory were not staged. Set --snpeff_jar and --snpeff_data." >&2
        exit 1
    fi

    # Nextflow stages these as symlinks. SnpEff hangs on a symlink -dataDir
    # and runs when given the canonical host path (Docker mounts that path).
    data_arg=()
    if [[ -d snpeff_data ]]; then
        data_arg=(-dataDir "\$(readlink -f snpeff_data)")
    fi

    if [[ -s snpEff.jar ]]; then
        ${java} -Xmx${mem} -jar "\$(readlink -f snpEff.jar)" \\
            "\${data_arg[@]}" -v \\
            -csvStats ${meta.id}_snpeff.csv \\
            ${db} ${meta.id}.input.vcf > ${meta.id}.snpeff.vcf
    else
        snpEff -Xmx${mem} "\${data_arg[@]}" -v \\
            -csvStats ${meta.id}_snpeff.csv \\
            ${db} ${meta.id}.input.vcf > ${meta.id}.snpeff.vcf
    fi
    rm -f ${meta.id}.input.vcf

    for g in snpEff_genes.txt ${meta.id}_snpeff.genes.txt genes.txt; do
        if [[ -f "\$g" && "\$g" != "${meta.id}_snpeff.genes.txt" ]]; then
            mv -f "\$g" ${meta.id}_snpeff.genes.txt || true
        fi
    done
    touch ${meta.id}_snpeff.genes.txt 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$({ if [[ -s snpEff.jar ]]; then ${java} -jar "\$(readlink -f snpEff.jar)" -version 2>&1 | head -1; else snpEff -version 2>&1 | head -1; fi; } || echo unknown)
    END_VERSIONS
    """
}

process BCFTOOLS_PACK_SNPEFF {
    tag "$meta.id"
    label 'process_low'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/annotation/snpeff" }, mode: params.publish_dir_mode,
        pattern: '*.{vcf.gz,tbi}'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.snpeff.vcf.gz"), path("${meta.id}.snpeff.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    !params.skip_snpeff && !params.skip_annotation

    script:
    """
    bcftools view -Oz -o ${meta.id}.snpeff.vcf.gz ${vcf}
    bcftools index -t ${meta.id}.snpeff.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """
}

process CLINICAL_REPORTS {
    tag "$meta.id"
    label 'process_medium'

    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/annotation/reports" }, mode: params.publish_dir_mode,
        pattern: '*.tsv'
    publishDir path: { "${params.outdir}/${meta.id}/${meta.lane}/annotation/logs" }, mode: params.publish_dir_mode,
        pattern: '*.log'

    input:
    tuple val(meta), path(snpeff_vcf), path(snpeff_tbi), path(bcf_vcf), path(bcf_tbi)
    path  clinvar   // optional NO_FILE
    path  dbsnp     // optional NO_FILE
    path  omim

    output:
    path "${meta.id}_clinical_report.tsv", emit: clinical
    path "${meta.id}_significant_findings.tsv", emit: significant
    path "${meta.id}_acmg_curation.tsv", emit: acmg
    path "${meta.id}_annotate.log", emit: log
    path "versions.yml", emit: versions

    when:
    !params.skip_annotation

    script:
    def clinvar_arg = clinvar.name != 'NO_FILE' ? "--clinvar-ref ${clinvar}" : ''
    def dbsnp_arg = dbsnp.name != 'NO_FILE' ? "--dbsnp-vcf ${dbsnp}" : ''
    """
    # Ensure report scripts can import clinical_disorder from bin/
    export PATH="\${PATH}"
    BINDIR=\$(dirname \$(command -v vcf_to_clinical_report.py || echo .))
    export PYTHONPATH="\${BINDIR}:\${PYTHONPATH:-}"

    {
      echo "Clinical reports for ${meta.id} lane=${meta.lane}"
      echo "SnpEff VCF: ${snpeff_vcf}"
      echo "bcftools annotated VCF: ${bcf_vcf}"
      echo "PYTHONPATH=\$PYTHONPATH"

      FINAL_VCF="${snpeff_vcf}"
      if [[ ! -s "\${FINAL_VCF}" ]]; then FINAL_VCF="${bcf_vcf}"; fi

      # Clinical report style (results-bcftools / results-gatk)
      python3 \${BINDIR}/vcf_to_clinical_report.py \\
          --vcf "\${FINAL_VCF}" --bcftools-vcf ${bcf_vcf} \\
          ${clinvar_arg} --bcftools bcftools \\
          --omim-table ${omim} \\
          --sample ${meta.id} --out ${meta.id}_clinical_report.tsv

      python3 \${BINDIR}/vcf_to_significant_report.py \\
          --vcf "\${FINAL_VCF}" --bcftools-vcf ${bcf_vcf} \\
          ${clinvar_arg} ${dbsnp_arg} --bcftools bcftools \\
          --omim-table ${omim} \\
          --sample ${meta.id} --out ${meta.id}_significant_findings.tsv

      python3 \${BINDIR}/vcf_to_acmg_template.py \\
          --vcf "\${FINAL_VCF}" --bcftools-vcf ${bcf_vcf} \\
          --sample ${meta.id} --out ${meta.id}_acmg_curation.tsv

      echo "DONE reports"
      ls -lh ${meta.id}_*.tsv
    } 2>&1 | tee ${meta.id}_annotate.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}

process COMPARE_PIPELINES {
    tag "$meta.id"
    label 'process_low'

    publishDir path: { "${params.outdir}/${meta.id}" }, mode: params.publish_dir_mode, pattern: 'pipeline_comparison.tsv'

    input:
    tuple val(meta), path(bcf_vcf), path(gatk_vcf)

    output:
    path "pipeline_comparison.tsv", emit: tsv
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p results-bcftools/variants results-gatk/variants \\
             results-bcftools/annotation/reports results-gatk/annotation/reports
    cp -f ${bcf_vcf} results-bcftools/variants/${meta.id}.filtered.vcf.gz
    cp -f ${gatk_vcf} results-gatk/variants/${meta.id}.pass.vcf.gz
    # compare script expects sample_dir with results-bcftools/ and results-gatk/
    compare_pipeline_results.py . ${meta.id}
    # script writes pipeline_comparison.tsv in cwd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
