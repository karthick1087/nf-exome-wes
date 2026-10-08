// Annotation + clinical reports matching results-*/annotation/reports/
include {
    BCFTOOLS_ANNOTATE_DBSNP
    BCFTOOLS_ANNOTATE_CLINVAR
    BCFTOOLS_ANNOTATE_GNOMAD
    BCFTOOLS_FINALIZE_ANNOT
    SNPEFF_ANNOTATE
    CLINICAL_REPORTS
} from '../../modules/local/annotate'

workflow ANNOTATE_REPORTS {
    take:
    vcf        // [meta, vcf, tbi]
    dbsnp
    clinvar
    gnomad
    omim

    main:
    ch_versions = Channel.empty()
    ch_work = vcf

    // Step 1 dbSNP (optional)
    if (params.dbsnp_vcf) {
        BCFTOOLS_ANNOTATE_DBSNP(ch_work, dbsnp)
        ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE_DBSNP.out.versions)
        ch_work = BCFTOOLS_ANNOTATE_DBSNP.out.vcf
    }

    // Step 2 ClinVar (optional)
    if (params.clinvar_vcf) {
        BCFTOOLS_ANNOTATE_CLINVAR(ch_work, clinvar)
        ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE_CLINVAR.out.versions)
        ch_work = BCFTOOLS_ANNOTATE_CLINVAR.out.vcf
    }

    // Step 3 gnomAD (optional)
    if (params.gnomad_vcf) {
        BCFTOOLS_ANNOTATE_GNOMAD(ch_work, gnomad)
        ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE_GNOMAD.out.versions)
        ch_work = BCFTOOLS_ANNOTATE_GNOMAD.out.vcf
    }

    BCFTOOLS_FINALIZE_ANNOT(ch_work)
    ch_versions = ch_versions.mix(BCFTOOLS_FINALIZE_ANNOT.out.versions)
    ch_annot = BCFTOOLS_FINALIZE_ANNOT.out.vcf

    // SnpEff
    SNPEFF_ANNOTATE(ch_annot)
    ch_versions = ch_versions.mix(SNPEFF_ANNOTATE.out.versions)

    // Join snpeff + fully_annotated for reports
    ch_snpeff = SNPEFF_ANNOTATE.out.vcf
        .ifEmpty { ch_annot.map { m, v, t -> [m, v, t] } }

    // Prefer snpeff when available; always pass bcftools annotated
    ch_for_reports = ch_annot
        .map { m, v, t -> [m.id, m, v, t] }
        .join(
            SNPEFF_ANNOTATE.out.vcf.map { m, v, t -> [m.id, v, t] },
            remainder: true
        )
        .map { id, m, bcf_v, bcf_t, snp_v, snp_t ->
            def sv = snp_v ?: bcf_v
            def st = snp_t ?: bcf_t
            [m, sv, st, bcf_v, bcf_t]
        }

    CLINICAL_REPORTS(ch_for_reports, clinvar, dbsnp, omim)
    ch_versions = ch_versions.mix(CLINICAL_REPORTS.out.versions)

    emit:
    clinical    = CLINICAL_REPORTS.out.clinical
    significant = CLINICAL_REPORTS.out.significant
    acmg        = CLINICAL_REPORTS.out.acmg
    snpeff_vcf  = SNPEFF_ANNOTATE.out.vcf
    annotated   = ch_annot
    versions    = ch_versions
}
