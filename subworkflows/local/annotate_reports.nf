// Annotation + clinical reports matching results-*/annotation/reports/
include {
    BCFTOOLS_ANNOTATE_DBSNP
    BCFTOOLS_ANNOTATE_CLINVAR
    BCFTOOLS_ANNOTATE_GNOMAD
    BCFTOOLS_FINALIZE_ANNOT
    SNPEFF_ANNOTATE
    BCFTOOLS_PACK_SNPEFF
    CLINICAL_REPORTS
} from '../../modules/local/annotate'

workflow ANNOTATE_REPORTS {
    take:
    vcf        // [meta, vcf, tbi]
    dbsnp
    clinvar
    gnomad
    omim
    snpeff_jar
    snpeff_data

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

    // SnpEff writes a plain VCF. Pack it in the bcftools image (bgzip + tabix).
    SNPEFF_ANNOTATE(ch_annot, snpeff_jar, snpeff_data)
    ch_versions = ch_versions.mix(SNPEFF_ANNOTATE.out.versions)
    BCFTOOLS_PACK_SNPEFF(SNPEFF_ANNOTATE.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_PACK_SNPEFF.out.versions)

    // Prefer snpeff when available; always pass bcftools annotated
    ch_for_reports = ch_annot
        .map { m, v, t -> [m.id, m, v, t] }
        .join(
            BCFTOOLS_PACK_SNPEFF.out.vcf.map { m, v, t -> [m.id, v, t] },
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
    snpeff_vcf  = BCFTOOLS_PACK_SNPEFF.out.vcf
    annotated   = ch_annot
    versions    = ch_versions
}
