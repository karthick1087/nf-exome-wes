// GATK clinical lane → results-gatk/
include { GATK_MARKDUPLICATES    } from '../../modules/local/gatk_markduplicates'
include { GATK_BASERECALIBRATOR  } from '../../modules/local/gatk_bqsr'
include { GATK_APPLYBQSR         } from '../../modules/local/gatk_bqsr'
include { GATK_HAPLOTYPECALLER   } from '../../modules/local/gatk_haplotypecaller'
include { GATK_VARIANTFILTRATION } from '../../modules/local/gatk_haplotypecaller'
include { GATK_COLLECTHSMETRICS  } from '../../modules/local/gatk_haplotypecaller'

workflow CALL_GATK {
    take:
    bam
    fasta
    fai
    dict
    bed
    dbsnp_bqsr
    mills

    main:
    ch_versions = Channel.empty()

    GATK_MARKDUPLICATES(bam)
    ch_versions = ch_versions.mix(GATK_MARKDUPLICATES.out.versions)
    ch_dedup = GATK_MARKDUPLICATES.out.bam

    // Always attempt BQSR process; module can no-op if dbsnp is NO_FILE via when:
    GATK_BASERECALIBRATOR(ch_dedup, fasta, fai, dict, dbsnp_bqsr, mills)
    ch_versions = ch_versions.mix(GATK_BASERECALIBRATOR.out.versions)

    // Join recal table when present; otherwise pass dedup BAM through
    ch_with_table = ch_dedup
        .map { m, b, i -> [m.id, m, b, i] }
        .join(
            GATK_BASERECALIBRATOR.out.table.map { m, t -> [m.id, t] },
            remainder: true
        )

    ch_apply = ch_with_table
        .filter { id, m, b, i, t -> t != null }
        .map { id, m, b, i, t -> [m, b, i, t] }

    ch_skip_bqsr = ch_with_table
        .filter { id, m, b, i, t -> t == null }
        .map { id, m, b, i, t -> [m, b, i] }

    GATK_APPLYBQSR(ch_apply, fasta, fai, dict)
    ch_versions = ch_versions.mix(GATK_APPLYBQSR.out.versions)

    ch_for_hc = GATK_APPLYBQSR.out.bam.mix(ch_skip_bqsr)

    GATK_HAPLOTYPECALLER(ch_for_hc, fasta, fai, dict, bed)
    ch_versions = ch_versions.mix(GATK_HAPLOTYPECALLER.out.versions)

    GATK_VARIANTFILTRATION(GATK_HAPLOTYPECALLER.out.vcf, fasta, fai, dict)
    ch_versions = ch_versions.mix(GATK_VARIANTFILTRATION.out.versions)

    GATK_COLLECTHSMETRICS(ch_for_hc, fasta, fai, dict, bed)
    ch_versions = ch_versions.mix(GATK_COLLECTHSMETRICS.out.versions)

    emit:
    vcf       = GATK_VARIANTFILTRATION.out.pass
    filtered  = GATK_VARIANTFILTRATION.out.filtered
    summary   = GATK_VARIANTFILTRATION.out.summary
    bam_recal = ch_for_hc
    versions  = ch_versions
}
