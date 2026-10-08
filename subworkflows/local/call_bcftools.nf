// BCFtools germline calling lane → results-bcftools/
include { BCFTOOLS_MPILEUP_CALL } from '../../modules/local/bcftools_call'
include { BCFTOOLS_VIEW_EXOME   } from '../../modules/local/bcftools_call'

workflow CALL_BCFTOOLS {
    take:
    bam     // [meta, bam, bai] with meta.lane = 'results-bcftools'
    fasta
    fai
    bed     // path or NO_FILE

    main:
    ch_versions = Channel.empty()

    BCFTOOLS_MPILEUP_CALL(bam, fasta, fai, bed)
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP_CALL.out.versions)

    // optional exome subset (non-blocking for annotation which uses filtered)
    BCFTOOLS_VIEW_EXOME(BCFTOOLS_MPILEUP_CALL.out.filtered, bed)
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW_EXOME.out.versions)

    emit:
    vcf      = BCFTOOLS_MPILEUP_CALL.out.filtered  // [meta, vcf, tbi]
    raw      = BCFTOOLS_MPILEUP_CALL.out.raw
    summary  = BCFTOOLS_MPILEUP_CALL.out.summary
    versions = ch_versions
}
