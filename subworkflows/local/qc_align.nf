// fastp → BWA-MEM → sorted BAM
include { FASTP   } from '../../modules/local/fastp'
include { BWA_MEM } from '../../modules/local/bwa_mem'

workflow QC_ALIGN {
    take:
    reads     // [meta, r1, r2]  meta.lane = results-bcftools | results-gatk
    fasta     // unused (kept for API compat); indices carry fasta
    index     // genome file collection
    fai

    main:
    ch_versions = Channel.empty()

    FASTP(reads)
    ch_versions = ch_versions.mix(FASTP.out.versions)

    BWA_MEM(
        FASTP.out.reads,
        index
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    emit:
    bam      = BWA_MEM.out.bam
    flagstat = BWA_MEM.out.flagstat
    stats    = BWA_MEM.out.stats
    html     = FASTP.out.html
    json     = FASTP.out.json
    versions = ch_versions
}
