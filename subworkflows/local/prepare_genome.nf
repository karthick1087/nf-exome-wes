// Prepare reference: reuse host indices via PREPARE_REF
include { PREPARE_REF } from '../../modules/local/prepare_ref'

workflow PREPARE_GENOME {
    take:
    fasta_path   // val string: absolute host path to fasta

    main:
    ch_versions = Channel.empty()

    PREPARE_REF(fasta_path)
    ch_versions = ch_versions.mix(PREPARE_REF.out.versions)

    // Split staged genome files into fasta / fai / dict / full set
    ch_files = PREPARE_REF.out.genome_files.collect()

    ch_fasta = ch_files.map { files ->
        def list = files instanceof List ? files : [files]
        def f = list.find { it.name.endsWith('.fasta') || it.name.endsWith('.fa') }
        if (!f) error "No FASTA in prepared genome files: ${list}"
        f
    }

    ch_fai = ch_files.map { files ->
        def list = files instanceof List ? files : [files]
        def f = list.find { it.name.endsWith('.fai') }
        if (!f) error "No .fai in prepared genome files"
        f
    }

    ch_dict = ch_files.map { files ->
        def list = files instanceof List ? files : [files]
        // prefer hg38.dict over hg38.fasta.dict
        def f = list.find { it.name.endsWith('.dict') && !it.name.endsWith('.fasta.dict') && !it.name.endsWith('.fa.dict') }
        if (!f) f = list.find { it.name.endsWith('.dict') }
        if (!f) error "No .dict in prepared genome files"
        f
    }

    emit:
    fasta    = ch_fasta
    fai      = ch_fai
    index    = ch_files
    dict     = ch_dict
    versions = ch_versions
}
