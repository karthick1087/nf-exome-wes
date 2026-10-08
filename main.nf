#!/usr/bin/env nextflow
/*
 * =============================================================================
 * nf-exome-wes  —  Whole-Exome Sequencing germline pipeline (Nextflow DSL2)
 *
 * Output layout (per sample):
 *   results/{sample}/results-bcftools/{qc,align,variants,annotation/reports}
 *   results/{sample}/results-gatk/{qc,align,variants,annotation/reports}
 *   results/{sample}/pipeline_comparison.tsv
 *
 * Reports:
 *   {id}_clinical_report.tsv
 *   {id}_significant_findings.tsv
 *   {id}_acmg_curation.tsv
 * =============================================================================
 */

nextflow.enable.dsl = 2

include { PREPARE_GENOME    } from './subworkflows/local/prepare_genome'
include { QC_ALIGN          } from './subworkflows/local/qc_align'
include { CALL_BCFTOOLS     } from './subworkflows/local/call_bcftools'
include { CALL_GATK         } from './subworkflows/local/call_gatk'
include { ANNOTATE_REPORTS  } from './subworkflows/local/annotate_reports'
include { COMPARE_PIPELINES } from './modules/local/annotate'

def helpMessage() {
    return """
    ================================================================================
    nf-exome-wes  v${workflow.manifest.version}
    ================================================================================

    Usage:
      nextflow run main.nf -profile local,site \\
          --input samplesheet.csv --ref /path/hg38.fasta --pipeline both

      nextflow run main.nf -profile docker \\
          --sample_id SAMPLE001 --r1 R1.fq.gz --r2 R2.fq.gz --ref hg38.fasta

    Input:
      --input           CSV: sample_id,sample_name,r1,r2
      --sample_id --r1 --r2 --sample_name

    Mode (matches GUI):
      --pipeline        bcftools | gatk | both   [default: bcftools]

    Reference:
      --ref             hg38 FASTA (required)
      --exome_bed       capture BED

    Filters / compute:
      --threads 16 --min_qual 30 --min_depth 10 --outdir results

    GATK:
      --clinical_res --dbsnp_bqsr --mills_indels --known_sites_dir

    Annotation:
      --dbsnp_vcf --clinvar_vcf --gnomad_vcf
      --snpeff_jar --snpeff_data --snpeff_db hg38
      --omim_table      optional OMIM genemap2 TSV (not bundled)
      --skip_annotation --skip_snpeff

    Output layout (mirrors results-bcftools & results-gatk):
      {outdir}/{sample}/results-bcftools/{qc,align,variants,annotation/reports}
      {outdir}/{sample}/results-gatk/{qc,align,variants,annotation/reports}
      {outdir}/{sample}/pipeline_comparison.tsv

    Profiles: local lab site conda mamba docker singularity slurm test
    """.stripIndent()
}

def buildSampleChannel() {
    if (params.input) {
        return Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def id   = row.sample_id?.toString()?.trim()
                def name = (row.sample_name ?: id)?.toString()?.trim()
                def r1   = file(row.r1?.toString()?.trim(), checkIfExists: true)
                def r2   = file(row.r2?.toString()?.trim(), checkIfExists: true)
                if (!id) error "samplesheet row missing sample_id: ${row}"
                [ [id: id, name: name], r1, r2 ]
            }
    }
    if (params.r1 && params.r2 && params.sample_id) {
        def id   = params.sample_id.toString()
        def name = (params.sample_name ?: id).toString()
        return Channel.of([
            [id: id, name: name],
            file(params.r1, checkIfExists: true),
            file(params.r2, checkIfExists: true)
        ])
    }
    error "Provide --input samplesheet.csv OR --sample_id + --r1 + --r2"
}

def noFile() {
    def f = file("${projectDir}/assets/NO_FILE")
    if (!f.exists()) {
        f.text = ''
    }
    return f
}

def firstExisting(List candidates) {
    candidates.find { it && file(it).exists() }
}

// A CLI path must exist. An empty param may use the first conventional file that is present.
def optionalFile(explicit, List fallbacks) {
    if (explicit) {
        return file(explicit, checkIfExists: true)
    }
    def hit = firstExisting(fallbacks)
    return hit ? file(hit) : noFile()
}

workflow {

    if (params.help) {
        log.info helpMessage()
        return
    }

    log.info """
    ================================================================================
     nf-exome-wes  ·  pipeline=${params.pipeline}  threads=${params.threads}
     outdir=${params.outdir}  min_qual=${params.min_qual}  min_depth=${params.min_depth}
    ================================================================================
    """.stripIndent()

    if (!params.ref) {
        error "Required: --ref /path/to/hg38.fasta\n" + helpMessage()
    }

    noFile()

    def home = System.getenv('HOME') ?: ''
    def clinical = home ? "${home}/reference/clinical_wes" : ''
    def known = clinical ? "${clinical}/known_sites" : ''

    // Use absolute host path so PREPARE_REF can symlink prebuilt BWA indices
    def refAbs = file(params.ref, checkIfExists: true).toAbsolutePath().toString()
    ch_bed = optionalFile(params.exome_bed, [
        clinical ? "${clinical}/xgen_exome_v2_hg38_nochr.bed" : null
    ])
    ch_dbsnp_annot = optionalFile(params.dbsnp_vcf, [
        known ? "${known}/Homo_sapiens_assembly38.dbsnp138.nochr.vcf.gz" : null
    ])
    ch_clinvar = optionalFile(params.clinvar_vcf, [
        home ? "${home}/reference/annovar/humandb/clinvar.vcf.gz" : null
    ])
    ch_gnomad = params.gnomad_vcf ? file(params.gnomad_vcf, checkIfExists: true) : noFile()

    // OMIM genemap2 is optional. A missing table leaves disorder lookup to ClinVar.
    def omimDefault = "${projectDir}/assets/resources/use_omim_table.txt"
    def omimPath = params.omim_table ?: omimDefault
    ch_omim = file(omimPath).exists() ? file(omimPath) : noFile()

    if (!params.snpeff_jar) {
        def jar = firstExisting([
            home ? "${home}/anaconda3/share/snpeff-4.3.1t-0/snpEff.jar" : null,
            home ? "${home}/reference/snpeff/snpEff.jar" : null,
        ])
        if (jar) params.snpeff_jar = jar
    }
    if (!params.snpeff_data) {
        def dataDir = firstExisting([
            home ? "${home}/anaconda3/share/snpeff-4.3.1t-0/data" : null,
            home ? "${home}/reference/snpeff/data" : null,
        ])
        if (dataDir) params.snpeff_data = dataDir
    }

    def dbsnp_bqsr_path = params.dbsnp_bqsr
    if (!dbsnp_bqsr_path && params.known_sites_dir) {
        def cand = file("${params.known_sites_dir}/Homo_sapiens_assembly38.dbsnp138.nochr.vcf.gz")
        if (cand.exists()) {
            dbsnp_bqsr_path = cand.toString()
        }
    }
    if (!dbsnp_bqsr_path && known) {
        def cand = file("${known}/Homo_sapiens_assembly38.dbsnp138.nochr.vcf.gz")
        if (cand.exists()) {
            dbsnp_bqsr_path = cand.toString()
        }
    }
    ch_dbsnp_bqsr = dbsnp_bqsr_path ? file(dbsnp_bqsr_path, checkIfExists: true) : noFile()
    ch_mills = params.mills_indels ? file(params.mills_indels, checkIfExists: true) : noFile()

    PREPARE_GENOME(Channel.value(refAbs))
    ch_versions = PREPARE_GENOME.out.versions

    ch_samples = buildSampleChannel()

    def run_bcf  = params.pipeline in ['bcftools', 'both']
    def run_gatk = params.pipeline in ['gatk', 'both']

    ch_reads = Channel.empty()
    if (run_bcf) {
        ch_reads = ch_reads.mix(
            ch_samples.map { meta, r1, r2 -> [ meta + [lane: 'results-bcftools'], r1, r2 ] }
        )
    }
    if (run_gatk) {
        ch_reads = ch_reads.mix(
            ch_samples.map { meta, r1, r2 -> [ meta + [lane: 'results-gatk'], r1, r2 ] }
        )
    }

    QC_ALIGN(
        ch_reads,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.index,
        PREPARE_GENOME.out.fai
    )
    ch_versions = ch_versions.mix(QC_ALIGN.out.versions)

    ch_bam_bcf  = QC_ALIGN.out.bam.filter { meta, bam, bai -> meta.lane == 'results-bcftools' }
    ch_bam_gatk = QC_ALIGN.out.bam.filter { meta, bam, bai -> meta.lane == 'results-gatk' }

    ch_bcf_vcf  = Channel.empty()
    ch_gatk_vcf = Channel.empty()
    ch_annot_in = Channel.empty()

    if (run_bcf) {
        CALL_BCFTOOLS(
            ch_bam_bcf,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.fai,
            ch_bed
        )
        ch_versions = ch_versions.mix(CALL_BCFTOOLS.out.versions)
        ch_bcf_vcf  = CALL_BCFTOOLS.out.vcf
        ch_annot_in = ch_annot_in.mix(CALL_BCFTOOLS.out.vcf)
    }

    if (run_gatk) {
        CALL_GATK(
            ch_bam_gatk,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.fai,
            PREPARE_GENOME.out.dict,
            ch_bed,
            ch_dbsnp_bqsr,
            ch_mills
        )
        ch_versions = ch_versions.mix(CALL_GATK.out.versions)
        ch_gatk_vcf = CALL_GATK.out.vcf
        ch_annot_in = ch_annot_in.mix(CALL_GATK.out.vcf)
    }

    if (!params.skip_annotation) {
        ANNOTATE_REPORTS(
            ch_annot_in,
            ch_dbsnp_annot,
            ch_clinvar,
            ch_gnomad,
            ch_omim
        )
        ch_versions = ch_versions.mix(ANNOTATE_REPORTS.out.versions)
    }

    if (params.pipeline == 'both') {
        ch_cmp = ch_bcf_vcf
            .map { m, v, t -> [ m.id, m, v ] }
            .join( ch_gatk_vcf.map { m, v, t -> [ m.id, v ] } )
            .map { id, m, bcf_v, gatk_v -> [ m, bcf_v, gatk_v ] }

        COMPARE_PIPELINES(ch_cmp)
    }

    ch_versions
        .unique()
        .collectFile(
            name: 'pipeline_software_versions.yml',
            storeDir: "${params.outdir}/pipeline_info",
            newLine: true
        )
}


