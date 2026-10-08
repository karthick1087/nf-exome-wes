# Output layout

Standard per-sample folders:

- `results-bcftools/` — bcftools germline lane  
- `results-gatk/` — GATK clinical lane  

Root: `{params.outdir}/{sample_id}/`

## results-bcftools

| Path | Description |
|------|-------------|
| `qc/{id}_R1.clean.fastq.gz` | fastp trimmed R1 |
| `qc/{id}_R2.clean.fastq.gz` | fastp trimmed R2 |
| `qc/{id}_fastp.html` | fastp HTML |
| `qc/{id}_fastp.json` | fastp JSON |
| `align/{id}.sorted.bam[.bai]` | BWA-MEM BAM |
| `align/{id}.flagstat.txt` | samtools flagstat |
| `align/{id}.stats.txt` | samtools stats |
| `variants/{id}.raw.vcf.gz` | bcftools call |
| `variants/{id}.filtered.vcf.gz` | QUAL/DP filter |
| `variants/{id}.exome.vcf.gz` | optional BED subset |
| `variants/{id}_summary.txt` | text summary |
| `annotation/bcftools/{id}.dbsnp.vcf.gz` | rsIDs |
| `annotation/bcftools/{id}.clinvar.vcf.gz` | ClinVar INFO |
| `annotation/bcftools/{id}.fully_annotated.vcf.gz` | merged |
| `annotation/snpeff/{id}.snpeff.vcf.gz` | SnpEff ANN |
| `annotation/reports/{id}_clinical_report.tsv` | full clinical table |
| `annotation/reports/{id}_significant_findings.tsv` | prioritised hits |
| `annotation/reports/{id}_acmg_curation.tsv` | ACMG worksheet |
| `annotation/logs/{id}_annotate.log` | report log |

## results-gatk

Same `qc/` / `annotation/` idea, plus:

| Path | Description |
|------|-------------|
| `align/{id}.dedup.bam` | MarkDuplicates |
| `align/{id}.recal.bam` | after BQSR |
| `align/{id}.recal.table` | BQSR table |
| `align/{id}.dedup_metrics.txt` | duplicate metrics |
| `variants/{id}.gatk.raw.vcf.gz` | HaplotypeCaller |
| `variants/{id}.filtered.vcf.gz` | VariantFiltration |
| `variants/{id}.pass.vcf.gz` | PASS only → annotation |
| `variants/{id}_clinical_summary.txt` | summary |
| `qc/{id}.hs_metrics.txt` | CollectHsMetrics if BED |

## Comparison

`{sample}/pipeline_comparison.tsv` — columns  
`section`, `field`, `bcftools`, `gatk`, `notes`
