# Analysis guide — nf-exome-wes

End-to-end germline WES analysis with two calling lanes and clinical-style reports.

## What is analysed

| Stage | Tools | Purpose |
|-------|--------|---------|
| QC / trim | fastp | Adapters, quality, pair correction |
| Align | BWA-MEM + samtools | Illumina PE → sorted BAM + flagstat/stats |
| Call (bcftools) | mpileup + call + filter | Germline VCF (QUAL / DP) |
| Call (GATK) | MarkDuplicates → BQSR → HaplotypeCaller → VariantFiltration | Clinical best-practice VCF |
| Annotate | bcftools (dbSNP, ClinVar, gnomAD) + SnpEff | rsIDs, clinical significance, effects |
| Report | Python suite | Clinical / significant / ACMG TSVs |
| Compare | compare_pipeline_results.py | Side-by-side bcftools vs GATK (`--pipeline both`) |

## Modes

```bash
--pipeline bcftools   # results-bcftools/ only
--pipeline gatk       # results-gatk/ only
--pipeline both       # both + pipeline_comparison.tsv
```

## Inputs you must provide

1. **Paired FASTQ** (R1/R2) — samplesheet or CLI  
2. **Reference FASTA** (`--ref`, e.g. GRCh38/hg38) with or without prebuilt BWA index  
3. Optional: **exome BED**, **ClinVar**, **dbSNP**, **SnpEff DB**, **BQSR known sites**, **OMIM genemap2** (`--omim_table`; not included in the repo)

Paths are always **your** filesystem (or container mounts). The pipeline does not download genomes.

### Examples

Samplesheet:

```csv
sample_id,sample_name,r1,r2
SAMPLE001,case_one,/data/SAMPLE001_R1.fastq.gz,/data/SAMPLE001_R2.fastq.gz
```

Params file (`assets/params.example.yml`):

```yaml
pipeline: both
input: assets/samplesheet.csv
ref: /refs/hg38.fasta
exome_bed: /refs/exome_targets.bed
dbsnp_vcf: /refs/dbsnp.vcf.gz
clinvar_vcf: /refs/clinvar.vcf.gz
snpeff_db: hg38
```

## Output tree (per sample)

```
results/{sample_id}/
├── results-bcftools/
│   ├── qc/
│   ├── align/
│   ├── variants/
│   └── annotation/{bcftools,snpeff,reports,logs}/
├── results-gatk/
│   ├── qc/
│   ├── align/
│   ├── variants/
│   └── annotation/{bcftools,snpeff,reports,logs}/
└── pipeline_comparison.tsv    # only if --pipeline both
```

## Report columns (clinical)

- `sample`, `gene`, `transcript`, `variant_type`, `hgvs_c`, `hgvs_p`
- `chrom`, `pos`, `rsid`, `ref`, `alt`, `zygosity`, `depth`, `vaf_pct`
- `impact`, `clinvar_sig`, `clinvar_disease`, `disorder`, `inheritance`
- `gnomad_af`, `classification`

**Significant findings** — prioritised pathogenic / high-impact subset.  
**ACMG template** — criteria columns for manual curation.

## Typical full analysis

```bash
nextflow run main.nf -profile docker \
  --input assets/samplesheet.csv \
  --ref /refs/hg38.fasta \
  --exome_bed /refs/exome.bed \
  --pipeline both \
  --threads 16 \
  --min_qual 30 \
  --min_depth 10 \
  --dbsnp_vcf /refs/dbsnp.vcf.gz \
  --clinvar_vcf /refs/clinvar.vcf.gz \
  --snpeff_jar /opt/snpEff/snpEff.jar \
  --snpeff_data /opt/snpEff/data \
  --snpeff_db hg38 \
  --dbsnp_bqsr /refs/dbsnp.vcf.gz \
  --outdir results
```

Resume:

```bash
nextflow run main.nf -profile docker -params-file myparams.yml -resume
```

Calling only:

```bash
nextflow run main.nf -profile conda \
  --input assets/samplesheet.csv --ref hg38.fasta \
  --pipeline bcftools --skip_annotation true
```

## Interpreting comparison TSV

| section | meaning |
|---------|---------|
| `summary` / `variant_count` | Total variants per lane |
| `target` / gene loci | Optional key loci (configurable in the compare script) |

Labels: `MATCH`, `BCF_ONLY`, `GATK_ONLY`, `NOT_FOUND`.

## Quality checks after a run

```bash
bcftools view -H results/*/results-bcftools/variants/*.filtered.vcf.gz | wc -l
bcftools view -H results/*/results-gatk/variants/*.pass.vcf.gz | wc -l
wc -l results/*/results-*/annotation/reports/*_clinical_report.tsv
wc -l results/*/results-*/annotation/reports/*_significant_findings.tsv
head results/*/results-bcftools/align/*.flagstat.txt
```

## Profiles

| Profile | When |
|---------|------|
| `local` | Tools already on PATH |
| `lab` | `export NXF_EXTRA_PATH=...` for extra tool dirs |
| `site` | Optional defaults under `$HOME/reference/` |
| `conda` / `mamba` | `environment.yml` |
| `docker` / `singularity` | Containers |
| `slurm` | HPC |
| `test` | CI |

## Resource notes

- Prebuilt BWA index next to `--ref` is reused (no multi-hour reindex)  
- Known-site / ClinVar VCFs should have `.tbi` beside them (pipeline can copy or build if missing)  
- WES PE full both-lane runs typically take hours on 16 cores  
