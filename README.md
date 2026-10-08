# nf-exome-wes

**GitHub-ready Nextflow DSL2 pipeline** for whole-exome sequencing (WES) germline analysis.

Runs **bcftools** and/or **GATK** clinical lanes and writes a standard per-sample layout:

```
results/{sample}/
├── results-bcftools/   qc → align → variants → annotation/reports
├── results-gatk/       qc → align (dedup/BQSR) → variants → annotation/reports
└── pipeline_comparison.tsv          # when --pipeline both
```

Reports (per lane):

- `{id}_clinical_report.tsv`
- `{id}_significant_findings.tsv`
- `{id}_acmg_curation.tsv`

## Features

- **Modes:** `--pipeline bcftools | gatk | both`
- **Filters:** `--threads`, `--min_qual`, `--min_depth`
- **Annotation:** dbSNP, ClinVar, optional gnomAD, SnpEff, clinical TSVs
- **Portable:** conda / Docker / Singularity / SLURM
- **Small repo** — no FASTQs, BAMs, or genomes committed

## Requirements

| Component | Notes |
|-----------|--------|
| [Nextflow](https://www.nextflow.io/) ≥ 23.04 | `curl -s https://get.nextflow.io \| bash` |
| One of | host tools, **conda**, **Docker**, or **Singularity** |
| Data | R1/R2 FASTQ + reference FASTA (+ optional BED, ClinVar, dbSNP, SnpEff DB) |

```bash
mamba env create -f environment.yml
mamba activate nf-exome-wes
```

## Quick start

### 1. Clone

```bash
git clone https://github.com/karthick1087/nf-exome-wes.git
cd nf-exome-wes
```

### 2. Samplesheet

```csv
sample_id,sample_name,r1,r2
SAMPLE001,case1,/data/SAMPLE001_R1.fastq.gz,/data/SAMPLE001_R2.fastq.gz
```

### 3. Run full analysis (both callers)

```bash
nextflow run main.nf -profile docker \
  --input assets/samplesheet.csv \
  --ref /refs/hg38.fasta \
  --pipeline both \
  --threads 16 \
  --min_qual 30 \
  --min_depth 10 \
  --outdir results
```

Or use a params file:

```bash
cp assets/params.example.yml myparams.yml
# edit paths...
nextflow run main.nf -profile conda -params-file myparams.yml
```

### 4. Single sample without CSV

```bash
nextflow run main.nf -profile local \
  --sample_id SAMPLE001 --sample_name SAMPLE001 \
  --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz \
  --ref hg38.fasta --pipeline bcftools
```

### 5. Resume

```bash
nextflow run main.nf -profile docker -params-file myparams.yml -resume
```

## How reference and other files are supplied

All large inputs are **paths you pass in** (CLI or params file). Nothing is baked into the repo.

| Input | How to set |
|-------|------------|
| FASTQ R1/R2 | Samplesheet columns `r1`,`r2` **or** `--r1` / `--r2` |
| Reference FASTA | `--ref /path/hg38.fasta` (required) |
| BWA / fai / dict | Auto-reused if next to the FASTA (`*.bwt`, `*.fai`, `*.dict`) |
| Exome BED | `--exome_bed /path/targets.bed` (optional) |
| dbSNP / ClinVar / gnomAD | `--dbsnp_vcf` `--clinvar_vcf` `--gnomad_vcf` |
| SnpEff | `--snpeff_jar` `--snpeff_data` `--snpeff_db` |
| BQSR known sites | `--dbsnp_bqsr` `--mills_indels` or `--known_sites_dir` |
| OMIM genemap2 | Optional `--omim_table`, or `$HOME/reference/clinical_wes/use_omim_table.txt` when that file exists. Not bundled. |

See [docs/ANALYSIS.md](docs/ANALYSIS.md) for a full checklist.

## Parameters

| Param | Default | Description |
|-------|---------|-------------|
| `--pipeline` | `bcftools` | `bcftools` \| `gatk` \| `both` |
| `--input` | — | Samplesheet CSV |
| `--sample_id` / `--r1` / `--r2` | — | Single-sample mode |
| `--ref` | — | Reference FASTA (**required**) |
| `--exome_bed` | — | Capture regions BED |
| `--threads` | `16` | CPUs for heavy steps |
| `--min_qual` | `30` | bcftools QUAL filter |
| `--min_depth` | `10` | bcftools FORMAT/DP filter |
| `--outdir` | `results` | Publish root |
| `--dbsnp_vcf` / `--clinvar_vcf` / `--gnomad_vcf` | — | Annotation DBs |
| `--snpeff_jar` / `--snpeff_data` / `--snpeff_db` | `hg38` | SnpEff |
| `--dbsnp_bqsr` / `--mills_indels` | — | GATK BQSR |
| `--omim_table` | — | Optional OMIM genemap2 TSV |
| `--skip_annotation` | `false` | Stop after VCF |
| `--skip_snpeff` | `false` | Skip SnpEff only |

Help: `nextflow run main.nf --help`

## Profiles

| Profile | Use |
|---------|-----|
| `local` | Tools on PATH |
| `lab` | Optional `NXF_EXTRA_PATH` / `GATK_LOCAL_JAR` env |
| `site` | 16 threads. Empty reference params use `$HOME/reference/...` only when that file exists |
| `conda` / `mamba` | `environment.yml` |
| `docker` / `singularity` | Containers |
| `slurm` | Cluster |
| `test` | CI |

Example: `-profile local,site` or `-profile docker`

## Repository layout

```
nf-exome-wes/
├── main.nf
├── nextflow.config
├── modules/local/
├── subworkflows/local/
├── bin/                    # clinical report scripts
├── assets/                 # samplesheets, params example, resources
├── conf/
├── docs/
├── environment.yml
├── .github/workflows/ci.yml
├── LICENSE
└── README.md
```

## Documentation

- [docs/ANALYSIS.md](docs/ANALYSIS.md) — analysis workflow & QC  
- [docs/output_layout.md](docs/output_layout.md) — published paths  


```bash
nextflow run karthick1087/nf-exome-wes -r main -profile docker \
  --input samplesheet.csv --ref hg38.fasta --pipeline both
```

## License

MIT — see [LICENSE](LICENSE).
