# nf-exome-wes

**Nextflow DSL2 pipeline** for whole-exome sequencing (WES) germline analysis.

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

## Install

Java 11 or 17 is required. Reads and the reference FASTA stay on your disk.

```bash
git clone https://github.com/karthick1087/nf-exome-wes.git
cd nf-exome-wes
./install.sh
conda activate nf-exome-wes
```

`./install.sh` installs Nextflow when it is missing and creates the `nf-exome-wes` conda env from `environment.yml` the first time. Run it again only to confirm the install. It does not rebuild an env that already exists.

Use `-profile conda` or `-profile mamba` if you do not want to activate the env. Use `-profile docker` when Docker is installed instead of conda.

## Run

Edit `assets/samplesheet.csv` so `r1` and `r2` are your FASTQ paths. Edit `assets/params.example.yml` so `ref` is your reference FASTA, or pass those paths on the command line.

Both callers, tools on `PATH`:

```bash
nextflow run . -profile local \
  --input assets/samplesheet.csv \
  --ref /path/to/hg38.fasta \
  --pipeline both \
  --threads 16 \
  --min_qual 30 \
  --min_depth 10 \
  --outdir results
```

One sample, bcftools only:

```bash
nextflow run . -profile local \
  --sample_id SAMPLE001 \
  --sample_name SAMPLE001 \
  --r1 /path/to/SAMPLE001_R1.fastq.gz \
  --r2 /path/to/SAMPLE001_R2.fastq.gz \
  --ref /path/to/hg38.fasta \
  --pipeline bcftools
```

Conda, using a params file:

```bash
cp assets/params.example.yml myparams.yml
nextflow run . -profile conda -params-file myparams.yml
```

Docker:

```bash
nextflow run . -profile docker -params-file myparams.yml
```

From GitHub, without a local clone:

```bash
nextflow run karthick1087/nf-exome-wes -r main -profile docker \
  --input samplesheet.csv \
  --ref /path/to/hg38.fasta \
  --pipeline both
```

Continue a stopped run with the same command plus `-resume`.

Every run prints `Requirements OK (<profile>)` before alignment. If a tool is missing it stops with `Missing requirements:` and does no analysis.

| Profile | Must be available before the run |
|---------|----------------------------------|
| `local`, `lab`, `site` | `fastp`, `bwa`, `samtools`, `bcftools`, `bgzip`, `tabix`, `python3`. Also `gatk` when `--pipeline` is `gatk` or `both`. SnpEff unless you pass `--skip_snpeff` or `--skip_annotation`. |
| `conda`, `mamba` | `conda` or `mamba`. The first run builds `environment.yml`. |
| `docker` | `docker`. Images download on first use. SnpEff uses `--snpeff_jar` and `--snpeff_data` inside the image when those paths exist. The bcftools image then compresses that VCF. |
| `singularity` | `singularity` or `apptainer`. |
| `slurm` | Same tools as `local`, plus a SLURM cluster. |
| `test` | Same tools as `local`. Uses 2 CPUs and skips annotation. |

`site` sets 16 threads. Empty reference parameters use `$HOME/reference/...` only when that file exists. `lab` forwards `NXF_EXTRA_PATH` and `GATK_LOCAL_JAR` into each task.

## Test

These checks do not need reads or a reference:

```bash
nextflow run . --help
nextflow config -profile local
nextflow config -profile conda
nextflow config -profile docker
```

`--help` must print `nf-exome-wes`. `nextflow config` must exit 0.

To test a real start, use your FASTQs and FASTA. The run is valid when the log contains `Requirements OK` and then task names such as `FASTP`. Stop it with Ctrl+C if you only wanted the startup check, then repeat the same command with `-resume` to continue.

```bash
nextflow run . -profile local \
  --sample_id SAMPLE001 \
  --sample_name SAMPLE001 \
  --r1 /path/to/SAMPLE001_R1.fastq.gz \
  --r2 /path/to/SAMPLE001_R2.fastq.gz \
  --ref /path/to/hg38.fasta \
  --pipeline bcftools
```

Same start from the GitHub revision:

```bash
nextflow run karthick1087/nf-exome-wes -r main --help
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
├── install.sh
├── .github/workflows/ci.yml
├── LICENSE
└── README.md
```

## Documentation

- [docs/ANALYSIS.md](docs/ANALYSIS.md) — analysis workflow and QC
- [docs/output_layout.md](docs/output_layout.md) — published paths

## License

MIT — see [LICENSE](LICENSE).
