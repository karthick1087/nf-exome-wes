# Changelog

## 1.0.0

- Initial public release of **nf-exome-wes**
- Dual lanes: `results-bcftools` and `results-gatk`
- Options: `--pipeline`, `--threads`, `--min_qual`, `--min_depth`
- Clinical reports: clinical / significant findings / ACMG TSVs
- Reuses prebuilt BWA indexes via symlink staging
- Robust known-site / ClinVar indexing when `.tbi` is not auto-staged
- Profiles: local, lab, site, conda, mamba, docker, singularity, slurm, test
- `site` sets `$HOME/reference` paths only when those files exist
- OMIM genemap2 is not bundled; pass `--omim_table` if you have a licensed copy
- GitHub CI: Nextflow config + Python syntax + conda env
- Docker BWA image tag points at a published Quay build (`bwa` + `samtools`)
- Docker bcftools calling no longer runs `samtools` (that image has no samtools). Alignment flagstat stays with the BAM
- Docker SnpEff runs `--snpeff_jar` with `--snpeff_data` mounted into the image, then compresses the VCF in the bcftools image. The SnpEff 5.2 image cannot read a 4.3 database and has no `bgzip` or `tabix`
- Discovered SnpEff, dbSNP, and ClinVar paths are kept in local variables. Nextflow ignores a second assignment to `params`, which had staged an empty file and made the SnpEff image try to download `hg38`
- GATK `4.6.1.0` image already includes `bcftools`, `samtools`, `tabix`, and `bgzip`
- `./install.sh` installs Nextflow and the conda env on first install
- Every run checks the profile's tools before any analysis step
