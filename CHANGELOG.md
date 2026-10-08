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
- `./install.sh` installs Nextflow and the conda env on first install
- Every run checks the profile's tools before any analysis step
