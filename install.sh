#!/usr/bin/env bash
# First install for nf-exome-wes. Safe to re-run: existing tools and the
# conda env are left in place. Every nextflow run checks requirements again.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

if ! command -v java >/dev/null 2>&1; then
    echo "Java is required to run Nextflow. Install a JDK (11 or 17) and re-run ./install.sh"
    exit 1
fi

if ! command -v nextflow >/dev/null 2>&1; then
    echo "Installing Nextflow into \$HOME/.local/bin"
    mkdir -p "${HOME}/.local/bin"
    curl -fsSL https://get.nextflow.io | bash
    mv -f nextflow "${HOME}/.local/bin/nextflow"
    case ":${PATH}:" in
        *":${HOME}/.local/bin:"*) ;;
        *) echo "Add ${HOME}/.local/bin to PATH, then open a new shell." ;;
    esac
else
    echo "Nextflow: $(command -v nextflow)"
fi

if command -v mamba >/dev/null 2>&1; then
    CONDA_BIN="$(command -v mamba)"
elif command -v micromamba >/dev/null 2>&1; then
    CONDA_BIN="$(command -v micromamba)"
elif command -v conda >/dev/null 2>&1; then
    CONDA_BIN="$(command -v conda)"
else
    CONDA_BIN=""
fi

if [[ -n "${CONDA_BIN}" ]]; then
    if "${CONDA_BIN}" env list | awk '{print $1}' | grep -qx 'nf-exome-wes'; then
        echo "Conda env nf-exome-wes already exists"
    else
        echo "Creating conda env nf-exome-wes from environment.yml"
        "${CONDA_BIN}" env create -f environment.yml
    fi
    echo "Activate it with: conda activate nf-exome-wes"
    echo "Or skip activation and run: nextflow run . -profile conda"
else
    echo "conda/mamba not found."
    echo "Install Miniforge, re-run ./install.sh, or use -profile docker (Docker must be installed)."
fi

GATK_WRAP="${HOME}/reference/tools/gatk-4.6.1.0/gatk"
if [[ -x "${GATK_WRAP}" ]] && ! command -v gatk >/dev/null 2>&1; then
    mkdir -p "${HOME}/.local/bin"
    ln -sfn "${GATK_WRAP}" "${HOME}/.local/bin/gatk"
    echo "Linked gatk into ${HOME}/.local/bin"
    echo "Export GATK_LOCAL_JAR=${HOME}/reference/tools/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar"
fi

echo "Install step finished. The next nextflow run checks these tools again before any work."
