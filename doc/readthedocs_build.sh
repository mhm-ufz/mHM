#!/usr/bin/env bash
set -euo pipefail

# Run from the repository root.
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_dir}"

conda_env_name="mhm-docs"
output_dir="${READTHEDOCS_OUTPUT:-_readthedocs}"

conda env create --quiet \
  --name "${conda_env_name}" \
  --file doc/rtd-environment.yml

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${conda_env_name}"

# Select and export the compiler for doc/pre-proc.sh.
FC="$(command -v gfortran)"
export FC
"${FC}" --version

PROJECT_NUMBER="$(cat version.txt)"
export PROJECT_NUMBER

# Compiler activation supplies FC for doc/pre-proc.sh.
"${FC:?Conda activation did not configure a Fortran compiler}" --version
doxygen --version
dot -V

# Generate HTML and LaTeX.
doxygen doc/doxygen.config

# Compile the PDF from the directory expected by the LaTeX header.
(
  cd latex
  tectonic -r 10 -Z continue-on-error refman.tex
)

# Stage the documentation for Read the Docs.
mkdir -p "${output_dir}/html" "${output_dir}/pdf"
cp -a html/. "${output_dir}/html/"
cp latex/refman.pdf "${output_dir}/html/mhm_doc.pdf"
cp latex/refman.pdf "${output_dir}/pdf/mhm.pdf"
