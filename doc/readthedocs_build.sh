#!/usr/bin/env bash
set -euo pipefail

# Run from the repository root.
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_dir}"

conda_env_name="mhm-docs"
output_dir="${READTHEDOCS_OUTPUT:-_readthedocs}"

# Create and activate the documentation environment.
conda env create --quiet \
  --name "${conda_env_name}" \
  --file doc/rtd-environment.yml

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${conda_env_name}"

# Supply version metadata and the compiler for doc/pre-proc.sh.
PROJECT_NUMBER="$(cat version.txt)"
export PROJECT_NUMBER

FC="$(command -v gfortran)"
export FC

"${FC}" --version
doxygen --version
dot -V
pdflatex --version

# Generate HTML and LaTeX, then build the PDF with system TeX Live.
doxygen doc/doxygen.config
pdf_log="doxygen_latex_make.log"
printf 'Building PDF (log: %s)...\n' "${pdf_log}"
if make -C latex \
    LATEX_CMD="pdflatex -interaction=nonstopmode -halt-on-error -file-line-error" \
    >"${pdf_log}" 2>&1; then
  printf 'PDF build complete.\n'
else
  printf 'PDF build failed. Last 80 lines of %s:\n' "${pdf_log}" >&2
  tail -n 80 "${pdf_log}" >&2
  exit 1
fi

# Stage HTML and PDF for Read the Docs.
mkdir -p "${output_dir}/html" "${output_dir}/pdf"
cp -a html/. "${output_dir}/html/"
cp latex/refman.pdf "${output_dir}/html/mhm_doc.pdf"
cp latex/refman.pdf "${output_dir}/pdf/mhm.pdf"
