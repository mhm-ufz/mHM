#!/usr/bin/env bash
set -euo pipefail

# Run from the repository root. With no argument, run all build phases.
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_dir}"

phase="${1:-all}"
conda_env_name="${MHM_DOCS_ENV:-mhm-docs}"
output_dir="${READTHEDOCS_OUTPUT:-_readthedocs}"
case "${phase}" in
  all|environment|doxygen|pdf|publish) ;;
  *) printf 'Unknown build phase: %s\n' "${phase}" >&2; exit 2 ;;
esac

# Separate RTD commands expose the current phase even if a build times out.
if [[ "${phase}" == all || "${phase}" == environment ]]; then
  printf 'Creating documentation environment...\n'
  conda env create --quiet \
    --name "${conda_env_name}" \
    --file doc/rtd-environment.yml
fi

if [[ "${phase}" == all || "${phase}" == doxygen || "${phase}" == pdf ]]; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "${conda_env_name}"
fi

if [[ "${phase}" == all || "${phase}" == doxygen ]]; then
  # Supply version metadata and the compiler for doc/pre-proc.sh.
  PROJECT_NUMBER="$(cat version.txt)"
  export PROJECT_NUMBER
  FC="$(command -v gfortran)"
  export FC

  "${FC}" --version
  doxygen --version
  dot -V
  printf 'Generating HTML and LaTeX...\n'
  doxygen doc/doxygen-rtd.config
  # Doxygen 1.18's external-search script also runs on ordinary pages, which
  # have no searchresults element. Guard its empty-query message there.
  sed -i 's/^    results.innerHTML =/    if (results) results.innerHTML =/' html/search/search.js
fi

if [[ "${phase}" == all || "${phase}" == pdf ]]; then
  # Keep detailed compiler output out of the RTD command log.
  pdflatex --version
  pdf_log="doxygen_latex_make.log"
  pdf_started=${SECONDS}
  printf 'Building PDF (log: %s)...\n' "${pdf_log}"
  if latexmk -norc -C -bibtex -cd latex/refman.tex >"${pdf_log}" 2>&1 &&
      latexmk -norc -pdf -cd \
      -interaction=nonstopmode -halt-on-error -file-line-error latex/refman.tex \
      >>"${pdf_log}" 2>&1; then
    printf 'PDF build complete in %s seconds.\n' "$((SECONDS - pdf_started))"
  else
    printf 'PDF build failed. Last 80 lines of %s:\n' "${pdf_log}" >&2
    tail -n 80 "${pdf_log}" >&2
    exit 1
  fi
fi

if [[ "${phase}" == all || "${phase}" == publish ]]; then
  mkdir -p "${output_dir}/html" "${output_dir}/pdf"
  cp -a html/. "${output_dir}/html/"
  cp latex/refman.pdf "${output_dir}/pdf/mhm.pdf"
  printf 'HTML and PDF staged for Read the Docs.\n'
fi
