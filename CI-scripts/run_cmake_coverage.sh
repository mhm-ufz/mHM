#!/usr/bin/env bash
set -euo pipefail
build_dir=$1
source_dir=$2
ctest --test-dir "${build_dir}" --output-on-failure
if [[ -x "${build_dir}/app/mhm" ]]; then
  bash "${source_dir}/CI-scripts/run_v6_namelists.sh" --log-dir "${build_dir}/coverage-runs" "${build_dir}/app/mhm"
fi
