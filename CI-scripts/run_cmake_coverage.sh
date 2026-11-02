#!/usr/bin/env bash
set -euo pipefail

build_dir="$(cd "${1:-.}" && pwd)"
source_dir="$(cd "${2:-$(dirname "${BASH_SOURCE[0]}")/..}" && pwd)"

ctest --test-dir "${build_dir}" --output-on-failure

# Stage writable cases separately; domain inputs and reference results stay intact.
mkdir -p "${build_dir}/coverage-runs"
run_dir="$(mktemp -d "${build_dir}/coverage-runs/mhm.XXXXXX")"
cp -a "${source_dir}/check" "${run_dir}/check"
ln -s "${source_dir}/test_domain" "${run_dir}/test_domain"
ln -s "${source_dir}/test_domain_2" "${run_dir}/test_domain_2"
printf 'Logs and outputs: %s\n' "${run_dir}"
cd "${run_dir}/check"
python -u run_mhm_checks.py --exe "${build_dir}/mhm" --log_path "${run_dir}"
