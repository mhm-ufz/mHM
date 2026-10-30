#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Run curated v6 smoke namelists with one or more mHM v6 driver executables.

Usage:
  run_v6_namelists.sh [--log-dir DIR] [--threads N] EXE [EXE ...]

Options:
  --log-dir DIR  Directory for per-run logs. Default: logs_v6
  --threads N    Set OMP_NUM_THREADS=N for each run
EOF
}

log_dir="logs_v6"
threads=0
executables=()

while (($# > 0)); do
  case "$1" in
    --log-dir)
      log_dir="$2"
      shift 2
      ;;
    --threads)
      threads="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      executables+=("$@")
      break
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
    *)
      executables+=("$1")
      shift
      ;;
  esac
done

if ((${#executables[@]} == 0)); then
  echo "At least one executable must be provided." >&2
  usage >&2
  exit 2
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "${script_dir}/.." && pwd)
cd "${repo_root}"

mkdir -p "${log_dir}"
log_dir=$(cd "${log_dir}" && pwd)

# shellcheck source=CI-scripts/v6_namelist_helpers.sh
source "${script_dir}/v6_namelist_helpers.sh"

run_namelist() {
  local exe="$1"
  local nml="$2"
  local exe_name nml_name log_path param_file out_file expected_failure status normalized_output normalized_expected
  local -a cmd=()

  exe_name=$(basename "${exe}")
  nml_name=$(basename "${nml}" .nml)
  log_path="${run_dir}/${nml_name}.log"
  param_file=$(parameter_file_for_namelist "${nml_name}")
  out_file=$(output_file_for_namelist "${nml_name}")

  cmd+=("${exe}" -n "${nml}" -p "${param_file}" -o "${out_file}")
  expected_failure=$(expected_failure_diagnostic "${nml_name}")

  echo "Running ${exe_name} with ${nml}"
  if ((threads > 0)); then
    if env OMP_NUM_THREADS="${threads}" "${cmd[@]}" >"${log_path}" 2>&1; then
      status=0
    else
      status=$?
    fi
  else
    if "${cmd[@]}" >"${log_path}" 2>&1; then
      status=0
    else
      status=$?
    fi
  fi

  if [[ -n "${expected_failure}" ]]; then
    if ((status == 0)); then
      echo "Expected failure succeeded: ${exe_name} ${nml}" >&2
      exit 1
    fi
    normalized_output=$(tr -d '[:space:]' <"${log_path}")
    normalized_expected=$(printf '%s' "${expected_failure}" | tr -d '[:space:]')
    if [[ "${normalized_output}" != *"${normalized_expected}"* ]]; then
      echo "Expected failure diagnostic missing: ${exe_name} ${nml}" >&2
      tail -n 50 "${log_path}" >&2 || true
      exit 1
    fi
    echo "Expected failure confirmed for ${exe_name} with ${nml}"
  elif ((status != 0)); then
      echo "Run failed: ${exe_name} ${nml}" >&2
      tail -n 50 "${log_path}" >&2 || true
      exit 1
  fi

  echo "Completed ${exe_name} with ${nml}"
  post_run_check "${nml_name}"
}

namelists=(
  test_nml/mpr_pet_aspect_minimal.nml
  test_nml/mpr_temporal_lc_pet_lai_minimal.nml
  test_nml/mpr_pet_hargreaves_minimal.nml
  test_nml/mpr_pet_priestley_taylor_minimal.nml
  test_nml/mpr_pet_penman_minimal.nml
  test_nml/mpr_pet_aspect_weights_minimal.nml
  test_nml/mpr_snow_pet_lai_minimal.nml
  test_nml/mpr_runoff_baseflow_minimal.nml
  test_nml/mpr_neutrons_desilets_minimal.nml
  test_nml/mpr_neutrons_cosmic_minimal.nml
  test_nml/mhm_interception_snow_passthrough_minimal.nml
  test_nml/mhm_runoff_baseflow_restart_write.nml
  test_nml/mhm_runoff_baseflow_restart_read.nml
  test_nml/mhm_output_minimal.nml
  test_nml/mrm_rivertemp_meteo_minimal.nml
  test_nml/mrm_coarse_routing_output_minimal.nml
  test_nml/mrm_minimal.nml
  test_nml/mrm_minimal1.nml
  test_nml/mrm_minimal2.nml
  test_nml/mrm_lake_map_minimal.nml
  test_nml/mrm_lake_restart_write.nml
  test_nml/mrm_lake_restart_reference.nml
  test_nml/mrm_lake_restart_only.nml
  test_nml/mlm_lake_case1_pt.nml
  test_nml/mlm_lake_case1_cfpt.nml
  test_nml/mlm_lake_case1_penman.nml
  test_nml/mlm_lake_case1_albedo.nml
  test_nml/mlm_lake_case1_restart_write.nml
  test_nml/mlm_lake_case1_restart_read.nml
  test_nml/mlm_lake_case1_restart_only.nml
  test_nml/mlm_lake_case1_topology_write.nml
  test_nml/mlm_lake_case1_topology_only.nml
  test_nml/mlm_lake_case_minus2_restart_write.nml
  test_nml/mlm_lake_case_minus2_restart_only.nml
  test_nml/mlm_lake_case_minus2_topology_write.nml
  test_nml/mlm_lake_case_minus2_topology_only.nml
)

for exe in "${executables[@]}"; do
  cd "${repo_root}"
  if [[ ! -x "${exe}" ]]; then
    echo "Executable not found or not executable: ${exe}" >&2
    exit 1
  fi

  exe=$(cd "$(dirname "${exe}")" && pwd)/$(basename "${exe}")
  run_dir=$(mktemp -d "${log_dir}/$(basename "${exe}").XXXXXX")
  cmake -DSOURCE_DIR="${repo_root}" -DRUN_DIR="${run_dir}" -P "${script_dir}/stage_run.cmake"
  echo "Logs and outputs: ${run_dir}"
  cd "${run_dir}"
  for nml in "${namelists[@]}"; do
    if [[ ! -f "${nml}" ]]; then
      echo "Smoke namelist not found: ${nml}" >&2
      exit 1
    fi
    run_namelist "${exe}" "${nml}"
  done
done
