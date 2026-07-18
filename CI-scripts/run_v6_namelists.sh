#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Run curated v6 smoke namelists with one or more mHM v6 driver executables.

Usage:
  run_v6_namelists.sh [--log-dir DIR] [--threads N] [--mpi N] EXE [EXE ...]

Options:
  --log-dir DIR  Directory for per-run logs. Default: logs_v6
  --threads N    Set OMP_NUM_THREADS=N for each non-MPI run
  --mpi N        Run each executable via "mpirun -n N"
EOF
}

log_dir="logs_v6"
threads=0
mpi_ranks=0
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
    --mpi)
      mpi_ranks="$2"
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

if ((mpi_ranks > 0 && threads > 0)); then
  echo "Use either --mpi or --threads, not both." >&2
  exit 2
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "${script_dir}/.." && pwd)
cd "${repo_root}"

mkdir -p "${log_dir}"

cleanup_generated_outputs() {
  rm -f ./mpr*.nc ./mrm*.nc ./mhm*.nc
}

parameter_file_for_namelist() {
  local nml_name="$1"

  case "${nml_name}" in
    mpr_pet_aspect_*|mpr_temporal_lc_pet_lai_*|mpr_pet_hargreaves_*|mpr_pet_priestley_taylor_*|mpr_pet_penman_*|mrm_rivertemp_*)
      printf '%s\n' "test_nml/mhm_parameter_v6_meteo.nml"
      ;;
    mpr_snow_pet_lai_*|mpr_runoff_baseflow_*|mpr_soil_moisture_*|mpr_neutrons_desilets_*|mpr_neutrons_cosmic_*|mhm_runoff_baseflow_restart_*|mhm_output_*|mhm_interception_snow_passthrough_*)
      printf '%s\n' "test_nml/mhm_parameter_v6_hydrology.nml"
      ;;
    mpr_*)
      printf '%s\n' "mhm-para-template.nml"
      ;;
    mrm_minimal*)
      printf '%s\n' "test_nml/mhm_parameter_v6_routing.nml"
      ;;
    *)
      echo "No v6 parameter fixture configured for ${nml_name}." >&2
      return 1
      ;;
  esac
}

output_file_for_namelist() {
  local nml_name="$1"

  case "${nml_name}" in
    mhm_output_*)
      printf '%s\n' "test_nml/mhm_output_v6_smoke.nml"
      ;;
    *)
      printf '%s\n' "mhm-output-template.nml"
      ;;
  esac
}

run_namelist() {
  local exe="$1"
  local nml="$2"
  local exe_name nml_name log_path param_file out_file
  local -a cmd=()

  exe_name=$(basename "${exe}")
  nml_name=$(basename "${nml}" .nml)
  log_path="${log_dir}/${exe_name}__${nml_name}.log"
  param_file=$(parameter_file_for_namelist "${nml_name}")
  out_file=$(output_file_for_namelist "${nml_name}")

  cmd+=("${exe}" -n "${nml}" -p "${param_file}" -o "${out_file}")

  if ((mpi_ranks > 0)); then
    cmd=(mpirun -n "${mpi_ranks}" "${cmd[@]}")
  fi

  echo "Running ${exe_name} with ${nml}"
  if ((threads > 0)); then
    if ! env OMP_NUM_THREADS="${threads}" "${cmd[@]}" >"${log_path}" 2>&1; then
      echo "Run failed: ${exe_name} ${nml}" >&2
      tail -n 50 "${log_path}" >&2 || true
      exit 1
    fi
  else
    if ! "${cmd[@]}" >"${log_path}" 2>&1; then
      echo "Run failed: ${exe_name} ${nml}" >&2
      tail -n 50 "${log_path}" >&2 || true
      exit 1
    fi
  fi

  echo "Completed ${exe_name} with ${nml}"
}

namelists=(
  test_nml/mpr_pet_aspect_minimal.nml
  test_nml/mpr_temporal_lc_pet_lai_minimal.nml
  test_nml/mpr_pet_hargreaves_minimal.nml
  test_nml/mpr_pet_priestley_taylor_minimal.nml
  test_nml/mpr_pet_penman_minimal.nml
  test_nml/mpr_pet_aspect_hourly6h_minimal.nml
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
  test_nml/mrm_minimal.nml
  test_nml/mrm_minimal1.nml
  test_nml/mrm_minimal2.nml
)

for exe in "${executables[@]}"; do
  if [[ ! -x "${exe}" ]]; then
    echo "Executable not found or not executable: ${exe}" >&2
    exit 1
  fi

  cleanup_generated_outputs
  for nml in "${namelists[@]}"; do
    if [[ ! -f "${nml}" ]]; then
      echo "Smoke namelist not found: ${nml}" >&2
      exit 1
    fi
    run_namelist "${exe}" "${nml}"
  done
done
