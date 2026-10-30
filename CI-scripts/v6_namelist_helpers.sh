#!/usr/bin/env bash
# Shared file selections and checks for curated v6 scenarios.

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
    mrm_minimal*|mrm_lake_restart_*)
      printf '%s\n' "test_nml/mhm_parameter_v6_routing.nml"
      ;;
    mlm_lake_case1_*|mlm_lake_case_minus2_*)
      printf '%s\n' "test_nml/mhm_parameter_v6_lake.nml"
      ;;
    mrm_lake_map_minimal)
      printf '%s\n' "test_nml/mhm_parameter_v6_routing.nml"
      ;;
    mrm_coarse_routing_output_minimal)
      printf '%s\n' "test_nml/mhm_parameter_v6_routing_slow.nml"
      ;;
    *)
      echo "No v6 parameter fixture configured for ${nml_name}." >&2
      return 1
      ;;
  esac
}

post_run_check() {
  local nml_name="$1"

  case "${nml_name}" in
    mlm_lake_case1_restart_only|mlm_lake_case1_topology_only)
      cmp mlm_lake_case1_restart_read.nc "${nml_name}.nc"
      cmp mlm_lake_case1_restart_read_mlm.nc "${nml_name}_mlm.nc"
      ;;
    mrm_lake_restart_only)
      if ! cmp -s mrm_lake_restart_reference.nc mrm_lake_restart_only.nc; then
        echo "Lake restart-only continuation differs from the Input-backed continuation" >&2
        exit 1
      fi
      if ! cmp -s mrm_lake_restart_reference_mlm.nc mrm_lake_restart_only_mlm.nc; then
        echo "mLM restart-only output differs from the Input-backed continuation" >&2
        exit 1
      fi
      ;;
  esac
}

expected_failure_diagnostic() {
  local nml_name="$1"

  case "${nml_name}" in
    mrm_lake_map_minimal)
      printf '%s\n' "mRM: lake_outflow has no provider"
      ;;
  esac
}

output_file_for_namelist() {
  local nml_name="$1"

  case "${nml_name}" in
    mhm_runoff_baseflow_restart_*)
      printf '%s\n' "test_nml/mhm_output_restart_daily.nml"
      ;;
    mhm_output_*)
      printf '%s\n' "test_nml/mhm_output_v6_smoke.nml"
      ;;
    mrm_coarse_routing_output_minimal)
      printf '%s\n' "test_nml/mrm_output_hourly.nml"
      ;;
    mlm_lake_case1_*|mlm_lake_case_minus2_*)
      printf '%s\n' "test_nml/mlm_output_lake_full.nml"
      ;;
    mrm_lake_restart_*)
      printf '%s\n' "test_nml/mlm_output_daily.nml"
      ;;
    *)
      printf '%s\n' "mhm-output-template.nml"
      ;;
  esac
}
