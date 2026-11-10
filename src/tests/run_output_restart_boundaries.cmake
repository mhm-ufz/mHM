# Output-only requirements and restart cadence apply only to active output.
function(run_case name namelist parameters outputs diagnostic)
  execute_process(
    COMMAND "${TEST_EXECUTABLE}" -n "${namelist}" -p "${parameters}" -o "${outputs}"
    WORKING_DIRECTORY "${RUN_DIRECTORY}"
    RESULT_VARIABLE status OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr)
  set(log "${stdout}\n${stderr}")
  file(WRITE "${RUN_DIRECTORY}/${name}.log" "${log}")
  if(diagnostic STREQUAL "")
    if(NOT "${status}" STREQUAL "0")
      message(FATAL_ERROR "${name} failed: ${log}")
    endif()
  else()
    string(REGEX REPLACE "[ \t\r\n]+" "" normalized_log "${log}")
    string(REGEX REPLACE "[ \t\r\n]+" "" normalized_diagnostic "${diagnostic}")
    string(FIND "${normalized_log}" "${normalized_diagnostic}" found)
    if("${status}" STREQUAL "0" OR found EQUAL -1)
      message(FATAL_ERROR "${name} did not fail with '${diagnostic}': ${log}")
    endif()
  endif()
endfunction()

file(READ "${RUN_DIRECTORY}/test_nml/mpr_pet_aspect_minimal.nml" pet_config)
string(REPLACE "&config_mhm" "&config_mhm\n  write_restart(1) = .true.\n  restart_output_path(1) = 'disabled_output_restart.nc'"
  disabled_config "${pet_config}")
file(WRITE "${RUN_DIRECTORY}/disabled_output.nml" "${disabled_config}")
file(WRITE "${RUN_DIRECTORY}/output_requests.nml"
  "&output_mhm\n output_frequency=72\n out_pet=.true.\n out_sm=.true.\n out_sm_all=.true.\n out_qd=.true.\n out_aet_all=.true.\n/\n")

# A PET-only run has no soil-saturation provider and ends between 72-hour output points.
file(REMOVE "${RUN_DIRECTORY}/disabled_output_restart.nc")
run_case(disabled_output disabled_output.nml test_nml/mhm_parameter_v6_meteo.nml output_requests.nml "")
if(NOT EXISTS "${RUN_DIRECTORY}/disabled_output_restart.nc")
  message(FATAL_ERROR "Disabled output prevented writing an mHM restart")
endif()

# The same cadence must still reject a restart between active output points.
string(REPLACE "&config_mhm" "&config_mhm\n  output_path(1) = 'active_output.nc'" active_config "${disabled_config}")
file(WRITE "${RUN_DIRECTORY}/active_output.nml" "${active_config}")
run_case(active_output active_output.nml test_nml/mhm_parameter_v6_meteo.nml output_requests.nml
  "mHM restart output time must coincide with an mHM output boundary")

# Lake and routing restarts follow the same rule when their namelists are present.
file(READ "${RUN_DIRECTORY}/test_nml/mlm_lake_case1_restart_write.nml" lake_config)
string(REGEX REPLACE "\n[ \t]*output_path\\(1\\)[^\n]*" "" lake_config "${lake_config}")
file(WRITE "${RUN_DIRECTORY}/lake_disabled_output.nml" "${lake_config}")
file(WRITE "${RUN_DIRECTORY}/lake_output_requests.nml"
  "&output_mlm\n output_frequency=72\n/\n&output_mrm\n output_frequency=72\n/\n")
file(REMOVE "${RUN_DIRECTORY}/mlm_lake_case1_restart_mlm.nc" "${RUN_DIRECTORY}/mlm_lake_case1_restart_mrm.nc")
run_case(lake_disabled_output lake_disabled_output.nml test_nml/mhm_parameter_v6_lake.nml lake_output_requests.nml "")
foreach(component mlm mrm)
  if(NOT EXISTS "${RUN_DIRECTORY}/mlm_lake_case1_restart_${component}.nc")
    message(FATAL_ERROR "Disabled output prevented writing an ${component} restart")
  endif()
  string(REPLACE "&config_${component}" "&config_${component}\n output_path(1) = 'active_${component}.nc'"
    active_config "${lake_config}")
  file(WRITE "${RUN_DIRECTORY}/active_${component}.nml" "${active_config}")
  if(component STREQUAL "mlm")
    set(label mLM)
  else()
    set(label mRM)
  endif()
  run_case(active_${component} active_${component}.nml test_nml/mhm_parameter_v6_lake.nml lake_output_requests.nml
    "${label} restart output time must coincide with an ${label} output boundary")
endforeach()

# Routing completion remains mandatory when the output file contains no groups.
file(READ "${RUN_DIRECTORY}/test_nml/mrm_coarse_routing_output_minimal.nml" routing_config)
string(REPLACE "1990-01-02 00:00" "1990-01-01 01:00" routing_config "${routing_config}")
file(WRITE "${RUN_DIRECTORY}/routing_boundary.nml" "${routing_config}")
file(WRITE "${RUN_DIRECTORY}/no_output.nml" "! Output intentionally omitted.\n")
run_case(routing_boundary routing_boundary.nml test_nml/mhm_parameter_v6_routing_slow.nml no_output.nml
  "mRM simulation end must be a routing boundary")

# Omitted groups and explicitly empty groups must have identical activation and
# restart cadence behavior. Runtime files stay in this test's staged run root.
file(WRITE "${RUN_DIRECTORY}/empty_output.nml" "&output_mhm\n/\n&output_mlm\n/\n&output_mrm\n/\n")
string(REPLACE "&config_mhm" "&config_mhm\n  output_path(1) = 'default_mhm.nc'"
  pet_defaults "${disabled_config}")
string(REPLACE "write_restart(1) = .true." "write_restart(1) = .false."
  pet_defaults "${pet_defaults}")
file(WRITE "${RUN_DIRECTORY}/pet_defaults.nml" "${pet_defaults}")
string(REPLACE "&config_mhm" "&config_mhm\n  output_path(1) = 'default_mhm.nc'"
  pet_restart "${disabled_config}")
file(WRITE "${RUN_DIRECTORY}/pet_defaults_restart.nml" "${pet_restart}")

set(lake_defaults "${lake_config}")
foreach(component mlm mrm)
  string(REPLACE "&config_${component}" "&config_${component}\n output_path(1) = 'default_${component}.nc'"
    lake_defaults "${lake_defaults}")
endforeach()
file(WRITE "${RUN_DIRECTORY}/lake_defaults.nml" "${lake_defaults}")
string(REPLACE "1990-01-02" "1990-01-01 12:00" lake_partial "${lake_config}")
file(WRITE "${RUN_DIRECTORY}/lake_no_path.nml" "${lake_partial}")

file(READ "${RUN_DIRECTORY}/test_nml/mrm_minimal_poi_output1.nml" routing_defaults)
string(REPLACE "1990-02-01" "1990-01-02" routing_defaults "${routing_defaults}")
string(REPLACE "&config_mrm" "&config_mrm\n output_path(1) = 'default_routing.nc'"
  routing_defaults "${routing_defaults}")
file(WRITE "${RUN_DIRECTORY}/routing_defaults.nml" "${routing_defaults}")

set(default_streams default_mhm.nc default_mlm.nc default_mrm.nc default_routing.nc
  mrm_node_output_test1.nc mrm_poi_output_test1.nc)
foreach(outputs no_output.nml empty_output.nml)
  foreach(output IN LISTS default_streams)
    file(REMOVE "${RUN_DIRECTORY}/${output}")
  endforeach()
  # With no output path, optional defaults must not constrain restart times.
  run_case(pet_no_path_${outputs} disabled_output.nml test_nml/mhm_parameter_v6_meteo.nml "${outputs}" "")
  run_case(lake_no_path_${outputs} lake_no_path.nml test_nml/mhm_parameter_v6_lake.nml "${outputs}" "")
  run_case(pet_defaults_${outputs} pet_defaults.nml test_nml/mhm_parameter_v6_meteo.nml "${outputs}" "")
  run_case(lake_defaults_${outputs} lake_defaults.nml test_nml/mhm_parameter_v6_lake.nml "${outputs}" "")
  run_case(routing_defaults_${outputs} routing_defaults.nml test_nml/mhm_parameter_v6_routing.nml "${outputs}" "")
  foreach(output IN LISTS default_streams)
    if(NOT EXISTS "${RUN_DIRECTORY}/${output}")
      message(FATAL_ERROR "Optional defaults did not create configured stream ${output}")
    endif()
    if(outputs STREQUAL "no_output.nml")
      configure_file("${RUN_DIRECTORY}/${output}" "${RUN_DIRECTORY}/${output}.reference" COPYONLY)
    else()
      execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${RUN_DIRECTORY}/${output}" "${RUN_DIRECTORY}/${output}.reference" RESULT_VARIABLE different)
      if(NOT "${different}" STREQUAL "0")
        message(FATAL_ERROR "Omitted and empty groups produced different ${output}")
      endif()
    endif()
  endforeach()
  run_case(pet_default_boundary_${outputs} pet_defaults_restart.nml test_nml/mhm_parameter_v6_meteo.nml "${outputs}"
    "mHM restart output time must coincide with an mHM output boundary")
  foreach(component mlm mrm)
    string(REPLACE "&config_${component}" "&config_${component}\n output_path(1) = 'partial_${component}.nc'"
      stream_config "${lake_partial}")
    file(WRITE "${RUN_DIRECTORY}/lake_${component}.nml" "${stream_config}")
    if(component STREQUAL "mlm")
      set(label mLM)
    else()
      set(label mRM)
    endif()
    run_case(lake_${component}_${outputs} lake_${component}.nml test_nml/mhm_parameter_v6_lake.nml "${outputs}"
      "${label} restart output time must coincide with an ${label} output boundary")
  endforeach()

  # Each routing stream independently activates the default daily cadence.
  string(REPLACE "1990-01-02" "1990-01-01 12:00" routing_partial "${routing_defaults}")
  string(REGEX REPLACE "\n[ \t]*output_(node_|poi_)?path\\(1\\)[^\n]*" "" routing_partial "${routing_partial}")
  file(WRITE "${RUN_DIRECTORY}/routing_no_path.nml" "${routing_partial}")
  run_case(routing_no_path_${outputs} routing_no_path.nml test_nml/mhm_parameter_v6_routing.nml "${outputs}" "")
  foreach(stream path node_path poi_path)
    string(REPLACE "&config_mrm" "&config_mrm\n output_${stream}(1) = 'partial_${stream}.nc'"
      stream_config "${routing_partial}")
    file(WRITE "${RUN_DIRECTORY}/routing_${stream}.nml" "${stream_config}")
    run_case(routing_${stream}_${outputs} routing_${stream}.nml test_nml/mhm_parameter_v6_routing.nml "${outputs}"
      "mRM restart output time must coincide with an mRM output boundary")
  endforeach()
endforeach()
