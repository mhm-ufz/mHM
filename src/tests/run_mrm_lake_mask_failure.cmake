execute_process(
  COMMAND
    "${TEST_EXECUTABLE}"
      -n test_nml/mrm_lake_map_minimal.nml
      -p test_nml/mhm_parameter_v6_routing.nml
      -o mhm-output-template.nml
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE stdout
  ERROR_VARIABLE stderr
)

set(output "${stdout}${stderr}")
if(result EQUAL 0)
  message(FATAL_ERROR "Lake-enabled mRM unexpectedly succeeded")
endif()
if(NOT output MATCHES "lake-aware level-3 topology is available, but mLM flux exchange and lake routing are not implemented")
  message(FATAL_ERROR "Lake-enabled mRM failed without the expected diagnostic:\n${output}")
endif()
