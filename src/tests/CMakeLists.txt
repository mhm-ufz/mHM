enable_testing()
include(${PROJECT_SOURCE_DIR}/CI-scripts/stage_run.cmake)

set(LIB_NAME mhm_lib)
find_package(PFUNIT QUIET)

if(PFUNIT_FOUND)
  message(STATUS "mHM: add unit tests")
  configure_file(test_parameters_expected.nml test_parameters_expected.nml COPYONLY)
  set(testfiles
    test_canopy_interc.pf
    test_exchange_variable.pf
    test_meteo_spatial_tools.pf
    test_meteo_temporal_tools.pf
    test_mlm_restart.pf
    test_mrm_output_timing.pf
    test_parameters.pf
    test_pet.pf
    test_river_router.pf
    test_river_tools.pf
    test_river_upscaler.pf
    test_runoff.pf
    test_snow_accum_melt.pf
    test_soil_moisture.pf
    test_time_config.pf
  )
  foreach(file ${testfiles})
    get_filename_component(testname ${file} NAME_WE)
    add_pfunit_ctest(${testname} TEST_SOURCES ${file} LINK_LIBRARIES ${LIB_NAME})
  endforeach()
else()
    message(STATUS "mHM: No test program found, install PFUNIT")
endif()

# Parameter validation is intentionally fatal inside mHM. These expected-error
# cases run as subprocesses because error stop terminates a pFUnit executable
# before pFUnit can record an assertion. The runner checks both termination and
# the expected diagnostic; ordinary parameter behavior remains in pFUnit.
add_executable(test_parameter_failures test_parameter_failures.f90)
target_link_libraries(test_parameter_failures PRIVATE ${LIB_NAME})
add_test(
  NAME test_parameter_failures
  COMMAND
    ${CMAKE_COMMAND}
      -DTEST_EXECUTABLE=$<TARGET_FILE:test_parameter_failures>
      -DPROJECT_SOURCE_DIR=${CMAKE_CURRENT_BINARY_DIR}/runs/test_parameter_failures
      -P ${CMAKE_CURRENT_SOURCE_DIR}/run_parameter_failures.cmake
)

add_executable(test_exchange_failures test_exchange_failures.f90)
target_link_libraries(test_exchange_failures PRIVATE ${LIB_NAME})
add_test(
  NAME test_exchange_failures
  COMMAND
    ${CMAKE_COMMAND}
      -DTEST_EXECUTABLE=$<TARGET_FILE:test_exchange_failures>
      -DPROJECT_SOURCE_DIR=${CMAKE_CURRENT_BINARY_DIR}/runs/test_exchange_failures
      -P ${CMAKE_CURRENT_SOURCE_DIR}/run_exchange_failures.cmake
)

add_executable(test_river_temperature test_river_temperature.f90)
target_link_libraries(test_river_temperature PRIVATE ${LIB_NAME})
add_test(NAME test_river_temperature COMMAND test_river_temperature)

add_executable(test_river_routing_optimization test_river_routing_optimization.f90)
target_link_libraries(test_river_routing_optimization PRIVATE ${LIB_NAME})
add_test(NAME test_river_routing_optimization COMMAND test_river_routing_optimization)

if(TARGET mhm)
  add_test(
    NAME test_mrm_river_temperature
    COMMAND
      $<TARGET_FILE:mhm>
        -n test_nml/mrm_rivertemp_meteo_minimal.nml
        -p test_nml/mhm_parameter_v6_meteo.nml
        -o test_nml/mrm_rivertemp_output.nml
  )
endif()

add_executable(test_parameter_reinitialize test_parameter_reinitialize.f90)
target_link_libraries(test_parameter_reinitialize PRIVATE ${LIB_NAME})
add_test(NAME test_parameter_reinitialize COMMAND test_parameter_reinitialize)

if(TARGET mhm)
  add_test(
    NAME test_mrm_intermediate_step
    COMMAND
      $<TARGET_FILE:mhm>
        -n test_nml/mrm_8h_minimal.nml
        -p test_nml/mhm_parameter_v6_routing.nml
        -o mhm-output-template.nml
  )
endif()

if(TARGET mhm)
  add_test(
    NAME test_mpr_lake_mask
    COMMAND
      $<TARGET_FILE:mhm>
        -n test_nml/mpr_lake_mask_minimal.nml
        -p test_nml/mhm_parameter_v6_hydrology.nml
        -o mhm-output-template.nml
  )
endif()

add_executable(test_mpr_multidomain_soil test_mpr_multidomain_soil.f90)
target_link_libraries(test_mpr_multidomain_soil PRIVATE ${LIB_NAME})
add_test(NAME test_mpr_multidomain_soil COMMAND test_mpr_multidomain_soil)

# No integration or unit test writes into the source checkout.
get_property(mhm_tests DIRECTORY PROPERTY TESTS)
foreach(test_name IN LISTS mhm_tests)
  set(run_dir "${CMAKE_CURRENT_BINARY_DIR}/runs/${test_name}")
  mhm_stage_run("${PROJECT_SOURCE_DIR}" "${run_dir}")
  set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY "${run_dir}")
endforeach()

if(TARGET mhm_coverage_CI)
  get_property(mhm_test_targets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
  add_dependencies(mhm_coverage_CI ${mhm_test_targets})
endif()
