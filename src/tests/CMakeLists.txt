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
    test_lake_bathymetry.pf
    test_lake_evaporation.pf
    test_lake_soil_lut.pf
    test_lake_spill.pf
    test_lake_spillway_length.pf
    test_lake_water_balance.pf
    test_meteo_temporal_tools.pf
    test_mlm_restart.pf
    test_mrm_output_timing.pf
    test_mpr_regionalization.pf
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

add_executable(test_forward_run test_forward_run.f90)
target_link_libraries(test_forward_run PRIVATE ${LIB_NAME})
add_test(NAME test_forward_run COMMAND test_forward_run)

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

add_executable(test_bulk_components test_bulk_components.f90)
target_link_libraries(test_bulk_components PRIVATE ${LIB_NAME})
add_test(NAME test_bulk_components COMMAND test_bulk_components)

add_executable(test_horizon_restart test_horizon_restart.f90)
target_link_libraries(test_horizon_restart PRIVATE ${LIB_NAME})
add_test(NAME test_horizon_restart COMMAND test_horizon_restart)

add_executable(test_layered_mhm_output test_layered_mhm_output.f90)
target_link_libraries(test_layered_mhm_output PRIVATE ${LIB_NAME})
add_test(NAME test_layered_mhm_output COMMAND test_layered_mhm_output)

add_executable(test_soil_horizons test_soil_horizons.f90)
target_link_libraries(test_soil_horizons PRIVATE ${LIB_NAME})
add_test(NAME test_soil_horizons COMMAND test_soil_horizons)
add_test(NAME test_horizon_failures COMMAND ${CMAKE_COMMAND}
  -DTEST_EXECUTABLE=$<TARGET_FILE:test_soil_horizons>
  -DRESTART_EXECUTABLE=$<TARGET_FILE:test_horizon_restart>
  -DRUN_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR}/runs/test_horizon_failures
  -P ${CMAKE_CURRENT_SOURCE_DIR}/run_horizon_failures.cmake)

set_source_files_properties(test_neutron_parallel.f90 PROPERTIES Fortran_PREPROCESS ON)
add_executable(test_neutron_parallel test_neutron_parallel.f90)
target_link_libraries(test_neutron_parallel PRIVATE ${LIB_NAME})
add_test(NAME test_neutron_parallel COMMAND test_neutron_parallel)
set_tests_properties(test_neutron_parallel PROPERTIES ENVIRONMENT "OMP_DYNAMIC=FALSE")

add_executable(test_lake_input_options test_lake_input_options.f90)
target_link_libraries(test_lake_input_options PRIVATE ${LIB_NAME})
add_test(NAME test_lake_input_options COMMAND test_lake_input_options)

add_executable(test_lake_borrowed_parameters test_lake_borrowed_parameters.f90)
target_link_libraries(test_lake_borrowed_parameters PRIVATE ${LIB_NAME})
add_test(NAME test_lake_borrowed_parameters COMMAND test_lake_borrowed_parameters)
add_test(NAME test_lake_borrowed_failures COMMAND ${CMAKE_COMMAND}
  -DTEST_EXECUTABLE=$<TARGET_FILE:test_lake_borrowed_parameters>
  -DRUN_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR}/runs/test_lake_borrowed_failures
  -P ${CMAKE_CURRENT_SOURCE_DIR}/run_lake_borrowed_failures.cmake)

add_executable(test_lake_continuation test_lake_continuation.f90)
target_link_libraries(test_lake_continuation PRIVATE ${LIB_NAME})
add_test(NAME test_lake_continuation COMMAND test_lake_continuation)
add_test(NAME test_mlm_restart_failures COMMAND ${CMAKE_COMMAND}
  -DTEST_EXECUTABLE=$<TARGET_FILE:test_lake_continuation>
  -DRUN_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR}/runs/test_mlm_restart_failures
  -P ${CMAKE_CURRENT_SOURCE_DIR}/run_mlm_restart_failures.cmake)

add_executable(test_mlm_parameters test_mlm_parameters.f90)
target_link_libraries(test_mlm_parameters PRIVATE ${LIB_NAME})
add_test(NAME test_mlm_parameters COMMAND test_mlm_parameters)

add_executable(test_lake_update test_lake_update.f90)
target_link_libraries(test_lake_update PRIVATE ${LIB_NAME})
add_test(NAME test_lake_update COMMAND test_lake_update)
set_tests_properties(test_lake_update PROPERTIES ENVIRONMENT "OMP_DYNAMIC=FALSE")

add_executable(test_meteo_lake_restart test_meteo_lake_restart.f90)
target_link_libraries(test_meteo_lake_restart PRIVATE ${LIB_NAME})
add_test(NAME test_meteo_lake_restart COMMAND test_meteo_lake_restart)
add_test(NAME test_meteo_lake_floatweights COMMAND test_meteo_lake_restart floatweights)
add_test(NAME test_meteo_lake_geometry COMMAND test_meteo_lake_restart geometry)
add_test(NAME test_meteo_lake_default_mode COMMAND test_meteo_lake_restart default_mode)
add_test(NAME test_meteo_lake_file_weighted COMMAND test_meteo_lake_restart file_weighted)
add_test(NAME test_meteo_lake_varying_weights COMMAND test_meteo_lake_restart varying_weights)
add_test(NAME test_meteo_lake_fractions COMMAND test_meteo_lake_restart fractions)
add_test(NAME test_meteo_lake_other_domain COMMAND test_meteo_lake_restart other_domain)
add_test(NAME test_meteo_lake_rebuild COMMAND test_meteo_lake_restart rebuild)
add_test(NAME test_meteo_lake_rebuild_weighted COMMAND test_meteo_lake_restart rebuild_weighted)
add_test(NAME test_lake_restart_failures COMMAND ${CMAKE_COMMAND}
  -DTEST_EXECUTABLE=$<TARGET_FILE:test_meteo_lake_restart>
  -DRUN_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR}/runs/test_lake_restart_failures
  -P ${CMAKE_CURRENT_SOURCE_DIR}/run_lake_restart_failures.cmake)

# No integration or unit test writes into the source checkout.
get_property(mhm_tests DIRECTORY PROPERTY TESTS)
foreach(test_name IN LISTS mhm_tests)
  set(run_dir "${CMAKE_CURRENT_BINARY_DIR}/runs/${test_name}")
  mhm_stage_run("${PROJECT_SOURCE_DIR}" "${run_dir}")
  if(test_name STREQUAL "test_lake_input_options")
    configure_file("${PROJECT_SOURCE_DIR}/examples/fixtures/lake_outlets.nc" "${run_dir}/mlm_lake_optional.nc" COPYONLY)
    configure_file("${PROJECT_SOURCE_DIR}/examples/fixtures/lake_outlets.nc" "${run_dir}/mlm_lake_topology_input.nc" COPYONLY)
  endif()
  if(test_name STREQUAL "test_forward_run")
    mhm_stage_run("${PROJECT_SOURCE_DIR}" "${run_dir}/alternate")
    foreach(api_run_dir "${run_dir}" "${run_dir}/alternate")
      configure_file("${api_run_dir}/test_nml/mrm_minimal.nml" "${api_run_dir}/mhm.nml" COPYONLY)
      configure_file("${api_run_dir}/test_nml/mhm_parameter_v6_routing.nml" "${api_run_dir}/parameters.nml" COPYONLY)
      configure_file("${api_run_dir}/mhm-output-template.nml" "${api_run_dir}/outputs.nml" COPYONLY)
    endforeach()
  endif()
  set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY "${run_dir}")
endforeach()

if(TARGET mhm_coverage_CI)
  get_property(mhm_test_targets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
  add_dependencies(mhm_coverage_CI ${mhm_test_targets})
endif()
