enable_testing()

set(LIB_NAME mhm_lib)
find_package(PFUNIT QUIET)

if(PFUNIT_FOUND)
  message(STATUS "mHM: add unit tests")
  configure_file(test_parameters_expected.nml test_parameters_expected.nml COPYONLY)
  file(GLOB testfiles *.pf)
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
      -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/run_parameter_failures.cmake
)

add_executable(test_parameter_reinitialize test_parameter_reinitialize.f90)
target_link_libraries(test_parameter_reinitialize PRIVATE ${LIB_NAME})
add_test(NAME test_parameter_reinitialize COMMAND test_parameter_reinitialize)
set_tests_properties(test_parameter_reinitialize PROPERTIES WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")

add_test(
  NAME test_mrm_intermediate_step
  COMMAND
    $<TARGET_FILE:mhm-driver-exe>
      -n test_nml/mrm_8h_minimal.nml
      -p test_nml/mhm_parameter_v6_routing.nml
      -o mhm-output-template.nml
)
set_tests_properties(test_mrm_intermediate_step PROPERTIES WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
