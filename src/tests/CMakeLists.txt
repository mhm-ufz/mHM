enable_testing()

set(LIB_NAME mhm_lib)
find_package(PFUNIT QUIET)

if(PFUNIT_FOUND)
  message(STATUS "mHM: add unit tests")
  file(GLOB testfiles *.pf)
  foreach(file ${testfiles})
    get_filename_component(testname ${file} NAME_WE)
    add_pfunit_ctest(${testname} TEST_SOURCES ${file} LINK_LIBRARIES ${LIB_NAME})
  endforeach()
else()
    message(STATUS "mHM: No test program found, install PFUNIT")
endif()

foreach(testname
    parameters_duplicate_fail
    parameters_group_duplicate_fail
    parameters_definition_bounds_fail
    parameters_initial_bounds_fail
    parameters_strict_bounds_fail
    parameters_set_bounds_fail
    parameters_write_bounds_fail
    parameters_write_nonfinite_fail
    parameters_write_size_fail
    parameter_input_missing_fail
)
  add_executable(test_${testname} test_${testname}.f90)
  target_link_libraries(test_${testname} PRIVATE ${LIB_NAME})
  add_test(NAME test_${testname} COMMAND test_${testname})
  set_tests_properties(test_${testname} PROPERTIES WILL_FAIL TRUE)
endforeach()

add_executable(test_parameter_reinitialize test_parameter_reinitialize.f90)
target_link_libraries(test_parameter_reinitialize PRIVATE ${LIB_NAME})
add_test(NAME test_parameter_reinitialize COMMAND test_parameter_reinitialize)
set_tests_properties(test_parameter_reinitialize PROPERTIES WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
