# pFUnit cannot catch a Fortran error stop. Run each fatal contract in a fresh
# process and require its specific diagnostic so an unrelated failure cannot
# satisfy the test.
function(expect_parameter_failure scenario diagnostic)
  execute_process(
    COMMAND "${TEST_EXECUTABLE}" "${scenario}"
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr
  )

  set(output "${stdout}\n${stderr}")
  if("${result}" STREQUAL "0")
    message(FATAL_ERROR "Parameter failure scenario '${scenario}' unexpectedly succeeded.\n${output}")
  endif()

  string(FIND "${output}" "${diagnostic}" diagnostic_position)
  if(diagnostic_position EQUAL -1)
    message(
      FATAL_ERROR
      "Parameter failure scenario '${scenario}' did not emit '${diagnostic}'.\n${output}"
    )
  endif()
endfunction()

expect_parameter_failure(
  duplicate-name
  "parameters%add_process: duplicate parameter 'factor' in process 'pet'."
)
expect_parameter_failure(
  duplicate-group
  "parameters%add_process: namelist group 'pet_m1' is already registered."
)
expect_parameter_failure(
  definition-bounds
  "parameters%seal: invalid bounds for 'pet.factor'."
)
expect_parameter_failure(
  initial-bounds
  "parameters%seal: initial value for 'pet.factor' is outside its bounds."
)
expect_parameter_failure(
  strict-bounds
  "parameters%seal: optimized parameter 'pet.factor' requires strict bounds."
)
expect_parameter_failure(
  set-bounds
  "parameters%set: value for 'pet.factor' is outside its bounds."
)
expect_parameter_failure(
  write-bounds
  "parameters%write_namelist: value for 'pet.factor' is outside its bounds."
)
expect_parameter_failure(
  write-nonfinite
  "parameters%write_namelist: value for 'pet.factor' is not finite."
)
expect_parameter_failure(
  write-size
  "parameters%write_namelist: expected"
)
expect_parameter_failure(
  missing-input
  "mRM: failed to validate parameter block 'routing_2': namelist not configured"
)
