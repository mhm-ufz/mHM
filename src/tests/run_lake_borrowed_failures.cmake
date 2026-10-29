foreach(mode frozen_ks frozen_kbf)
  execute_process(COMMAND "${TEST_EXECUTABLE}" "${mode}" WORKING_DIRECTORY "${RUN_DIRECTORY}"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(status EQUAL 0 OR NOT "${output}${error}" MATCHES "parameters generating restored")
    message(FATAL_ERROR "${mode}: expected failure and diagnostic; status=${status}\n${output}${error}")
  endif()
endforeach()
