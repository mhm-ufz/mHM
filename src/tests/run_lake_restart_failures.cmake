foreach(mode empty negative nan bounds duplicate normalization temporal offsets)
  execute_process(COMMAND "${TEST_EXECUTABLE}" "${mode}" WORKING_DIRECTORY "${RUN_DIRECTORY}"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(status EQUAL 0 OR NOT "${output}${error}" MATCHES "malformed lake support restart")
    message(FATAL_ERROR "${mode}: expected failure and diagnostic; status=${status}\n${output}${error}")
  endif()
endforeach()
foreach(mode origin resolution orientation mask ids mapping conflict explicit_false omitted_conflict missing_fraction fraction_conflict rebuild_weights weight_origin)
  execute_process(COMMAND "${TEST_EXECUTABLE}" "${mode}" WORKING_DIRECTORY "${RUN_DIRECTORY}"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(status EQUAL 0 OR NOT "${output}${error}" MATCHES "Meteo:.*(requires compatible|IDs do not match|disagrees|conflicts|rebuilding weighted|lake temporal weight file geometry|must be configured consistently)")
    message(FATAL_ERROR "${mode}: expected failure and diagnostic; status=${status}\n${output}${error}")
  endif()
endforeach()
