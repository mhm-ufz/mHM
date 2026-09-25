function(expect_exchange_failure scenario diagnostic)
  execute_process(
    COMMAND "${TEST_EXECUTABLE}" "${scenario}"
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr
  )

  set(output "${stdout}\n${stderr}")
  if("${result}" STREQUAL "0")
    message(FATAL_ERROR "Exchange failure scenario '${scenario}' unexpectedly succeeded.\n${output}")
  endif()

  string(REGEX REPLACE "[ \t\r\n]+" "" normalized_output "${output}")
  string(REGEX REPLACE "[ \t\r\n]+" "" normalized_diagnostic "${diagnostic}")
  string(FIND "${normalized_output}" "${normalized_diagnostic}" diagnostic_position)
  if(diagnostic_position EQUAL -1)
    message(FATAL_ERROR "Exchange failure scenario '${scenario}' did not emit '${diagnostic}'.\n${output}")
  endif()
endfunction()

expect_exchange_failure(duplicate-provider "Meteo: duplicate provider declaration for test; existing provider is Input.")
expect_exchange_failure(missing-provider "mHM: test has no provider.")
expect_exchange_failure(storage-shape-change "Input: storage shape or stepping changed for test.")
expect_exchange_failure(missing-data "mRM: <unnamed> data not connected.")
expect_exchange_failure(storage-copy-shape "Exchange storage shape mismatch for test.")
expect_exchange_failure(static-stepping "Input: static exchange field has non-static stepping: test.")
expect_exchange_failure(alias-unconnected "test: pass-through source is not connected for alias.")
expect_exchange_failure(storage-rank "Input: storage rank mismatch for test.")
expect_exchange_failure(data-rank "MPR: test has unexpected layered shape.")
