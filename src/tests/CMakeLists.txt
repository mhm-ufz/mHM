enable_testing()

set(LIB_NAME mhm_lib)
find_package(PFUNIT)

if(PFUNIT_FOUND)
  set(pfunit_suffix ".pf")
  file(GLOB testfiles ./*${pfunit_suffix})
  foreach(file ${testfiles})
    get_filename_component(testname ${file} NAME_WE)
    add_pfunit_ctest(${testname}
      TEST_SOURCES ${file}
      LINK_LIBRARIES ${LIB_NAME}
    )
  endforeach()
else()
    message(STATUS "No test program found, install PFUNIT")
endif()

