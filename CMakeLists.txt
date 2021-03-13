# mHM cmake script
cmake_minimum_required(VERSION 3.12)
# additional cmake-modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)
# get version and date from file
include(version)
get_version(MHM_VER MHM_VER_DEV MHM_DATE)
# create the project
project(mhm
  VERSION ${MHM_VER}
  DESCRIPTION "The mesoscale Hydrological Model"
  HOMEPAGE_URL "https://www.ufz.de/mhm"
  LANGUAGES Fortran
  )
# add full version and date to pre-processor flags (qoutes need in before hand)
add_compile_definitions(MHMVERSION='${MHM_VER_DEV}' MHMDATE='${MHM_DATE}')
# Compile
add_subdirectory(./src)  # the lib folder containing mhm_lib
add_executable(mhm ./src/mHM/mhm_driver.f90)
target_link_libraries(mhm PRIVATE mhm_lib)
# setup coverage with GNU
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU" AND CMAKE_WITH_COVERAGE)
  include(CodeCoverage)
  append_coverage_compiler_flags()
  SETUP_TARGET_FOR_COVERAGE_LCOV(
    NAME mhm_coverage_CI
    EXECUTABLE ../CI-scripts/run_cmake_coverage.sh
    DEPENDENCIES mhm mhm_lib
    EXCLUDE src/lib/*
    GENHTML_ARGS -t "mHM coverage" --html-prolog ../doc/html_files/cov_header.prolog
  )
endif()

# Usually that works fine, except, one is on a module system and tries to execute the executable
# in the end without having the modules loaded. A workaround is provided using the variable
# CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT
# if this variable is set to ON (do not set the variable inside of this cmake file), then the
# paths are added to the INSTALL_RPATH, and via the second command also to the build.
# It is a bit of a mess and workaround though.
if (CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT)
  set_target_properties(mhm
    PROPERTIES
    INSTALL_RPATH "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}"
    BUILD_WITH_INSTALL_RPATH ON
  )
endif()
