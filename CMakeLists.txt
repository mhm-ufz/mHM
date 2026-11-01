# common (default) options to cmake
# -DCMAKE_BUILD_TYPE=Release     - compile in debug, release, or relwithdebinfo mode
# -DBUILD_MHM_DRIVER=ON          - build mHM executable
# -DBUILD_MHM_PYBIND=ON          - build python bindings
# -DBUILD_MHM_LIB_SHARED=OFF     - build mhm_lib as a shared library
# -DBUILD_MHM_LIB_PIC=OFF        - build mhm_lib with position independent code
# -DBUILD_TESTING=ON             - build mHM pfUnit tests
# -DCMAKE_WITH_COVERAGE=ON       - build with GNU coverage support
# -DCMAKE_WITH_MPI=ON            - build with MPI support
# -DCMAKE_WITH_OpenMP=ON         - build with OpenMP support
# -DMHM_ENABLE_NATIVE=ON         - enable host-native release tuning for local/HPC builds
# -DCMAKE_WITH_GPROF=ON          - build with GNU gprof instrumentation
# -DCMAKE_VERBOSE_MAKEFILE=OFF   - see all the commands
# -DCMAKE_FIND_LIBRARY_CUSTOM_LIB_SUFFIX=64     - find extra library paths on some systems
cmake_minimum_required(VERSION 3.18)

# get version and date from files (version.txt and version_date.txt)
include(cmake/version.cmake)
get_version(MHM SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

# create the project
project(mhm
  VERSION ${MHM_VER}
  DESCRIPTION "The mesoscale Hydrological Model"
  HOMEPAGE_URL "https://www.ufz.de/mhm"
  LANGUAGES Fortran C
)

set(MHM_IS_TOP_LEVEL OFF)
if(DEFINED PROJECT_IS_TOP_LEVEL)
  set(MHM_IS_TOP_LEVEL ${PROJECT_IS_TOP_LEVEL})
elseif(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(MHM_IS_TOP_LEVEL ON)
endif()

message(STATUS "mHM: build type - ${CMAKE_BUILD_TYPE}")

option(BUILD_MHM_DRIVER "Build mHM with executable driver." ON)
option(BUILD_MHM_PYBIND "Build mHM python bindings." OFF)
option(BUILD_MHM_LIB_SHARED "Build mhm library as shared." OFF)
option(BUILD_MHM_LIB_PIC "Build mhm library as position independent." OFF)
option(BUILD_TESTING "Build with pfUnit tests." OFF)
option(CMAKE_WITH_COVERAGE "Build with GNU coverage support." OFF)
option(CMAKE_WITH_MPI "Build the module with MPI support." OFF)
set(CMAKE_WITH_OpenMP_DEFAULT OFF)
if(DEFINED ENV{MHM_BUILD_PARALLEL})
  set(CMAKE_WITH_OpenMP_DEFAULT "$ENV{MHM_BUILD_PARALLEL}")
endif()
option(CMAKE_WITH_OpenMP "Build the module with OpenMP support." ${CMAKE_WITH_OpenMP_DEFAULT})
option(MHM_ENABLE_NATIVE "Enable host-native tuning for Release builds." OFF)
option(CMAKE_WITH_GPROF "Build with GNU gprof instrumentation." OFF)

if(BUILD_MHM_PYBIND)
  set(BUILD_MHM_LIB_SHARED OFF)
  set(BUILD_MHM_LIB_PIC ON)
  set(BUILD_MHM_DRIVER OFF)
endif()

if(MHM_IS_TOP_LEVEL)
  include(CTest)
endif()

if(CMAKE_WITH_COVERAGE)
  if(NOT BUILD_TESTING OR NOT BUILD_MHM_DRIVER)
    message(FATAL_ERROR "mHM: CMAKE_WITH_COVERAGE requires BUILD_TESTING=ON and BUILD_MHM_DRIVER=ON.")
  endif()
  if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message(FATAL_ERROR "mHM: CMAKE_WITH_COVERAGE requires GNU Fortran.")
  endif()
endif()

add_subdirectory(src)

if(BUILD_MHM_DRIVER)
  add_executable(mhm src/mHM/mhm_driver.f90)
  target_link_libraries(mhm PRIVATE mhm_lib)
  install(TARGETS mhm DESTINATION bin)
endif()

if(BUILD_MHM_PYBIND)
  add_subdirectory(pybind/src)
endif()

# add full version and date to pre-processor flags (quotes need in before hand)
target_compile_definitions(mhm_lib PRIVATE
  MHMVERSION='${MHM_VER_DEV}'
  MHMDATE='${MHM_DATE}'
)

# setup coverage with GNU
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU" AND CMAKE_WITH_COVERAGE AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/CI-scripts" AND TARGET mhm)
  message(STATUS "mHM: add coverage")
  include(cmake/CodeCoverage.cmake)
  append_coverage_compiler_flags_to_target(mhm_lib)
  append_coverage_compiler_flags_to_target(mhm)
  SETUP_TARGET_FOR_COVERAGE_LCOV(
    NAME mhm_coverage_CI
    EXECUTABLE ${CMAKE_CURRENT_SOURCE_DIR}/CI-scripts/run_cmake_coverage.sh
    EXECUTABLE_ARGS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES mhm mhm_lib
    # Preserve overlapping GNU Fortran source-line counters with LCOV 2.2+.
    LCOV_ARGS --rc check_data_consistency=0
    GENHTML_ARGS --rc check_data_consistency=0 -t "mHM coverage"
      --prefix ${CMAKE_CURRENT_SOURCE_DIR}
      --html-prolog ${CMAKE_CURRENT_SOURCE_DIR}/doc/html_files/cov_header.prolog
  )
endif()

# add pfunit test folder
if(BUILD_TESTING AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/tests")
  add_subdirectory(src/tests)
endif()
