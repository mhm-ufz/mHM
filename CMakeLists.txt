# common (default) options to cmake
# -DCMAKE_BUILD_TYPE=Release     - compile in debug, release, or relwithdebinfo mode
# -DMHM_BUILD_DRIVER=ON          - build mHM executables
# -DMHM_BUILD_PYBIND=ON          - build python bindings
# -DMHM_BUILD_LIB_SHARED=OFF     - build mhm_lib as a shared library
# -DMHM_BUILD_LIB_PIC=OFF        - build mhm_lib with position independent code
# -DMHM_BUILD_TESTING=ON         - build mHM pfUnit tests
# -DMHM_WITH_COVERAGE=ON         - build with GNU coverage support
# -DMHM_WITH_MPI=ON              - build with MPI support
# -DMHM_WITH_OpenMP=ON           - build with OpenMP support
# -DMHM_ENABLE_NATIVE=ON         - enable host-native release tuning for local/HPC builds
# -DMHM_WITH_GPROF=ON            - build with GNU gprof instrumentation
# -DCMAKE_VERBOSE_MAKEFILE=OFF   - see all the commands
# -DCMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT=OFF   - build the library independent of the module system,
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

include(CMakeDependentOption)

set(MHM_WITH_OpenMP_DEFAULT OFF)
if(DEFINED ENV{MHM_BUILD_PARALLEL})
  set(MHM_WITH_OpenMP_DEFAULT "$ENV{MHM_BUILD_PARALLEL}")
endif()

option(MHM_BUILD_DRIVER "Build mHM executables." ON)
option(MHM_BUILD_PYBIND "Build mHM python bindings." OFF)
option(MHM_BUILD_LIB_SHARED "Build mhm_lib as a shared library." OFF)
option(MHM_BUILD_LIB_PIC "Build mhm_lib with position independent code." OFF)
option(MHM_BUILD_TESTING "Build mHM pfUnit tests." OFF)
cmake_dependent_option(MHM_WITH_COVERAGE "Build with GNU coverage support." OFF
  "MHM_BUILD_TESTING" OFF)
option(MHM_WITH_MPI "Build the module with MPI support." OFF)
option(MHM_WITH_OpenMP "Build the module with OpenMP support." ${MHM_WITH_OpenMP_DEFAULT})
option(MHM_ENABLE_NATIVE "Enable host-native tuning for Release builds." OFF)
option(MHM_WITH_GPROF "Build with GNU gprof instrumentation." OFF)

if(MHM_BUILD_PYBIND)
  set(MHM_BUILD_LIB_SHARED OFF)
  set(MHM_BUILD_LIB_PIC ON)
  set(MHM_BUILD_DRIVER OFF)
endif()

if(MHM_IS_TOP_LEVEL)
  set(BUILD_TESTING "${MHM_BUILD_TESTING}" CACHE BOOL "Enable testing via CTest.")
  include(CTest)
endif()

if(MHM_WITH_COVERAGE)
  if(NOT BUILD_TESTING)
    message(FATAL_ERROR "mHM: MHM_WITH_COVERAGE requires MHM_BUILD_TESTING=ON.")
  endif()
  if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message(FATAL_ERROR "mHM: MHM_WITH_COVERAGE requires GNU Fortran.")
  endif()
endif()

add_subdirectory(src)

if(MHM_BUILD_DRIVER)
  if(EXISTS app)
    message(STATUS "${PROJECT_NAME}: found app")
    add_subdirectory(app)
  endif()
  add_executable(mhm src/mHM/mhm_driver.f90)
  target_link_libraries(mhm PRIVATE mhm_lib)
  install(TARGETS mhm DESTINATION bin)
endif()

if(MHM_BUILD_PYBIND)
  add_subdirectory(pybind/src)
endif()

# add full version and date to pre-processor flags (qoutes need in before hand)
target_compile_definitions(mhm_lib PRIVATE
  MHMVERSION='${MHM_VER_DEV}'
  MHMDATE='${MHM_DATE}'
)

# setup coverage with GNU
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU" AND MHM_WITH_COVERAGE AND EXISTS "CI-scripts")
  message(STATUS "mHM: add coverage")
  include(cmake/CodeCoverage.cmake)
  append_coverage_compiler_flags_to_target(mhm_lib)
  set(mhm_coverage_dependencies mhm_lib)
  if(TARGET mhm)
    append_coverage_compiler_flags_to_target(mhm)
    list(APPEND mhm_coverage_dependencies mhm)
  endif()
  SETUP_TARGET_FOR_COVERAGE_LCOV(
    NAME mhm_coverage_CI
    EXECUTABLE ../CI-scripts/run_cmake_coverage.sh
    DEPENDENCIES ${mhm_coverage_dependencies}
    GENHTML_ARGS -t "mHM coverage" --html-prolog ../doc/html_files/cov_header.prolog
  )
endif()

if(MHM_BUILD_TESTING AND BUILD_TESTING AND EXISTS "src/tests")
  add_subdirectory(src/tests)
endif()

# With this, paths are added to the INSTALL_RPATH, and via the second command also to the build.
if(CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT AND TARGET mhm)
  set_target_properties(mhm
    PROPERTIES
    INSTALL_RPATH "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}"
    BUILD_WITH_INSTALL_RPATH ON
  )
endif()
