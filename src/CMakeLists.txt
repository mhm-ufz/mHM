list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/cmake-modules)

include(FetchContent)

set (LIB_NAME mhm_lib)
# use all mo_*.f90/F90 files for the library
file(GLOB sources_mhm ./mHM/*mo_*.*90)
file(GLOB sources_mrm ./mRM/*mo_*.*90)
file(GLOB sources_common ./common/*mo_*.*90)
list(APPEND sources ${sources_mhm} ${sources_mrm} ${sources_common})
option(BUILD_MHM_LIB_SHARED "Build mhm library as shared." OFF)
if(BUILD_MHM_LIB_SHARED)
  add_library(${LIB_NAME} SHARED ${sources})
else()
  add_library(${LIB_NAME} STATIC ${sources})
endif()
target_include_directories(${LIB_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

include(CPM)
CPMAddPackage(
  NAME            forces
  GIT_REPOSITORY  https://git.ufz.de/chs/forces.git
  GIT_TAG         develop
  OPTIONS
    "BUILD_TESTING OFF"
)
if(BUILD_MHM_LIB_SHARED)
  set_property(TARGET forces PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
target_link_libraries(${LIB_NAME} PUBLIC forces)

# Namelist files for MPR
set(GENERATE_TFS_CONFIG
  "${PROJECT_SOURCE_DIR}/mpr_mhm.nml"
  CACHE FILEPATH "path to config file for MPR"
)
set(GENERATE_TFS_PARAM
  "${PROJECT_SOURCE_DIR}/mpr_global_parameter_mhm.nml"
  CACHE FILEPATH "path to config file with extra parameters for MPR"
)
CPMAddPackage(
  NAME            mpr
  GIT_REPOSITORY  https://git.ufz.de/chs/mpr.git
  GIT_TAG         develop
  # SOURCE_SUBDIR   src
  OPTIONS
    "BUILD_TESTING OFF"
)
if(BUILD_MHM_LIB_SHARED)
  set_property(TARGET mpr PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
target_link_libraries(${LIB_NAME} PUBLIC mpr)

# add all compile options (MPI, OpenMP, Lapack, Coverage)
include(compileoptions)
if (CMAKE_WITH_MPI)
  target_compile_definitions(${LIB_NAME} PRIVATE MPI)
  target_link_libraries(${LIB_NAME} PRIVATE MPI::MPI_Fortran)
endif()
if (CMAKE_WITH_OpenMP)
  target_link_libraries(${LIB_NAME} PRIVATE OpenMP::OpenMP_Fortran)
endif()

# by setting compile options and definitions PUBLIC, they are also used by
# programms linking agains it (mhm exe in this case)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  target_compile_definitions(${LIB_NAME} PRIVATE GFORTRAN)
  target_compile_options(${LIB_NAME} PUBLIC
    -ffree-form -ffixed-line-length-132
    $<$<CONFIG:DEBUG>:-Og -g -Wall -Wextra -fimplicit-none -fbacktrace -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=snan -pedantic-errors>
    $<$<CONFIG:RELEASE>:-O3>
    $<$<BOOL:${CMAKE_WITH_COVERAGE}>:-g -fprofile-arcs -ftest-coverage>
    $<$<BOOL:${CMAKE_WITH_GPROF}>:-pg>
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # https://discourse.cmake.org/t/preserving-options-with-spaces-in-add-compile-options/1551/2
  target_compile_definitions(${LIB_NAME} PRIVATE INTEL)
  target_compile_options(${LIB_NAME} PUBLIC
    -nofixed "SHELL:-assume byterecl" "SHELL:-fp-model source" -m64 "SHELL:-assume realloc-lhs"
    $<$<CONFIG:DEBUG>:-g "SHELL:-warn all" "SHELL:-check all" -debug -traceback -fp-stack-check -O0>
    $<$<CONFIG:RELEASE>:-O3 -qoverride-limits>
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_compile_definitions(${LIB_NAME} PRIVATE NAG)
  target_compile_options(${LIB_NAME} PUBLIC
    -colour -unsharedf95 -ideclient
    $<$<CONFIG:DEBUG>:-gline -g -nan -O0 -C=all -strict95 -ieee=stop>
    $<$<CONFIG:RELEASE>:-O4 -ieee=full>
  )
endif()
