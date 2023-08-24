set (LIB_NAME mhm_lib)
# use all mo_*.f90/F90 files for the library
file(GLOB sources_mpr MPR/*mo_*.*90)
file(GLOB sources_mhm mHM/*mo_*.*90)
file(GLOB sources_mrm mRM/*mo_*.*90)
file(GLOB sources_common common/*mo_*.*90)
file(GLOB sources_meteo meteo/*mo_*.*90)
file(GLOB sources_coupling coupling/*mo_*.*90)
list(APPEND sources ${sources_mpr} ${sources_mhm} ${sources_mrm} ${sources_common} ${sources_meteo} ${sources_coupling})
option(BUILD_MHM_LIB_SHARED "Build mhm library as shared." OFF)
if(BUILD_MHM_LIB_SHARED)
  add_library(${LIB_NAME} SHARED ${sources})
else()
  add_library(${LIB_NAME} STATIC ${sources})
endif()
target_include_directories(${LIB_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# add FORCES
include(../cmake/CPM.cmake)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../forces")
  message(STATUS "mHM: found local forces directory")
  set(CPM_forces_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../forces" CACHE PATH "Local source path for FORCES.")
else()
  set(CPM_forces_SOURCE $ENV{MHM_BUILD_FORCES_PATH} CACHE PATH "Local source path for FORCES.")
endif()
# when building sdist download forces to top folder with CPM by setting "source cache"
if(${SKBUILD_STATE} MATCHES "sdist")
  set(FORCES_WITH_NETCDF OFF)  # we don't need netcdf to create sdist
  message(STATUS "mHM: downloading forces for sdist")
  set(CPM_SOURCE_CACHE "${CMAKE_CURRENT_SOURCE_DIR}/deps")
  set(CPM_forces_SOURCE "")
endif()
# check forces version file
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/../version_forces.txt" ver_forces LIMIT_COUNT 1)
# cmake directive to specify forces version: CPM_forces_VERSION
set(CPM_forces_VERSION "${ver_forces}" CACHE STRING "FORCES version to download with CPM.")
CPMAddPackage("https://git.ufz.de/chs/forces.git#${CPM_forces_VERSION}")
if(BUILD_MHM_LIB_SHARED)
  set_property(TARGET forces PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
target_link_libraries(${LIB_NAME} PUBLIC forces)

if(${SKBUILD_STATE} MATCHES "sdist")
    file(REMOVE_RECURSE "${CMAKE_CURRENT_SOURCE_DIR}/../forces")
    file(GLOB _fpath RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/deps/forces" "${CMAKE_CURRENT_SOURCE_DIR}/deps/forces/*")
    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/deps/forces/${_fpath}" DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/..")
    file(RENAME "${CMAKE_CURRENT_SOURCE_DIR}/../${_fpath}" "${CMAKE_CURRENT_SOURCE_DIR}/../forces")
    message(STATUS "mHM: copy forces to ./forces")
endif()

# check mhm specific env-vars
if(DEFINED ENV{MHM_BUILD_PARALLEL})
  set(CMAKE_WITH_OpenMP $ENV{MHM_BUILD_PARALLEL})
endif()
# add all compile options (MPI, OpenMP, Lapack, Coverage)
include(../cmake/compileoptions.cmake)
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
    $<$<NOT:$<BOOL:${BUILD_MHM_PYBIND}>>:-ffree-form -ffixed-line-length-132>
    $<$<CONFIG:DEBUG>:-Og -g -Wall -Wextra -fimplicit-none -fbacktrace -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=snan -pedantic>
    $<$<CONFIG:RELEASE>:-O3>
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
