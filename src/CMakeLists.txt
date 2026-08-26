set (LIB_NAME mhm_lib)
# use all mo_*.f90/F90 files for the library
file(GLOB sources_mpr MPR/*mo_*.*90)
file(GLOB sources_mhm mHM/*mo_*.*90)
file(GLOB sources_mrm mRM/*mo_*.*90)
file(GLOB sources_common common/*mo_*.*90)
file(GLOB sources_meteo meteo/*mo_*.*90)
file(GLOB sources_coupling coupling/*mo_*.*90)
file(GLOB sources_exchange exchange/*mo_*.*90)
file(GLOB sources_nml nml/*nml_*.*90)
list(APPEND sources ${sources_mpr} ${sources_mhm} ${sources_mrm} ${sources_common} ${sources_meteo} ${sources_coupling} ${sources_exchange} ${sources_nml})

if(MHM_BUILD_LIB_SHARED)
  add_library(${LIB_NAME} SHARED ${sources})
else()
  add_library(${LIB_NAME} STATIC ${sources})
endif()

if(MHM_BUILD_LIB_PIC)
  set_property(TARGET ${LIB_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

target_include_directories(${LIB_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# add FORCES
# when building sdist download forces to top folder with CPM by setting "source cache"
if(DEFINED SKBUILD_STATE AND "${SKBUILD_STATE}" MATCHES "sdist")
  set(FORCES_WITH_NETCDF OFF)  # we don't need netcdf to create sdist
  message(STATUS "mHM: downloading forces for sdist")
  file(REMOVE_RECURSE "${CMAKE_CURRENT_SOURCE_DIR}/../forces")
endif()
# use CPM
include(../cmake/CPM.cmake)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../forces")
  message(STATUS "mHM: found local forces directory")
  set(CPM_forces_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../forces" CACHE PATH "Local source path for FORCES.")
else()
  set(CPM_forces_SOURCE $ENV{MHM_BUILD_FORCES_PATH} CACHE PATH "Local source path for FORCES.")
endif()
# check forces version file
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/../version_forces.txt" ver_forces LIMIT_COUNT 1)
# cmake directive to specify forces version: CPM_forces_VERSION
set(CPM_forces_VERSION "${ver_forces}" CACHE STRING "FORCES version to download with CPM.")
set(FORCES_WITH_MPI ${MHM_WITH_MPI} CACHE BOOL
  "Build FORCES with MPI support (follows MHM_WITH_MPI)." FORCE
)
set(FORCES_WITH_OpenMP ${MHM_WITH_OpenMP} CACHE BOOL
  "Build FORCES with OpenMP support (follows MHM_WITH_OpenMP)." FORCE
)
set(FORCES_OpenMP_LINK_PUBLIC ${MHM_WITH_OpenMP} CACHE BOOL
  "Build FORCES with public linkage to OpenMP (follows MHM_WITH_OpenMP)." FORCE
)
set(FORCES_ENABLE_NATIVE ${MHM_ENABLE_NATIVE} CACHE BOOL
  "Enable host-native tuning in FORCES (follows MHM_ENABLE_NATIVE)." FORCE
)
CPMAddPackage(
  NAME forces
  VERSION "${CPM_forces_VERSION}"
  URL
    "https://github.com/mhm-ufz/FORCES/archive/${CPM_forces_VERSION}.tar.gz"
    "https://git.ufz.de/chs/forces/-/archive/${CPM_forces_VERSION}/forces-${CPM_forces_VERSION}.tar.gz"
  EXCLUDE_FROM_ALL YES
  SYSTEM YES
)
if(DEFINED SKBUILD_METADATA_DIR)
  install(
    FILES "${FORCES_SOURCE_DIR}/LICENSE.md"
    DESTINATION "${SKBUILD_METADATA_DIR}/licenses"
    RENAME "FORCES_LICENSE.md"
  )
  install(
    FILES "${FORCES_SOURCE_DIR}/AUTHORS"
    DESTINATION "${SKBUILD_METADATA_DIR}/licenses"
    RENAME "FORCES_AUTHORS"
  )
endif()
if(MHM_BUILD_LIB_SHARED OR MHM_BUILD_LIB_PIC)
  set_property(TARGET forces PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
target_link_libraries(${LIB_NAME} PUBLIC forces)

if(DEFINED SKBUILD_STATE AND "${SKBUILD_STATE}" MATCHES "sdist")
    file(COPY "${FORCES_SOURCE_DIR}" DESTINATION "${CMAKE_CURRENT_SOURCE_DIR}/..")
    get_filename_component(_fpath "${FORCES_SOURCE_DIR}" NAME)
    file(RENAME "${CMAKE_CURRENT_SOURCE_DIR}/../${_fpath}" "${CMAKE_CURRENT_SOURCE_DIR}/../forces")
    message(STATUS "mHM: copied forces from '${FORCES_SOURCE_DIR}' to './forces'.")
endif()

# add all compile options (MPI, OpenMP, Coverage)
if(MHM_WITH_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  if(NOT MPI_Fortran_HAVE_F08_MODULE)
    message(FATAL_ERROR "mHM: MPI found, but the 'mpi_f08' module is not available.")
  endif()
  message(STATUS "mHM: use MPI")
  target_compile_definitions(${LIB_NAME} PRIVATE MPI)
  target_link_libraries(${LIB_NAME} PRIVATE MPI::MPI_Fortran)
endif()

if(MHM_ENABLE_NATIVE)
  message(STATUS "mHM: enable host-native Release tuning (opt-in local/HPC mode)")
else()
  message(STATUS "mHM: use portable Release tuning")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  target_compile_definitions(${LIB_NAME} PRIVATE "GFORTRAN")
  target_compile_options(${LIB_NAME} PUBLIC
    "$<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>"
  )
  target_compile_options(${LIB_NAME} PRIVATE
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:DEBUG>>:-Og;-Wall;-Wextra;-Wimplicit-interface;-Wsurprising;-fimplicit-none;-fbacktrace;-fcheck=all,no-array-temps;-ffpe-trap=zero,overflow,underflow;-finit-real=snan;-pedantic-errors>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELWITHDEBINFO>>:-fbacktrace>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<BOOL:${MHM_WITH_GPROF}>>:-pg>"
  )
  target_link_options(${LIB_NAME} PUBLIC
    "$<$<AND:$<LINK_LANGUAGE:Fortran>,$<BOOL:${MHM_WITH_GPROF}>>:-pg>"
  )
  if(MHM_ENABLE_NATIVE)
    target_compile_options(${LIB_NAME} PRIVATE
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELEASE>>:-march=native>"
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELWITHDEBINFO>>:-march=native>"
    )
  endif()
  # https://www.scivision.dev/gfortran-15-external-argument-mismatch
  # AND: passed procedures with pointer attribute should be called with positional arguments
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "15.0" AND
     CMAKE_Fortran_COMPILER_VERSION VERSION_LESS          "15.2")
    target_compile_options(${LIB_NAME} PRIVATE
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:DEBUG>>:-Wno-external-argument-mismatch>"
    )
  endif()
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  target_compile_definitions(${LIB_NAME} PRIVATE "INTEL")
  target_compile_options(${LIB_NAME} PRIVATE
    "$<$<COMPILE_LANGUAGE:Fortran>:-nofixed;SHELL:-assume byterecl;-standard-realloc-lhs>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:DEBUG>>:-fp-model=source;SHELL:-warn all;SHELL:-check all;SHELL:-check noarg_temp_created;-debug;-traceback;-fp-stack-check;-O0>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELEASE>>:-qoverride-limits>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELWITHDEBINFO>>:-traceback>"
  )
  if(MHM_ENABLE_NATIVE)
    target_compile_options(${LIB_NAME} PRIVATE
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELEASE>>:-xHost>"
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELWITHDEBINFO>>:-xHost>"
    )
  endif()
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  target_compile_definitions(${LIB_NAME} PRIVATE "INTEL")
  if(WIN32)
    target_compile_options(${LIB_NAME} PRIVATE
      "$<$<COMPILE_LANGUAGE:Fortran>:/nofixed;/assume:byterecl;/standard-realloc-lhs;/heap-arrays>"
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:DEBUG>>:/fp:source;/warn:all;/check:all;/check:noarg_temp_created;/traceback;/Od>"
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELEASE>>:/Qoverride-limits>"
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELWITHDEBINFO>>:/traceback>"
    )
    if(MHM_ENABLE_NATIVE)
      target_compile_options(${LIB_NAME} PRIVATE
        "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELEASE>>:/QxHost>"
        "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELWITHDEBINFO>>:/QxHost>"
      )
    endif()
  else()
    target_compile_options(${LIB_NAME} PRIVATE
      "$<$<COMPILE_LANGUAGE:Fortran>:-nofixed;SHELL:-assume byterecl;-standard-realloc-lhs>"
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:DEBUG>>:-fp-model=source;SHELL:-warn all;SHELL:-check all;SHELL:-check noarg_temp_created;SHELL:-check nouninit;-traceback;-O0>"
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELEASE>>:-qoverride-limits>"
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELWITHDEBINFO>>:-traceback>"
    )
    if(MHM_ENABLE_NATIVE)
      target_compile_options(${LIB_NAME} PRIVATE
        "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELEASE>>:-xHost>"
        "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELWITHDEBINFO>>:-xHost>"
      )
    endif()
  endif()
endif()

# flang (LLVMFlang)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
  target_compile_options(${LIB_NAME} PRIVATE
    "$<$<COMPILE_LANGUAGE:Fortran>:-frealloc-lhs>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:DEBUG>>:-fimplicit-none;-pedantic;-Wall>"
  )
  if(MHM_ENABLE_NATIVE)
    target_compile_options(${LIB_NAME} PRIVATE
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELEASE>>:-march=native>"
      "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELWITHDEBINFO>>:-march=native>"
    )
  endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  if(MHM_ENABLE_NATIVE)
    message(STATUS "mHM: host-native Release tuning has no effect for NAG")
  endif()
  target_compile_definitions(${LIB_NAME} PRIVATE "NAG")
  if(MHM_WITH_OpenMP)
    # NAG's -gline traceback instrumentation destabilizes explicit OpenMP calls.
    set(NAG_DEBUG_FLAGS "-nan;-O0;-C;-C=alias;-C=dangling;-strict95;-ieee=full")
    set(NAG_RELWITHDEBINFO_FLAGS "-ieee=full")
  else()
    set(NAG_DEBUG_FLAGS "-gline;-nan;-O0;-C;-C=alias;-C=dangling;-strict95;-ieee=full")
    set(NAG_RELWITHDEBINFO_FLAGS "-gline;-ieee=full")
  endif()
  target_compile_options(${LIB_NAME} PRIVATE
    "$<$<COMPILE_LANGUAGE:Fortran>:-fpp;-colour;-unsharedf95;-ideclient>"
    # "-C=all" is not set, only "-C -C=alias -C=dangling" and "-ieee=full" instead of "-ieee=stop" because
    # this effectively omits the -C=intovf flag which checks for integer overflow
    # we need to exclude that as the random number generator relies on that technique
    # -ieee=full is needed for mo_utils (is_nan, is_finite etc. fails with -ieee=stop)
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:DEBUG>>:${NAG_DEBUG_FLAGS}>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELEASE>>:-O4;-ieee=full>"
    "$<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:RELWITHDEBINFO>>:${NAG_RELWITHDEBINFO_FLAGS}>"
  )
endif()
