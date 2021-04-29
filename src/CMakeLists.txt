list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/cmake-modules)

include(FetchContent)

set (LIB_NAME mhm_lib)
# use all mo_*.f90 files for the library
file(GLOB sources_mhm ./mHM/*mo_*.f90)
file(GLOB sources_mrm ./mRM/*mo_*.f90)
file(GLOB sources_common ./common/*mo_*.f90)
file(GLOB sources_common_mHM_mRM ./common_mHM_mRM/*mo_*.f90)
file(GLOB sources_mpr ./MPR/*mo_*.f90)
list(APPEND sources ${sources_mhm} ${sources_mrm} ${sources_common} ${sources_common_mHM_mRM} ${sources_mpr})
add_library(${LIB_NAME} ${sources})

set (FORCES_NAME lightweight_fortran_lib)
find_package(${FORCES_NAME})

if(NOT ${FORCES_NAME}_FOUND)
  FetchContent_Declare(
    ${FORCES_NAME}
    GIT_REPOSITORY git@git.ufz.de:ottor/lightweight_fortran_lib.git
    GIT_TAG        temporary_test
    SOURCE_SUBDIR  src
  )
  set (${FORCES_NAME}_BUILD_TESTING OFF)
  FetchContent_MakeAvailable(${FORCES_NAME})
endif()
target_link_libraries(${LIB_NAME} PUBLIC ${FORCES_NAME})

set (MPR_NAME MPR)
find_package(${MPR_NAME})
if(NOT ${MPR_NAME}_FOUND)
  FetchContent_Declare(
    ${MPR_NAME}
    GIT_REPOSITORY git@git.ufz.de:chs/mpr.git
    GIT_TAG        1.0.0
  )
  set (${MPR_NAME}_BUILD_TESTING OFF)
  FetchContent_MakeAvailable(${MPR_NAME})
endif()
# target_link_libraries(${LIB_NAME} PUBLIC ${MPR_NAME})

# get pre-processor flag for current compiler (-fpp or -cpp) and 'cpp_definitions' function
include(checkfortranpreprocessor)
get_preproc_flag(XPP_FLAG)
# Add definitions. These should later be set via the cache or "cmake -DCMAKE_ABSOFT=ON .." for example.
cpp_definitions(pgiFortran OFF "Code exchange for pgi compiler dependent issues")
cpp_definitions(MPR_STANDALONE OFF "If set to ON, only MPR is compiled")
cpp_definitions(ABSOFT OFF "Documentation to be added.")
# add all compile options (MPI, OpenMP, Lapack, Coverage)
include(compileoptions)
# add MPI pre-processor directive (also for mhm exe)
if (CMAKE_WITH_MPI)
  target_compile_definitions(${LIB_NAME} PUBLIC MPI)
endif()
# all compile and link options are PUBLIC in order to be forwared to 'mhm' exe
# target_link_options only available in cmake 3.13 (but NETCDF_LDFLAGS_OTHER not working)
set_property(TARGET ${LIB_NAME} PROPERTY LINK_FLAGS
  "${MPI_Fortran_LINK_FLAGS} ${OpenMP_Fortran_FLAGS} ${LAPACK_LINKER_FLAGS}"
)
target_compile_options(${LIB_NAME} PUBLIC
  ${XPP_FLAG}
  ${MPI_Fortran_COMPILE_FLAGS}
  ${OpenMP_Fortran_FLAGS}
)
target_include_directories(${LIB_NAME} PUBLIC
  ${MPI_Fortran_INCLUDE_PATH}
  ${OpenMP_Fortran_LIBRARY}
  ${CMAKE_CURRENT_BINARY_DIR} # needed for mhm exe
)
target_link_libraries(${LIB_NAME} PUBLIC
  ${MPI_Fortran_LIBRARIES}
  ${OpenMP_Fortran_LIBRARIES}
  ${LAPACK_LIBRARIES}
)
# by setting compile options and definitions PUBLIC, they are also used by
# programms linking agains it (mhm exe in this case)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  target_compile_definitions(${LIB_NAME} PUBLIC GFORTRAN)
  target_compile_options(${LIB_NAME} PUBLIC
    -ffree-form -ffixed-line-length-132
    $<$<CONFIG:DEBUG>:-pedantic-errors -Wall -W -O -g -Wno-maybe-uninitialized>
    $<$<CONFIG:RELEASE>:-O3>
    $<$<BOOL:${CMAKE_WITH_COVERAGE}>:-g -fprofile-arcs -ftest-coverage>
    $<$<BOOL:${CMAKE_WITH_GPROF}>:-pg>
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # https://discourse.cmake.org/t/preserving-options-with-spaces-in-add-compile-options/1551/2
  target_compile_definitions(${LIB_NAME} PUBLIC INTEL)
  target_compile_options(${LIB_NAME} PUBLIC
    -nofixed "SHELL:-assume byterecl" "SHELL:-fp-model source" -m64 "SHELL:-assume realloc-lhs"
    $<$<CONFIG:DEBUG>:-g "SHELL:-warn all" "SHELL:-check all" -debug -traceback -fp-stack-check -O0>
    $<$<CONFIG:RELEASE>:-O3 -qoverride-limits>
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_compile_definitions(${LIB_NAME} PUBLIC NAG)
  target_compile_options(${LIB_NAME} PUBLIC
    -colour -unsharedf95 -ideclient
    $<$<CONFIG:DEBUG>:-g -nan -O0 -C=all -strict95 -ieee=stop>
    $<$<CONFIG:RELEASE>:-O4 -ieee=full>
  )
endif()
