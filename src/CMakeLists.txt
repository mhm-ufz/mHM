list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/cmake-modules)

set(LIB_NAME mhm_lib)
# use all mo_*.f90 files for the library
file(GLOB sources_mhm ./mHM/*mo_*.f90)
file(GLOB sources_mrm ./mRM/*mo_*.f90)
file(GLOB sources_common ./common/*mo_*.f90)
list(APPEND sources ${sources_mhm} ${sources_mrm} ${sources_common})
add_library(${LIB_NAME} ${sources})

include(CPM)
CPMAddPackage(
  NAME            forces
  GIT_REPOSITORY  https://git.ufz.de/chs/forces.git
  GIT_TAG         fix_error_message
  OPTIONS
    "BUILD_TESTING OFF"
)
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
  GIT_TAG         cmake_python_option
  # SOURCE_SUBDIR   src
  OPTIONS
    "BUILD_TESTING OFF"
)
target_link_libraries(${LIB_NAME} PUBLIC mpr)

# get pre-processor flag for current compiler (-fpp or -cpp) and 'cpp_definitions' function
include(fortranpreprocessor)
get_preproc_flag(XPP_FLAG)
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
    $<$<CONFIG:DEBUG>:-Og -g -Wall -Wextra -fimplicit-none -fbacktrace -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=snan -pedantic-errors>
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
