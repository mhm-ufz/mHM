list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake-modules)
set (LIB_NAME mhm_lib)
# use all mo_*.f90 files for the library
file(GLOB_RECURSE sources ./*mo_*.f90)
add_library(${LIB_NAME} ${sources})
# add all compile options (MPI, OpenMP, Lapack, Coverage)
include(compileoptions)

# the FindNetCDFF.cmake automatically used by find_package(NetCDFF ...)
set(NETCDF_F90 "YES")  # required interface
find_package(NetCDFF REQUIRED)
# from that module we gain the following variables:
# NETCDF_INCLUDES	: the include directory
# NETCDF_LINK_LIBRARIES : the absolute path to and with the libraries
# NETCDF_CFLAGS_OTHER	: additional compilation flags
# NETCDF_LDFLAGS_OTHER	: additional linking flags

# get pre-processor flag for current compiler (either -fpp or -cpp)
# and 'cpp_definitions' function
include(checkfortranpreprocessor)
get_preproc_flag(XPP_FLAG)
# Add definitions. These should later be set via the cache line file and only have a default value here.
cpp_definitions("-DpgiFortran" "CMAKE_pgiFortran" "OFF" "Code exchange for pgi compiler dependent issues")
cpp_definitions("-DMPR_STANDALONE" "CMAKE_MPR_STANDALONE" "OFF" "If set to ON, only MPR is compiled")
cpp_definitions("-DABSOFT" "CMAKE_ABSOFT" "OFF" "Documentation to be added. If you you are developer, you might edit this string in CMakeLists.txt")

# all compile and link options are PUBLIC in order to be forwared to 'mhm' exe
# target_link_options only available in cmake 3.13
if(CMAKE_VERSION VERSION_LESS 3.13)
  set_property(TARGET ${LIB_NAME} PROPERTY LINK_FLAGS
    "${NETCDF_LDFLAGS_OTHER} ${MPI_Fortran_LINK_FLAGS} ${OpenMP_Fortran_FLAGS} ${LAPACK_LINKER_FLAGS}"
  )
else()
  target_link_options(${LIB_NAME} PUBLIC
    ${NETCDF_LDFLAGS_OTHER}
    ${MPI_Fortran_LINK_FLAGS}
    ${OpenMP_Fortran_FLAGS}
    ${LAPACK_LINKER_FLAGS}
  )
endif()
target_compile_options(${LIB_NAME} PUBLIC
  ${XPP_FLAG}
  ${NETCDF_CFLAGS_OTHER}
  ${MPI_Fortran_COMPILE_FLAGS}
  ${OpenMP_Fortran_FLAGS}
)
target_include_directories(${LIB_NAME} PUBLIC
  ${NETCDF_INCLUDES}
  ${MPI_Fortran_INCLUDE_PATH}
  ${OpenMP_Fortran_LIBRARY}
  ${CMAKE_CURRENT_BINARY_DIR} # needed for mhm exe
)
target_link_libraries(${LIB_NAME} PUBLIC
  ${NETCDF_LINK_LIBRARIES}
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
  )
  if(CMAKE_WITH_COVERAGE)
    target_compile_options(${LIB_NAME} PUBLIC -g -fprofile-arcs -ftest-coverage)
  endif()
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

if (CMAKE_WITH_MPI)
  target_compile_definitions(${LIB_NAME} PUBLIC MPI)
endif()

# Show compile options
get_target_property(MAIN_CFLAGS ${LIB_NAME} COMPILE_OPTIONS)
message(STATUS "mhm_lib compiler flags are: ${MAIN_CFLAGS}")
# Show compile options
get_target_property(MAIN_CDEFS ${LIB_NAME} COMPILE_DEFINITIONS)
message(STATUS "mhm_lib compiler defs are: ${MAIN_CDEFS}")

# automatically enable testing (OFF by default)
option(BUILD_TESTING "" OFF)
include(CTest)

# add pfunit test folder
if(BUILD_TESTING)
  add_subdirectory(./tests)
endif()
