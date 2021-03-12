set (LIB_NAME mhm_lib)
# use all mo_*.f90 files for the library
file(GLOB_RECURSE sources ./*mo_*.f90)
add_library(${LIB_NAME} ${sources})
# target_include_directories(${LIB_NAME} PRIVAT ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake-modules)
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

target_link_options(${LIB_NAME} PUBLIC ${NETCDF_LDFLAGS_OTHER})
target_include_directories(${LIB_NAME} PUBLIC
  ${NETCDF_INCLUDES}
  ${MPI_Fortran_INCLUDE_PATH}
  ${OpenMP_Fortran_LIBRARY}
  )
target_compile_options(${LIB_NAME} PUBLIC ${NETCDF_CFLAGS_OTHER} ${XPP_FLAG})
target_link_libraries(${LIB_NAME} ${NETCDF_LINK_LIBRARIES} ${MPI_Fortran_LIBRARIES})
target_include_directories(${LIB_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
# target_include_directories(${LIB_NAME} PUBLIC ${NETCDF_INCLUDES})

# by setting compile options and definitions PUBLIC, they are also used by
# programms linking agains it (mhm exe in this case)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  target_compile_options(${LIB_NAME} PUBLIC
    -ffree-form -ffixed-line-length-132
    $<$<CONFIG:DEBUG>:-pedantic-errors -Wall -W -O -g -Wno-maybe-uninitialized>
    $<$<CONFIG:RELEASE>:-O3>
    )
  target_compile_definitions(${LIB_NAME} PUBLIC GFORTRAN)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  target_compile_options(${LIB_NAME} PUBLIC
  -nofixed -assume byterecl -fp-model source -m64 -assume realloc-lhs
    $<$<CONFIG:DEBUG>:-warn all -g -debug -traceback -fp-stack-check -O0 -debug -check all>
    $<$<CONFIG:RELEASE>:-O3 -qoverride-limits>
  )
  target_compile_definitions(${LIB_NAME} PUBLIC INTEL)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_compile_options(${LIB_NAME} PUBLIC
    -ieee=full-colour -unsharedf95 -ideclient
    $<$<CONFIG:DEBUG>:-g -nan -O0 -C=all -strict95 -ieee=stop>
    $<$<CONFIG:RELEASE>:-O4 -ieee=full>
  )
  target_compile_definitions(${LIB_NAME} PUBLIC NAG)
endif()

# add pfunit test folder
if(BUILD_TESTING)
  add_subdirectory(./tests)
endif()
