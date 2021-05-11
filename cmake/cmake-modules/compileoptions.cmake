# The variable "CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT" can be set before executing cmake via a cache command:
# $cmake -DCMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT=ON ..
# or cache file:
# $cmake -C ../CMakeCacheFiles/eve ..
# or after executing CMake editing the CMakeCache.txt, preferably with a corresponding cmake editor i.e ccmake
option(CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT "build the module INDEPENDENT of the module system, so the build in the build tree works even after a module purge")
message(STATUS "build INDEPENDENT of module system ${CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT}")

# The variable "CMAKE_WITH_MPI" can be set before executing cmake via a cache command:
# $cmake -DCMAKE_WITH_MPI=ON ..
# or in a cache file:
# $cmake -C ../CMakeCacheFiles/example
# or after executing CMake editing the CMakeCache.txt,
# preferably with a corresponding cmake editor i.e. ccmake
# same with OpenMP, lapack, coverage (all OFF by default)
option(CMAKE_WITH_MPI "build the module with MPI, so it can be executed using mpirun")
option(CMAKE_WITH_OpenMP "build the module with OpenMP parallelization")
option(CMAKE_WITH_LAPACK "build the module with lapack library")
option(CMAKE_WITH_COVERAGE "build the module with gcov coverage support")
option(CMAKE_WITH_GPROF "enable generation of profiling information with gnus gprof utility tool")

# if cmake provides a findLIBRARY module, this gets invoked via find_package(LIBRARY)
if (CMAKE_WITH_MPI)
  # find if there is an MPI setup on the system and if so, set corresponding variables
  find_package(MPI)
  if (NOT ${MPI_Fortran_FOUND})
    message(FATAL_ERROR "MPI required but not found")
  else()
    message(STATUS "found MPI_Fortran_COMPILER ${MPI_Fortran_COMPILER}")
  endif()
endif()
if (CMAKE_WITH_OpenMP)
  # find if there is an OpenMP setup on the system and if so, set corresponding variables
  find_package(OpenMP)
  if (NOT ${OpenMP_Fortran_FOUND})
    message(FATAL_ERROR "OpenMP required but not found")
  endif()
endif()
if (CMAKE_WITH_LAPACK)
  # find if there is an LAPACK library on the system and if so, set corresponding variables
  find_package(LAPACK)
  if (NOT ${LAPACK_FOUND})
    message(FATAL_ERROR "lapack required but not found")
  endif()
endif()
