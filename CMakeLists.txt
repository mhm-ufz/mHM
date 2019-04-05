# This cmake file is not meant to be edited for a special setup.
# For special setups use cache line files or command line options, as described a few
# lines ahead concerning module independend builds
cmake_minimum_required(VERSION 3.5)

project(mhm Fortran)

# The variable "CMAKE_BUILD_MODULE_SYSTEM_INDEPENDEND" can be set before executing cmake via a cache command:
# $cmake -DCMAKE_BUILD_MODULE_SYSTEM_INDEPENDEND:STRING=ON ..
# or cache file: 
# $cmake -C ../CMakeCacheFiles/eve ..
# or after executing CMake editing the CMakeCache.txt, preferably with a corresponding cmake editor i.e ccmake
set(CMAKE_BUILD_MODULE_SYSTEM_INDEPENDEND OFF CACHE STRING "build the module independend of the module system, so the build in the build tree works even after a module purge")
message(STATUS "build independend of module system ${CMAKE_BUILD_MODULE_SYSTEM_INDEPENDEND}")

# set specific place where to search for the netCDF directory
set(CMAKE_NETCDF_DIR " " CACHE STRING "set set specific place where to search for the netCDF directory")
message(STATUS "search in additional directory ${CMAKE_NETCDF_DIR} for netCDF")

# The variable "CMAKE_WITH_MPI" can be set before executing cmake via a cache command:
# $cmake -DCMAKE_WITH_MPI:STRING=ON ..
# or in a cache file:
# $cmake -C ../CMakeCacheFiles/example
# or after executing CMake editing the CMakeCache.txt, preferably with a corresponding cmake editor i.e. ccmake
set(CMAKE_WITH_MPI OFF CACHE STRING "build the module with MPI, so it can be executed using mpirun")

# additional cmake-modules created for the purpose of finding netCDF or other libraries ly in the source_directory in
# a folder named cmake-modules. This command tells cmake to search there for Find<module>.cmake files
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

set (NETCDF_F90 "YES")
# the FindNetCDFF.cmake file can be found after we added the cmake-modules folder to the CMAKE_MODULE_PATH
# the build fails, if it is not present
find_package(NetCDFF REQUIRED)
# from that module we gain the following variables:
# NETCDF_INCLUDES	: the include directory
# NETCDF_LINK_LIBRARIES : the absolute path to and with the libraries
# NETCDF_CFLAGS_OTHER	: additional compilation flags
# NETCDF_LDFLAGS_OTHER	: additional linking flags
if (CMAKE_WITH_MPI)
        # find if there is an MPI setup on the system and if so, set corresponding variables
        find_package(MPI)
        if (${MPI_FOUND})
                message(FATAL_ERROR "MPI required but not found")
        else()
                message(STATUS "found MPI_Fortran_COMPILER ${MPI_Fortran_COMPILER}")
        endif()
endif()

include_directories(${NETCDF_INCLUDES} ${MPI_Fortran_INCLUDE_PATH})

# ifort and gfortran need the flag -cpp to interpret definitions like -DMRM2MHM
# the nag compiler is not able to interpret the flag -cpp but can interpret these definitions anyway
# so we check whether the compiler is able to use the flag -cpp
# for that we need the module CheckFortranCompilerFlag
include(CheckFortranCompilerFlag)
CHECK_Fortran_COMPILER_FLAG("-cpp" CPP_FLAG)
# if the flag exists, we add it to the compilation flags
if (CPP_FLAG)
	set(ADDITIONAL_GCC_FLAGS "-cpp")
endif()

# copied from another cmake file, maybe interesting later
#if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
#    set(bounds "-fbounds-check")
#endif()
#if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
#    set(dialect "-stand f08 -free -implicitnone")
#    set(bounds "-check bounds")
#endif()
#
#set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

# Add definitions. These should later be set via the cache line file and only
# have a defauld value here. These part only concerns user build definitions like MRM2MHM
add_definitions(-DMRM2MHM=.true.)
#add_definitions(-DpgiFortran=.false.)
#add_definitions(-DMPR_STANDALONE=.false.)
#add_definitions(-DABSOFT=.false.)
#add_definitions(-D=.false.)

# Compile.
#
# This is not recommended but wished by some members of our working group. With this command
# all files in src with the ending f90 or h are written into sources, and therefore linked
# together to the executable, if relevant or not
file(GLOB_RECURSE sources  src/*.f90 src/*.h)
# this command is able to create dependencies, compile and add the sources in the right order
add_executable(mhm ${sources})

# the libraries are added to the executable by the linker, using the full path and using the
# rpath option, except the libraries are located in ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}.
target_link_libraries(mhm ${NETCDF_LINK_LIBRARIES} ${MPI_Fortran_LIBRARIES})
set_property(TARGET mhm PROPERTY COMPILE_FLAGS "${CMAKE_Fortran_FLAGS} ${ADDITIONAL_GCC_FLAGS} ${NETCDF_CFLAGS_OTHER} ${MPI_Fortran_COMPILE_FLAGS}")
set_property(TARGET mhm PROPERTY LINK_FLAGS "${NETCDF_LDFLAGS_OTHER} ${MPI_Fortran_LINK_FLAGS}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -pedantic-errors -Wall -W -O -g -Wno-maybe-uninitialized")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_DEBUG} -O3")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form -ffixed-line-length-132")
message(STATUS "the following debug flags will be used: ${CMAKE_Fortran_FLAGS_DEBUG}")
# Usually that works fine, except, one is on a module system and tries to execute the executable
# in the end without having the modules loaded. A workaround is provided using the variable
# CMAKE_BUILD_MODULE_SYSTEM_INDEPENDEND
# if this variable is set to ON (do not set the variable inside of this cmake file), then the
# paths are added to the INSTALL_RPATH, and via the second command also to the build.
# It is a bit of a mess and workaround though.
if (CMAKE_BUILD_MODULE_SYSTEM_INDEPENDEND)
        set_target_properties(mhm PROPERTIES INSTALL_RPATH "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
        set_target_properties(mhm PROPERTIES BUILD_WITH_INSTALL_RPATH ON)
endif()
#set_target_properties(mhm PROPERTIES SKIP_BUILD_RPATH OFF)
#set_target_properties(mhm PROPERTIES INSTALL_RPATH_USE_LINKPATH ON)

# possibility to create symlinks to the build or to the source directory, so a copy of mhm does not have to be done. On the other
# hand, only for the test basins the executable would not be copied or linked to a predictable location. So it may be a good
# idea to not copy it with the cmake setup
#add_custom_target(link_namelists ALL "${CMAKE_COMMAND}" -E create_symlink "${CMAKE_CURRENT_SOURCE_DIR}/mhm.nml" "${CMAKE_CURRENT_BINARY_DIR}/mhm.nml")
