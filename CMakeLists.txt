cmake_minimum_required(VERSION 3.5)

project(mhm Fortran)

#set(BUILD_WITH_INSTALL_RPATH true)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)
#sometimes PkgConfig is not on systems like mac
#there is a way to make it conditional, if it is there use this
#find_package(PkgConfig REQUIRED)
# Find the native NetCDF includes and library
#include(cmake-modules/FindNetCDF.cmake)
# PkgConfig spezifisch
#pkg_check_modules(NetCDF REQUIRED)

set (NETCDF_F90 "YES")
find_package(NetCDFF REQUIRED)
include_directories(${NETCDF_INCLUDES})

#include(CHECK_Fortran_COMPILER_FLAG)
#include(CheckFortranSourceCompiles)
#include(CMakeCheckCompilerFlagCommonPatterns)
include(CheckFortranCompilerFlag)

CHECK_Fortran_COMPILER_FLAG("-cpp" CPP_FLAG)
if (${CPP_FLAG})
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
set(ADDITIONAL_GCC_FLAGS "-cpp")
endif(${CPP_FLAG})

# Options
#add_definitions(${ADDITIONAL_GCC_FLAGS})
add_definitions(-DMRM2MHM=.true.)
#add_definitions(${NETCDF_CFLAGS_OTHER})
#add_definitions(-DpgiFortran=.false.)
#add_definitions(-DMPR_STANDALONE=.false.)
#add_definitions(-DABSOFT=.false.)
#add_definitions(-D=.false.)

#set(pgiFortran OFF)
#configure_file(config.h.in config.h)

#include_directories(${HDF5_INCLUDE_DIRS})
#include_directories(${NETCDF_LINK_LIBRARIES})
#
# Compile.
#
file(GLOB_RECURSE sources  src/*.f90 src/*.h)
add_executable(mhm ${sources})

target_link_libraries(mhm ${NETCDF_LINK_LIBRARIES})
set_property(TARGET mhm PROPERTY COMPILE_FLAGS "${ADDITIONAL_GCC_FLAGS} ${NETCDF_CFLAGS_OTHER}")
set_property(TARGET mhm PROPERTY LINK_FLAGS "${NETCDF_LDFLAGS_OTHER}")
set_target_properties(mhm PROPERTIES INSTALL_RPATH "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
set_target_properties(mhm PROPERTIES BUILD_WITH_INSTALL_RPATH ON)
#set_target_properties(mhm PROPERTIES INSTALL_RPATH $ENV{LD_LIBRARY_PATH})
#message(STATUS "rpath variable $ENV{LD_LIBRARY_PATH}")
#set_target_properties(mhm PROPERTIES SKIP_BUILD_RPATH OFF)
#set_target_properties(mhm PROPERTIES INSTALL_RPATH_USE_LINKPATH ON)
#set_target_properties(mhm PROPERTIES BUILT_WITH_INSTALL_RPATH ON)
#target_link_libraries (mhm ${NETCDF_LIBRARIES_C})

add_custom_target(link_namelists ALL "${CMAKE_COMMAND}" -E create_symlink "${CMAKE_CURRENT_SOURCE_DIR}/mhm.nml" "${CMAKE_CURRENT_BINARY_DIR}/mhm.nml")
