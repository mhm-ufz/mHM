# This cmake file is not meant to be edited for a special setup.
# For special setups use cache line files or command line options, as described a few
# lines ahead concerning module INDEPENDENT builds
cmake_minimum_required(VERSION 3.5)
# check version file
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/version.txt")
	file(STRINGS "version.txt" MHMVERSION LIMIT_COUNT 1)
else()
	set(MHMVERSION "0.0.0-dev0") # default version
endif()

# version should be of the form (semver.org):
# - 1.2.3-dev0 (development with number)
# - 1.2.3-rc1  (release candidate with number)
# - 1.2.3      (release)
# remove possible "v" prefix and find major.minor.patch version
string(REGEX MATCH "^v?([0-9]+)" _ ${MHMVERSION})
set(ver_major ${CMAKE_MATCH_1})
string(REGEX MATCH "^v?[0-9]+\.([0-9]+)" _ ${MHMVERSION})
set(ver_minor ${CMAKE_MATCH_1})
string(REGEX MATCH "^v?[0-9]+\.[0-9]+\.([0-9]+)" _ ${MHMVERSION})
set(ver_patch ${CMAKE_MATCH_1})
# find pre-release tag
string(REGEX MATCH ".*-(.+)" _ ${MHMVERSION})
set(ver_pre ${CMAKE_MATCH_1})

# create the version string for cmake (fill up with 0)
if ("${ver_major}" STREQUAL "")
	set(ver_major 0) # default version
endif()
if ("${ver_minor}" STREQUAL "")
	set(ver_minor 0) # default version
endif()
if ("${ver_patch}" STREQUAL "")
	set(ver_patch 0) # default version
endif()
set(ver_final ${ver_major}.${ver_minor}.${ver_patch}) # full version

# check date file (if not a development version)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/version_date.txt" AND (NOT (${ver_pre} MATCHES "^dev.*")))
	file(STRINGS "version_date.txt" MHMDATE LIMIT_COUNT 1)
else()
	string(TIMESTAMP MHMDATE "%Y-%m-%d") # current date
endif()

# create the project
project(mhm
        VERSION ${ver_final}
        DESCRIPTION "The mesoscale Hydrological Model"
        HOMEPAGE_URL "https://www.ufz.de/mhm"
        LANGUAGES Fortran)
message("mhm VERSION: ${mhm_VERSION} (from ${MHMVERSION})")
message("mhm DATE:    ${MHMDATE}")
# add version to pre-processor flags (qoutes need in before hand)
add_compile_definitions(MHMVERSION=\"${MHMVERSION}\")
# add date to pre-processor flags (qoutes need in before hand)
add_compile_definitions(MHMDATE=\"${MHMDATE}\")

# The variable "CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT" can be set before executing cmake via a cache command:
# $cmake -DCMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT:STRING=ON ..
# or cache file:
# $cmake -C ../CMakeCacheFiles/eve ..
# or after executing CMake editing the CMakeCache.txt, preferably with a corresponding cmake editor i.e ccmake
set(CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT OFF CACHE STRING "build the module INDEPENDENT of the module system, so the build in the build tree works even after a module purge")
message(STATUS "build INDEPENDENT of module system ${CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT}")

# set specific place where to search for the netCDF directory
set(CMAKE_NETCDF_DIR " " CACHE STRING "set set specific place where to search for the netCDF directory")
message(STATUS "search in additional directory ${CMAKE_NETCDF_DIR} for netCDF")

# The variable "CMAKE_WITH_MPI" can be set before executing cmake via a cache command:
# $cmake -DCMAKE_WITH_MPI:STRING=ON ..
# or in a cache file:
# $cmake -C ../CMakeCacheFiles/example
# or after executing CMake editing the CMakeCache.txt, preferably with a corresponding cmake editor i.e. ccmake
set(CMAKE_WITH_MPI OFF CACHE STRING "build the module with MPI, so it can be executed using mpirun")
# same with OpenMP
set(CMAKE_WITH_OpenMP OFF CACHE STRING "build the module with OpenMP parallelization")
# same with lapack
set(CMAKE_WITH_LAPACK OFF CACHE STRING "build the module with lapack library")
# same with coverage
set(CMAKE_WITH_COVERAGE OFF CACHE STRING "build the module with gcov coverage support")

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

# if cmake provides a findLIBRARY module, this gets invoked via find_package(LIBRARY)
if (CMAKE_WITH_MPI)
        # find if there is an MPI setup on the system and if so, set corresponding variables
        find_package(MPI)
        if (NOT ${MPI_Fortran_FOUND})
                message(FATAL_ERROR "MPI required but not found")
        else()
                message(STATUS "found MPI_Fortran_COMPILER ${MPI_Fortran_COMPILER}")
        endif()
	add_definitions("-DMPI")
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

include_directories(${NETCDF_INCLUDES} ${MPI_Fortran_INCLUDE_PATH} ${OpenMP_Fortran_LIBRARY})

# ifort and gfortran need the flag -cpp to interpret pre-processor directives
# the nag compiler is not able to interpret the flag -cpp but can interpret these definitions anyway
# so we check whether the compiler is able to use the flag -cpp
# for that we need the module CheckFortranCompilerFlag
include(CheckFortranCompilerFlag)
CHECK_Fortran_COMPILER_FLAG("-cpp" CPP_FLAG)
# if the flag exists, we add it to the compilation flags
if (CPP_FLAG)
	set(ADDITIONAL_GCC_FLAGS "-cpp")
endif()

# this function adds definitions but also creates a corresponding CMAKE variable with CACHE STRING
# i.e.:
# The variable "${defCMakeName}" can be set before executing cmake via a cache command cmake -D...
# or in a cache file:
# $cmake -C ../CMakeCacheFiles/example
# or after executing CMake editing the CMakeCache.txt, preferably with a corresponding cmake editor i.e. ccmake
# cmake ..
function(cpp_definitions defName defCMakeName value cacheString)
	set(${defCMakeName} "${value}" CACHE STRING "${cacheString}")
	if (${defCMakeName})
		add_definitions("${defName}")
	endif()
endfunction()

# add pre-processor -fpp for NAG
CHECK_Fortran_COMPILER_FLAG("-fpp" FPP_FLAG)
# if the flag exists, we add it to the compilation flags
if (FPP_FLAG)
	set(ADDITIONAL_GCC_FLAGS "-fpp")
endif()

# Add definitions. These should later be set via the cache line file and only have a default value here.
cpp_definitions("-DpgiFortran" "CMAKE_pgiFortran" "OFF" "Code exchange for pgi compiler dependent issues")
cpp_definitions("-DMPR_STANDALONE" "CMAKE_MPR_STANDALONE" "OFF" "If set to ON, only MPR is compiled")
cpp_definitions("-DABSOFT" "CMAKE_ABSOFT" "OFF" "Documentation to be added. If you you are developer, you might edit this string in CMakeLists.txt")

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
target_link_libraries(mhm ${NETCDF_LINK_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${OpenMP_Fortran_LIBRARIES} ${LAPACK_LIBRARIES})
set_property(TARGET mhm PROPERTY COMPILE_FLAGS "${CMAKE_Fortran_FLAGS} ${ADDITIONAL_GCC_FLAGS} ${NETCDF_CFLAGS_OTHER} ${MPI_Fortran_COMPILE_FLAGS} ${OpenMP_Fortran_FLAGS}")
set_property(TARGET mhm PROPERTY LINK_FLAGS "${NETCDF_LDFLAGS_OTHER} ${MPI_Fortran_LINK_FLAGS} ${OpenMP_Fortran_FLAGS} ${LAPACK_LINKER_FLAGS}")
# set compiling flags for debug and realese version
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form -ffixed-line-length-132")
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -pedantic-errors -Wall -W -O -g -Wno-maybe-uninitialized")
	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3")
	cpp_definitions("-DGFORTRAN" "CMAKE_GFORTRAN" "ON" "Code exchange for gfortran compiler dependent issues")
	if (CMAKE_WITH_COVERAGE)
		include(CodeCoverage)
		APPEND_COVERAGE_COMPILER_FLAGS()
		SETUP_TARGET_FOR_COVERAGE_LCOV(NAME mhm_coverage_CI
			EXECUTABLE ../CI-scripts/run_cmake_coverage.sh
			DEPENDENCIES mhm
			EXCLUDE src/lib/*
			GENHTML_ARGS -t "mHM coverage" --html-prolog ../doc/html_files/cov_header.prolog)
	endif()
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -nofixed -assume byterecl -fp-model source -m64 -assume realloc-lhs ") # precise -> source: suppress warning, computation identical
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -warn all -g -debug -traceback -fp-stack-check -O0 -debug -check all")
	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -qoverride-limits")
	cpp_definitions("-DINTEL" "CMAKE_INTEL" "ON" "Code exchange for intel compiler dependent issues")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -colour -unsharedf95 -ideclient")
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g -nan -O0 -C=all -strict95 -ieee=stop")
	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O4 -ieee=full")
	cpp_definitions("-DNAG" "CMAKE_NAG" "ON" "Code exchange for NAG compiler dependent issues")
endif()
message(STATUS "the following debug flags will be used: ${CMAKE_Fortran_FLAGS_DEBUG}")
# Usually that works fine, except, one is on a module system and tries to execute the executable
# in the end without having the modules loaded. A workaround is provided using the variable
# CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT
# if this variable is set to ON (do not set the variable inside of this cmake file), then the
# paths are added to the INSTALL_RPATH, and via the second command also to the build.
# It is a bit of a mess and workaround though.
if (CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT)
        set_target_properties(mhm PROPERTIES INSTALL_RPATH "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
        set_target_properties(mhm PROPERTIES BUILD_WITH_INSTALL_RPATH ON)
endif()
#set_target_properties(mhm PROPERTIES SKIP_BUILD_RPATH OFF)
#set_target_properties(mhm PROPERTIES INSTALL_RPATH_USE_LINKPATH ON)

# possibility to create symlinks to the build or to the source directory, so a copy of mhm does not have to be done. On the other
# hand, only for the test domains the executable would not be copied or linked to a predictable location. So it may be a good
# idea to not copy it with the cmake setup
#add_custom_target(link_namelists ALL "${CMAKE_COMMAND}" -E create_symlink "${CMAKE_CURRENT_SOURCE_DIR}/mhm.nml" "${CMAKE_CURRENT_BINARY_DIR}/mhm.nml")
