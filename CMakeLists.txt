# mHM cmake script
cmake_minimum_required(VERSION 3.12)
# additional cmake-modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

# get version and date from file
include(version)
get_version(MHM_VER MHM_VER_DEV MHM_DATE)

# create the project
project(mhm
	VERSION ${MHM_VER}
	DESCRIPTION "The mesoscale Hydrological Model"
	HOMEPAGE_URL "https://www.ufz.de/mhm"
	LANGUAGES Fortran
	)
# add full version to pre-processor flags (qoutes need in before hand)
add_compile_definitions(MHMVERSION='${MHM_VER_DEV}')
# add date to pre-processor flags (qoutes need in before hand)
add_compile_definitions(MHMDATE='${MHM_DATE}')

# add all compile options (MPI, OpenMP, Lapack, Coverage)
include(compileoptions)

# the FindNetCDFF.cmake file can be found after we added the cmake-modules folder to the CMAKE_MODULE_PATH
# the build fails, if it is not present
set(NETCDF_F90 "YES")
find_package(NetCDFF REQUIRED)
# from that module we gain the following variables:
# NETCDF_INCLUDES	: the include directory
# NETCDF_LINK_LIBRARIES : the absolute path to and with the libraries
# NETCDF_CFLAGS_OTHER	: additional compilation flags
# NETCDF_LDFLAGS_OTHER	: additional linking flags
include_directories(${NETCDF_INCLUDES} ${MPI_Fortran_INCLUDE_PATH} ${OpenMP_Fortran_LIBRARY})

# get pre-processor flag for current compiler (either -fpp or -cpp)
include(checkfortranpreprocessor)
get_preproc_flag(XPP_FLAG)

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
set_property(TARGET mhm PROPERTY COMPILE_FLAGS "${CMAKE_Fortran_FLAGS} ${XPP_FLAG} ${NETCDF_CFLAGS_OTHER} ${MPI_Fortran_COMPILE_FLAGS} ${OpenMP_Fortran_FLAGS}")
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
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -nofixed -assume byterecl -fp-model source -m64 -assume realloc-lhs ") # precise -> source: suppress warning, computation identical
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -warn all -g -debug -traceback -fp-stack-check -O0 -debug -check all")
	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -qoverride-limits")
	cpp_definitions("-DINTEL" "CMAKE_INTEL" "ON" "Code exchange for intel compiler dependent issues")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -colour -unsharedf95 -ideclient")
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g -nan -O0 -C=all -strict95 -ieee=stop")
	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O4 -ieee=full")
	cpp_definitions("-DNAG" "CMAKE_NAG" "ON" "Code exchange for NAG compiler dependent issues")
else()
	message(FATAL_ERROR "Unsupported Fortran Compiler: ${CMAKE_Fortran_COMPILER_ID}")
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
