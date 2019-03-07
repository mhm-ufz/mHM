# The cmake config approach of finding libraries
# This is the cleanest way. If the package itself provides a CMake config and
# defines the way of building the package, it is a good way to use the hints.
# On some systems that config file does not exist or cannot be found (strange setup module
# systems for instance).  In that case we fall back to another method. PkgConfig would be a
# feasible approach, see below. But we then go the nf-config way.
#****************************************************************

# The PkgConfig approach of finding libraries
# PkgConfig does not exist on many systems, so we do not use this approach
# Anyway, one can see quite well here how the variables may be called and what we need
#*************************************************************************************
#find_package(PkgConfig REQUIRED)
#pkg_check_modules(_NETCDF REQUIRED netcdf-fortran)
##find_path(NETCDF_INCLUDE_DIR HINTS ${_NETCDF_INCLUDE_DIRS})
#find_library(NETCDF_LIBRARIES NAMES ${_NETCDF_LIBRARIES} HINTS ${_NETCDF_LIBRARY_DIRS})
#set(NETCDF_CFLAGS_OTHER "${_NETCDF_CFLAGS_OTHER}" CACHE STRING "Additional compiler flags for NetCDF")
#set(NETCDF_LDFLAGS_OTHER "${_NETCDF_LDFLAGS_OTHER}" CACHE STRING "Additional link flags for NetCDF")
#set(NETCDF_INCLUDE_DIR "${_NETCDF_INCLUDE_DIRS}" CACHE STRING "Include directories for NetCDF")
##message(STATUS "NETCDF VARIABLES: libraries: ${NETCDF_LIBRARIES}, flags: ${NETCDF_FLAGS}, dirs: ${NETCDF_INCLUDE_DIR}")
#find_package_handle_standard_args(NETCDF REQUIRED_VARS NETCDF_LIBRARIES)

# The nf-config approach
# nf-config is a readable file and a program existent on any system, although, on MacOS everything
# is commented except an "echo "nf-config is not implemented yet""...
# executing "nf-config" writes out a description of every parameter of the file and "nf-config --all" prints how the
# parameters are set
# with these parameters we are able to set the library, include directory and flag variables accordingly.
#********************************************************************************************************

find_program(NETCDF_CONFIG nf-config
  HINTS NETCDF_DIR ENV NETCDF_DIR)
message(STATUS "found ${NETCDF_CONFIG}")
execute_process(COMMAND ${NETCDF_CONFIG} --includedir OUTPUT_VARIABLE NETCDF_INCLUDES OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "netcdff includes ${NETCDF_INCLUDES}")
execute_process(COMMAND ${NETCDF_CONFIG} --fflags OUTPUT_VARIABLE NETCDF_CFLAGS_OTHER OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "netcdff netcdf link flags ${NETCDF_CFLAGS_OTHER}")
execute_process(COMMAND ${NETCDF_CONFIG} --flibs OUTPUT_VARIABLE NETCDF_LDFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "netcdff netcdf library link flags ${NETCDF_LDFLAGS}")

string(REPLACE " " ";" NETCDF_LDFLAGS_LIST ${NETCDF_LDFLAGS})
foreach(flag ${NETCDF_LDFLAGS_LIST})
	if (flag MATCHES "^-L(.*)")
		list(APPEND _search_paths ${CMAKE_MATCH_1})
		continue()
	endif()
	if (flag MATCHES "^-l(.*)")
		set(_pkg_search "${CMAKE_MATCH_1}")
	else()
		list(APPEND _link_flags "${flag}")
		continue()
	endif()

	if(_search_paths)
		# Firstly search in -L paths
		find_library(pkgcfg_lib_NETCDF_${_pkg_search}
			NAMES ${_pkg_search}
			HINTS ${_search_paths} NO_DEFAULT_PATH)
	endif()
	find_library(pkgcfg_lib_NETCDF_${_pkg_search}
		NAMES ${_pkg_search}
		${_find_opts})
	list(APPEND _libs "${pkgcfg_lib_NETCDF_${_pkg_search}}")
endforeach()

set(NETCDF_LINK_LIBRARIES "${_libs}")
message(STATUS "found netcdf libraries ${NETCDF_LINK_LIBRARIES}")
set(NETCDF_LDFLAGS_OTHER "${_link_flags}")
message(STATUS "found netcdf other flags ${NETCDF_LDFLAGS_OTHER}")
