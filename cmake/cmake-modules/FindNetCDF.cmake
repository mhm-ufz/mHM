# taken directly from https://github.com/WRF-CMake/wrf/blob/wrf-cmake/cmake/FindNetCDF.cmake on 2020/05/06
# and added "$ENV{..}/{..}" to HINTS args in multiple find_* calls
# https://github.com/Kitware/VTK/blob/master/CMake/FindNetCDF.cmake
# =========================================================================
#
# Copyright (c) 1993-2015 Ken Martin, Will Schroeder, Bill Lorensen
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
#  * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
#    of any contributors may be used to endorse or promote products derived
#    from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# =========================================================================
#
# - Find NetCDF
# Find the native NetCDF includes and library
#
#  NETCDF_INCLUDE_DIR  - user modifiable choice of where netcdf headers are
#  NETCDF_LIBRARY      - user modifiable choice of where netcdf libraries are
#
# Your package can require certain interfaces to be FOUND by setting these
#
#  NETCDF_F77         - require the F77 interface and link the fortran library
#  NETCDF_F90         - require the F90 interface and link the fortran library
#
# Or equivalently by calling FindNetCDF with a COMPONENTS argument containing one or
# more of "CXX;F77;F90".
#
# When interfaces are requested the user has access to interface specific hints:
#
#  NETCDF_${LANG}_INCLUDE_DIR - where to search for interface header files
#  NETCDF_${LANG}_LIBRARY     - where to search for interface libraries
#
# This module returns these variables for the rest of the project to use.
#
#  NETCDF_FOUND          - True if NetCDF found including required interfaces (see below)
#  NETCDF_LIBRARIES      - All netcdf related libraries.
#  NETCDF_INCLUDE_DIRS   - All directories to include.
#  NETCDF_HAS_INTERFACES - Whether requested interfaces were found or not.
#  NETCDF_${LANG}_INCLUDE_DIRS/NETCDF_${LANG}_LIBRARIES - C/F70/F90 only interface
#
# Normal usage would be:
#  set (NETCDF_F90 "YES")
#  find_package (NetCDF REQUIRED)
#  target_link_libraries (uses_everthing ${NETCDF_LIBRARIES})
#  target_link_libraries (only_uses_f90 ${NETCDF_F90_LIBRARIES})

#search starting from user editable cache var
if (NETCDF_INCLUDE_DIR AND NETCDF_LIBRARY)
  # Already in cache, be silent
  set (NETCDF_FIND_QUIETLY TRUE)
endif ()

set(USE_DEFAULT_PATHS "")
if(NETCDF_NO_DEFAULT_PATH)
	set(USE_DEFAULT_PATHS "NO_DEFAULT_PATH")
endif()

find_path (NETCDF_INCLUDE_DIR netcdf.h
  HINTS "${NETCDF_DIR}/include" "$ENV{NETCDF_DIR}/include")
mark_as_advanced (NETCDF_INCLUDE_DIR)
set (NETCDF_C_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR})

find_library (NETCDF_LIBRARY NAMES netcdf
  HINTS "${NETCDF_DIR}/lib" "$ENV{NETCDF_DIR}/lib"
  PATH_SUFFIXES "x86_64-linux-gnu"
  )
mark_as_advanced (NETCDF_LIBRARY)

set (NETCDF_C_LIBRARIES ${NETCDF_LIBRARY})

#start finding requested language components
set (NetCDF_libs "")
set (NetCDF_includes "${NETCDF_INCLUDE_DIR}")

get_filename_component (NetCDF_lib_dirs "${NETCDF_LIBRARY}" PATH)
set (NETCDF_HAS_INTERFACES "YES") # will be set to NO if we're missing any interfaces

macro (NetCDF_check_interface lang header libs)
  if (NETCDF_${lang})
    #search starting from user modifiable cache var
    find_path (NETCDF_${lang}_INCLUDE_DIR NAMES ${header}
      HINTS "${NETCDF_INCLUDE_DIR}"
            "${NETCDF_FORTRAN_DIR}/include"
            "$ENV{NETCDF_FORTRAN_DIR}/include"
      ${USE_DEFAULT_PATHS})

    find_library (NETCDF_${lang}_LIBRARY NAMES ${libs}
      HINTS "${NetCDF_lib_dirs}"
            "${NETCDF_FORTRAN_DIR}/lib"
            "$ENV{NETCDF_FORTRAN_DIR}/lib"
      PATH_SUFFIXES "x86_64-linux-gnu"
      ${USE_DEFAULT_PATHS})

    mark_as_advanced (NETCDF_${lang}_INCLUDE_DIR NETCDF_${lang}_LIBRARY)

    #export to internal varS that rest of project can use directly
    set (NETCDF_${lang}_LIBRARIES ${NETCDF_${lang}_LIBRARY})
    set (NETCDF_${lang}_INCLUDE_DIRS ${NETCDF_${lang}_INCLUDE_DIR})

    if (NETCDF_${lang}_INCLUDE_DIR AND NETCDF_${lang}_LIBRARY)
      list (APPEND NetCDF_libs ${NETCDF_${lang}_LIBRARY})
      list (APPEND NetCDF_includes ${NETCDF_${lang}_INCLUDE_DIR})
    else ()
      if ( ${lang} STREQUAL "F90")
      	find_program(NETCDFF_CONFIG nf-config
          HINTS ${CMAKE_NETCDF_DIR})
        message(STATUS "found ${NETCDFF_CONFIG}")
        execute_process(COMMAND ${NETCDFF_CONFIG} --includedir OUTPUT_VARIABLE NETCDF_${lang}_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "netcdff includes ${NETCDF_INCLUDES}")
        execute_process(COMMAND ${NETCDFF_CONFIG} --fflags OUTPUT_VARIABLE NETCDF_CFLAGS_OTHER OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "netcdff netcdf link flags ${NETCDF_CFLAGS_OTHER}")
        execute_process(COMMAND ${NETCDFF_CONFIG} --flibs OUTPUT_VARIABLE NETCDF_LDFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "netcdff netcdf library link flags ${NETCDF_LDFLAGS}")

        if (CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT)
          find_program(NETCDF_CONFIG nc-config
            HINTS ${CMAKE_NETCDF_DIR})
          execute_process(COMMAND ${NETCDF_CONFIG} --libs OUTPUT_VARIABLE NETCDF_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
          message(STATUS "netcdf library link flags ${NETCDF_LIBS}")
        endif()

        # In a clean cmake setup the libraries are included via the target_link_libraries and not
        # via flags. Cmake creates system dependend flags and rpaths using the libraries itself.
        # nf-config on the other hand gives us a list of flags, linking with -l and -L.
        # we cut the flag string into seperated flags and create libraries and other flags from it
        string(REPLACE " " ";" NETCDF_LDFLAGS_LIST "${NETCDF_LDFLAGS} ${NETCDF_LIBS}")
        foreach(flag ${NETCDF_LDFLAGS_LIST})
          # message(STATUS "${flag}")
          if (flag MATCHES "^-L(.*)")
            list(APPEND _search_paths ${CMAKE_MATCH_1})
            continue()
          endif()
          if (flag MATCHES "^-l(.*)")
            set(_pkg_search "${CMAKE_MATCH_1}")
          else()
            string(CONCAT _link_flags "${_link_flags}" " " "${flag}")
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
            HINTS ENV LD_LIBRARY_PATH)
          message(STATUS "found ${pkgcfg_lib_NETCDF_${_pkg_search}}")
          list(APPEND _libs "${pkgcfg_lib_NETCDF_${_pkg_search}}")
        endforeach()

        # those are set one level up
        set(NETCDF_LDFLAGS_OTHER "${_link_flags}")
        message(STATUS "found netcdf other flags ${NETCDF_LDFLAGS_OTHER}")

        list (APPEND NETCDF_LINK_LIBRARIES "${_libs}")
        message(STATUS "found netcdf libraries ${NETCDF_LINK_LIBRARIES}")
        list (APPEND NetCDF_includes ${NETCDF_${lang}_INCLUDE_DIR})
      else ()
        set (NETCDF_HAS_INTERFACES "NO")
        message (STATUS "Failed to find NetCDF interface for ${lang}")
      endif ()
    endif ()
  endif ()
endmacro ()

list (FIND NetCDF_FIND_COMPONENTS "F77" _nextcomp)
if (_nextcomp GREATER -1)
  set (NETCDF_F77 1)
endif ()
list (FIND NetCDF_FIND_COMPONENTS "F90" _nextcomp)
if (_nextcomp GREATER -1)
  set (NETCDF_F90 1)
endif ()
NetCDF_check_interface (F77 netcdf.inc  netcdff)
NetCDF_check_interface (F90 netcdf.mod  netcdff)

#export accumulated results to internal varS that rest of project can depend on
list (APPEND NetCDF_libs "${NETCDF_C_LIBRARIES}")
set (NETCDF_LIBRARIES ${NetCDF_libs})
set (NETCDF_INCLUDE_DIRS ${NetCDF_includes})

# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NetCDF
  DEFAULT_MSG NETCDF_LIBRARIES NETCDF_INCLUDE_DIRS NETCDF_HAS_INTERFACES)