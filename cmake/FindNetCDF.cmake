# FindNetCDF.cmake - find NetCDF libraries
# ========================================
# safely find netcdf c/fortran
# To specify a particular NetCDF library, use
#
#     cmake -DNetCDF_ROOT=/path/to/netcdff -B build
#
# or set environment variable NetCDF_ROOT=/path/to/netcdff
#
# MIT License
# -----------
#[[
  Copyright (c) 2020 - 2022 CHS Developers

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
]]

# function to search for netcdf-C
function(search_netcdf_c)

  if(PkgConfig_FOUND AND NOT NetCDF_C_LIBRARY)
    pkg_search_module(pkg_conf_nc QUIET netcdf)
  endif()

  # include dirs
  find_path(NetCDF_C_INCLUDE_DIR
    NAMES netcdf.h
    HINTS ${pkg_conf_nc_INCLUDE_DIRS}
  )
  if(NOT NetCDF_C_INCLUDE_DIR)
    return()
  endif()

  # library
  find_library(NetCDF_C_LIBRARY
    NAMES netcdf
    HINTS ${pkg_conf_nc_LIBRARY_DIRS} ${pkg_conf_nc_LIBDIR}
  )
  if(NOT NetCDF_C_LIBRARY)
    return()
  endif()

  # set required variables in parent scope
  set(NetCDF_C_FOUND true PARENT_SCOPE)
  set(NetCDF_C_INCLUDE_DIR ${NetCDF_C_INCLUDE_DIR} PARENT_SCOPE)
  set(NetCDF_C_LIBRARY ${NetCDF_C_LIBRARY} PARENT_SCOPE)

endfunction(search_netcdf_c)

# function to search for netcdf-Fortran
function(search_netcdf_fortran)

  if(PkgConfig_FOUND AND NOT NetCDF_Fortran_LIBRARY)
    pkg_search_module(pkg_conf_nf QUIET netcdf-fortran netcdf)
  endif()

  # include dirs
  find_path(NetCDF_Fortran_INCLUDE_DIR
    names netcdf.mod
    HINTS ${pkg_conf_nf_INCLUDE_DIRS}
  )
  if(NOT NetCDF_Fortran_INCLUDE_DIR)
    return()
  endif()

  # library
  find_library(NetCDF_Fortran_LIBRARY
    NAMES netcdff
    HINTS ${pkg_conf_nf_LIBRARY_DIRS} ${pkg_conf_nf_LIBDIR}
  )
  if(NOT NetCDF_Fortran_LIBRARY)
    return()
  endif()

  # set required variables in parent scope
  set(NetCDF_Fortran_FOUND true PARENT_SCOPE)
  set(NetCDF_Fortran_INCLUDE_DIR ${NetCDF_Fortran_INCLUDE_DIR} PARENT_SCOPE)
  set(NetCDF_Fortran_LIBRARY ${NetCDF_Fortran_LIBRARY} PARENT_SCOPE)

endfunction(search_netcdf_fortran)

# search for PkgConfig
find_package(PkgConfig QUIET)

# try CMake built-in finder for NetCDF-c.
find_package(netCDF CONFIG QUIET)
if(netCDF_FOUND)
  set(NetCDF_C_FOUND "${netCDF_FOUND}")
  set(NetCDF_C_INCLUDE_DIR "${netCDF_INCLUDE_DIR}")
  set(NetCDF_C_LIBRARY "${netCDF_LIBRARIES}")
  set(NetCDF_VERSION "${NetCDFVersion}")
else()
  search_netcdf_c()
endif()
set(_nc_req_vars ${NetCDF_C_LIBRARY})

# search for netcdf-fortran if wanted
if(Fortran IN_LIST NetCDF_FIND_COMPONENTS)
  search_netcdf_fortran()
  list(APPEND _nc_req_vars ${NetCDF_Fortran_LIBRARY})
endif()

# hide the cached variables from cmake GUI
mark_as_advanced(
  NetCDF_C_INCLUDE_DIR
  NetCDF_Fortran_INCLUDE_DIR
  NetCDF_C_LIBRARY
  NetCDF_Fortran_LIBRARY
)

# check the requirements
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF
  REQUIRED_VARS _nc_req_vars
  HANDLE_COMPONENTS
)

# define imported targets
if(NetCDF_FOUND)
  set(NetCDF_C_INCLUDE_DIRS ${NetCDF_C_INCLUDE_DIR})
  set(NetCDF_C_LIBRARIES ${NetCDF_C_LIBRARY})
  add_library(NetCDF::NetCDF_C UNKNOWN IMPORTED)
  set_target_properties(NetCDF::NetCDF_C PROPERTIES
    IMPORTED_LOCATION "${NetCDF_C_LIBRARY}"
    INTERFACE_LINK_LIBRARIES "${NetCDF_C_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_C_INCLUDE_DIR}"
  )
  if(NetCDF_Fortran_FOUND)
    set(NetCDF_Fortran_INCLUDE_DIRS ${NetCDF_Fortran_INCLUDE_DIR})
    set(NetCDF_Fortran_LIBRARIES ${NetCDF_Fortran_LIBRARY})
    add_library(NetCDF::NetCDF_Fortran UNKNOWN IMPORTED)
    set_target_properties(NetCDF::NetCDF_Fortran PROPERTIES
      IMPORTED_LOCATION "${NetCDF_Fortran_LIBRARY}"
      INTERFACE_LINK_LIBRARIES "${NetCDF_Fortran_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_Fortran_INCLUDE_DIR}"
    )
  endif()
endif()
