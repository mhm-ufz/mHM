# FindYAC.cmake - find YAC coupler libraries
# ==========================================
# safely find YAC
# To specify a particular YAC library, use
#
#     -DYAC_ROOT=/path/to/yac
#
# or set environment variable YAC_ROOT=/path/to/yac
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

# function to search for yac
function(search_yac)

  if(PkgConfig_FOUND AND NOT YAC_LIB AND NOT YAC_INCLUDE_DIR)
    set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${YAC_ROOT}/lib/pkgconfig")
    pkg_check_modules(YAC REQUIRED yac)
  else()
    return()
  endif()

  # library
  if(CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT)
    execute_process(
      COMMAND bash "-c" "pkg-config --libs yac"
      OUTPUT_VARIABLE YAC_LIB_TEMP1
    )
    string(STRIP ${YAC_LIB_TEMP1} YAC_LIB_TEMP1)
    execute_process(
      COMMAND bash "-c" "pkg-config --libs-only-L yac | sed 's/-L/-Wl,-rpath,/g'"
      OUTPUT_VARIABLE YAC_LIB_TEMP2
    )
    string(STRIP ${YAC_LIB_TEMP2} YAC_LIB_TEMP2)
    string(CONCAT YAC_LIBRARY "${YAC_LIB_TEMP1}" " ${YAC_LIB_TEMP2}")
  else()
    set(YAC_LIBRARY ${YAC_LDFLAGS})
  endif()

  # set required variables in parent scope
  set(YAC_FOUND true PARENT_SCOPE)
  set(YAC_INCLUDE_DIR ${YAC_INCLUDE_DIRS} PARENT_SCOPE)
  set(YAC_LIB ${YAC_LIBRARY} PARENT_SCOPE)

endfunction(search_yac)

find_package(PkgConfig REQUIRED)
search_yac()
set(_yac_req_vars ${YAC_LIB})

mark_as_advanced (YAC_LIB YAC_INCLUDE_DIR)

# check the requirements
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(YAC
  REQUIRED_VARS _yac_req_vars
  HANDLE_COMPONENTS
)
