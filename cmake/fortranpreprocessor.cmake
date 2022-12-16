# fortranpreprocessor.cmake - pre-processor helpers
# =================================================
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

# get_preproc_flag : get correct fortran pre-processor flag
#
# get_preproc_flag(FLAG)
#
# will store pre-processor flag in "FLAG" (-fpp or -cpp)
function(get_preproc_flag preprog_flag)
  include(CheckFortranCompilerFlag)
  # check pre-processor -fpp for NAG or Intel
  CHECK_Fortran_COMPILER_FLAG("-fpp" FPP_FLAG)
  # check -cpp for other compilers
  CHECK_Fortran_COMPILER_FLAG("-cpp" CPP_FLAG)
  # if the flag exists, we add it to the compilation flags
  if(FPP_FLAG)
    set(preprog_flag_ "-fpp")
  elseif(CPP_FLAG)
    set(preprog_flag_ "-cpp")
  else()
    message(FATAL_ERROR "Compiler does not support -fpp or -cpp")
  endif()
  set(${preprog_flag} ${preprog_flag_} PARENT_SCOPE)
  message(STATUS "Pre-Processor Flag: '${preprog_flag_}'")
endfunction()

# this function adds definitions but also creates a corresponding CMAKE_* variable with CACHE STRING
# $cmake -C ../CMakeCacheFiles/example
# adding cpp_definitions(DEFNAME OFF "some info")
# will enable the pre-processor directive DEFNAME and can be controlled by
# cmake -DCMAKE_DEFNAME=ON ..
function(cpp_definitions defName value cacheString)
  option(CMAKE_${defName} ${cacheString} ${value})
  if (CMAKE_${defName})
    add_compile_definitions(${defName})
  endif()
endfunction()

# this code is based on from https://github.com/fortran-lang/stdlib/blob/master
# MIT License
#
# Copyright (c) 2019 Fortran stdlib developers
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Preprocesses a list of files with given preprocessor and preprocessor options
#
# Args:
#     preproc [in]: Preprocessor program
#     preprocopts [in]: Preprocessor options
#     srcext [in]: File extension of the source files
#     trgext [in]: File extension of the target files
#     srcfiles [in]: List of the source files
#     trgfiles [out]: Contains the list of the preprocessed files on exit
#
function(preprocess preproc preprocopts srcext trgext srcfiles trgfiles)

  set(_trgfiles)
  foreach(srcfile IN LISTS srcfiles)
    string(REGEX REPLACE "\\.${srcext}$" ".${trgext}" trgfile ${srcfile})
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${trgfile}
      COMMAND ${preproc} ${preprocopts} ${CMAKE_CURRENT_SOURCE_DIR}/${srcfile} ${CMAKE_CURRENT_BINARY_DIR}/${trgfile}
      MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${srcfile})
    list(APPEND _trgfiles ${CMAKE_CURRENT_BINARY_DIR}/${trgfile})
  endforeach()
  set(${trgfiles} ${_trgfiles} PARENT_SCOPE)

endfunction()

# Preprocesses fortran files with fypp.
#
# It assumes that source files have the ".fypp" extension. Target files will be
# created the extension ".f90". The FYPP variable must contain the path to the
# fypp-preprocessor.
#
# Args:
#     fyppopts [in]: Options to pass to fypp.
#     fyppfiles [in]: Files to be processed by fypp
#     f90files [out]: List of created f90 files on exit
#
function (fypp_f90 fyppopts fyppfiles f90files)
  preprocess("${FYPP}" "${fyppopts}" "fypp" "f90" "${fyppfiles}" _f90files)
  set(${f90files} ${_f90files} PARENT_SCOPE)
endfunction()