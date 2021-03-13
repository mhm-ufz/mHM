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
    add_compile_definitions("${defName}")
  endif()
endfunction()
