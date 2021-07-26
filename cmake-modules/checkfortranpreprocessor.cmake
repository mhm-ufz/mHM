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
