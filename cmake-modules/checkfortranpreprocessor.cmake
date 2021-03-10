# get_preproc_flag : get correct fortran pre-processor flag
#
# get_preproc_flag(FLAG)
#
# will store pre-processor flag in "FLAG" (-fpp or -cpp)
function(get_preproc_flag preprog_flag)
  include(CheckFortranCompilerFlag)
  # add pre-processor -fpp for NAG or Intel
  CHECK_Fortran_COMPILER_FLAG("-fpp" FPP_FLAG)
  # if the flag exists, we add it to the compilation flags
  if (FPP_FLAG)
    set(${preprog_flag} "-fpp" PARENT_SCOPE)
  else()
    # check -cpp for other compilers
    CHECK_Fortran_COMPILER_FLAG("-cpp" CPP_FLAG)
    # if the flag exists, we add it to the compilation flags
    if (CPP_FLAG)
      set(${preprog_flag} "-cpp" PARENT_SCOPE)
    else()
      message(FATAL_ERROR "Compiler does not support -fpp or -cpp")
    endif()
  endif()
  message(STATUS "Pre-Processor Flag: '${preprog_flag}'")
endfunction()
