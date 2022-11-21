enable_language(C)
find_package(PythonExtensions REQUIRED)
find_package(NumPy REQUIRED)
find_package(F2PY REQUIRED)
# rpath should be set to look locally for the needed .so file
if (APPLE)
  set(rpath "@loader_path")
  set_target_properties(mhm
    PROPERTIES
    INSTALL_RPATH "@loader_path"
    BUILD_RPATH "@loader_path"
  )
else()
  set(rpath "$$ORIGIN")  # need double $ to not expand the variable
  set_target_properties(mhm
    PROPERTIES
    BUILD_RPATH_USE_ORIGIN ON
    INSTALL_RPATH "$ORIGIN"
    BUILD_RPATH "$ORIGIN"
  )
endif()

install(TARGETS mhm_lib LIBRARY DESTINATION .)
install(TARGETS mhm RUNTIME DESTINATION .)

get_target_property(LIB_INCLUDES mhm_lib INCLUDE_DIRECTORIES)
message("includes: ${LIB_INCLUDES}")

set(f2py_module_name "wrapper")
set(fortran_src_file "${CMAKE_CURRENT_SOURCE_DIR}/${f2py_module_name}.f90")
set(generated_module_file ${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
add_custom_target(${f2py_module_name} ALL DEPENDS ${generated_module_file} mhm_lib)

add_custom_command(
  OUTPUT ${generated_module_file}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  # use Command Line Tool Mode of CMake to set LDFLAGS for f2py (rpath)
  COMMAND ${CMAKE_COMMAND} -E env 'LDFLAGS=-Wl,-rpath,${rpath}' env 'NPY_DISTUTILS_APPEND_FLAGS=1'
    ${F2PY_EXECUTABLE}
    -m ${f2py_module_name}
    -c ${fortran_src_file}
    -I${LIB_INCLUDES}
    -L${LIB_INCLUDES}
    -lmhm_lib
    --f90exec=${CMAKE_Fortran_COMPILER}
)
install(FILES ${generated_module_file} DESTINATION .)
