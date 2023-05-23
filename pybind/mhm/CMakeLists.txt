cmake_minimum_required(VERSION 3.17.2...3.24)
enable_language(C)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

install(TARGETS mhm_lib LIBRARY DESTINATION mhm)
if(BUILD_MHM_DRIVER)
  # link to mhm-lib next to executable
  if (APPLE)
    set_target_properties(mhm
      PROPERTIES
      INSTALL_RPATH "@loader_path"
      BUILD_RPATH "@loader_path"
    )
  else()
    set_target_properties(mhm
      PROPERTIES
      BUILD_RPATH_USE_ORIGIN ON
      INSTALL_RPATH "$ORIGIN"
      BUILD_RPATH "$ORIGIN"
    )
  endif()
  install(TARGETS mhm RUNTIME DESTINATION mhm)
endif()

# Grab the variables from a local Python installation F2PY headers
execute_process(
  COMMAND "${Python_EXECUTABLE}" -c
    "import os; from numpy import f2py; print(f2py.get_include() if hasattr(f2py, 'get_include') else os.path.join(os.path.dirname(f2py.__file__), 'src'))"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(f2py_module_name "wrapper")
set(f2py_module_name_c "${f2py_module_name}module.c")
set(f2py_module_name_f "${f2py_module_name}-f2pywrappers2.f90")
set(f2py_fortranobject "${F2PY_INCLUDE_DIR}/fortranobject.c")
set(fortran_src_file "${CMAKE_CURRENT_SOURCE_DIR}/${f2py_module_name}.f90")
add_custom_target(genpyf ALL DEPENDS ${fortran_src_file})

# numpy 1.17.3 (oldest-supported-numpy for py3.8) has no --f2cmap option in f2py
# so we copy the .f2py_f2cmap file to the binary directory where f2py is executed
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/.f2py_f2cmap" "${CMAKE_CURRENT_BINARY_DIR}/.f2py_f2cmap" COPYONLY)

# only use f2py to create c-wrapper for fortran module
add_custom_command(
  OUTPUT ${f2py_module_name_c} ${f2py_module_name_f}
  COMMAND ${Python_EXECUTABLE} -m "numpy.f2py"
    -m ${f2py_module_name}
    ${fortran_src_file}
    --lower
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${fortran_src_file}
)

# add python module for wrapper
python_add_library(${f2py_module_name}
  MODULE WITH_SOABI
  ${f2py_module_name_c}
  ${f2py_fortranobject}
  ${f2py_module_name_f}
  ${fortran_src_file}
)

# link mhm to wrapper
target_include_directories(${f2py_module_name} PUBLIC ${F2PY_INCLUDE_DIR})
target_link_libraries(${f2py_module_name} PUBLIC mhm_lib Python::NumPy)

# linker flags (rpath)
if (UNIX)
  if (APPLE)
    set_target_properties(${f2py_module_name}
      PROPERTIES
      LINK_FLAGS  "-Wl,-dylib,-undefined,dynamic_lookup"
      INSTALL_RPATH "@loader_path"
      BUILD_RPATH "@loader_path"
    )
  else()
    set_target_properties(${f2py_module_name}
      PROPERTIES
      LINK_FLAGS  "-Wl,--allow-shlib-undefined"
      BUILD_RPATH_USE_ORIGIN ON
      INSTALL_RPATH "$ORIGIN"
      BUILD_RPATH "$ORIGIN"
    )
  endif()
endif()
add_dependencies(${f2py_module_name} genpyf)

# install
install(TARGETS ${f2py_module_name} DESTINATION mhm)
