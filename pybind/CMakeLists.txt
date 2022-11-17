cmake_minimum_required(VERSION 3.14.0)

project(mhm_pybind)
find_package(PythonExtensions REQUIRED)
find_package(NumPy REQUIRED)
find_package(F2PY REQUIRED)

add_subdirectory(src/mhm)
