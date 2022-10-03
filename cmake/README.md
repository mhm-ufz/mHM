# Cmake Fortran Scripts

Cmake scripts for Fortran Projects at CHS

## Cmake Modules

Modules under `cmake-modules` provide different functions and macros
to be used in `CMakeLists.txt` files.

### `fortranpreprocessor.cmake`
Based on: https://github.com/fortran-lang/stdlib/blob/master/cmake/stdlib.cmake
- handles `fypp` pre-processor

Additional functions:
- `get_preproc_flag`: checks for the pre-processor flag of the current compiler (`-fpp` or `-cpp`)
- `cpp_definitions`: adds a compile definition and a corresponding cmake option

Can be included with:
```cmake
include(fortranpreprocessor)
```

### `CodeCoverage.cmake`
Taken from: https://github.com/bilke/cmake-modules/blob/master/CodeCoverage.cmake

Can be included and set up with:
```cmake
include(CodeCoverage)
append_coverage_compiler_flags_to_target(PROJECT)
SETUP_TARGET_FOR_COVERAGE_LCOV(
  NAME PROJECT_coverage_CI
  EXECUTABLE PROJECT
  DEPENDENCIES PROJECT
  GENHTML_ARGS -t "PROJECT coverage"
)
```

### `compileoptions.cmake`
Adds cmake compile options
- `CMAKE_BUILD_MODULE_SYSTEM_INDEPENDENT`: setting r_path
- `CMAKE_WITH_MPI`: use MPI
- `CMAKE_WITH_OpenMP`: use OpenMP
- `CMAKE_WITH_LAPACK`: use LAPACK bindings
- `CMAKE_WITH_COVERAGE`: Coverage calculation
- `CMAKE_WITH_GPROF`: gprof profiling information (see [here](https://www.thegeekstuff.com/2012/08/gprof-tutorial/) for tutorial)
- will also search for `MPI`, `OpenMP` and `LAPACK` if the option are `ON`

Can be included with:
```cmake
include(compileoptions)
```

### `FindNetCDF.cmake`

Can be used like:
```cmake
find_package(NetCDF COMPONENTS Fortran)
target_link_libraries(<target> PUBLIC NetCDF::NetCDF_Fortran)
```

To specify a particular NetCDF library, use

    cmake -DNetCDF_ROOT=/path/to/netcdff -B build

or set environment variable `NetCDF_ROOT=/path/to/netcdff`.

### `version.cmake`

Provides a function to read version and date from given files `version.txt` and (optional) `version_date.txt`.

Can be included and used with (`PROJECT_VER_DEV` will hold the develop version string):
```cmake
include(version)
get_version(PROJECT_VER PROJECT_VER_DEV PROJECT_DATE)
```

If a development version was found (e.g. `1.0.0-dev0`), it will try to determine the latest git commit short hash and add it (e.g. `1.0.0-dev0+02d7e12`).
If there is no `version_date.txt` file or a development version was found, the script will try to determine the latest git commit date or use the current date.
`PROJECT_VER` will hold the simple version to be used with cmakes `project` command (e.g. `1.0.0`).
A prefixed `v` in the version will be ignored.

### `CPM.cmake` (v0.35.6)

CPM.cmake is a CMake script that adds dependency management capabilities to CMake.
It's built as a thin wrapper around CMake's [FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html)
module that adds version control, caching, a simple API [and more](#comparison-to-pure-fetchcontent--externalproject).

Copied from: https://github.com/cpm-cmake/CPM.cmake

Can be used like:
```cmake
include(CPM)
CPMAddPackage("https://git.ufz.de/chs/forces.git@0.3.2")
target_link_libraries(<project> PUBLIC forces)
```


## Cmake Cache Files

Cache files provide templates for specific cmake setups.
Cmake can be controlled with such a file (`exampleFile`) in the following way
```bash
cmake -C exampleFile ..
```

Have a look at the example files in `cmake-cache-files` taken from `mHM`.
