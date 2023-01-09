# Cmake Fortran Scripts

Cmake scripts for Fortran Projects at CHS

## Cmake Modules

These modules provide different functions and macros to be used in `CMakeLists.txt` files.

To use a module, you either need to update the cmake module path (assuming these files are in a top-level folder called `cmake/`):

    list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

and then include the module by name, e.g.:

    include(CPM)

or include the module directly by path (relative to the `CMakeLists.txt` file)

    include(cmake/CPM.cmake)

For package finders (e.g. YAC or NetCDF) you **need** to update the cmake module path.

### `fortranpreprocessor.cmake`

Draws inspiration from https://github.com/fortran-lang/stdlib/blob/master/cmake/stdlib.cmake (MIT license)
to handle `fypp` pre-processor.

Additional functions:
- `get_preproc_flag`: checks for the pre-processor flag of the current compiler (`-fpp` or `-cpp`)
- `cpp_definitions`: adds a compile definition and a corresponding cmake option

Can be included with:
```cmake
include(fortranpreprocessor)
```

### `CodeCoverage.cmake`

Copied from: https://github.com/bilke/cmake-modules/blob/master/CodeCoverage.cmake (Copyright Lars Bilke, BSL 1.0 license, see file)

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
- will also search for `MPI`, `OpenMP` and `LAPACK` if the respective option is `ON`

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

### `FindYAC.cmake`

Can be used like:
```cmake
find_package(YAC)
target_include_directories(<target> PUBLIC ${YAC_INCLUDE_DIR})
target_link_libraries(<target> PUBLIC ${YAC_LIB})
```

To specify a particular NetCDF library, use

    cmake -DYAC_ROOT=/path/to/yac -B build

or set environment variable `YAC_ROOT=/path/to/yac`.

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

### `CPM.cmake` (v0.36.0)

Copied from: https://github.com/cpm-cmake/CPM.cmake (Copyright Lars Melchior and contributors, MIT license, see file)

CPM.cmake is a CMake script that adds dependency management capabilities to CMake.
It's built as a thin wrapper around CMake's [FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html)
module that adds version control, caching, a simple API [and more](#comparison-to-pure-fetchcontent--externalproject).

Can be used like:
```cmake
include(CPM)
CPMAddPackage("https://git.ufz.de/chs/forces.git@0.3.2")
target_link_libraries(<project> PUBLIC forces)
```

## Adding to your project

In order to add these files to your project, you can either just copy the desired files or use
[git-subrepo](https://github.com/ingydotnet/git-subrepo) to add this repository as a folder:
```bash
git subrepo clone --branch=main --force https://git.ufz.de/chs/cmake-fortran-scripts.git cmake
```
You should replace `main` with a specific tag (like `v1.8`) to use a released version.

## License

MIT License (MIT)

Copyright (c) 2020 - 2022 CHS Developers
