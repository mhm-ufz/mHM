# mHM - Python bindings

[TOC]

Python bindings to control mHM.

The wrapper (`mhm/wrapper.f90`) is just a small layer on top of the
interfaces provided by mHM to be compatible with [f2py](https://numpy.org/doc/stable/f2py/index.html).


## Installation

There is a [PyPI package](https://pypi.org/project/mhm) to install the latest release:

```bash
pip install mhm
```

Installing the mHM Python package will provide the `mhm` command to execute mHM the traditional way.

In order to compile the Python bindings from scratch you need:
1. [Python](https://www.python.org/) with version at least v3.8 and [pip](https://pip.pypa.io/)
2. a Fortran, a C and a C++ compiler (set the environment variables `FC` (and `F77`), `CC` and `CXX` accordingly).
    In case of gcc, this could look like:
    ```bash
    export FC="gfortran"
    export F77="gfortran"
    export CC="gcc"
    export CXX="g++"
    ```
3. [NetCDF-Fortran](https://github.com/Unidata/netcdf-fortran) installed in your system path

See the [Compilation](../doc/INSTALL.md) instructions for these dependencies.

You can also use a conda environment (set up with [miniforge](https://mhm-ufz.org/guides/install-unix/) for example)
to get everything:
```bash
conda install -y pip netcdf-fortran fortran-compiler c-compiler cxx-compiler
```

To compile everything after cloning/downloading, you can use pip:

```bash
pip install -v .
```

To install it directly from the git repository you can type:

```bash
pip install -v git+https://git.ufz.de/mhm/mhm.git
```

### Environment variables

The following environment variables can be used to control the compilation and installation of the python bindings for mHM:

- `SKBUILD_CMAKE_BUILD_TYPE=[Release|Debug]`: build type for the mhm library (default: `Release`)
- `MHM_BUILD_FORCES_PATH=<path>`: custom path to forces source dir (default: None)
- `MHM_BUILD_PARALLEL=[0|1]`: whether to use OpenMP with mHM (default: `0`)


## Test domain download tool

Together with the Python bindings comes a command line script to download the test domains:
```bash
mhm-download --verbose --branch develop --domain 1 --path mhm_domain/
```

You can then run mHM on this test domain with:
```bash
mhm mhm_domain/
```

You can get help on how to use this script with `mhm-download -h`:
```
$ mhm-download -h
usage: mhm-download [-h] [-V] [-v] [-b BRANCH] [-d {1,2}] [-p PATH]

Download tool to retrieve the test domains for mHM.

optional arguments:
  -h, --help            show this help message and exit
  -V, --version         display version information
  -v, --verbose         be verbose (default: False)
  -b BRANCH, --branch BRANCH
                        branch, tag, or commit of the mHM repository to take the test domain from,
                        by default tag determined from the mHM version (default: None)
  -d {1,2}, --domain {1,2}
                        test domain '1' or '2' (default: 1)
  -p PATH, --path PATH  destination path for the downloaded folder,
                        by default the original folder name in the current directory (default: None)
```

Within python scripts, you can use this tool with `mhm.download_test`. See below for examples.


## Documentation

See `mhm.tools` and `wrapper.f90` for further information on the provided routines.


## Examples

If you have cloned the repository, you can do the following to simply run mhm without optimization:

```python
import mhm

# download test domain 1
mhm.download_test(path="example_domain")
# run the downloaded example
mhm.model.init(cwd="example_domain")
mhm.model.run()
mhm.model.finalize()
```

Or you can do the following to control each timestep:
```python
import mhm

# assuming to run from the mhm repo root
mhm.model.init()
mhm.run.prepare()
ndomians = mhm.run.get_ndomains()
for i in range(1, ndomians + 1):
    mhm.run.prepare_domain(domain=i) # 0 by default
    while not mhm.run.finished():
        mhm.run.do_time_step()
        mhm.run.write_output()
    mhm.run.finalize_domain()
mhm.run.finalize()
mhm.model.finalize()
```

See also the `examples` folder.


## License

LGPLv3 (c) 2005-2023 mHM-Developers
