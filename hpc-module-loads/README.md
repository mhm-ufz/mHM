# HPC Fortran Module Loads

Module load scripts on HPC Clusters for Fortran Projects at CHS.

## Toolchains at EVE

All these scripts will load:

- the respective compilers and set `FC`, `F77`, `CC` and `CXX` env-var (optional MPI support)
- netCDF-Fortran
- CMake
- the MPR Python Environment (_except chs-conda environment_)
- pFUnit - Fortran unit testing framework (_not available for GNU 6.4_)

### Usage
- Conda environment with gfortran:
  ```bash
  source eve.chs-conda01 # or
  source eve.chs-conda02
  ```
- GNU 6.4 compiler (`foss/2018a` Toolchain):
  ```bash
  source eve.gfortran64 # or
  source eve.gfortran64MPI
  ```
- GNU 7.3 compiler (`foss/2018b` Toolchain):
  ```bash
  source eve.gfortran73 # or
  source eve.gfortran73MPI
  ```
- GNU 8.3 compiler (`foss/2019b` Toolchain):
  ```bash
  source eve.gfortran83 # or
  source eve.gfortran83MPI
  ```
- GNU 10.2 compiler (`foss/2020b` Toolchain):
  ```bash
  source eve.gfortran102 # or
  source eve.gfortran102MPI
  ```
- GNU 12.2 compiler (`foss/2022b` Toolchain):
  ```bash
  source eve.gfortran122 # or
  source eve.gfortran122MPI
  ```
- Intel 18 compiler (`iomkl/2018b` Toolchain):
  ```bash
  source eve.intel18 # or
  source eve.intel18MPI
  ```
- Intel 19 compiler (`iomkl/2020a` Toolchain):
  ```bash
  source eve.intel19 # or
  source eve.intel19MPI
  ```
- Intel 19.1.3 compiler (`iomkl/2020b` Toolchain):
  ```bash
  source eve.intel2020b # or
  source eve.intel2020bMPI
  ```
- NAG 6.2 compiler:
  ```bash
  source eve.nagfor62
  ```
- PGI 19 compiler:
  ```bash
  source eve.pgi19
  ```

## Toolchains at Juwels

All these scripts will load:

- the respective fortran compiler and set `FC` env-var
- netCDF-Fortran
- CMake

### Usage

- GNU compiler with MPI
  ```bash
  source juwels.gfortranMPI
  ```
- Intel compiler with MPI
  ```bash
  source juwels.intelMPI
  ```

## Toolchains on Atos

All these scripts will load:

- the respective fortran compiler and set `FC` env-var
- MPI
- netCDF-Fortran parallel
- CMake

## Toolchain on Levante

This script will load the following modules on [Levante](https://docs.dkrz.de/doc/levante/index.html) at [DKRZ](https://www.dkrz.de):

```bash
source levante.gfortran112
```

- git
- gfortran 11.2 compiler and set `FC` env-var
- netCDF-Fortran 4.5.3
- CMake (build tools)

## Toolchain on LUMI

This script will load the following modules on [LUMI](https://www.lumi-supercomputer.eu/):

```bash
source lumi.gfortran112
```

- gfortran 11.2 compiler and set `FC` env-var
- netCDF-Fortran
- CMake (build tools)

## License

MIT License (MIT)

Copyright (c) 2020 - 2023 CHS Developers
