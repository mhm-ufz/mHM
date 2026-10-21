# mHM v6: Build and run

The v6 implementation runs independent domains through Input, MPR, meteo, mHM,
mLM, and mRM containers. Each domain owns its exchange and component state.
The execution entrypoint is `app/driver.F90`.

The CLI delegates to `run_forward` in `mo_mhm_runner`. Fortran callers can import it
with `use mo_mhm_runner, only: run_forward`. Its `main_file`, `para_file`, `out_file`,
and `cwd` arguments are optional, defaulting to `mhm.nml`, `parameters.nml`,
`outputs.nml`, and `.`. These defaults are defined in `mo_common_constants` and
shared with the CLI. Callers configure logging; `cwd` resolves relative paths
without changing the process working directory.

Run the commands below from the repository root.

## Build

Requirements: CMake 3.18+, a Fortran and C compiler, NetCDF support, and FORCES.
CMake uses a local `forces/` checkout when available; otherwise it downloads the
version selected in `version_forces.txt`.

```bash
cmake -S . -B debug -DCMAKE_BUILD_TYPE=Debug -DMHM_BUILD_TESTING=ON
cmake --build debug --parallel 4
ctest --test-dir debug --output-on-failure
```

pFUnit is optional for configuration, but required to run the `.pf` unit tests.
The standalone regression and failure-path tests run without pFUnit.

For a release build:

```bash
cmake -S . -B release -DCMAKE_BUILD_TYPE=Release
cmake --build release --parallel 4
cmake --install release --prefix ./install
```

The executable is `debug/app/mhm` or `release/app/mhm`, installed as `bin/mhm`.
Enable OpenMP with `-DMHM_WITH_OpenMP=ON`. MPI execution is currently unavailable.
Use `-DMHM_BUILD_DRIVER=OFF` for a library-only build.

## Run and validate

```bash
./debug/app/mhm examples/domain_01
./debug/app/mhm examples/domain_02
./debug/app/mhm examples/
./debug/app/mhm -n combine.nml examples/
./debug/app/mhm -n test_nml/mrm_minimal.nml \
  -p test_nml/mhm_parameter_v6_routing.nml -o mhm-output-template.nml
bash CI-scripts/run_v6_namelists.sh --log-dir logs_v6 ./debug/app/mhm
```

The positional directory is the run root for relative configuration, input,
output, and restart paths. Default files within it are `mhm.nml`,
`parameters.nml`, and `outputs.nml`. Use `mhm --help` for CLI options.

The curated smoke runner stages inputs in fresh run directories under its log
directory. CTest similarly uses separate build-tree run directories.

## Configuration and Python

`nml-schemas/*.yml` and `nml-config.toml` define the namelist interfaces.
Generated modules, templates, and reference documentation are updated with:

```bash
nml-tools generate --config nml-config.toml
```

Namelist reference documents remain in `doc/nml/`. The Doxygen documentation
includes the scientific overview, source API, installation guidance, and
release history. Generate HTML and LaTeX from the repository root with:

```bash
PROJECT_NUMBER="$(cat version.txt)" doxygen doc/doxygen.config
```

Open `html/index.html` to browse it. The PDF build uses
`cd latex && tectonic -r 10 refman.tex`; documentation CI publishes HTML and PDF.
The theme, headers, references, and figures remain in `doc/`.

The Python package provides `mhm.run_forward`, `mhm.__version__`, `mhm.f_version`,
and `mhm.download_examples`. Its `mhm` console script and `python -m mhm` use the
same forward runner and default paths as the native driver. Use `mhm-download`
or `python -m mhm.download` to retrieve example domains; automatic downloads try
GitHub first and then UFZ GitLab. Stepwise model-control bindings and calibration
interfaces will be implemented against the v6 architecture later.
See [Python package guidance](../pybind/README.md).

mHM is distributed under LGPLv3+; see `LICENSE.md`, `COPYING`,
`COPYING.LESSER`, and `AUTHORS`.

## Soil horizons and layered output

The soil-depth namelist still specifies layer **bottoms** in millimetres. NetCDF
restart and output coordinates now specify layer **tops**, with bounds,
`units="mm"`, and `positive="down"`. In fixed-depth mode (`soil_db_mode=1`), all
interfaces are fixed. In mode 0, the terminal bound repeats the final layer top.
For example, layer tops `[0, 200, 500]` have vertices `[0, 200, 500, 500]`.
This is an **mHM-specific convention**: coincident final bounds indicate that mHM
does not provide the spatially variable bottom depth. They do not imply zero
physical thickness. Generic CF applications may interpret this marker as a zero-sized
cell, so these bounds must not be used to calculate bottom-layer thickness or volume.
The coordinate's `comment` explains the convention. No boundary-provenance attribute
or L1 `soil_depth` field is written. Capacities and root fractions continue to use
soil-specific depths. Both neutron formulations continue to exclude the last layer.

Existing v6 restarts with bottom-referenced coordinates are read through their
bounds. In mode 0, old nominal terminal depths are normalized to the repeated final
top; old `bottom_boundary_type` attributes are accepted and checked against the
configured soil mode. Layer counts and internal interfaces must match; fixed-depth
terminal bounds must also match. Mode-1 COSMIC restarts with invalid L3 must be regenerated.

The following outputs now share a `z` dimension rather than numbering variables:

| Previous names | Current field | Meaning |
| --- | --- | --- |
| `SWC_L01`, `SWC_L02`, … | `SWC` | Soil water storage [mm] |
| `SM_L01`, `SM_L02`, … | `SM` | Relative saturation [1] |
| `soil_infil_L01`, … | `soil_infil` | Accumulated infiltration [mm] |
| `aET_L01`, … | `aET_soil` | Accumulated soil evapotranspiration [mm] |

Select layer 1 with `SWC.isel(z=0)` in xarray, for example. CDL dimension order is
`(time,z,y,x)` or `(time,z,lat,lon)`, including `z=1` for one-layer models. Total
`aET`, `SM_Lall`, and other ordinary fields keep their horizontal dimensions.
Output switches, time bounds, precision, compression, state averaging, flux sums,
and pervious-area factors are unchanged.
