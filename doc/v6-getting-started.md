# mHM v6: Build and run {#v6_getting_started}

The v6 implementation runs independent domains through Input, MPR, meteo, mHM,
mLM, and mRM containers. Each domain owns its exchange and component state.
The execution entrypoint is `app/driver.F90`.

For the architecture and a practical workflow for extending components, see the
[v6 developer guide](@ref v6_development).

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
python -m pip install 'nml-tools==0.6.0'
nml-tools generate --config nml-config.toml
nml-tools check --config nml-config.toml --diff
```

GitLab installs nml-tools using `required-version` from `nml-config.toml` and runs
the generated-file check before compilation. All-optional output groups
may be omitted from an existing output file; omitted groups use the same defaults
as empty groups. Configured output paths still determine which streams are active.

Namelist reference documents remain in `doc/nml/`. The Doxygen documentation
includes the scientific overview, source API, installation guidance, and
release history. Generate HTML and LaTeX from the repository root with:

```bash
PROJECT_NUMBER="$(cat version.txt)" doxygen doc/doxygen.config
```

Open `html/index.html` to browse it. The documentation environment in
`doc/rtd-environment.yml` pins Doxygen 1.18.0. The vendored theme is
[doxygen-awesome v2.5.0](https://github.com/jothepro/doxygen-awesome-css/releases/tag/v2.5.0);
its styles and extensions are copied unchanged into `doc/doxygen-awesome-css/`.
The header shows Doxygen's standard project version from `version.txt`;
On Read the Docs, the navigation dropdown offers version selection and a PDF
download for the current version. Disable the floating flyout under
`Settings > Addons > Flyout menu` in RTD; the dropdown still receives its Addons data.

Build the PDF with system TeX Live using `make -C latex`. The extra stylesheet
`doc/mhm_api_tables.sty` uses `longtable` for API parameter tables and
`xltabular` for simple documentation tables to reduce PDF build time. Merged
cells, captions, and labels use Doxygen's standard table renderer. Install the
corresponding TeX Live packages for local PDF builds, including `newunicodechar`
for symbols in source listings. The Unicode mapping supports pdfLaTeX and XeTeX
(including Tectonic). GitLab compiles with Tectonic using `-Z continue-on-errors`,
runs TeX Live's `makeindex`, then compiles again to include the index. Its Fortran
compiler comes from the `chs-fortran-conda` module's `FC` setting. Each pipeline
documents its checked-out commit and retains `html/`, `mhm_doc.pdf`, and the
Doxygen/PDF logs as artifacts.
GitLab Pages publishes only the `main` branch build at its root; Read the Docs
provides versioned documentation. Tectonic chooses its own reruns;
`-r 10` would force ten additional passes rather than set a maximum.
Read the Docs uses `latexmk -pdf` to schedule compiler, bibliography, and index
passes from changed files, and retains the complete manual with source listings.
Environment creation, Doxygen, PDF compilation, and publishing
run as separate commands so the build log identifies the current phase.
`doc/readthedocs_build.sh` runs all phases when called without arguments, logs
the PDF build to `doxygen_latex_make.log`, and prints the last 80 log lines on
failure. Set `MHM_DOCS_ENV` to use a different Conda environment name.
The theme, headers, references, and figures remain in `doc/`.

### Coverage reports

Coverage uses a serial GNU Debug build with pFUnit installed:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DMHM_BUILD_DRIVER=ON \
  -DMHM_BUILD_TESTING=ON -DBUILD_TESTING=ON -DMHM_WITH_COVERAGE=ON -DMHM_WITH_OpenMP=OFF
cmake --build build --parallel 4
cmake --build build --target mhm_coverage_CI
mkdir -p coverage
cp -a build/mhm_coverage_CI/. coverage/
lcov --rc check_data_consistency=0 --add-tracefile build/mhm_coverage_CI.info \
  --substitute "s#^(\Q$PWD\E|\Q$(pwd -P)\E)/##" --output-file coverage/coverage.info
lcov_cobertura coverage/coverage.info --base-dir . --output coverage/coverage.xml
```

Install LCOV 2.0 or newer (`lcov` and `genhtml`) and Python's `lcov_cobertura`
before running these commands. Use the `gcov` executable from the compiler's GCC version;
set CMake's `GCOV_PATH` if the default differs. Coverage builds require pFUnit
so that the unit suites cannot be silently skipped. The coverage target runs
CTest and the curated v6 smoke scenarios, retaining unexecuted code in its baseline.

LCOV strips both logical and physical checkout prefixes, preserving Fortran filename
case and line counters. The commands stage HTML, LCOV, and Cobertura XML in `coverage/`.
The HTML report uses the repository root as its prefix, showing `app/` and `src/` paths.
Function/line consistency checks are disabled because GNU Fortran can map distinct
functions to overlapping source lines; the original counters are retained.
GitLab coverage uses `chs-fortran-conda/0.6.0`, including LCOV 2.5, and takes its
compiler from the module settings. GitLab reads the Cobertura XML for MR diff
annotations and exposes the HTML report as a CI artifact; coverage is not copied to Pages.
The GitHub Coverage workflow builds pFUnit 4.12.0 with GNU Fortran 13 and uploads
the same report format through the Coveralls GitHub Action. It uses the built-in
GitHub token, so no additional Coveralls repository secret is needed. GitHub
retains coverage artifacts for 14 days. GitLab uses the instance's default artifact
expiry and the project's setting for retaining the latest successful artifacts.

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
