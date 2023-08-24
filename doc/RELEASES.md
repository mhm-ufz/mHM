# mHM Release Notes

[TOC]

## mHM v5.13.1 (Aug 2023)

### Enhancements
- Meteo Coupling ([!161](https://git.ufz.de/mhm/mhm/-/merge_requests/161))
  - added coupling capabilities for meteo data to mHM
  - new namelist `coupling` in `mhm.nml` with the following entries:
    - `case`: 0 or 1 to dis-/enable coupling, 0 by default
    - `meteo_timestep`: 1 or 24 to set the coupling meteo time-step in hours, 0 by default
    - `meteo_time_ref_endpoint`: bool, whether the time-stamp of the coupled meteo data refers to the endpoint of the time interval, .false. by default
    - `meteo_expect_pre`: bool, whether precipitation is expected from the coupling, .false. by default
    - `meteo_expect_temp`: bool, whether temperature is expected from the coupling, .false. by default
    - `meteo_expect_pet`: bool, whether PET is expected from the coupling, .false. by default
    - `meteo_expect_tmin`: bool, whether daily min. temperature is expected from the coupling, .false. by default
    - `meteo_expect_tmax`: bool, whether daily max. temperature is expected from the coupling, .false. by default
    - `meteo_expect_netrad`: bool, whether net radiation is expected from the coupling, .false. by default
    - `meteo_expect_absvappress`: bool, whether absolute vapour pressure is expected from the coupling, .false. by default
    - `meteo_expect_windspeed`: bool, whether windspeed is expected from the coupling, .false. by default
    - `meteo_expect_ssrd`: bool, whether short wave radiation is expected from the coupling, .false. by default
    - `meteo_expect_strd`: bool, whether long wave radiation is expected from the coupling, .false. by default
    - `meteo_expect_tann`: bool, whether annual mean air temperature is expected from the coupling, .false. by default
  - new namelist entry in `directories_mHM`:
    - `dir_meteo_header`: folder containing the ascii header for level2, to be used when precipitation is coupled and not present, by default `<dir_Precipitation>`

- Python bindings ([!161](https://git.ufz.de/mhm/mhm/-/merge_requests/161))
  - added unit-tests for coupling meteo data
  - added `get_mask` to get domain mask for selected level
  - added `compressed` argument to `get_variable` to get flat array
  - added `set_meteo` routine to set meteo-data for current time-step (uses new `mhm.set.meteo` low level function)
  - added `mhm.model.config_coupling` to set coupling configuration as an alternative to the `coupling` namelist
  - added `mhm.get.number_of_horizons` to get the number of horizons of the current setup

- Python setup ([!163](https://git.ufz.de/mhm/mhm/-/merge_requests/163))
  - switching to [scikit-build-core](https://scikit-build-core.readthedocs.io) als build backend to be future-prove
  - force using the [Ninja build system](https://ninja-build.org), fixing the [conda feedstock ninja patch](https://github.com/conda-forge/mhm-feedstock/blob/a77b4c0dae25bfcc3e77e87f6c5cdec8dff6339f/recipe/conda_force_ninja.patch#L7)
  - mhm console script now always included
  - slightly restructured python package to separate python files from wrapper for clean wheels
  - editable installs more stable now
  - cleaned up env-var configs for python installation

### Bugfixes

- reset verbosity level when finalizing model
- Cmake: fixed compiling issue with Python bindings on MacOS using clang ([see conda feedstock patch](https://github.com/conda-forge/mhm-feedstock/blob/a77b4c0dae25bfcc3e77e87f6c5cdec8dff6339f/recipe/conda_clang_fixes.patch#L7))
- Documentation: now built with doxygen 1.9.7 for better type documentations and some minor fixes ([!167](https://git.ufz.de/mhm/mhm/-/merge_requests/167))

## mHM v5.13.0 (May 2023)

### Enhancements

- added simple namelists in each test-domain folder for single domain run ([!144](https://git.ufz.de/mhm/mhm/-/merge_requests/144))
  - can be run by `mhm ./test_domain`
- Python bindings and wheel distribution for mHM ([!138](https://git.ufz.de/mhm/mhm/-/merge_requests/138))
- updated FORCES to v0.5 ([!141](https://git.ufz.de/mhm/mhm/-/merge_requests/141), [!157](https://git.ufz.de/mhm/mhm/-/merge_requests/157))
- added silent mode (`--quiet / -q`) ([!141](https://git.ufz.de/mhm/mhm/-/merge_requests/141))
- Output: unification, selectable time reference point and bounds for all coordinates ([!139](https://git.ufz.de/mhm/mhm/-/merge_requests/139))
  - added option `output_time_reference[_mrm]` to `m[h|r]m_outputs.nml` where users can now select the reference time stamp for the current time interval:
    ```
    location of reference time point in outputs:
    0: start of the time interval (i.e. 1990-01-01 00:00 for daily output on 1990-01-01)
    1: center of the time interval (i.e. 1990-01-01 12:00 for daily output on 1990-01-01)
    2: end of the time interval (i.e. 1990-01-02 00:00 for daily output on 1990-01-01)
    ```
    - This will be `0` by default. Before it was always the end of the time interval -1h (that is why we always had 23h as timestamp)
    - When using the start of the time interval we actually match the time-stamps of the gauge time-series (that was shifted in `discharge.nc`)
    - When using center of time interval output unit may change to minutes (to get center of an hour)
  - added `time_bnds` to all netcdf outputs to explicitly describe the time-span for the given time-stamp (following cf-conventions)
  - addad `axis` attribute to all coordinates in all netcdf outputs (following cf-conventions)
  - added `bounds` for all spatial coordinate axes (either easting/northing or lon/lat)
  - unified output handling in `mo_nc_output`
  - use new output handler for groundwater coupling
- test-domain download helper command added with python-bindings `mhm-download` ([!149](https://git.ufz.de/mhm/mhm/-/merge_requests/149))
- added a dependency install script `CI-scripts/install-deps` ([!152](https://git.ufz.de/mhm/mhm/-/merge_requests/152))
- added demonstration script to optimize mHM with [Spotpy](https://github.com/thouska/spotpy) ([!145](https://git.ufz.de/mhm/mhm/-/merge_requests/145))

### Changes

- all errors are now printed to std-err and raise exit code 1 ([!141](https://git.ufz.de/mhm/mhm/-/merge_requests/141))
- changed command line option short name for version from `-v` to `-V` (`-v` reserved for verbosity) ([!141](https://git.ufz.de/mhm/mhm/-/merge_requests/141))

### Bugfixes

- fixed indexing in groundwater coupling calculations ([!139](https://git.ufz.de/mhm/mhm/-/merge_requests/139))
- fixed reading restart without original input files ([!154](https://git.ufz.de/mhm/mhm/-/merge_requests/154), [!156](https://git.ufz.de/mhm/mhm/-/merge_requests/156))
- output: fixed single precision output writing ([!159](https://git.ufz.de/mhm/mhm/-/merge_requests/159))
- several internal cleanups


## mHM v5.12.0 (Oct 2022)

### Enhancements

- model can now run with hourly forcing ([!39](https://git.ufz.de/mhm/mhm/-/merge_requests/39))
- new neutron module ([!92](https://git.ufz.de/mhm/mhm/-/merge_requests/92))
- add [FORCES](https://git.ufz.de/chs/forces/) as a dependency ([!103](https://git.ufz.de/mhm/mhm/-/merge_requests/103), [!124](https://git.ufz.de/mhm/mhm/-/merge_requests/124))
- new command line interface to pass used paths to mhm ([!105](https://git.ufz.de/mhm/mhm/-/merge_requests/105))
- new target function (`opti_function=34`) for baseflow index optimization ([!116](https://git.ufz.de/mhm/mhm/-/merge_requests/116))
- add L1_temp_calc and L1_prec_calc as global variables for coupling ([!120](https://git.ufz.de/mhm/mhm/-/merge_requests/120))
- new sub-daily discharge output ([!114](https://git.ufz.de/mhm/mhm/-/merge_requests/114))
- updated install instructions for all systems ([!124](https://git.ufz.de/mhm/mhm/-/merge_requests/124))
- better error messages when reading LUTs or time-series ([!128](https://git.ufz.de/mhm/mhm/-/merge_requests/128))
- new post-proc R scripts ([!102](https://git.ufz.de/mhm/mhm/-/merge_requests/102), [!113](https://git.ufz.de/mhm/mhm/-/merge_requests/113))
- Output:
  - add '_FillValue' alongside 'missing_value' in NetCDF output ([!127](https://git.ufz.de/mhm/mhm/-/merge_requests/127))
  - add snow melt flux output ([!126](https://git.ufz.de/mhm/mhm/-/merge_requests/126))
  - add line breaks in ConfigFile.log for simulation periods ([!125](https://git.ufz.de/mhm/mhm/-/merge_requests/125))
  - discharge.nc is now written with modern netcdf routines ([!135](https://git.ufz.de/mhm/mhm/-/merge_requests/135))

### Changes

- refactoring of mhm_eval to use interfaces ([!108](https://git.ufz.de/mhm/mhm/-/merge_requests/108))
- PGI support was dropped ([!108](https://git.ufz.de/mhm/mhm/-/merge_requests/108))
- change units of easting/northing in output to "m" ([!125](https://git.ufz.de/mhm/mhm/-/merge_requests/125))

### Bugfixes

- mo_grid now uses `nint` to calculate nrows and ncols in order to be stable when cell size factor is close to 1 ([!56](https://git.ufz.de/mhm/mhm/-/merge_requests/56))
- bug fix of mo_restart on reading the correct dimension of land-cover ([!121](https://git.ufz.de/mhm/mhm/-/merge_requests/121))
- cmake doesn't alter `CMAKE_MODULE_PATH` anymore ([!122](https://git.ufz.de/mhm/mhm/-/merge_requests/122))
- slope_tmp now allocated in L11_calc_celerity ([!110](https://git.ufz.de/mhm/mhm/-/merge_requests/110))

## mHM v5.11.2 (Jul 2021)

### Enhancements

- documentation modernized with [doxygen-awesome-css](https://github.com/jothepro/doxygen-awesome-css) (!86)
- cmake update to be able to install mHM (`cmake --install`) (!85)
- added pFUnit tests thanks to Nicola Döring (!76)
- link to a new [YouTube tutorial](https://youtu.be/FGJOcYEzbP4) for compiling mHM with cygwin by Mehmet Cüneyd Demirel added to the documentation (!74)
- NetCDF output: add deflate and precision option to namelists (!73)
- refactor cmake workflow (!72)

### Bugfixes

- fixed: `mrm` tried to write output even if routing was switched off (!82)
- unreachable `else` branch in `feddes_et_reduction` removed (!77)
- unnecessary `inout` variable intent in `soil_moisture` removed (!77)


## mHM v5.11.1 (Mar 2021)

### Enhancements

- added compile information for cygwin (!68)

### Bugfixes

- removed note about mHM 5.10 from the README
- smhorizon: tmp_rootfraccoef was corrected directly if it is not between 0 and 1, but actually  FCnorm should always be between 0 and 1 (!67)

## mHM v5.11 (Feb 2021)

### Experimental Features

- river temperature routing was implemented in an alpha version 0.1 (!37)
  - this feature is in an experimental stage and should not be considered stable!

### Enhancements

- introduced central version files `version.txt` and `version_data.txt` (!51)
- added Feddes and global FC dependency on root fraction coef. at SM process(3)=4 (!43)
- Online documentation generated with doxygen: https://mhm.pages.ufz.de/mhm/develop/ (!44)
- CI/CD with GitLab Runner (!11, !13, !14, !28, !32, !48, !50)
  - building on EVE for multiple compiler (GNU 7.3/8.3, Intel 18/19, NAG 6.2)
  - building debug/release serial/parallel
  - memcheck with valgrind
  - running all check-cases with all compiled versions
  - calculation of coverage
  - new checking script `run_mhm_checks.py`
- the domain loop is now parallelized with MPI
- objective function for boxcox-transformed streamflow
- post processing script for probabilistic forecasts
- different module load scripts for EVE
- Objective function from separate mhm calls (!7)
- new data type for simulated gridded optidata (!10)
- new datatype datetimeinfo (!16)
- added module `mo_os` to check files and directories (!41, !57)

### Changes

- internal: "basin" renamed to "domain"
- TWS input file changed from ascii to netCDF (!9)
- Switched to cell wise kge of et and tws in opti_function 33 (!12)
- restart files are now given by name (!34)
- removed mRM standalone and statically integrate mRM into mHM (!53)
- removed the old makefile and legacy checking scripts (!55)
- minimal Cmake version is now `3.12` (!58)

### Bugfixes

- Finalparam.nml is now written with specific format (Intel/GNU compatibility) (#40)
- FinalParam.nml routing section bug fixed (#49, !25)
- dirEvapotranspiration is now allocated before writing
- cmake: netcdf link flags where separated by ";"
- sharing of L0 domain now working
- added L1_jarvis_thresh_c1 to restart file for process id 2 AND 3 (#29, !15)
- allowing higher routing resolution than hydrology (!21)
- domainID not set correctly for mRM if restart is activated (!30)
- mHM states_fluxes netCDF output was curvilinear even if coordinate system is set to regular latlon (#98, !31)
- missmatch in messages about written mhm fluxes (!42)
- Fixing wrongly matched IDs from L1 to L11 when routing resolution (L11) is finer than L1 resolution (!45)
- The length in net_startup was only cut in case there are less then 2 lengths (!46)
- corrected unit attributes for lat lon variables (!47)
- Allow run mHM and mRM without any observed gauge for processCase(8) = 2 / 3 (#27, !52)

## mHM v5.10 (June 2019)

### New Features:

- New routing process introduced `processCase(8) = 3` [see Thober et al, 2019, GMDD, in press](https://doi.org/10.5194/gmd-2019-13) for more details
- mRM is decoupled from mHM and mRM now resides in `deps/mrm` as an independent submodule [more information on handling submodules](https://git.ufz.de/howto/gitlab)
- New option to compile mHM with cmake is provided, see more details under [cmake manual][1].
- Visualization/animation R script (producing PDF and GIF) of mHM netcdf files included under post-proc [animate1.R][2]
- New option for coupling of mRM to a groundwater model (`gw_coupling = .true.`). The river head can be computed based on the Manning equation.

[1]: https://git.ufz.de/mhm/mhm/blob/develop/doc/INSTALL.md
[2]: https://git.ufz.de/mhm/mhm/blob/develop/post-proc/animate1.R

### Bugs resolved from release 5.9:

- Enable use of i8 for time_data in common/mo_read_forcing_nc.f90,
  otherwise netcdf time stamps with initial dates prior to 1900 were
  wrong (due to overflow)

### Known bugs:

1. Adaptive routing does not allow to run without at least 1 gauge specification
2. Incompatibility of Finalparam.nml format between Intel and GNU
3. If ProcessOption(3) is set to 3 and optimization is activated, the
  created FinalParam.nml misses the header for the namelist of the
  soil moisture parameters.
4. Land cover scenes cannot be changed between the run generating the
  restart file and the run using the restart file. This causes
  unpredictable behaviour by the model.
5. Simulation period must span overall land cover scenes specified in the namelist.
6. Cut-off for link length is calculated with missing values, but those
should be neglected.
7. Using a higher routing resolution than hydrology resolution may
  cause segmentation faults because mapping from L1id on L11 is not
  working correctly
8. If ProcessOption(5) is set to -1 and optimization is activated, the
  created FinalParam.nml misses the header for the namelist of the
  PET process.

### Restrictions:

- For `gfortran` compilers mHM supports only v4.8 and higher.
- If you wish to use features connected to ground albedo neutrons (`processCase(9)`), please contact [Martin Schrön](mailto:martin.schroen@ufz.de).

## mHM v5.9 (July 2018)

### New Features:

- Major restructuralization of the mHM code.
- MPR is now executed before mHM is run.
- MPR can be compiled as a standalone tool.
- mHM without mRM can be compiled.
- The code in general is strictly reorganized into modules that belong
to their respective processes (MPR, mHM, mRM). This is not only done
for constants, global variables and so on, but also for every module
and the reading of input as well as the namelists.
- This leads to the creation of those folders:
 - `./common` (code shared by MPR, mHM and mRM) included for every compilation option
 - `./lib` (code shared by MPR, mHM and mRM) included for every compilation option
 - `./MPR` (code for MPR), not included for mRM standalone
 - `./common_mHM_mRM` (code shared by mHM and mRM), not included for MPR standalone
 - `./mHM` (code for mHM), not included for MPR and mRM standalone
 - `./mRM` (code for mRM), not included for MPR standalone and mHM without routing
- Code is reformatted (indentation=2, spacing unified)
- Removed many duplicate code parts (e.g. shared between mHM and mRM)
- Check cases are minimized (reduced output, shorter time periods (<=2 years), less basins)
- Check cases can now be set up more easily (by use of model_wrapper for automatic creation of new nml files)
- Check cases now run a python script for output comparison, advantage: tolerance now also allowed for ascii-output and support of 4-D netCDF files without time dimension
- fSealed is now an effective parameter
- mHM effective parameters now all have three dimensions internally (nCells_L1, [iHorizon, LAI-Time], nLCoverScenes)
- Introduction of new derived types:
 - `Grid` (merge of basin_info, basin_info_mrm, gridGeoRef, nCells, longitude, latitude, Id) used for each level (0, 1, 11, 2) individually
 - `GridRemapper` (merge of lower_bound, upper_bound, etc.)
- mhm_eval and mrm_eval now have common procedure interface (needed for fully flexible optimisation)
- A new post-processor can check and adapt some fields for doxygen generation

- Changes that are not backwards-compatible:
 - Restart files are restructured (now contain only the minimum required for restart):
 - MPR: effective parameters at L1 + grid information
 - mHM: effective parameters, states and fluxes at L1 + grid information
 - mRM: routing-specific parameters, states, fluxes, configs at L1/L11 + grid information
 - mhm.nml is restructured and not backwards compatible mainly due to the reorganization of modules in dependency of their processes
 - gridded LAI values are now used for all effective parameters (no fallback to LAIclasses anymore)
 - canopy height used for aerodynamic resistance is now scaled with the actual LAI timeseries and not with a dummy timeseries of intensive orchard

- Considerable improvements and reduced redundancy in estimation of an empirical distribution of slopes (sort function in the mo_startup).
Example for the Australian domain (180 million L0 cells), it reduced time for sorting slope from 32hours to only 1 minute.
- mtCLIM preprocessor in pre-proc/mtCLIM, based on [mtCLIM v4.3](http://www.ntsg.umt.edu/project/mt-clim.php). This code
is able to estimate humidity (vapore pressure or vapore pressure deficit) and incoming shortwave radiation based on meteorological variables
(minimum and maximum air temperature, precipitation) and morphological characteristics of the underlying terrain (digital elevation model, slope, aspect).
- mHM2OGS preprocessor in pre-proc/GIS2FEM3. This code converts the
triangular-wise or quadrilateral-wise recharge data from mHM
into the nodal source terms of a three dimensional finite element model for [OGS](https://www.ufz.de/index.php?en=38011).
- Updated documentation for CYGWIN installations under Windows 7
and 10.
- Added new objective function number 31: weighted NSE (NSE is
weighted with observed discharge)
- Removal of bin files and related code (only nc and ascii files are used)

### Bugs resolved from release 5.8:

- Enable use of i8 for time_data in common/mo_read_forcing_nc.f90,
  otherwise netcdf time stamps with initial dates prior to 1900 were
  wrong (due to overflow)

### Known bugs:

None.

### Restrictions:

- For `gfortran` compilers mHM supports only v4.8 and higher.
- If you wish to use features connected to ground albedo neutrons (`processCase(9)`), please contact [Martin Schrön](mailto:martin.schroen@ufz.de).
- If you wish to use the multi-scale Routing Model as stand-alone version, please contact [Stephan Thober](mailto:stephan.thober@ufz.de).





## mHM v5.8 (Dec 2017)

### New Features:

- Implementation of a new process for PET correction based on LAI at PET `process(5)=-1` (Cuneyd Demirel + *GEUS* colleagues);
- Pre-processor code for SOILGRIDS data as used for the EDgE project (Rohini Kumar);
- Reduced computational time of the neutron forward model COSMIC by factors of 30--100 (Maren Kaluza);
- Compression of the netCDF output files (David Schaefer);
- Optional project description added into the [mhm.nml](mhm.nml)

### Bugs resolved from release 5.7:

- `processCase(3)=3` did not work when compiled with openMP.
- openMP declarations missing in [mo_mpr_smhorizons.f90](src/mHM/mo_mpr_smhorizons.f90) for the case of `iFlag_soil=1`

### Known bugs:

None.

### Restrictions:

- For `gfortran` compilers mHM supports only v4.8 and higher.
- If you wish to use a special process description of evapotranspiration (`processCase(4)`) please contact [Matthias Zink](mailto:matthias.zink@ufz.de).
- If you wish to use features connected to ground albedo neutrons (`processCase(9)`), please contact [Martin Schrön](mailto:martin.schroen@ufz.de).
- If you wish to use the multi-scale Routing Model as stand-alone version, please contact [Stephan Thober](mailto:stephan.thober@ufz.de).


## mHM v5.7 (Jun 2017)

### New Features:

- New process descriptions for the soil evapotranspiration module:
    - Field capacity dependency to root fraction coefficient (`processCase(3)=2`) -- implemented by *GEUS*.
    - Jarvis (1989, J. Hydrol.) evapotranspiration reduction (`processCase(3)=3`) -- implemented by *GEUS*.
- Use local, monthly LAI climatology instead of look-up table, i.e., [LAI_classdefinition.txt](test_basin/input/morph/LAI_classdefinition.txt) (`timeStep_LAI_input=1`).
- New objective functions for model calibration
    - Calibration of mHM using catchment average evapotransipration (`opti_function=27`).
    - Calibration of mHM using soil moisture and streamflow simultaneously (`opti_function=28`).

### Bugs resolved:

- Calibration using `processCase(8)=2` (routing with adaptive timestep) does properly work now.
- Streamflow output is now properly written to the NetCDF.

### Known bugs:

None

### Restrictions:

- For `gfortran` compilers mHM supports only v4.8 and higher.
- If you wish to use a special process description of evapotranspiration (`processCase(4)`) please contact [Matthias Zink](mailto:matthias.zink@ufz.de).
- If you wish to use features connected to ground albedo neutrons (`processCase(9)`), please contact [Martin Schrön](mailto:martin.schroen@ufz.de).
- If you wish to use the multi-scale Routing Model as stand-alone version, please contact [Stephan Thober](mailto:stephan.thober@ufz.de).


## mHM v5.6 (Dec 2016)

### New Features:

- **Routing extended:** Implementation of a new parametrization for the routing model (`processCase(8)=2`). This routing option is based on an adaptive time step to improve the scalability and transferability of the model as well as a significant reduction in run time. The adaptive time step is calculated as ratio of routing resolution and celerity, the latter can be given as parameter in [mhm_parameter.nml](mhm_parameter.nml).

### Bugs resolved:

- Any model time step from 1 h to 24 h can be chosen (in releases v5.4 and v5.5 only 1 h worked properly).
- Estimation of the Hargreaves-Samani PET for high altitudes works properly know (there have been numerical issues for high latitude values).
- Reading catchment outlets from the *restart* file works now (bug appeared in v5.5).

### Known bugs:

- Calibration using `processCase(8)=2` (adaptive timestep) does not work, please use `processCase(8)=1`.

### Restrictions:

- For `gfortran` compilers mHM supports only v4.8 and higher.
- If you wish to use a special process description of evapotranspiration (`processCase(4)`) please contact [Matthias Zink](mailto:matthias.zink@ufz.de).
- If you wish to use features connected to ground albedo neutrons (`processCase(9)`), please contact [Martin Schrön](mailto:martin.schroen@ufz.de).
- If you wish to use the multi-scale Routing Model as stand-alone version, please contact [Stephan Thober](mailto:stephan.thober@ufz.de).


## mHM v5.5 (Jun 2016)

### New Features:

- Routing works on domains with multiple outlets (e.g., continental level).
- New option for providing soil data. They can be provided as predefined layers (one map per layer).
- Speed up of mHM for big domains, due to reformulations in the model start up.
- Pre-processing: new tools for i) cutting out a catchment from a existing dataset, ii) estimation of Hargreaves-Samani evapotranspiration, and iii) enlarging the grids of the input files.

### Bugs resolved:

- Assigning routing parameters is done properly now.

### Known bugs:

- Specifying a model time step of 24h (in [mhm.nml](mhm.nml)) does not work, please stick with the default time step (1h)

### Restrictions:

- For `gfortran` compilers mHM supports only v4.8 and higher.
- If you wish to use a special process description of evapotranspiration (`processCase(4)`) please contact [Matthias Zink](mailto:matthias.zink@ufz.de).
- If you wish to use features connected to ground albedo neutrons (`processCase(9)`), please contact [Martin Schrön](mailto:martin.schroen@ufz.de).
- If you wish to use the multi-scale Routing Model as stand-alone version, please contact [Stephan Thober](mailto:stephan.thober@ufz.de).


## mHM v5.4 (Dec 2015)

### New Features:

- The routing of mHM can be used as a stand alone version/independent software called *multiscale Routing Model mRM*, e.g. for coupling to other environmental models.
- A new output, i.e. fields of routed discharge, is now available. They are stored in `mRM_fluxes_and_states.nc` and are controlled by a new namelist contained in the [mrm_outputs.nml](mrm_outputs.nml) file, which is an optional file.
- New calibration objectives have been incorporated. It is now possible, additionally to the former objectives, to calibrate mHM against additional input data:
    - total water storage (e.g. GRACE) and discharge simultaneously (`opti_function=15`), and/or
    - cosmic ray neutron counts (`opti_function=17`).
- New post-processing: a mHM python class for reading all inputs and outputs of a model run can be found in [post-proc/](post-proc/).
- Reorganization of the NetCDF writing in mHM to simply future implementations of additional outputs.

### Bugs resolved:

- Calibration with catchment average soil moisture (`opti_function=10`) works properly now.
- Discharge output for multi basin runs with different time periods for each basin works properly now.
- Hargreaves-Samani PET calculation (`processCase(5)=1`) is valid on southern hemisphere too now.

### Known bugs:

- None.

### Restrictions:

- For `gfortran` compilers mHM supports only v4.8 and higher.
- If you wish to use a special process description of evapotranspiration (`processCase(4)`) please contact [Matthias Zink](mailto:matthias.zink@ufz.de).
- If you wish to use features connected to ground albedo neutrons (`processCase(9)`), please contact [Martin Schrön](mailto:martin.schroen@ufz.de).
- If you wish to use the multi-scale Routing Model as stand-alone version, please contact [Stephan Thober](mailto:stephan.thober@ufz.de).


## mHM v5.3 (Jun 2015)

### New Features:

- Simulation period and warming days can be now given per basin (see `time_periods` in [mhm.nml](mhm.nml))
- Enabling use of MPI (set `mpi=true` in [Makefile](Makefile) and use `#ifdef MPI` for MPI specific code)
- Optional input data can be loaded for example to calibrate against soil moisture (see `optional_data` in [mhm.nml](mhm.nml))
- Generation of ground albedo cosmic-ray neutrons (see `processCase(10)` in [mhm.nml](mhm.nml)); these calculations are based on the [COSMIC](http://cosmos.hwr.arizona.edu/Software/cosmic.html) code, which was originally written by Rafael Rosolem. Please contact [Martin Schrön](mailto:martin.schroen@ufz.de) if you like to use this new feature.
- Several new objective functions, e.g. calibrating the Kling-Gupta efficiency of catchment's average soil moisture (`opti_function=10`) or calibrating multiple basins regarding Kling-Gupta efficiency of discharge (`opti_function=14`) among others; calibration against soil moisture is still purpose to research (`opti_function=10-14`). The interested user may contact [Matthias Zink](mailto:matthias.zink@ufz.de) for further details.

### Bugs resolved:

- Calibration using potential evapotranspiration from input file (i.e. `processCase(4)=0`) is now working properly

### Known bugs:

- Compiling mHM with the recent Cygwin version under Windows is leading to an error message indicating circular dependencies. The reason for this is Unicode characters in some source code files. Please contact [mhm-admin@ufz.de](mhm-admin@ufz.de), if you get this error message. We will provide you the files with cleaned characters..

### Restrictions:

- If you wish to use a special process description of evapotranspiration (i.e. Hargreaves-samani, Priestley-Taylor, or Penman-Monteith) please contact [Matthias Zink](mailto:matthias.zink@ufz.de). The special cases are set in [mhm.nml](mhm.nml) (see `processCase(4)`).
- If you wish to use the new feature of calculating neutron counts please contact [Martin Schrön](mailto:martin.schroen@ufz.de). The feature can be enabled in [mhm.nml](mhm.nml) (see `processCase(8)`).


## mHM 5.2 (Dec 2014)

### New Features:

- Chunk-wise reading of input data (see `timestep_model_inputs`)
- Complete revision of writing netCDF files
- Possibility to discard multi-scale parameter regionalization (MPR) calculations (see `perform_mpr`)
- Several process descriptions of evapotranspiration implemented (see `processCase(4)`): Read PET, Hargreaves-Samani, Priestley-Taylor, Penman-Monteith. Please contact [Matthias Zink](mailto:matthias.zink@ufz.de) if you use one of the last three options, since the code is not under GNU Public license up to now.
- Adding routines for signature calculations of time series (see `mo_signatures`)
- New objective function for calibrating discharge with Kling-Gupta efficiency measure (KGE, see `opti_function`)
- New output variables (see [mhm_outputs.nml](mhm_outputs.nml))
- Sorting algorithm changed to public available library orderpack (see `mo_orderpack`)

### Bugs resolved:

- Some variables in restart file where not assigned correctly
- Variables not initialized correctly

### Known bugs:

- Calibration using PET values read from the input file (`processCase(5)=0`) is running, but yields wrong results due to a wrong initialization of variables. **The bug is resolved and will be released with version 5.3.**


## mHM 5.1 (Jun 2014)

### New Features:

- OpenMP handling of routines such as the multi-scale parameter regionalization (MPR)
- Multi-scale implementation, i.e. running mHM simultaneously in several basins with different resolutions
- Automatic check case framework, i.e. testing new implementations on their validity and back-compatibility
- Implementation of inflow gauges, i.e. feeding discharge time series from upstream areas at catchment boundaries
- File [gaugeinfo.txt](https://git.ufz.de/mhm/mhm/blob/da4702f3d91bee97edee5bc5726c4617a89f44a5/gaugeinfo.txt) specifying gauging stations is now part of namelist [mhm.nml](mhm.nml)
- Code is now free of Numerical Recipes proprietary code
- Can now run on a single cell (no routing performed) Hydrological modelling resolution (L1) equal to morphological input data resolution (L0) possible
- Windows compatible (with Cygwin)
- Support of regular geographic coordinate systems (e.g lat-lon) in addition to equal-area coordinate systems (UTM)

### Bugs resolved:

- Initialization of states was not correct when running mHM in calibration mode.
- Calculated parameter values ([mhm_parameters.nml](mhm_parameters.nml)) not necessarily in bound (check added).
- Aggregation/Disaggregation of meteorological data corrected.
- Forecast with mHM did not work because modelling period was restricted to discharge data period.
- Wrong mapping of evaluation discharge gauges for runs involving multiple gauges.

### Known bugs:

- Print out of River network in Config File is wrong for Multi-Basin setup, i.e., the River network is always properly written for the first basin, but not properly for subsequent basins when these are either different ones or the same one with a different Hydrology or Routing resolution.
- mHM does not abort if x-axis of L0 (morphological data) and L2 (meteorological data) do not span over exactly the same range.



## mHM 5.0 (Dec 2013)

### New Features:

- Full modular version
- Automatic documentation by doxygen
- Running mHM for multiple basin simultaneously
- Definition of 8 major processes:
    - interception,
    - snow,
    - soil moisture,
    - direct runoff,
    - evapotranspiration,
    - interflow,
    - percolation,
    - routing
- Choice of different descriptions of processes possible
- Input in binary `*.bin` or netcdf `*.nc` format
- Various calibration routines and objective functions
- Consistent numerical precision handling of variables

### Known bugs:

- None.
