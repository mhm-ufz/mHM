# Namelist Reference <!--! {#nml} -->

<!--! [TOC] -->

This section provides an overview of the namelists used in mHM for model configuration.
Each namelist is documented on a separate page, detailing its parameters and their usage.

<a name="nml_config"></a>

## Configuration Namelists

The configuration namelists define the model setup and runtime behavior.
They are all listed below, one page per namelist:

- [Project configuration](config_project.md): Configuration for the overall project setup `&config_project`.
- [Domain configuration](config_domain.md): Configuration for model domains `&config_domain`.
- [Processes configuration](config_processes.md): Configuration for process case selection `&config_processes`.
- [Grid resolution configuration](config_resolution.md): Configuration for model grid resolutions `&config_resolution`.
- [Optimization configuration](config_optimize.md): Configuration for optimization settings `&config_optimize`.
- [Time configuration](config_time.md): Configuration for time settings `&config_time`.
- [Look-up table configuration](config_lut.md): Lookup-table paths for soil, geology and LAI `&config_lut`.
- [Input configuration](config_input.md): Configuration for input data `&config_input`.
- [Coupling configuration](config_coupling.md): Configuration for coupling settings `&config_coupling`.
- [Meteorological configuration](config_meteo.md): Configuration for meteorological component `&config_meteo`.
- [MPR configuration](config_mpr.md): Configuration for the multiscale parameter regionalization (MPR) `&config_mpr`.
- [mHM model configuration](config_mhm.md): Configuration for the hydrological component `&config_mhm`.
- [mLM configuration](config_mlm.md): Configuration for the lake component (mLM) `&config_mlm`.
- [mRM configuration](config_mrm.md): Configuration for the routing component (mRM) `&config_mrm`.
- [Observations configuration](config_observations.md): Configuration for reference observations `&config_observations`.

<a name="nml_parameters"></a>

## Parameter Namelists

[Convert legacy parameter namelists to the v6 format](../parameter-migration-v6.md).

Parameter namelists use the buffer-compatible `parameter_t` derived type in this fixed order:

`value, optimize, min, max`

These fields are used as follows:

- `value`: input value and initial calibration value.
- `optimize`: calibration switch (`.true.` = calibrate, `.false.` = keep fixed).
- `min`: minimum admissible parameter value.
- `max`: maximum admissible parameter value.

Bounds are supplied by each parameter schema, so fixed parameters can use `parameter = value`
and calibration with default bounds can use `parameter = value, .true.`. The complete form
`parameter = value, optimize, min, max` overrides those defaults.

They are all listed below, one page per namelist:

- [Interception - Case 1](parameter_interception_1.md): Parameters for interception (case 1) `&interception_1`.
- [Snow - Case 1](parameter_snow_1.md): Parameters for snow (case 1) `&snow_1`.
- [Soil moisture - Case 1](parameter_soil_moisture_1.md): Parameters for soil moisture (case 1) `&soil_moisture_1`.
- [Soil moisture - Case 2](parameter_soil_moisture_2.md): Parameters for soil moisture (case 2) `&soil_moisture_2`.
- [Soil moisture - Case 3](parameter_soil_moisture_3.md): Parameters for soil moisture (case 3) `&soil_moisture_3`.
- [Soil moisture - Case 4](parameter_soil_moisture_4.md): Parameters for soil moisture (case 4) `&soil_moisture_4`.
- [Direct runoff - Case 1](parameter_direct_runoff_1.md): Parameters for direct runoff (case 1) `&direct_runoff_1`.
- [PET - Case -2](parameter_pet_m2.md): Parameters for PET (case -2 - aspect correction) `&pet_m2`.
- [PET - Case -1](parameter_pet_m1.md): Parameters for PET (case -1 - LAI correction) `&pet_m1`.
- [PET - Case 1](parameter_pet_1.md): Parameters for PET (case 1 - Hargreaves-Samani) `&pet_1`.
- [PET - Case 2](parameter_pet_2.md): Parameters for PET (case 2 - Priestley-Taylor) `&pet_2`.
- [PET - Case 3](parameter_pet_3.md): Parameters for PET (case 3 - Penman-Monteith) `&pet_3`.
- [Interflow - Case 1](parameter_interflow_1.md): Parameters for interflow (case 1) `&interflow_1`.
- [Percolation - Case 1](parameter_percolation_1.md): Parameters for percolation (case 1) `&percolation_1`.
- [Routing - Case 1](parameter_routing_1.md): Parameters for routing (case 1 - Muskingum) `&routing_1`.
- [Routing - Case 2](parameter_routing_2.md): Parameters for routing (case 2 - adaptive timestep) `&routing_2`.
- [Routing - Case 3](parameter_routing_3.md): Parameters for routing (case 3 - varying celerity) `&routing_3`.
- [Neutrons - Case 1](parameter_neutrons_1.md): Parameters for neutrons (case 1 - Desilets) `&neutrons_1`.
- [Neutrons - Case 2](parameter_neutrons_2.md): Parameters for neutrons (case 2 - COSMIC) `&neutrons_2`.
- [Baseflow - Case 1](parameter_baseflow_1.md): Parameters for geological baseflow recession (case 1) `&baseflow_1`.
- [Lakes - Case 1](parameter_lakes_1.md): Parameters for natural lakes (case 1) `&lakes_1`.
- [River temperature - Case 1](parameter_river_temperature_1.md): Parameters for river temperature (case 1) `&river_temperature_1`.

<a name="nml_outputs"></a>

## Output Namelists

The output namelists define which model outputs are generated during simulation.
They are all listed below, one page per namelist:

- [mHM output configuration](output_mhm.md): Output configuration for mHM `&output_mhm`.
- [mLM output configuration](output_mlm.md): Output configuration for mLM `&output_mlm`.
- [mRM output configuration](output_mrm.md): Output configuration for mRM `&output_mrm`.


<!--!
\htmlonly
<div hidden>
\endhtmlonly
\latexonly
\iffalse
\endlatexonly
\subpage config_project
\subpage config_domain
\subpage config_processes
\subpage config_resolution
\subpage config_optimize
\subpage config_time
\subpage config_lut
\subpage config_input
\subpage config_coupling
\subpage config_meteo
\subpage config_mpr
\subpage config_mhm
\subpage config_mlm
\subpage config_mrm
\subpage config_observations
\subpage parameter_migration_v6
\subpage interception_1
\subpage snow_1
\subpage soil_moisture_1
\subpage soil_moisture_2
\subpage soil_moisture_3
\subpage soil_moisture_4
\subpage direct_runoff_1
\subpage pet_m2
\subpage pet_m1
\subpage pet_1
\subpage pet_2
\subpage pet_3
\subpage interflow_1
\subpage percolation_1
\subpage routing_1
\subpage routing_2
\subpage routing_3
\subpage neutrons_1
\subpage neutrons_2
\subpage baseflow_1
\subpage lakes_1
\subpage river_temperature_1
\subpage output_mhm
\subpage output_mlm
\subpage output_mrm
\latexonly
\fi
\endlatexonly
\htmlonly
</div>
\endhtmlonly
-->
