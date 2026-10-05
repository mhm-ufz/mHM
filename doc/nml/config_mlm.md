# mLM configuration {#config_mlm}

[TOC]

Configuration for the v6 lake component (mLM).

**Namelist**: `config_mlm`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [input_path](#input_path) | string array | no | no | mLM lake input path |
| [output_path](#output_path) | string array | no | no | mLM lake-point output path |
| [read_restart](#read_restart) | logical array | no | no | Read mLM restart |
| [restart_input_path](#restart_input_path) | string array | no | no | mLM restart input path |
| [write_restart](#write_restart) | logical array | no | no | Write mLM restart |
| [restart_output_path](#restart_output_path) | string array | no | no | mLM restart output path |
| [force_initialize_lake_at_restart](#force_initialize_lake_at_restart) | logical | no | no | Force lake level initialization at restart |
| [force_initialize_lake_at_a_date](#force_initialize_lake_at_a_date) | logical | no | no | Force lake level initialization at a date |
| [evaporation_method](#evaporation_method) | integer | no | no | Lake evaporation method |
| [input_windspeed_height](#input_windspeed_height) | real | no | no | Input windspeed height [m] |
| [varying_albedo](#varying_albedo) | logical | no | no | Variable water albedo |
| [varying_lambda](#varying_lambda) | logical | no | no | Variable latent heat of vaporisation |
| [varying_psychro](#varying_psychro) | logical | no | no | Variable psychrometric constant |
| [albedo_method](#albedo_method) | integer | no | no | Water albedo calculation method |
| [hanasaki_2006](#hanasaki_2006) | logical | no | no | Hanasaki et al. (2006) release target |
| [point_demand](#point_demand) | logical | no | no | Point demand input |
| [bathymetry_type](#bathymetry_type) | integer array | no | no | Lake hAV curve type |
| [natural](#natural) | logical array | no | no | Natural lake flag |
| [use_watersupply](#use_watersupply) | logical array | no | no | Water-supply use flag |
| [use_irrigation](#use_irrigation) | logical array | no | no | Irrigation use flag |
| [service_distance](#service_distance) | real array | no | no | Downstream service distance [km] |
| [bathymetry_filename](#bathymetry_filename) | string array | no | no | Lake bathymetry file |
| [release_target_filename](#release_target_filename) | string array | no | no | Release target file |
| [watersupply_demand_filename](#watersupply_demand_filename) | string array | no | no | Water-supply demand file |
| [irrigation_demand_filename](#irrigation_demand_filename) | string array | no | no | Irrigation demand file |
| [ibt_filename](#ibt_filename) | string array | no | no | Inter-basin transfer file |
| [level_filename](#level_filename) | string array | no | no | Observed water-level file |
| [area_filename](#area_filename) | string array | no | no | Observed surface-area file |
| [spill_filename](#spill_filename) | string array | no | no | Observed spill file |
| [print_lake_maps](#print_lake_maps) | logical | no | no | Print lake maps |
| [print_links_table](#print_links_table) | logical | no | no | Print links table |

## Field details

### input_path

mLM lake input path `input_path`

Directory where the lake input data for each domain is located.

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no

### output_path

mLM lake-point output path `output_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no

### read_restart

Read mLM restart `read_restart`

Summary:
- Type: `logical, dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `.false.`

### restart_input_path

mLM restart input path `restart_input_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no

### write_restart

Write mLM restart `write_restart`

Summary:
- Type: `logical, dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `.false.`

### restart_output_path

mLM restart output path `restart_output_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no

### force_initialize_lake_at_restart

Force lake level initialization at restart `force_initialize_lake_at_restart`

Flag for lake water level initialization for a simulation started from a restart file.
Only used if the mLM restart is read (read_restart = true).
- if true : forces lake_init_level as the initial lake level (use during warm-up based restart).
- if false : initializes the lake level from the restart file (use for chronology runs,
  forecasting, or restart tests during development).

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### force_initialize_lake_at_a_date

Force lake level initialization at a date `force_initialize_lake_at_a_date`

Flag for lake water level initialization at the start of lake_init_date.
A valid lake_init_date must be provided when this is enabled.
- if true : forces initial_lake_level as the initial lake level at the start of lake_init_date.
- if false : no forced initialization.

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### evaporation_method

Lake evaporation method `evaporation_method`

Method used to estimate open-water evaporation from the lake surface.
- 1 : Penman (1948). Needs temperature, downward short- and long-wave radiation,
  windspeed and dew-point temperature.
- 2 : Priestley-Taylor (1972). Needs temperature and downward short- and long-wave
  radiation only (no windspeed or dew-point temperature).

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `1`
- Allowed values: `1`, `2`

### input_windspeed_height

Input windspeed height [m] `input_windspeed_height`

mLM requires windspeed data at 2 m height. If the input windspeed data is provided at a
height other than 2 m, mLM applies the Hellmann exponential law to derive the 2 m windspeed.
Enter the height for which the input windspeed data is provided, in metres.

Summary:
- Type: `real(dp)`
- Declared required: no
- Input required: no
- Default: `2.0`
- Minimum: `>= 0.0`
- Examples: `10.0`

### varying_albedo

Variable water albedo `varying_albedo`

Use a variable water albedo in the evaporation calculation. The varying albedo depends on
the solar geometry and therefore needs the geographic latitude/longitude of each lake. On a
spherical (lon/lat) grid these come from the lake point coordinates directly; on a projected
(cartesian X/Y) grid they are taken from the level-0 geographic auxiliary coordinates, which
must be supplied through the latlon file (config_input: latlon_path, with morph_lat_var /
morph_lon_var, default lat_l0 / lon_l0). The run stops if varying_albedo is enabled but the
lake geographic coordinates are unavailable.

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.true.`

### varying_lambda

Variable latent heat of vaporisation `varying_lambda`

Use a variable latent heat of vaporisation of water in the evaporation calculation.

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.true.`

### varying_psychro

Variable psychrometric constant `varying_psychro`

Use a variable psychrometric constant in the evaporation calculation.

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.true.`

### albedo_method

Water albedo calculation method `albedo_method`

Method of water albedo calculation. Only used if varying_albedo is true.
- 1 : Fresnel equation.
- 2 : Grishchenko (1959). See Appendix C, Figure C1 in Shrestha et al. (2024),
  Water Resources Research, for the varying albedo calculation.

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `2`
- Allowed values: `1`, `2`

### hanasaki_2006

Hanasaki et al. (2006) release target `hanasaki_2006`

Use the temperature release target after Hanasaki et al. (2006).
If enabled, the per-lake reservoir capacity and mean annual inflow must be provided.
Reference: Hanasaki, N., Kanae, S., & Oki, T. (2006). A reservoir operation scheme
for global river routing models. Journal of Hydrology, 327(1-2), 22-41.
doi:10.1016/j.jhydrol.2005.11.011

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.true.`

### point_demand

Point demand input `point_demand`

Flag deciding whether demand data is provided as point time series or as a space-time
field (netCDF).
- if true : lake_filename_Demand_<demandtype> must be provided for each activated use flag.
- if false : demand field netCDF files must be present in the demand input folder.

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### bathymetry_type

Lake hAV curve type `bathymetry_type`

Type of height-area-volume (hAV) relation used per lake.
- 0 : input bathymetry curves (bathymetry_filename must be provided).
- 1 : linear h~V relation with constant area.
- 2 : Liebe (2005), one half of a pyramid with square base cut diagonally at the base.
  Reference: Liebe, J. (2005). Estimation of water storage capacity and evaporation
  losses of small reservoirs in the Upper East Region of Ghana. MSc thesis, Rheinische
  Friedrich-Wilhelms-Universitaet Bonn.

Summary:
- Type: `integer(i4), dimension(n_domains, n_lakes)`
- Declared required: no
- Input required: no
- Default: `0`
- Allowed values: `0`, `1`, `2`

### natural

Natural lake flag `natural`

Whether the lake is a naturally existing lake (true) or an artificial reservoir (false).
A natural lake cannot have a commission date and its water elevation may rise above the
dam crest level to produce spill; a reservoir's elevation is capped at the dam crest level.

Summary:
- Type: `logical, dimension(n_domains, n_lakes)`
- Declared required: no
- Input required: no
- Default: `.true.`

### use_watersupply

Water-supply use flag `use_watersupply`

Enable consumptive water-supply use for the lake.

Summary:
- Type: `logical, dimension(n_domains, n_lakes)`
- Declared required: no
- Input required: no
- Default: `.false.`

### use_irrigation

Irrigation use flag `use_irrigation`

Enable consumptive irrigation use for the lake.

Summary:
- Type: `logical, dimension(n_domains, n_lakes)`
- Declared required: no
- Input required: no
- Default: `.false.`

### service_distance

Downstream service distance [km] `service_distance`

Downstream distance over which the reservoir provides irrigation and water-supply service.
If not defined or set to -9999, a default value of 600 km is used.

Summary:
- Type: `real(dp), dimension(n_domains, n_lakes)`
- Declared required: no
- Input required: no
- Default: `600.0`

### bathymetry_filename

Lake bathymetry file `bathymetry_filename`

Name of the file with the lake bathymetry table. Required when bathymetry_type is 0.
The table must have non-repeating values in ascending order (lower to higher elevations).

Summary:
- Type: `character(len=buf), dimension(n_domains, n_lakes)`
- Item format: `file-path`
- Declared required: no
- Input required: no

### release_target_filename

Release target file `release_target_filename`

Name of the file with the downstream lake release-target discharge time series.

Summary:
- Type: `character(len=buf), dimension(n_domains, n_lakes)`
- Item format: `file-path`
- Declared required: no
- Input required: no

### watersupply_demand_filename

Water-supply demand file `watersupply_demand_filename`

Name of the file with the water-supply demand discharge time series.

Summary:
- Type: `character(len=buf), dimension(n_domains, n_lakes)`
- Item format: `file-path`
- Declared required: no
- Input required: no

### irrigation_demand_filename

Irrigation demand file `irrigation_demand_filename`

Name of the file with the irrigation demand discharge time series.

Summary:
- Type: `character(len=buf), dimension(n_domains, n_lakes)`
- Item format: `file-path`
- Declared required: no
- Input required: no

### ibt_filename

Inter-basin transfer file `ibt_filename`

Name of the file with the optional inter-basin transfer inflow time series.

Summary:
- Type: `character(len=buf), dimension(n_domains, n_lakes)`
- Item format: `file-path`
- Declared required: no
- Input required: no

### level_filename

Observed water-level file `level_filename`

Name of the file with the observed lake water-level time series.

Summary:
- Type: `character(len=buf), dimension(n_domains, n_lakes)`
- Item format: `file-path`
- Declared required: no
- Input required: no

### area_filename

Observed surface-area file `area_filename`

Name of the file with the observed lake surface-area time series.

Summary:
- Type: `character(len=buf), dimension(n_domains, n_lakes)`
- Item format: `file-path`
- Declared required: no
- Input required: no

### spill_filename

Observed spill file `spill_filename`

Name of the file with the observed lake spill time series.

Summary:
- Type: `character(len=buf), dimension(n_domains, n_lakes)`
- Item format: `file-path`
- Declared required: no
- Input required: no

### print_lake_maps

Print lake maps `print_lake_maps`

Print maps of the lake, lake catchment, lake demand catchment, and lake circular catchment.
Optional output used for model verification.

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.true.`

### print_links_table

Print links table `print_links_table`

Print the table of link (river-lake, river-node) network variables.
Optional output used for model verification.

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.true.`

## Example

```fortran
&config_mlm
  input_path(:) = ""
  output_path(:) = ""
  read_restart(:) = .false.
  restart_input_path(:) = ""
  write_restart(:) = .false.
  restart_output_path(:) = ""
  force_initialize_lake_at_restart = .false.
  force_initialize_lake_at_a_date = .false.
  evaporation_method = 1
  input_windspeed_height = 10.0
  varying_albedo = .true.
  varying_lambda = .true.
  varying_psychro = .true.
  albedo_method = 2
  hanasaki_2006 = .true.
  point_demand = .false.
  bathymetry_type(:,:) = 0
  natural(:,:) = .true.
  use_watersupply(:,:) = .false.
  use_irrigation(:,:) = .false.
  service_distance(:,:) = 600.0
  bathymetry_filename(:,:) = ""
  release_target_filename(:,:) = ""
  watersupply_demand_filename(:,:) = ""
  irrigation_demand_filename(:,:) = ""
  ibt_filename(:,:) = ""
  level_filename(:,:) = ""
  area_filename(:,:) = ""
  spill_filename(:,:) = ""
  print_lake_maps = .true.
  print_links_table = .true.
/
```

