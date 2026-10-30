# mLM configuration {#config_mlm}

[TOC]

Configuration for the v6 lake component (mLM).

**Namelist**: `config_mlm`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [n_lakes](#n_lakes) | integer array | no | no | Number of lakes |
| [output_path](#output_path) | string array | no | no | mLM lake-point output file |
| [read_restart](#read_restart) | logical array | no | no | Read mLM restart |
| [restart_input_path](#restart_input_path) | string array | no | no | mLM restart input file |
| [write_restart](#write_restart) | logical array | no | no | Write mLM restart |
| [restart_save_topology](#restart_save_topology) | logical array | no | no | Include level-0 lake topology in restart |
| [restart_output_path](#restart_output_path) | string array | no | no | mLM restart output file |
| [force_initialize_lake_at_restart](#force_initialize_lake_at_restart) | logical | no | no | Force lake level initialization at restart |
| [force_initialize_lake_at_a_date](#force_initialize_lake_at_a_date) | logical | no | no | Force lake level initialization at a date |
| [evaporation_method](#evaporation_method) | integer | no | no | Lake evaporation method |
| [input_windspeed_height](#input_windspeed_height) | real | no | no | Input windspeed height [m] |
| [climatological_windspeed](#climatological_windspeed) | real | no | no | Climatological 2 m windspeed [m s-1] |
| [varying_albedo](#varying_albedo) | logical | no | no | Variable water albedo |
| [varying_lambda](#varying_lambda) | logical | no | no | Variable latent heat of vaporisation |
| [varying_psychro](#varying_psychro) | logical | no | no | Variable psychrometric constant |
| [albedo_method](#albedo_method) | integer | no | no | Water albedo calculation method |
| [bathymetry_type](#bathymetry_type) | integer array | no | no | Lake hAV curve type |
| [natural](#natural) | logical array | no | no | Natural lake flag |

## Field details

### n_lakes

Number of lakes `n_lakes`

Number of lakes for each configured domain.

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Declared required: no
- Input required: no
- Minimum: `>= 1`
- Examples: `[1]`

### output_path

mLM lake-point output file `output_path`

Lake-point output file for each domain, given with its full path.

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

mLM restart input file `restart_input_path`

mLM restart input file for each domain, given with its full path.

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

### restart_save_topology

Include level-0 lake topology in restart `restart_save_topology`

Include the lake footprint grid and map when writing a restart. Reading detects
topology independently. Compact restarts retain per-lake footprint areas;
meteo restart support permits their use without the original level-0 geometry.

Summary:
- Type: `logical, dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `.false.`

### restart_output_path

mLM restart output file `restart_output_path`

mLM restart output file for each domain, given with its full path.

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
- 3 : Calibration-free Priestley-Taylor, after Su and Singh (2023), adapted
  to the open-water Penman formulation. The constant Priestley-Taylor coefficient is
  replaced by an analytical coefficient that varies with available energy, vapour-pressure
  deficit and temperature, evaluated with the climatological windspeed (climatological_windspeed)
  rather than the instantaneous wind so that noisy per-step wind does not swing the evaporation.
  Needs the same forcing as Penman: temperature, downward short- and long-wave radiation,
  windspeed and dew-point temperature.

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `1`
- Allowed values: `1`, `2`, `3`

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

### climatological_windspeed

Climatological 2 m windspeed [m s-1] `climatological_windspeed`

Long-term mean (climatological) windspeed at 2 m, used only by the calibration-free Priestley-Taylor
method (evaporation_method = 3) inside its coefficient, in place of the instantaneous wind.
Decoupling the coefficient from the per-step wind keeps the evaporation robust to noisy
instantaneous windspeed. Ignored by the Penman (1) and Priestley-Taylor (2) methods.

Summary:
- Type: `real(dp)`
- Declared required: no
- Input required: no
- Default: `2.0`
- Minimum: `>= 0.0`
- Examples: `1.5`, `3.0`

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

### bathymetry_type

Lake hAV curve type `bathymetry_type`

Type of height-area-volume (hAV) relation used per lake.
- 0 : input bathymetry curves.
- 1 : linear h~V relation with constant area.
- 2 : Liebe (2005), one half of a pyramid with square base cut diagonally at the base.
  Reference: Liebe, J. (2005). Estimation of water storage capacity and evaporation
  losses of small reservoirs in the Upper East Region of Ghana. MSc thesis, Rheinische
  Friedrich-Wilhelms-Universitaet Bonn.

Summary:
- Type: `integer(i4), dimension(n_domains, max_lakes)`
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
- Type: `logical, dimension(n_domains, max_lakes)`
- Declared required: no
- Input required: no
- Default: `.true.`

## Example

```fortran
&config_mlm
  n_lakes(:) = 1
  output_path(:) = ""
  read_restart(:) = .false.
  restart_input_path(:) = ""
  write_restart(:) = .false.
  restart_save_topology(:) = .false.
  restart_output_path(:) = ""
  force_initialize_lake_at_restart = .false.
  force_initialize_lake_at_a_date = .false.
  evaporation_method = 1
  input_windspeed_height = 10.0
  climatological_windspeed = 1.5
  varying_albedo = .true.
  varying_lambda = .true.
  varying_psychro = .true.
  albedo_method = 2
  bathymetry_type(:,:) = 0
  natural(:,:) = .true.
/
```

