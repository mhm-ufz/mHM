# MPR configuration {#config_mpr}

[TOC]

Configuration for the multiscale parameter regionalization in mHM.

**Namelist**: `config_mpr`

## Fields

| Name | Type | Required | Info |
| --- | --- | --- | --- |
| [soil_db_mode](#soil_db_mode) | integer array | no | Soil database mode |
| [tillage_depth](#tillage_depth) | integer array | no | Tillage depth |
| [n_layers](#n_layers) | integer array | no | Number of soil layers |
| [soil_depth](#soil_depth) | integer array | no | Soil horizon depth |
| [fracSealed_cityArea](#fracsealed_cityarea) | real array | no | Sealed fraction of city area |
| [land_cover_path](#land_cover_path) | string array | no | Land cover path |
| [land_cover_var](#land_cover_var) | string array | no | Land cover variable |
| [lai_time_step](#lai_time_step) | integer array | no | LAI time step |
| [lai_path](#lai_path) | string array | no | LAI path |
| [lai_var](#lai_var) | string array | no | LAI variable |
| [soil_lut_path](#soil_lut_path) | string array | no | Soil LUT path |
| [geo_lut_path](#geo_lut_path) | string array | no | Geology LUT path |
| [lai_lut_path](#lai_lut_path) | string array | no | LAI LUT path |
| [read_restart](#read_restart) | logical array | no | Read restart |
| [restart_input_path](#restart_input_path) | string array | no | Restart input path |
| [write_restart](#write_restart) | logical array | no | Write restart |
| [restart_output_path](#restart_output_path) | string array | no | Restart output path |

## Field details

### soil_db_mode

Soil database mode `soil_db_mode`

Flag to handle multiple soil database types; valid for all domains.
- 0: classical mHM soil database with soil_class.asc and soil_classdefinition.txt.
  Soil horizons follow nSoilHorizons_mHM and soil_Depth for n-1 layers; last layer depth comes from the LUT.
- 1: harmonised horizon-specific soil database with soil_class_horizon_XX.asc and soil_classdefinition_iFlag_soilDB_1.txt.
  Soil depth is provided for all horizons; the last horizon depth is fixed and uniform;
  tillage depth must match a horizon depth.

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Default: `0`
- Allowed values: `0`, `1`

### tillage_depth

Tillage depth `tillage_depth`

Soil depth down to which organic matter is possible [mm].

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no

### n_layers

Number of soil layers `n_layers`

Number of mHM soil layers for each configured domain.

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Minimum: `>= 1`

### soil_depth

Soil horizon depth `soil_depth`

Bottom depth of soil horizons with respect to the ground surface [mm], positive downwards.
- if soil_db_mode = 0: provide depths for horizons 1..n-1; last horizon depth comes from the LUT.
- if soil_db_mode = 1: provide depths for all horizons 1..n; a soil_class_horizon_XX.asc is required for each.
Tillage depth should match one of the specified horizon depths.

Summary:
- Type: `integer(i4), dimension(max_layers, n_domains)`
- Required: no
- Default: `0`

### fracSealed_cityArea

Sealed fraction of city area `fracSealed_cityArea`

Fraction of area within city assumed to be fully sealed [0.0-1.0].

Summary:
- Type: `real(dp), dimension(n_domains)`
- Required: no

### land_cover_path

Land cover path `land_cover_path`

NetCDF land-cover dataset path.

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### land_cover_var

Land cover variable `land_cover_var`

Land-cover variable name in the NetCDF dataset.

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"land_cover"`

### lai_time_step

LAI time step `lai_time_step`

Time step for LAI input data [days]:
- if = 1: annual cycle monthly gridded LAI values;
- if = 0: no LAI input data used, static LAI from LUT;
- if = -1 : daily gridded LAI input data;
- if = -2 : monthly gridded LAI input data;
- if = -3 : yearly gridded LAI input data;

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Allowed values: `-3`, `-2`, `-1`, `0`, `1`

### lai_path

LAI path `lai_path`

Gridded LAI dataset path (if `lai_time_step` < 0 or = 1).

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### lai_var

LAI variable `lai_var`

LAI variable name in the gridded NetCDF dataset.

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"lai"`

### soil_lut_path

Soil LUT path `soil_lut_path`

Soil look-up table path.

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### geo_lut_path

Geology LUT path `geo_lut_path`

Geology look-up table path.

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### lai_lut_path

LAI LUT path `lai_lut_path`

LAI look-up table path.

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### read_restart

Read restart `read_restart`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### restart_input_path

Restart input path `restart_input_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Examples: `["mpr_restart_in.nc"]`

### write_restart

Write restart `write_restart`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### restart_output_path

Restart output path `restart_output_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Examples: `["mpr_restart_out.nc"]`

## Example

```fortran
&config_mpr
  soil_db_mode(:) = 0
  tillage_depth(:) = 0
  n_layers(:) = 0
  soil_depth(:,:) = 0
  fracSealed_cityArea(:) = 0.0
  land_cover_path(:) = ""
  land_cover_var(:) = "land_cover"
  lai_time_step(:) = -3
  lai_path(:) = ""
  lai_var(:) = "lai"
  soil_lut_path(:) = ""
  geo_lut_path(:) = ""
  lai_lut_path(:) = ""
  read_restart(:) = .false.
  restart_input_path(:) = "mpr_restart_in.nc"
  write_restart(:) = .false.
  restart_output_path(:) = "mpr_restart_out.nc"
/
```

