# Coupling configuration {#config_coupling}

[TOC]

Coupling grid and variable settings for meteorological, hydrological, and morphology components.
Arrays are indexed by domain (dimension 1).
Grid parameters define resolutions and coordinate systems for coupled runs.
Coupled flags indicate whether inputs are provided by an external model.

**Namelist**: `config_coupling`

## Fields

| Name | Type | Required | Info |
| --- | --- | --- | --- |
| [meteo_grid_nx](#meteo_grid_nx) | integer array | no | Meteo grid size in x-direction |
| [meteo_grid_ny](#meteo_grid_ny) | integer array | no | Meteo grid size in y-direction |
| [meteo_grid_xll](#meteo_grid_xll) | real array | no | Meteo grid x origin |
| [meteo_grid_yll](#meteo_grid_yll) | real array | no | Meteo grid y origin |
| [meteo_grid_cellsize](#meteo_grid_cellsize) | real array | no | Meteo grid cell size |
| [meteo_grid_ydir](#meteo_grid_ydir) | integer array | no | Meteo grid y direction |
| [meteo_grid_coordsys](#meteo_grid_coordsys) | integer array | no | Meteo grid coordinate system |
| [hydro_grid_nx](#hydro_grid_nx) | integer array | no | Hydro grid size in x-direction |
| [hydro_grid_ny](#hydro_grid_ny) | integer array | no | Hydro grid size in y-direction |
| [hydro_grid_xll](#hydro_grid_xll) | real array | no | Hydro grid x origin |
| [hydro_grid_yll](#hydro_grid_yll) | real array | no | Hydro grid y origin |
| [hydro_grid_cellsize](#hydro_grid_cellsize) | real array | no | Hydro grid cell size |
| [hydro_grid_ydir](#hydro_grid_ydir) | integer array | no | Hydro grid y direction |
| [hydro_grid_coordsys](#hydro_grid_coordsys) | integer array | no | Hydro grid coordinate system |
| [morph_grid_nx](#morph_grid_nx) | integer array | no | Morph grid size in x-direction |
| [morph_grid_ny](#morph_grid_ny) | integer array | no | Morph grid size in y-direction |
| [morph_grid_xll](#morph_grid_xll) | real array | no | Morph grid x origin |
| [morph_grid_yll](#morph_grid_yll) | real array | no | Morph grid y origin |
| [morph_grid_cellsize](#morph_grid_cellsize) | real array | no | Morph grid cell size |
| [morph_grid_ydir](#morph_grid_ydir) | integer array | no | Morph grid y direction |
| [morph_grid_coordsys](#morph_grid_coordsys) | integer array | no | Morph grid coordinate system |
| [pre_coupled](#pre_coupled) | logical array | no | Precipitation coupled |
| [pet_coupled](#pet_coupled) | logical array | no | Potential evapotranspiration coupled |
| [temp_coupled](#temp_coupled) | logical array | no | Air temperature coupled |
| [tann_coupled](#tann_coupled) | logical array | no | Air temperature annual mean coupled |
| [tmin_coupled](#tmin_coupled) | logical array | no | Air temperature daily minimum coupled |
| [tmax_coupled](#tmax_coupled) | logical array | no | Air temperature daily maximum coupled |
| [ssrd_coupled](#ssrd_coupled) | logical array | no | Surface shortwave radiation coupled |
| [strd_coupled](#strd_coupled) | logical array | no | Surface thermal radiation coupled |
| [netrad_coupled](#netrad_coupled) | logical array | no | Net radiation coupled |
| [eabs_coupled](#eabs_coupled) | logical array | no | Vapor pressure coupled |
| [wind_coupled](#wind_coupled) | logical array | no | Wind speed coupled |
| [runoff_coupled](#runoff_coupled) | logical array | no | Runoff coupled |
| [runoff_sealed_coupled](#runoff_sealed_coupled) | logical array | no | Sealed runoff coupled |
| [interflow_fast_coupled](#interflow_fast_coupled) | logical array | no | Fast interflow coupled |
| [interflow_slow_coupled](#interflow_slow_coupled) | logical array | no | Slow interflow coupled |
| [baseflow_coupled](#baseflow_coupled) | logical array | no | Baseflow coupled |
| [dem_coupled](#dem_coupled) | logical array | no | DEM coupled |
| [slope_coupled](#slope_coupled) | logical array | no | Slope coupled |
| [aspect_coupled](#aspect_coupled) | logical array | no | Aspect coupled |
| [geo_class_coupled](#geo_class_coupled) | logical array | no | Geology class coupled |
| [soil_class_coupled](#soil_class_coupled) | logical array | no | Soil class coupled |
| [lai_class_coupled](#lai_class_coupled) | logical array | no | LAI class coupled |
| [river_width_coupled](#river_width_coupled) | logical array | no | River width coupled |
| [meteo_mask_coupled](#meteo_mask_coupled) | logical array | no | Meteorological mask coupled |
| [hydro_mask_coupled](#hydro_mask_coupled) | logical array | no | Hydrological mask coupled |
| [morph_mask_coupled](#morph_mask_coupled) | logical array | no | Morphology mask coupled |
| [hydro_latlon_coupled](#hydro_latlon_coupled) | logical array | no | Hydrological latlon coupled |
| [morph_latlon_coupled](#morph_latlon_coupled) | logical array | no | Morphological latlon coupled |
| [route_latlon_coupled](#route_latlon_coupled) | logical array | no | Routing latlon coupled |

## Field details

### meteo_grid_nx

Meteo grid size in x-direction `meteo_grid_nx`

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Minimum: `>= 1`
- Examples: `[10]`

### meteo_grid_ny

Meteo grid size in y-direction `meteo_grid_ny`

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Minimum: `>= 1`
- Examples: `[10]`

### meteo_grid_xll

Meteo grid x origin `meteo_grid_xll`

Summary:
- Type: `real(dp), dimension(n_domains)`
- Required: no

### meteo_grid_yll

Meteo grid y origin `meteo_grid_yll`

Summary:
- Type: `real(dp), dimension(n_domains)`
- Required: no

### meteo_grid_cellsize

Meteo grid cell size `meteo_grid_cellsize`

Summary:
- Type: `real(dp), dimension(n_domains)`
- Required: no
- Minimum: `> 0.0`
- Examples: `[1.0]`

### meteo_grid_ydir

Meteo grid y direction `meteo_grid_ydir`

Direction of the y-axis for the meteorological grid.
- 0 : top-down y-axis
- 1 : bottom-up y-axis

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Default: `0`
- Allowed values: `0`, `1`

### meteo_grid_coordsys

Meteo grid coordinate system `meteo_grid_coordsys`

Coordinate system for the meteorological grid.
- 0 : cartesian
- 1 : latlon

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Default: `0`
- Allowed values: `0`, `1`

### hydro_grid_nx

Hydro grid size in x-direction `hydro_grid_nx`

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Minimum: `>= 1`
- Examples: `[10]`

### hydro_grid_ny

Hydro grid size in y-direction `hydro_grid_ny`

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Minimum: `>= 1`
- Examples: `[10]`

### hydro_grid_xll

Hydro grid x origin `hydro_grid_xll`

Summary:
- Type: `real(dp), dimension(n_domains)`
- Required: no

### hydro_grid_yll

Hydro grid y origin `hydro_grid_yll`

Summary:
- Type: `real(dp), dimension(n_domains)`
- Required: no

### hydro_grid_cellsize

Hydro grid cell size `hydro_grid_cellsize`

Summary:
- Type: `real(dp), dimension(n_domains)`
- Required: no
- Minimum: `> 0.0`
- Examples: `[1.0]`

### hydro_grid_ydir

Hydro grid y direction `hydro_grid_ydir`

Direction of the y-axis for the hydrological grid.
- 0 : top-down y-axis
- 1 : bottom-up y-axis

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Default: `0`
- Allowed values: `0`, `1`

### hydro_grid_coordsys

Hydro grid coordinate system `hydro_grid_coordsys`

Coordinate system for the hydrological grid.
- 0 : cartesian
- 1 : latlon

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Default: `0`
- Allowed values: `0`, `1`

### morph_grid_nx

Morph grid size in x-direction `morph_grid_nx`

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Minimum: `>= 1`
- Examples: `[10]`

### morph_grid_ny

Morph grid size in y-direction `morph_grid_ny`

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Minimum: `>= 1`
- Examples: `[10]`

### morph_grid_xll

Morph grid x origin `morph_grid_xll`

Summary:
- Type: `real(dp), dimension(n_domains)`
- Required: no

### morph_grid_yll

Morph grid y origin `morph_grid_yll`

Summary:
- Type: `real(dp), dimension(n_domains)`
- Required: no

### morph_grid_cellsize

Morph grid cell size `morph_grid_cellsize`

Summary:
- Type: `real(dp), dimension(n_domains)`
- Required: no
- Minimum: `> 0.0`
- Examples: `[1.0]`

### morph_grid_ydir

Morph grid y direction `morph_grid_ydir`

Direction of the y-axis for the morphology grid.
- 0 : top-down y-axis
- 1 : bottom-up y-axis

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Default: `0`
- Allowed values: `0`, `1`

### morph_grid_coordsys

Morph grid coordinate system `morph_grid_coordsys`

Coordinate system for the morphology grid.
- 0 : cartesian
- 1 : latlon

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Default: `0`
- Allowed values: `0`, `1`

### pre_coupled

Precipitation coupled `pre_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### pet_coupled

Potential evapotranspiration coupled `pet_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### temp_coupled

Air temperature coupled `temp_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### tann_coupled

Air temperature annual mean coupled `tann_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### tmin_coupled

Air temperature daily minimum coupled `tmin_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### tmax_coupled

Air temperature daily maximum coupled `tmax_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### ssrd_coupled

Surface shortwave radiation coupled `ssrd_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### strd_coupled

Surface thermal radiation coupled `strd_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### netrad_coupled

Net radiation coupled `netrad_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### eabs_coupled

Vapor pressure coupled `eabs_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### wind_coupled

Wind speed coupled `wind_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### runoff_coupled

Runoff coupled `runoff_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### runoff_sealed_coupled

Sealed runoff coupled `runoff_sealed_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### interflow_fast_coupled

Fast interflow coupled `interflow_fast_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### interflow_slow_coupled

Slow interflow coupled `interflow_slow_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### baseflow_coupled

Baseflow coupled `baseflow_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### dem_coupled

DEM coupled `dem_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### slope_coupled

Slope coupled `slope_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### aspect_coupled

Aspect coupled `aspect_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### geo_class_coupled

Geology class coupled `geo_class_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### soil_class_coupled

Soil class coupled `soil_class_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### lai_class_coupled

LAI class coupled `lai_class_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### river_width_coupled

River width coupled `river_width_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### meteo_mask_coupled

Meteorological mask coupled `meteo_mask_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### hydro_mask_coupled

Hydrological mask coupled `hydro_mask_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### morph_mask_coupled

Morphology mask coupled `morph_mask_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### hydro_latlon_coupled

Hydrological latlon coupled `hydro_latlon_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### morph_latlon_coupled

Morphological latlon coupled `morph_latlon_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### route_latlon_coupled

Routing latlon coupled `route_latlon_coupled`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

## Example

```fortran
&config_coupling
  meteo_grid_nx(:) = 10
  meteo_grid_ny(:) = 10
  meteo_grid_xll(:) = 0.0
  meteo_grid_yll(:) = 0.0
  meteo_grid_cellsize(:) = 1.0
  meteo_grid_ydir(:) = 0
  meteo_grid_coordsys(:) = 0
  hydro_grid_nx(:) = 10
  hydro_grid_ny(:) = 10
  hydro_grid_xll(:) = 0.0
  hydro_grid_yll(:) = 0.0
  hydro_grid_cellsize(:) = 1.0
  hydro_grid_ydir(:) = 0
  hydro_grid_coordsys(:) = 0
  morph_grid_nx(:) = 10
  morph_grid_ny(:) = 10
  morph_grid_xll(:) = 0.0
  morph_grid_yll(:) = 0.0
  morph_grid_cellsize(:) = 1.0
  morph_grid_ydir(:) = 0
  morph_grid_coordsys(:) = 0
  pre_coupled(:) = .false.
  pet_coupled(:) = .false.
  temp_coupled(:) = .false.
  tann_coupled(:) = .false.
  tmin_coupled(:) = .false.
  tmax_coupled(:) = .false.
  ssrd_coupled(:) = .false.
  strd_coupled(:) = .false.
  netrad_coupled(:) = .false.
  eabs_coupled(:) = .false.
  wind_coupled(:) = .false.
  runoff_coupled(:) = .false.
  runoff_sealed_coupled(:) = .false.
  interflow_fast_coupled(:) = .false.
  interflow_slow_coupled(:) = .false.
  baseflow_coupled(:) = .false.
  dem_coupled(:) = .false.
  slope_coupled(:) = .false.
  aspect_coupled(:) = .false.
  geo_class_coupled(:) = .false.
  soil_class_coupled(:) = .false.
  lai_class_coupled(:) = .false.
  river_width_coupled(:) = .false.
  meteo_mask_coupled(:) = .false.
  hydro_mask_coupled(:) = .false.
  morph_mask_coupled(:) = .false.
  hydro_latlon_coupled(:) = .false.
  morph_latlon_coupled(:) = .false.
  route_latlon_coupled(:) = .false.
/
```

