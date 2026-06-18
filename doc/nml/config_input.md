# Input configuration {#config_input}

[TOC]

Paths and variable names for input data used by mHM.
Arrays are indexed by domain (dimension 1). Most paths are optional.
Variable name entries define the NetCDF variable names to read.

**Namelist**: `config_input`

## Fields

| Name | Type | Required | Info |
| --- | --- | --- | --- |
| [chunking](#chunking) | integer array | no | Chunking for input data |
| [time_stamp_location](#time_stamp_location) | integer array | no | NetCDF time-stamp location |
| [latlon_path](#latlon_path) | string array | no | Latlon specification file path |
| [morph_latlon](#morph_latlon) | logical array | no | DEM in latlon coordinates |
| [pre_path](#pre_path) | string array | no | Precipitation input |
| [pet_path](#pet_path) | string array | no | Potential evapotranspiration input |
| [temp_path](#temp_path) | string array | no | Air temperature input |
| [tann_path](#tann_path) | string array | no | Air temperature annual mean input |
| [tmin_path](#tmin_path) | string array | no | Air temperature daily minimum input |
| [tmax_path](#tmax_path) | string array | no | Air temperature daily maximum input |
| [ssrd_path](#ssrd_path) | string array | no | Surface shortwave radiation downwards input |
| [strd_path](#strd_path) | string array | no | Surface thermal radiation downwards input |
| [netrad_path](#netrad_path) | string array | no | Net radiation input |
| [eabs_path](#eabs_path) | string array | no | Vapor pressure input |
| [wind_path](#wind_path) | string array | no | Wind speed input |
| [meteo_mask_path](#meteo_mask_path) | string array | no | Meteorological mask file path |
| [runoff_path](#runoff_path) | string array | no | Runoff input |
| [runoff_sealed_path](#runoff_sealed_path) | string array | no | Sealed runoff input |
| [interflow_fast_path](#interflow_fast_path) | string array | no | Fast interflow input |
| [interflow_slow_path](#interflow_slow_path) | string array | no | Slow interflow input |
| [baseflow_path](#baseflow_path) | string array | no | Baseflow input |
| [hydro_mask_path](#hydro_mask_path) | string array | no | Hydrological mask file path |
| [dem_path](#dem_path) | string array | no | DEM input |
| [slope_path](#slope_path) | string array | no | Slope input |
| [aspect_path](#aspect_path) | string array | no | Aspect input |
| [fdir_path](#fdir_path) | string array | no | Flow direction input |
| [facc_path](#facc_path) | string array | no | Flow accumulation input |
| [geo_class_path](#geo_class_path) | string array | no | Geology class input |
| [soil_class_path](#soil_class_path) | string array | no | Soil class input |
| [soil_horizon_class_path](#soil_horizon_class_path) | string array | no | Soil horizon class input |
| [lai_class_path](#lai_class_path) | string array | no | LAI class input |
| [river_width_path](#river_width_path) | string array | no | River width input |
| [morph_mask_path](#morph_mask_path) | string array | no | Morphology mask file path |
| [pre_var](#pre_var) | string array | no | Precipitation variable name |
| [pet_var](#pet_var) | string array | no | Potential evapotranspiration variable name |
| [temp_var](#temp_var) | string array | no | Air temperature variable name |
| [tann_var](#tann_var) | string array | no | Air temperature annual mean variable name |
| [tmin_var](#tmin_var) | string array | no | Air temperature daily minimum variable name |
| [tmax_var](#tmax_var) | string array | no | Air temperature daily maximum variable name |
| [ssrd_var](#ssrd_var) | string array | no | Surface shortwave radiation variable name |
| [strd_var](#strd_var) | string array | no | Surface thermal radiation variable name |
| [netrad_var](#netrad_var) | string array | no | Net radiation variable name |
| [eabs_var](#eabs_var) | string array | no | Vapor pressure variable name |
| [wind_var](#wind_var) | string array | no | Wind speed variable name |
| [meteo_mask_var](#meteo_mask_var) | string array | no | Meteorological mask variable name |
| [runoff_var](#runoff_var) | string array | no | Runoff variable name |
| [runoff_sealed_var](#runoff_sealed_var) | string array | no | Sealed runoff variable name |
| [interflow_fast_var](#interflow_fast_var) | string array | no | Fast interflow variable name |
| [interflow_slow_var](#interflow_slow_var) | string array | no | Slow interflow variable name |
| [baseflow_var](#baseflow_var) | string array | no | Baseflow variable name |
| [hydro_mask_var](#hydro_mask_var) | string array | no | Hydrological mask variable name |
| [dem_var](#dem_var) | string array | no | DEM variable name |
| [slope_var](#slope_var) | string array | no | Slope variable name |
| [aspect_var](#aspect_var) | string array | no | Aspect variable name |
| [fdir_var](#fdir_var) | string array | no | Flow direction variable name |
| [facc_var](#facc_var) | string array | no | Flow accumulation variable name |
| [geo_class_var](#geo_class_var) | string array | no | Geology class variable name |
| [soil_class_var](#soil_class_var) | string array | no | Soil class variable name |
| [lai_class_var](#lai_class_var) | string array | no | LAI class variable name |
| [river_width_var](#river_width_var) | string array | no | River width variable name |
| [morph_mask_var](#morph_mask_var) | string array | no | Morphology mask variable name |
| [hydro_lat_var](#hydro_lat_var) | string array | no | Hydrological latitude variable name |
| [hydro_lon_var](#hydro_lon_var) | string array | no | Hydrological longitude variable name |
| [morph_lat_var](#morph_lat_var) | string array | no | Morphology latitude variable name |
| [morph_lon_var](#morph_lon_var) | string array | no | Morphology longitude variable name |
| [route_lat_var](#route_lat_var) | string array | no | Routing latitude variable name |
| [route_lon_var](#route_lon_var) | string array | no | Routing longitude variable name |

## Field details

### chunking

Chunking for input data `chunking`

Chunking configuration for reading input data (array dimension 1: domain)
- n>0 : after each n time-steps
- 0 : only at beginning of the run
- -1 : daily
- -2 : monthly
- -3 : yearly

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Default: `0`
- Minimum: `>= -3`

### time_stamp_location

NetCDF time-stamp location `time_stamp_location`

NetCDF time-stamp location: when no time-bounds are given in the input data.
- 0 : beginning of time step (default)
- 1 : center of time step
- 2 : end of time step

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Default: `0`
- Allowed values: `0`, `1`, `2`

### latlon_path

Latlon specification file path `latlon_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### morph_latlon

DEM in latlon coordinates `morph_latlon`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### pre_path

Precipitation input `pre_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### pet_path

Potential evapotranspiration input `pet_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### temp_path

Air temperature input `temp_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### tann_path

Air temperature annual mean input `tann_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### tmin_path

Air temperature daily minimum input `tmin_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### tmax_path

Air temperature daily maximum input `tmax_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### ssrd_path

Surface shortwave radiation downwards input `ssrd_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### strd_path

Surface thermal radiation downwards input `strd_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### netrad_path

Net radiation input `netrad_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### eabs_path

Vapor pressure input `eabs_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### wind_path

Wind speed input `wind_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### meteo_mask_path

Meteorological mask file path `meteo_mask_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### runoff_path

Runoff input `runoff_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### runoff_sealed_path

Sealed runoff input `runoff_sealed_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### interflow_fast_path

Fast interflow input `interflow_fast_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### interflow_slow_path

Slow interflow input `interflow_slow_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### baseflow_path

Baseflow input `baseflow_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### hydro_mask_path

Hydrological mask file path `hydro_mask_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### dem_path

DEM input `dem_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### slope_path

Slope input `slope_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### aspect_path

Aspect input `aspect_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### fdir_path

Flow direction input `fdir_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### facc_path

Flow accumulation input `facc_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### geo_class_path

Geology class input `geo_class_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### soil_class_path

Soil class input `soil_class_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### soil_horizon_class_path

Soil horizon class input `soil_horizon_class_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### lai_class_path

LAI class input `lai_class_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### river_width_path

River width input `river_width_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### morph_mask_path

Morphology mask file path `morph_mask_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### pre_var

Precipitation variable name `pre_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"pre"`

### pet_var

Potential evapotranspiration variable name `pet_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"pet"`

### temp_var

Air temperature variable name `temp_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"tavg"`

### tann_var

Air temperature annual mean variable name `tann_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"tann"`

### tmin_var

Air temperature daily minimum variable name `tmin_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"tmin"`

### tmax_var

Air temperature daily maximum variable name `tmax_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"tmax"`

### ssrd_var

Surface shortwave radiation variable name `ssrd_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"ssrd"`

### strd_var

Surface thermal radiation variable name `strd_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"strd"`

### netrad_var

Net radiation variable name `netrad_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"net_rad"`

### eabs_var

Vapor pressure variable name `eabs_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"eabs"`

### wind_var

Wind speed variable name `wind_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"windspeed"`

### meteo_mask_var

Meteorological mask variable name `meteo_mask_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"mask"`

### runoff_var

Runoff variable name `runoff_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"runoff"`

### runoff_sealed_var

Sealed runoff variable name `runoff_sealed_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"runoff_sealed"`

### interflow_fast_var

Fast interflow variable name `interflow_fast_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"interflow_fast"`

### interflow_slow_var

Slow interflow variable name `interflow_slow_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"interflow_slow"`

### baseflow_var

Baseflow variable name `baseflow_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"baseflow"`

### hydro_mask_var

Hydrological mask variable name `hydro_mask_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"mask"`

### dem_var

DEM variable name `dem_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"dem"`

### slope_var

Slope variable name `slope_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"slope"`

### aspect_var

Aspect variable name `aspect_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"aspect"`

### fdir_var

Flow direction variable name `fdir_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"fdir"`

### facc_var

Flow accumulation variable name `facc_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"facc"`

### geo_class_var

Geology class variable name `geo_class_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"geology_class"`

### soil_class_var

Soil class variable name `soil_class_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"soil_class"`

### lai_class_var

LAI class variable name `lai_class_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"LAI_class"`

### river_width_var

River width variable name `river_width_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"P_bkfl"`

### morph_mask_var

Morphology mask variable name `morph_mask_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"mask"`

### hydro_lat_var

Hydrological latitude variable name `hydro_lat_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"lat"`

### hydro_lon_var

Hydrological longitude variable name `hydro_lon_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"lon"`

### morph_lat_var

Morphology latitude variable name `morph_lat_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"lat_l0"`

### morph_lon_var

Morphology longitude variable name `morph_lon_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"lon_l0"`

### route_lat_var

Routing latitude variable name `route_lat_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"lat_l11"`

### route_lon_var

Routing longitude variable name `route_lon_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"lon_l11"`

## Example

```fortran
&config_input
  chunking(:) = 0
  time_stamp_location(:) = 0
  latlon_path(:) = ""
  morph_latlon(:) = .false.
  pre_path(:) = ""
  pet_path(:) = ""
  temp_path(:) = ""
  tann_path(:) = ""
  tmin_path(:) = ""
  tmax_path(:) = ""
  ssrd_path(:) = ""
  strd_path(:) = ""
  netrad_path(:) = ""
  eabs_path(:) = ""
  wind_path(:) = ""
  meteo_mask_path(:) = ""
  runoff_path(:) = ""
  runoff_sealed_path(:) = ""
  interflow_fast_path(:) = ""
  interflow_slow_path(:) = ""
  baseflow_path(:) = ""
  hydro_mask_path(:) = ""
  dem_path(:) = ""
  slope_path(:) = ""
  aspect_path(:) = ""
  fdir_path(:) = ""
  facc_path(:) = ""
  geo_class_path(:) = ""
  soil_class_path(:) = ""
  soil_horizon_class_path(:) = ""
  lai_class_path(:) = ""
  river_width_path(:) = ""
  morph_mask_path(:) = ""
  pre_var(:) = "pre"
  pet_var(:) = "pet"
  temp_var(:) = "tavg"
  tann_var(:) = "tann"
  tmin_var(:) = "tmin"
  tmax_var(:) = "tmax"
  ssrd_var(:) = "ssrd"
  strd_var(:) = "strd"
  netrad_var(:) = "net_rad"
  eabs_var(:) = "eabs"
  wind_var(:) = "windspeed"
  meteo_mask_var(:) = "mask"
  runoff_var(:) = "runoff"
  runoff_sealed_var(:) = "runoff_sealed"
  interflow_fast_var(:) = "interflow_fast"
  interflow_slow_var(:) = "interflow_slow"
  baseflow_var(:) = "baseflow"
  hydro_mask_var(:) = "mask"
  dem_var(:) = "dem"
  slope_var(:) = "slope"
  aspect_var(:) = "aspect"
  fdir_var(:) = "fdir"
  facc_var(:) = "facc"
  geo_class_var(:) = "geology_class"
  soil_class_var(:) = "soil_class"
  lai_class_var(:) = "LAI_class"
  river_width_var(:) = "P_bkfl"
  morph_mask_var(:) = "mask"
  hydro_lat_var(:) = "lat"
  hydro_lon_var(:) = "lon"
  morph_lat_var(:) = "lat_l0"
  morph_lon_var(:) = "lon_l0"
  route_lat_var(:) = "lat_l11"
  route_lon_var(:) = "lon_l11"
/
```

