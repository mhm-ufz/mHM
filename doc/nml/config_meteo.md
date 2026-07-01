# Meteorological configuration {#config_meteo}

[TOC]

Configuration for meteorological input data handling in mHM.
Meteorological weights can be used to disaggregate daily data to hourly values.

**Namelist**: `config_meteo`

## Fields

| Name | Type | Required | Info |
| --- | --- | --- | --- |
| [read_meteo_weights](#read_meteo_weights) | logical array | no | Read meteorological weights |
| [pre_weights_path](#pre_weights_path) | string array | no | Precipitation weights path |
| [pet_weights_path](#pet_weights_path) | string array | no | Potential evapotranspiration weights path |
| [temp_weights_path](#temp_weights_path) | string array | no | Surface downward shortwave radiation weights path |
| [ssrd_weights_path](#ssrd_weights_path) | string array | no | Surface downward shortwave radiation weights path |
| [strd_weights_path](#strd_weights_path) | string array | no | Surface downward longwave radiation weights path |
| [pre_weights_var](#pre_weights_var) | string array | no | Precipitation weights variable name |
| [pet_weights_var](#pet_weights_var) | string array | no | Potential evapotranspiration weights variable name |
| [temp_weights_var](#temp_weights_var) | string array | no | Average temperature weights variable name |
| [ssrd_weights_var](#ssrd_weights_var) | string array | no | Surface downward shortwave radiation weights variable name |
| [strd_weights_var](#strd_weights_var) | string array | no | Surface downward longwave radiation weights variable name |
| [frac_night_pre](#frac_night_pre) | real array | no | Fraction of nightly precipitation |
| [frac_night_pet](#frac_night_pet) | real array | no | Fraction of nightly potential evapotranspiration |
| [frac_night_temp](#frac_night_temp) | real array | no | Fraction of nightly temperature |
| [frac_night_ssrd](#frac_night_ssrd) | real array | no | Fraction of nightly surface downward shortwave radiation |
| [frac_night_strd](#frac_night_strd) | real array | no | Fraction of nightly surface downward longwave radiation |
| [share_frac](#share_frac) | logical | no | Share fractions between domains |

## Field details

### read_meteo_weights

Read meteorological weights `read_meteo_weights`

Read meteorological weights from file instead of using fixed night ratios.
The dimension for the weights are in FORTRAN-notation (rows, colums, months=12, hours=24)
and in C-notation like in NetCDF: (hours=24, months=12, colums, rows).

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### pre_weights_path

Precipitation weights path `pre_weights_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Required: no

### pet_weights_path

Potential evapotranspiration weights path `pet_weights_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Required: no

### temp_weights_path

Surface downward shortwave radiation weights path `temp_weights_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Required: no

### ssrd_weights_path

Surface downward shortwave radiation weights path `ssrd_weights_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Required: no

### strd_weights_path

Surface downward longwave radiation weights path `strd_weights_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Required: no

### pre_weights_var

Precipitation weights variable name `pre_weights_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"pre_weights"`

### pet_weights_var

Potential evapotranspiration weights variable name `pet_weights_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"pet_weights"`

### temp_weights_var

Average temperature weights variable name `temp_weights_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"tavg_weights"`

### ssrd_weights_var

Surface downward shortwave radiation weights variable name `ssrd_weights_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"ssrd_weights"`

### strd_weights_var

Surface downward longwave radiation weights variable name `strd_weights_var`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Default: `"strd_weights"`

### frac_night_pre

Fraction of nightly precipitation `frac_night_pre`

Night ratio for precipitation (only if read_meteo_weights = .false.).
Only night values (18h-5h) required because day values are the complement to 1.
(array dimension 1: hour, dimension 2: domain)

Summary:
- Type: `real(dp), dimension(12, n_domains)`
- Required: no
- Examples: `[0.46, 0.5, 0.52, 0.51, 0.48, 0.5, 0.49, 0.48, 0.52, 0.56, 0.5, 0.47]`

### frac_night_pet

Fraction of nightly potential evapotranspiration `frac_night_pet`

Night ratio for potential evapotranspiration (only if read_meteo_weights = .false.).
Only night values (18h-5h) required because day values are the complement to 1.
(array dimension 1: hour, dimension 2: domain)

Summary:
- Type: `real(dp), dimension(12, n_domains)`
- Required: no
- Examples: `[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]`

### frac_night_temp

Fraction of nightly temperature `frac_night_temp`

Additive night deviation for temperature (only if read_meteo_weights = .false.).
Only night values (18h-5h) required because day values are the complement to 1.
(array dimension 1: hour, dimension 2: domain)

Summary:
- Type: `real(dp), dimension(12, n_domains)`
- Required: no
- Examples: `[-0.76, -1.3, -1.88, -2.38, -2.72, -2.75, -2.74, -3.04, -2.44, -1.6, -0.94, -0.53]`

### frac_night_ssrd

Fraction of nightly surface downward shortwave radiation `frac_night_ssrd`

Night ratio for surface downward shortwave radiation (only if read_meteo_weights = .false.).
Only night values (18h-5h) required because day values are the complement to 1.
(array dimension 1: hour, dimension 2: domain)

Summary:
- Type: `real(dp), dimension(12, n_domains)`
- Required: no
- Examples: `[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]`

### frac_night_strd

Fraction of nightly surface downward longwave radiation `frac_night_strd`

Night ratio for surface downward longwave radiation (only if read_meteo_weights = .false.).
Only night values (18h-5h) required because day values are the complement to 1.
(array dimension 1: hour, dimension 2: domain)

Summary:
- Type: `real(dp), dimension(12, n_domains)`
- Required: no
- Examples: `[0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45]`

### share_frac

Share fractions between domains `share_frac`

Share the same fraction values between all domains taking the values from the first domain.
This only applies if `read_meteo_weights` is `.false.`.
If set to false, different fraction values can be set for each domain.

Summary:
- Type: `logical`
- Required: no
- Default: `.true.`

## Example

```fortran
&config_meteo
  read_meteo_weights(:) = .false.
  pre_weights_path(:) = ""
  pet_weights_path(:) = ""
  temp_weights_path(:) = ""
  ssrd_weights_path(:) = ""
  strd_weights_path(:) = ""
  pre_weights_var(:) = "pre_weights"
  pet_weights_var(:) = "pet_weights"
  temp_weights_var(:) = "tavg_weights"
  ssrd_weights_var(:) = "ssrd_weights"
  strd_weights_var(:) = "strd_weights"
  frac_night_pre(:,1) = 0.46, 0.5, 0.52, 0.51, 0.48, 0.5, 0.49, 0.48, 0.52, 0.56, 0.5, 0.47
  frac_night_pet(:,1) = 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1
  frac_night_temp(:,1) = -0.76, -1.3, -1.88, -2.38, -2.72, -2.75, -2.74, -3.04, -2.44, -1.6, -0.94, -0.53
  frac_night_ssrd(:,1) = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  frac_night_strd(:,1) = 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45
  share_frac = .true.
/
```

