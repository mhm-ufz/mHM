# Observations configuration {#config_observations}

[TOC]

Configuration for observation input data in mHM.

**Namelist**: `config_observations`

## Fields

| Name | Type | Required | Info |
| --- | --- | --- | --- |
| [sm_path](#sm_path) | string array | no | Soil moisture data path |
| [neutrons_path](#neutrons_path) | string array | no | Neutron data path |
| [et_path](#et_path) | string array | no | Evapotranspiration data path |
| [tws_path](#tws_path) | string array | no | Domain average TWS path |
| [BFI_obs](#bfi_obs) | real array | no | Baseflow index per domain |
| [sm_horizons](#sm_horizons) | integer | no | Number of mHM soil moisture horizons |
| [sm_time_step](#sm_time_step) | integer | no | Time step of soil moisture |
| [et_time_step](#et_time_step) | integer | no | Time step of evapotranspiration |
| [tws_time_step](#tws_time_step) | integer | no | Time step of total water storage |

## Field details

### sm_path

Soil moisture data path `sm_path`

Path to soil moisture input data.
Expected variable name: "sm".

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### neutrons_path

Neutron data path `neutrons_path`

Path to neutron input data.
Expected variable name: "neutrons".

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### et_path

Evapotranspiration data path `et_path`

Path to evapotranspiration input data.
Expected variable name: "et".

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### tws_path

Domain average TWS path `tws_path`

Path to total water storage input data.
Expected variable name: "twsa".

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no

### BFI_obs

Baseflow index per domain `BFI_obs`

Baseflow index per domain. Only needed if not calculated (BFI_calc = .false.).
You can overwrite single BFI values to not calculate them internally (if BFI_calc = .true.).

Summary:
- Type: `real(dp), dimension(n_domains)`
- Required: no

### sm_horizons

Number of mHM soil moisture horizons `sm_horizons`

Number of mHM soil moisture horizons for which the soil moisture input is representative for (counted top to down).

Summary:
- Type: `integer(i4)`
- Required: no

### sm_time_step

Time step of soil moisture `sm_time_step`

Time step of soil moisture input data.
- -1 : daily
- -2 : monthly
- -3 : yearly

Summary:
- Type: `integer(i4)`
- Required: no

### et_time_step

Time step of evapotranspiration `et_time_step`

Time step of evapotranspiration input data.
- -1 : daily
- -2 : monthly
- -3 : yearly

Summary:
- Type: `integer(i4)`
- Required: no

### tws_time_step

Time step of total water storage `tws_time_step`

Time step of total water storage input data.
- -1 : daily
- -2 : monthly
- -3 : yearly

Summary:
- Type: `integer(i4)`
- Required: no

## Example

```fortran
&config_observations
  sm_path(:) = ""
  neutrons_path(:) = ""
  et_path(:) = ""
  tws_path(:) = ""
  BFI_obs(:) = 0.0
  sm_horizons = 0
  sm_time_step = 0
  et_time_step = 0
  tws_time_step = 0
/
```

