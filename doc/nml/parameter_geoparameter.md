# Geological parameters {#geoparameter}

[TOC]

Parameters for geoparameter.

**Namelist**: `geoparameter`

## Fields

| Name | Type | Required | Info |
| --- | --- | --- | --- |
| [GeoParam](#geoparam) | real array | yes | Geological parameters |

## Field details

### GeoParam

Geological parameters `GeoParam`

Geological parameters (ordering according to class-definition file)
These parameters are NOT REGIONALIZED yet, i.e. these are <beta> and not <gamma>

Summary:
- Type: `real(dp), dimension(5, max_geo_units)`
- Flexible tail dims: 1
- Required: yes
- Examples: `[1.0, 1000.0, 100.0, 1.0, 1.0]`

## Example

```fortran
&geoparameter
  GeoParam(:,1) = 1.0, 1000.0, 100.0, 1.0, 1.0
/
```

