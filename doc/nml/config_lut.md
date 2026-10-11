# Look-up table configuration {#config_lut}

[TOC]

Domain-indexed look-up table (LUT) paths provided as metadata independent of the
process containers. Attached to the exchange of each domain so that any container
(land-surface MPR, Input, mLM) can resolve the shared class-to-property tables
without duplicating the path in its own process configuration.

**Namelist**: `config_lut`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [soil_lut_path](#soil_lut_path) | string array | no | no | Soil look-up table path |
| [geo_lut_path](#geo_lut_path) | string array | no | no | Geology look-up table path |

## Field details

### soil_lut_path

Soil look-up table path `soil_lut_path`

Soil look-up table (class -> texture) path. Resolved on the exchange and read by the
land-surface MPR (soil texture per soil class) and by the Input container (top-horizon
sand and clay of the lake footprint soil classes for the detached lake multiscale
parameter regionalization of the lakebed conductivity).

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no

### geo_lut_path

Geology look-up table path `geo_lut_path`

Geology look-up table (class -> GeoParam index + geological unit) path. Resolved on the
exchange and read by the land-surface MPR and by the mLM container (per-geological-unit
baseflow recession values over the lake footprint for the detached lake multiscale
parameter regionalization of the lake baseflow recession coefficient).

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no

## Example

```fortran
&config_lut
  soil_lut_path(:) = ""
  geo_lut_path(:) = ""
/
```

