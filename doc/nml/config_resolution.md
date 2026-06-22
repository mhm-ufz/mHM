# Grid resolution configuration {#config_resolution}

[TOC]

Domain-indexed target resolutions for internally derived grids.

**Namelist**: `config_resolution`

## Fields

| Name | Type | Required | Info |
| --- | --- | --- | --- |
| [hydro](#hydro) | real array | no | Hydrological grid resolution |
| [route](#route) | real array | no | Routing grid resolution |

## Field details

### hydro

Hydrological grid resolution `hydro`

Target resolution for deriving the hydrological grid from the morphology grid.

Summary:
- Type: `real(dp), dimension(n_domains)`
- Required: no
- Minimum: `> 0.0`

### route

Routing grid resolution `route`

Target resolution for deriving the routing grid from the morphology grid.

Summary:
- Type: `real(dp), dimension(n_domains)`
- Required: no
- Minimum: `> 0.0`

## Example

```fortran
&config_resolution
  hydro(:) = 0.0
  route(:) = 0.0
/
```

