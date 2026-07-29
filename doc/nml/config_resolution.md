# Grid resolution configuration {#config_resolution}

[TOC]

Domain-indexed target resolutions for internally derived grids.

**Namelist**: `config_resolution`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [hydro](#hydro) | real array | no | no | Hydrological grid resolution |
| [route](#route) | real array | no | no | Routing grid resolution |

## Field details

### hydro

Hydrological grid resolution `hydro`

Target resolution for deriving the hydrological grid from the morphological grid.
A value of 0 requires an existing grid and disables internal grid derivation.
Unit is degrees for geographic coordinate systems and meters otherwise, as determined from the morphological grid.

Summary:
- Type: `real(dp), dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `0.0`
- Minimum: `>= 0.0`
- Examples: `[0.1]`

### route

Routing grid resolution `route`

Target resolution for deriving the routing grid from the morphological grid.
A value of 0 requires an existing grid and disables internal grid derivation.
Unit is degrees for geographic coordinate systems and meters otherwise, as determined from the morphological grid.

Summary:
- Type: `real(dp), dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `0.0`
- Minimum: `>= 0.0`
- Examples: `[0.5]`

## Example

```fortran
&config_resolution
  hydro(:) = 0.1
  route(:) = 0.5
/
```

