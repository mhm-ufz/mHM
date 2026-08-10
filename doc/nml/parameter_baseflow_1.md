# Baseflow - Case 1 {#baseflow_1}

[TOC]

Geological baseflow recession parameters.

**Namelist**: `baseflow_1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [baseflow_recession](#baseflow_recession) | type(parameter_t) array | yes | yes | Baseflow recession time for each geological unit [d] |

## Field details

### baseflow_recession

Baseflow recession time for each geological unit [d] `baseflow_recession`

Values are addressed by the first column of each geology class-definition file.

Summary:
- Type: `type(parameter_t), dimension(n_geo_units)`
- Declared required: yes
- Input required: yes
- Default: `{min: 1.0, max: 1000.0}` (broadcast item default)

Components:
- `baseflow_recession%value`: `real(dp)`; declared required yes; input required yes
- `baseflow_recession%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `baseflow_recession%min`: `real(dp)`; declared required yes; input required no; default `1.0` (item default)
- `baseflow_recession%max`: `real(dp)`; declared required yes; input required no; default `1000.0` (item default)

## Derived types

### `parameter_t`

Calibration parameter

A model parameter with optional calibration metadata.

- Ownership: imported from `mo_parameter_types`
- Buffer-compatible: yes
- Component order: value, optimize, min, max
- **Declaration-order contract:** the imported Fortran type must declare components in the resolved schema order shown above.
- `value`: `real(dp)`; declared required yes; input required yes
- `optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `min`: `real(dp)`; declared required yes; input required yes
- `max`: `real(dp)`; declared required yes; input required yes

## Example

```fortran
&baseflow_1
  baseflow_recession(:)%value = 100.0
  baseflow_recession(:)%optimize = .true.
  baseflow_recession(:)%min = 1.0
  baseflow_recession(:)%max = 1000.0
/
```

