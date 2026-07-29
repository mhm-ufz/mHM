# Routing - Case 3 {#routing_3}

[TOC]

Parameters for varying-celerity routing.

**Namelist**: `routing_3`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [slope_factor](#slope_factor) | type(parameter_t) | yes | yes | Slope factor |

## Field details

### slope_factor

Slope factor `slope_factor`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.1, upper_bound: 100.0}`

Components:
- `slope_factor%value`: `real(dp)`; declared required yes; input required yes
- `slope_factor%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `slope_factor%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.1` (object default)
- `slope_factor%upper_bound`: `real(dp)`; declared required yes; input required no; default `100.0` (object default)

## Derived types

### `parameter_t`

Calibration parameter

A model parameter with optional calibration metadata.

- Ownership: `nml_helper`
- Buffer-compatible: yes
- Component order: value, optimize, lower_bound, upper_bound
- `value`: `real(dp)`; declared required yes; input required yes
- `optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `lower_bound`: `real(dp)`; declared required yes; input required yes
- `upper_bound`: `real(dp)`; declared required yes; input required yes

## Example

```fortran
&routing_3
  slope_factor%value = 30.0
  slope_factor%optimize = .false.
  slope_factor%lower_bound = 0.1
  slope_factor%upper_bound = 100.0
/
```

