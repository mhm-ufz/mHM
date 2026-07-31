# Routing - Case 2 {#routing_2}

[TOC]

Parameters for constant-celerity routing.

**Namelist**: `routing_2`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [streamflow_celerity](#streamflow_celerity) | type(parameter_t) | yes | yes | Streamflow celerity |

## Field details

### streamflow_celerity

Streamflow celerity `streamflow_celerity`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.1, max: 15.0}`

Components:
- `streamflow_celerity%value`: `real(dp)`; declared required yes; input required yes
- `streamflow_celerity%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `streamflow_celerity%min`: `real(dp)`; declared required yes; input required no; default `0.1` (object default)
- `streamflow_celerity%max`: `real(dp)`; declared required yes; input required no; default `15.0` (object default)

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
&routing_2
  streamflow_celerity%value = 1.5
  streamflow_celerity%optimize = .false.
  streamflow_celerity%min = 0.1
  streamflow_celerity%max = 15.0
/
```

