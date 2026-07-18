# Interception - Case 1 {#interception_1}

[TOC]

Parameters for interception.

**Namelist**: `interception_1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [canopy_interception_factor](#canopy_interception_factor) | type(parameter_t) | yes | yes | Multiplier relating LAI to interception storage |

## Field details

### canopy_interception_factor

Multiplier relating LAI to interception storage `canopy_interception_factor`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.15, upper_bound: 0.4}`

Components:
- `canopy_interception_factor%value`: `real(dp)`; declared required yes; input required yes
- `canopy_interception_factor%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `canopy_interception_factor%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.15` (object default)
- `canopy_interception_factor%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.4` (object default)

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
&interception_1
  canopy_interception_factor%value = 0.15
  canopy_interception_factor%optimize = .true.
  canopy_interception_factor%lower_bound = 0.15
  canopy_interception_factor%upper_bound = 0.4
/
```

