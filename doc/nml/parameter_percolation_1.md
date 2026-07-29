# Percolation - Case 1 {#percolation_1}

[TOC]

Parameters for percolation and karst recharge.

**Namelist**: `percolation_1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [recharge_coefficient](#recharge_coefficient) | type(parameter_t) | yes | yes | Recharge coefficient |
| [karstic_recharge_factor](#karstic_recharge_factor) | type(parameter_t) | yes | yes | Karstic recharge factor |

## Field details

### recharge_coefficient

Recharge coefficient `recharge_coefficient`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.0, upper_bound: 50.0}`

Components:
- `recharge_coefficient%value`: `real(dp)`; declared required yes; input required yes
- `recharge_coefficient%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `recharge_coefficient%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `recharge_coefficient%upper_bound`: `real(dp)`; declared required yes; input required no; default `50.0` (object default)

### karstic_recharge_factor

Karstic recharge factor `karstic_recharge_factor`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: -5.0, upper_bound: 5.0}`

Components:
- `karstic_recharge_factor%value`: `real(dp)`; declared required yes; input required yes
- `karstic_recharge_factor%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `karstic_recharge_factor%lower_bound`: `real(dp)`; declared required yes; input required no; default `-5.0` (object default)
- `karstic_recharge_factor%upper_bound`: `real(dp)`; declared required yes; input required no; default `5.0` (object default)

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
&percolation_1
  recharge_coefficient%value = 35.0
  recharge_coefficient%optimize = .true.
  recharge_coefficient%lower_bound = 0.0
  recharge_coefficient%upper_bound = 50.0
  karstic_recharge_factor%value = -1.0
  karstic_recharge_factor%optimize = .true.
  karstic_recharge_factor%lower_bound = -5.0
  karstic_recharge_factor%upper_bound = 5.0
/
```

