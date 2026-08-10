# Percolation - Case 1 {#percolation_1}

[TOC]

Parameters for percolation and karst recharge.

**Namelist**: `percolation_1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [recharge_coefficient](#recharge_coefficient) | type(parameter_t) | yes | yes | Percolation and recharge time [d] |
| [karstic_recharge_factor](#karstic_recharge_factor) | type(parameter_t) | yes | yes | Karstic recharge factor [1] |

## Field details

### recharge_coefficient

Percolation and recharge time [d] `recharge_coefficient`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 50.0}`

Components:
- `recharge_coefficient%value`: `real(dp)`; declared required yes; input required yes
- `recharge_coefficient%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `recharge_coefficient%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `recharge_coefficient%max`: `real(dp)`; declared required yes; input required no; default `50.0` (object default)

### karstic_recharge_factor

Karstic recharge factor [1] `karstic_recharge_factor`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: -5.0, max: 5.0}`

Components:
- `karstic_recharge_factor%value`: `real(dp)`; declared required yes; input required yes
- `karstic_recharge_factor%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `karstic_recharge_factor%min`: `real(dp)`; declared required yes; input required no; default `-5.0` (object default)
- `karstic_recharge_factor%max`: `real(dp)`; declared required yes; input required no; default `5.0` (object default)

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
&percolation_1
  recharge_coefficient%value = 35.0
  recharge_coefficient%optimize = .true.
  recharge_coefficient%min = 0.0
  recharge_coefficient%max = 50.0
  karstic_recharge_factor%value = -1.0
  karstic_recharge_factor%optimize = .true.
  karstic_recharge_factor%min = -5.0
  karstic_recharge_factor%max = 5.0
/
```

