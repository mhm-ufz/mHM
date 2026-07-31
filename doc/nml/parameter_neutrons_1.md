# Neutrons - Case 1 {#neutrons_1}

[TOC]

Parameters for the experimental Desilets neutron formulation.

**Namelist**: `neutrons_1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [desilets_n0](#desilets_n0) | type(parameter_t) | yes | yes | Desilets dry neutron count N0 |
| [desilets_lw0](#desilets_lw0) | type(parameter_t) | yes | yes | Desilets lattice-water parameter LW0 |
| [desilets_lw1](#desilets_lw1) | type(parameter_t) | yes | yes | Desilets lattice-water parameter LW1 |

## Field details

### desilets_n0

Desilets dry neutron count N0 `desilets_n0`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 300.0, max: 2000.0}`

Components:
- `desilets_n0%value`: `real(dp)`; declared required yes; input required yes
- `desilets_n0%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `desilets_n0%min`: `real(dp)`; declared required yes; input required no; default `300.0` (object default)
- `desilets_n0%max`: `real(dp)`; declared required yes; input required no; default `2000.0` (object default)

### desilets_lw0

Desilets lattice-water parameter LW0 `desilets_lw0`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 0.2}`

Components:
- `desilets_lw0%value`: `real(dp)`; declared required yes; input required yes
- `desilets_lw0%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `desilets_lw0%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `desilets_lw0%max`: `real(dp)`; declared required yes; input required no; default `0.2` (object default)

### desilets_lw1

Desilets lattice-water parameter LW1 `desilets_lw1`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 0.05}`

Components:
- `desilets_lw1%value`: `real(dp)`; declared required yes; input required yes
- `desilets_lw1%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `desilets_lw1%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `desilets_lw1%max`: `real(dp)`; declared required yes; input required no; default `0.05` (object default)

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
&neutrons_1
  desilets_n0%value = 1500.0
  desilets_n0%optimize = .false.
  desilets_n0%min = 300.0
  desilets_n0%max = 2000.0
  desilets_lw0%value = 0.1783
  desilets_lw0%optimize = .false.
  desilets_lw0%min = 0.0
  desilets_lw0%max = 0.2
  desilets_lw1%value = 0.0
  desilets_lw1%optimize = .false.
  desilets_lw1%min = 0.0
  desilets_lw1%max = 0.05
/
```

