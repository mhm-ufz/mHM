# Neutrons - Case 2 {#neutrons_2}

[TOC]

Parameters for the experimental COSMIC neutron formulation.

**Namelist**: `neutrons_2`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [cosmic_n0](#cosmic_n0) | type(parameter_t) | yes | yes | COSMIC N0 parameter |
| [cosmic_n1](#cosmic_n1) | type(parameter_t) | yes | yes | COSMIC N1 parameter |
| [cosmic_n2](#cosmic_n2) | type(parameter_t) | yes | yes | COSMIC N2 parameter |
| [cosmic_alpha0](#cosmic_alpha0) | type(parameter_t) | yes | yes | COSMIC alpha0 parameter |
| [cosmic_alpha1](#cosmic_alpha1) | type(parameter_t) | yes | yes | COSMIC alpha1 parameter |
| [cosmic_l30](#cosmic_l30) | type(parameter_t) | yes | yes | COSMIC L30 parameter |
| [cosmic_l31](#cosmic_l31) | type(parameter_t) | yes | yes | COSMIC L31 parameter |
| [cosmic_lw0](#cosmic_lw0) | type(parameter_t) | yes | yes | COSMIC lattice-water parameter LW0 |
| [cosmic_lw1](#cosmic_lw1) | type(parameter_t) | yes | yes | COSMIC lattice-water parameter LW1 |

## Field details

### cosmic_n0

COSMIC N0 parameter `cosmic_n0`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 300.0, max: 2000.0}`

Components:
- `cosmic_n0%value`: `real(dp)`; declared required yes; input required yes
- `cosmic_n0%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `cosmic_n0%min`: `real(dp)`; declared required yes; input required no; default `300.0` (object default)
- `cosmic_n0%max`: `real(dp)`; declared required yes; input required no; default `2000.0` (object default)

### cosmic_n1

COSMIC N1 parameter `cosmic_n1`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.01, max: 10.0}`

Components:
- `cosmic_n1%value`: `real(dp)`; declared required yes; input required yes
- `cosmic_n1%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `cosmic_n1%min`: `real(dp)`; declared required yes; input required no; default `0.01` (object default)
- `cosmic_n1%max`: `real(dp)`; declared required yes; input required no; default `10.0` (object default)

### cosmic_n2

COSMIC N2 parameter `cosmic_n2`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.01, max: 10.0}`

Components:
- `cosmic_n2%value`: `real(dp)`; declared required yes; input required yes
- `cosmic_n2%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `cosmic_n2%min`: `real(dp)`; declared required yes; input required no; default `0.01` (object default)
- `cosmic_n2%max`: `real(dp)`; declared required yes; input required no; default `10.0` (object default)

### cosmic_alpha0

COSMIC alpha0 parameter `cosmic_alpha0`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.01, max: 10.0}`

Components:
- `cosmic_alpha0%value`: `real(dp)`; declared required yes; input required yes
- `cosmic_alpha0%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `cosmic_alpha0%min`: `real(dp)`; declared required yes; input required no; default `0.01` (object default)
- `cosmic_alpha0%max`: `real(dp)`; declared required yes; input required no; default `10.0` (object default)

### cosmic_alpha1

COSMIC alpha1 parameter `cosmic_alpha1`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.01, max: 10.0}`

Components:
- `cosmic_alpha1%value`: `real(dp)`; declared required yes; input required yes
- `cosmic_alpha1%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `cosmic_alpha1%min`: `real(dp)`; declared required yes; input required no; default `0.01` (object default)
- `cosmic_alpha1%max`: `real(dp)`; declared required yes; input required no; default `10.0` (object default)

### cosmic_l30

COSMIC L30 parameter `cosmic_l30`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 26.56, max: 424.78}`

Components:
- `cosmic_l30%value`: `real(dp)`; declared required yes; input required yes
- `cosmic_l30%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `cosmic_l30%min`: `real(dp)`; declared required yes; input required no; default `26.56` (object default)
- `cosmic_l30%max`: `real(dp)`; declared required yes; input required no; default `424.78` (object default)

### cosmic_l31

COSMIC L31 parameter `cosmic_l31`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: -118.3, max: 200.28}`

Components:
- `cosmic_l31%value`: `real(dp)`; declared required yes; input required yes
- `cosmic_l31%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `cosmic_l31%min`: `real(dp)`; declared required yes; input required no; default `-118.3` (object default)
- `cosmic_l31%max`: `real(dp)`; declared required yes; input required no; default `200.28` (object default)

### cosmic_lw0

COSMIC lattice-water parameter LW0 `cosmic_lw0`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 0.2}`

Components:
- `cosmic_lw0%value`: `real(dp)`; declared required yes; input required yes
- `cosmic_lw0%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `cosmic_lw0%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `cosmic_lw0%max`: `real(dp)`; declared required yes; input required no; default `0.2` (object default)

### cosmic_lw1

COSMIC lattice-water parameter LW1 `cosmic_lw1`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 0.05}`

Components:
- `cosmic_lw1%value`: `real(dp)`; declared required yes; input required yes
- `cosmic_lw1%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `cosmic_lw1%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `cosmic_lw1%max`: `real(dp)`; declared required yes; input required no; default `0.05` (object default)

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
&neutrons_2
  cosmic_n0%value = 1500.0
  cosmic_n0%optimize = .false.
  cosmic_n0%min = 300.0
  cosmic_n0%max = 2000.0
  cosmic_n1%value = 1.0
  cosmic_n1%optimize = .false.
  cosmic_n1%min = 0.01
  cosmic_n1%max = 10.0
  cosmic_n2%value = 1.0
  cosmic_n2%optimize = .false.
  cosmic_n2%min = 0.01
  cosmic_n2%max = 10.0
  cosmic_alpha0%value = 1.0
  cosmic_alpha0%optimize = .false.
  cosmic_alpha0%min = 0.01
  cosmic_alpha0%max = 10.0
  cosmic_alpha1%value = 1.0
  cosmic_alpha1%optimize = .false.
  cosmic_alpha1%min = 0.01
  cosmic_alpha1%max = 10.0
  cosmic_l30%value = 106.1942
  cosmic_l30%optimize = .false.
  cosmic_l30%min = 26.56
  cosmic_l30%max = 424.78
  cosmic_l31%value = 40.9879
  cosmic_l31%optimize = .false.
  cosmic_l31%min = -118.3
  cosmic_l31%max = 200.28
  cosmic_lw0%value = 0.1783
  cosmic_lw0%optimize = .false.
  cosmic_lw0%min = 0.0
  cosmic_lw0%max = 0.2
  cosmic_lw1%value = 0.0
  cosmic_lw1%optimize = .false.
  cosmic_lw1%min = 0.0
  cosmic_lw1%max = 0.05
/
```

