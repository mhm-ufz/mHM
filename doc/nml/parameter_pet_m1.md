# PET - Case -1 {#pet_m1}

[TOC]

Parameters for LAI correction of externally supplied PET.

**Namelist**: `pet_m1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [pet_a_forest](#pet_a_forest) | type(parameter_t) | yes | yes | PET correction coefficient a for forest |
| [pet_a_impervious](#pet_a_impervious) | type(parameter_t) | yes | yes | PET correction coefficient a for impervious areas |
| [pet_a_pervious](#pet_a_pervious) | type(parameter_t) | yes | yes | PET correction coefficient a for pervious areas |
| [pet_b](#pet_b) | type(parameter_t) | yes | yes | PET correction coefficient b |
| [pet_c](#pet_c) | type(parameter_t) | yes | yes | PET correction coefficient c |

## Field details

### pet_a_forest

PET correction coefficient a for forest `pet_a_forest`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.3, max: 1.3}`

Components:
- `pet_a_forest%value`: `real(dp)`; declared required yes; input required yes
- `pet_a_forest%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `pet_a_forest%min`: `real(dp)`; declared required yes; input required no; default `0.3` (object default)
- `pet_a_forest%max`: `real(dp)`; declared required yes; input required no; default `1.3` (object default)

### pet_a_impervious

PET correction coefficient a for impervious areas `pet_a_impervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.3, max: 1.3}`

Components:
- `pet_a_impervious%value`: `real(dp)`; declared required yes; input required yes
- `pet_a_impervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `pet_a_impervious%min`: `real(dp)`; declared required yes; input required no; default `0.3` (object default)
- `pet_a_impervious%max`: `real(dp)`; declared required yes; input required no; default `1.3` (object default)

### pet_a_pervious

PET correction coefficient a for pervious areas `pet_a_pervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.3, max: 1.3}`

Components:
- `pet_a_pervious%value`: `real(dp)`; declared required yes; input required yes
- `pet_a_pervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `pet_a_pervious%min`: `real(dp)`; declared required yes; input required no; default `0.3` (object default)
- `pet_a_pervious%max`: `real(dp)`; declared required yes; input required no; default `1.3` (object default)

### pet_b

PET correction coefficient b `pet_b`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 1.5}`

Components:
- `pet_b%value`: `real(dp)`; declared required yes; input required yes
- `pet_b%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `pet_b%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `pet_b%max`: `real(dp)`; declared required yes; input required no; default `1.5` (object default)

### pet_c

PET correction coefficient c `pet_c`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: -2.0, max: 0.0}`

Components:
- `pet_c%value`: `real(dp)`; declared required yes; input required yes
- `pet_c%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `pet_c%min`: `real(dp)`; declared required yes; input required no; default `-2.0` (object default)
- `pet_c%max`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)

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
&pet_m1
  pet_a_forest%value = 0.3
  pet_a_forest%optimize = .true.
  pet_a_forest%min = 0.3
  pet_a_forest%max = 1.3
  pet_a_impervious%value = 0.8
  pet_a_impervious%optimize = .true.
  pet_a_impervious%min = 0.3
  pet_a_impervious%max = 1.3
  pet_a_pervious%value = 1.3
  pet_a_pervious%optimize = .true.
  pet_a_pervious%min = 0.3
  pet_a_pervious%max = 1.3
  pet_b%value = 1.5
  pet_b%optimize = .true.
  pet_b%min = 0.0
  pet_b%max = 1.5
  pet_c%value = -0.7
  pet_c%optimize = .true.
  pet_c%min = -2.0
  pet_c%max = 0.0
/
```

