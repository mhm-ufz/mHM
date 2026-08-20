# PET - Case 1 {#pet_1}

[TOC]

Parameters for the Hargreaves-Samani PET method.

**Namelist**: `pet_1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [correction_factor_min](#correction_factor_min) | type(parameter_t) | yes | yes | Minimum PET aspect-correction factor |
| [correction_factor_delta](#correction_factor_delta) | type(parameter_t) | yes | yes | Delta added to the minimum PET aspect-correction factor |
| [aspect_threshold](#aspect_threshold) | type(parameter_t) | yes | yes | Aspect threshold |
| [hargreaves_samani_coefficient](#hargreaves_samani_coefficient) | type(parameter_t) | yes | yes | Hargreaves-Samani coefficient |

## Field details

### correction_factor_min

Minimum PET aspect-correction factor `correction_factor_min`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.7, max: 1.3}`

Components:
- `correction_factor_min%value`: `real(dp)`; declared required yes; input required yes
- `correction_factor_min%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `correction_factor_min%min`: `real(dp)`; declared required yes; input required no; default `0.7` (object default)
- `correction_factor_min%max`: `real(dp)`; declared required yes; input required no; default `1.3` (object default)

### correction_factor_delta

Delta added to the minimum PET aspect-correction factor `correction_factor_delta`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 0.2}`

Components:
- `correction_factor_delta%value`: `real(dp)`; declared required yes; input required yes
- `correction_factor_delta%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `correction_factor_delta%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `correction_factor_delta%max`: `real(dp)`; declared required yes; input required no; default `0.2` (object default)

### aspect_threshold

Aspect threshold `aspect_threshold`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 160.0, max: 200.0}`

Components:
- `aspect_threshold%value`: `real(dp)`; declared required yes; input required yes
- `aspect_threshold%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `aspect_threshold%min`: `real(dp)`; declared required yes; input required no; default `160.0` (object default)
- `aspect_threshold%max`: `real(dp)`; declared required yes; input required no; default `200.0` (object default)

### hargreaves_samani_coefficient

Hargreaves-Samani coefficient `hargreaves_samani_coefficient`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0016, max: 0.003}`

Components:
- `hargreaves_samani_coefficient%value`: `real(dp)`; declared required yes; input required yes
- `hargreaves_samani_coefficient%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `hargreaves_samani_coefficient%min`: `real(dp)`; declared required yes; input required no; default `0.0016` (object default)
- `hargreaves_samani_coefficient%max`: `real(dp)`; declared required yes; input required no; default `0.003` (object default)

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
&pet_1
  correction_factor_min%value = 0.93
  correction_factor_min%optimize = .true.
  correction_factor_min%min = 0.7
  correction_factor_min%max = 1.3
  correction_factor_delta%value = 0.19
  correction_factor_delta%optimize = .true.
  correction_factor_delta%min = 0.0
  correction_factor_delta%max = 0.2
  aspect_threshold%value = 171.0
  aspect_threshold%optimize = .true.
  aspect_threshold%min = 160.0
  aspect_threshold%max = 200.0
  hargreaves_samani_coefficient%value = 0.0023
  hargreaves_samani_coefficient%optimize = .true.
  hargreaves_samani_coefficient%min = 0.0016
  hargreaves_samani_coefficient%max = 0.003
/
```

