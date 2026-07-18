# PET - Case 1 {#pet_1}

[TOC]

Parameters for the Hargreaves-Samani PET method.

**Namelist**: `pet_1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [correction_factor_min](#correction_factor_min) | type(parameter_t) | yes | yes | Minimum PET aspect-correction factor |
| [correction_factor_max](#correction_factor_max) | type(parameter_t) | yes | yes | Maximum PET aspect-correction factor |
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
- Default: `{lower_bound: 0.7, upper_bound: 1.3}`

Components:
- `correction_factor_min%value`: `real(dp)`; declared required yes; input required yes
- `correction_factor_min%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `correction_factor_min%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.7` (object default)
- `correction_factor_min%upper_bound`: `real(dp)`; declared required yes; input required no; default `1.3` (object default)

### correction_factor_max

Maximum PET aspect-correction factor `correction_factor_max`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.0, upper_bound: 0.2}`

Components:
- `correction_factor_max%value`: `real(dp)`; declared required yes; input required yes
- `correction_factor_max%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `correction_factor_max%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `correction_factor_max%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.2` (object default)

### aspect_threshold

Aspect threshold `aspect_threshold`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 160.0, upper_bound: 200.0}`

Components:
- `aspect_threshold%value`: `real(dp)`; declared required yes; input required yes
- `aspect_threshold%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `aspect_threshold%lower_bound`: `real(dp)`; declared required yes; input required no; default `160.0` (object default)
- `aspect_threshold%upper_bound`: `real(dp)`; declared required yes; input required no; default `200.0` (object default)

### hargreaves_samani_coefficient

Hargreaves-Samani coefficient `hargreaves_samani_coefficient`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.0016, upper_bound: 0.003}`

Components:
- `hargreaves_samani_coefficient%value`: `real(dp)`; declared required yes; input required yes
- `hargreaves_samani_coefficient%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `hargreaves_samani_coefficient%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.0016` (object default)
- `hargreaves_samani_coefficient%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.003` (object default)

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
&pet_1
  correction_factor_min%value = 0.93
  correction_factor_min%optimize = .true.
  correction_factor_min%lower_bound = 0.7
  correction_factor_min%upper_bound = 1.3
  correction_factor_max%value = 0.19
  correction_factor_max%optimize = .true.
  correction_factor_max%lower_bound = 0.0
  correction_factor_max%upper_bound = 0.2
  aspect_threshold%value = 171.0
  aspect_threshold%optimize = .true.
  aspect_threshold%lower_bound = 160.0
  aspect_threshold%upper_bound = 200.0
  hargreaves_samani_coefficient%value = 0.0023
  hargreaves_samani_coefficient%optimize = .true.
  hargreaves_samani_coefficient%lower_bound = 0.0016
  hargreaves_samani_coefficient%upper_bound = 0.003
/
```

