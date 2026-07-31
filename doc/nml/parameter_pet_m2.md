# PET - Case -2 {#pet_m2}

[TOC]

Parameters for aspect correction of externally supplied PET.

**Namelist**: `pet_m2`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [correction_factor_min](#correction_factor_min) | type(parameter_t) | yes | yes | Minimum PET aspect-correction factor |
| [correction_factor_max](#correction_factor_max) | type(parameter_t) | yes | yes | Maximum PET aspect-correction factor |
| [aspect_threshold](#aspect_threshold) | type(parameter_t) | yes | yes | Aspect threshold |

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

### correction_factor_max

Maximum PET aspect-correction factor `correction_factor_max`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 0.2}`

Components:
- `correction_factor_max%value`: `real(dp)`; declared required yes; input required yes
- `correction_factor_max%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `correction_factor_max%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `correction_factor_max%max`: `real(dp)`; declared required yes; input required no; default `0.2` (object default)

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
&pet_m2
  correction_factor_min%value = 0.9
  correction_factor_min%optimize = .true.
  correction_factor_min%min = 0.7
  correction_factor_min%max = 1.3
  correction_factor_max%value = 0.1
  correction_factor_max%optimize = .true.
  correction_factor_max%min = 0.0
  correction_factor_max%max = 0.2
  aspect_threshold%value = 180.0
  aspect_threshold%optimize = .true.
  aspect_threshold%min = 160.0
  aspect_threshold%max = 200.0
/
```

