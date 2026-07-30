# PET - Case 2 {#pet_2}

[TOC]

Parameters for the Priestley-Taylor PET method.

**Namelist**: `pet_2`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [priestley_taylor_coefficient](#priestley_taylor_coefficient) | type(parameter_t) | yes | yes | Priestley-Taylor coefficient |
| [priestley_taylor_lai_correction](#priestley_taylor_lai_correction) | type(parameter_t) | yes | yes | Priestley-Taylor LAI correction factor |

## Field details

### priestley_taylor_coefficient

Priestley-Taylor coefficient `priestley_taylor_coefficient`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.75, max: 1.75}`

Components:
- `priestley_taylor_coefficient%value`: `real(dp)`; declared required yes; input required yes
- `priestley_taylor_coefficient%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `priestley_taylor_coefficient%min`: `real(dp)`; declared required yes; input required no; default `0.75` (object default)
- `priestley_taylor_coefficient%max`: `real(dp)`; declared required yes; input required no; default `1.75` (object default)

### priestley_taylor_lai_correction

Priestley-Taylor LAI correction factor `priestley_taylor_lai_correction`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: -0.5, max: 0.2}`

Components:
- `priestley_taylor_lai_correction%value`: `real(dp)`; declared required yes; input required yes
- `priestley_taylor_lai_correction%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `priestley_taylor_lai_correction%min`: `real(dp)`; declared required yes; input required no; default `-0.5` (object default)
- `priestley_taylor_lai_correction%max`: `real(dp)`; declared required yes; input required no; default `0.2` (object default)

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
&pet_2
  priestley_taylor_coefficient%value = 1.19
  priestley_taylor_coefficient%optimize = .true.
  priestley_taylor_coefficient%min = 0.75
  priestley_taylor_coefficient%max = 1.75
  priestley_taylor_lai_correction%value = 0.058
  priestley_taylor_lai_correction%optimize = .true.
  priestley_taylor_lai_correction%min = -0.5
  priestley_taylor_lai_correction%max = 0.2
/
```

