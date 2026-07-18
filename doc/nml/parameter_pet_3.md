# PET - Case 3 {#pet_3}

[TOC]

Parameters for the Penman-Monteith PET method.

**Namelist**: `pet_3`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [canopy_height_forest](#canopy_height_forest) | type(parameter_t) | yes | yes | Canopy height for forest |
| [canopy_height_impervious](#canopy_height_impervious) | type(parameter_t) | yes | yes | Canopy height for impervious areas |
| [canopy_height_pervious](#canopy_height_pervious) | type(parameter_t) | yes | yes | Canopy height for pervious areas |
| [displacement_height_coefficient](#displacement_height_coefficient) | type(parameter_t) | yes | yes | Displacement-height coefficient |
| [momentum_roughness_length_coefficient](#momentum_roughness_length_coefficient) | type(parameter_t) | yes | yes | Momentum roughness-length coefficient |
| [heat_roughness_length_coefficient](#heat_roughness_length_coefficient) | type(parameter_t) | yes | yes | Heat roughness-length coefficient |
| [stomatal_resistance](#stomatal_resistance) | type(parameter_t) | yes | yes | Stomatal resistance |

## Field details

### canopy_height_forest

Canopy height for forest `canopy_height_forest`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 15.0, upper_bound: 40.0}`

Components:
- `canopy_height_forest%value`: `real(dp)`; declared required yes; input required yes
- `canopy_height_forest%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `canopy_height_forest%lower_bound`: `real(dp)`; declared required yes; input required no; default `15.0` (object default)
- `canopy_height_forest%upper_bound`: `real(dp)`; declared required yes; input required no; default `40.0` (object default)

### canopy_height_impervious

Canopy height for impervious areas `canopy_height_impervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.01, upper_bound: 0.5}`

Components:
- `canopy_height_impervious%value`: `real(dp)`; declared required yes; input required yes
- `canopy_height_impervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `canopy_height_impervious%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.01` (object default)
- `canopy_height_impervious%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.5` (object default)

### canopy_height_pervious

Canopy height for pervious areas `canopy_height_pervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.1, upper_bound: 5.0}`

Components:
- `canopy_height_pervious%value`: `real(dp)`; declared required yes; input required yes
- `canopy_height_pervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `canopy_height_pervious%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.1` (object default)
- `canopy_height_pervious%upper_bound`: `real(dp)`; declared required yes; input required no; default `5.0` (object default)

### displacement_height_coefficient

Displacement-height coefficient `displacement_height_coefficient`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.5, upper_bound: 0.85}`

Components:
- `displacement_height_coefficient%value`: `real(dp)`; declared required yes; input required yes
- `displacement_height_coefficient%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `displacement_height_coefficient%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.5` (object default)
- `displacement_height_coefficient%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.85` (object default)

### momentum_roughness_length_coefficient

Momentum roughness-length coefficient `momentum_roughness_length_coefficient`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.09, upper_bound: 0.16}`

Components:
- `momentum_roughness_length_coefficient%value`: `real(dp)`; declared required yes; input required yes
- `momentum_roughness_length_coefficient%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `momentum_roughness_length_coefficient%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.09` (object default)
- `momentum_roughness_length_coefficient%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.16` (object default)

### heat_roughness_length_coefficient

Heat roughness-length coefficient `heat_roughness_length_coefficient`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.07, upper_bound: 0.13}`

Components:
- `heat_roughness_length_coefficient%value`: `real(dp)`; declared required yes; input required yes
- `heat_roughness_length_coefficient%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `heat_roughness_length_coefficient%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.07` (object default)
- `heat_roughness_length_coefficient%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.13` (object default)

### stomatal_resistance

Stomatal resistance `stomatal_resistance`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 10.0, upper_bound: 200.0}`

Components:
- `stomatal_resistance%value`: `real(dp)`; declared required yes; input required yes
- `stomatal_resistance%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `stomatal_resistance%lower_bound`: `real(dp)`; declared required yes; input required no; default `10.0` (object default)
- `stomatal_resistance%upper_bound`: `real(dp)`; declared required yes; input required no; default `200.0` (object default)

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
&pet_3
  canopy_height_forest%value = 15.0
  canopy_height_forest%optimize = .true.
  canopy_height_forest%lower_bound = 15.0
  canopy_height_forest%upper_bound = 40.0
  canopy_height_impervious%value = 0.02
  canopy_height_impervious%optimize = .true.
  canopy_height_impervious%lower_bound = 0.01
  canopy_height_impervious%upper_bound = 0.5
  canopy_height_pervious%value = 0.11
  canopy_height_pervious%optimize = .true.
  canopy_height_pervious%lower_bound = 0.1
  canopy_height_pervious%upper_bound = 5.0
  displacement_height_coefficient%value = 0.64
  displacement_height_coefficient%optimize = .true.
  displacement_height_coefficient%lower_bound = 0.5
  displacement_height_coefficient%upper_bound = 0.85
  momentum_roughness_length_coefficient%value = 0.095
  momentum_roughness_length_coefficient%optimize = .true.
  momentum_roughness_length_coefficient%lower_bound = 0.09
  momentum_roughness_length_coefficient%upper_bound = 0.16
  heat_roughness_length_coefficient%value = 0.075
  heat_roughness_length_coefficient%optimize = .true.
  heat_roughness_length_coefficient%lower_bound = 0.07
  heat_roughness_length_coefficient%upper_bound = 0.13
  stomatal_resistance%value = 56.0
  stomatal_resistance%optimize = .true.
  stomatal_resistance%lower_bound = 10.0
  stomatal_resistance%upper_bound = 200.0
/
```

