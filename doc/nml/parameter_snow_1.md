# Snow - Case 1 {#snow_1}

[TOC]

Parameters for the degree-day snow module.

**Namelist**: `snow_1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [snow_threshold_temperature](#snow_threshold_temperature) | type(parameter_t) | yes | yes | Threshold for rain and snow partitioning [degC] |
| [degree_day_factor_forest](#degree_day_factor_forest) | type(parameter_t) | yes | yes | Degree-day factor for forest [m degC-1] |
| [degree_day_factor_impervious](#degree_day_factor_impervious) | type(parameter_t) | yes | yes | Degree-day factor for impervious areas [m degC-1] |
| [degree_day_factor_pervious](#degree_day_factor_pervious) | type(parameter_t) | yes | yes | Degree-day factor for pervious areas [m degC-1] |
| [degree_day_factor_precipitation](#degree_day_factor_precipitation) | type(parameter_t) | yes | yes | Precipitation-dependent degree-day factor increase [degC-1] |
| [max_degree_day_factor_forest](#max_degree_day_factor_forest) | type(parameter_t) | yes | yes | Maximum degree-day factor for forest [m degC-1] |
| [max_degree_day_factor_impervious](#max_degree_day_factor_impervious) | type(parameter_t) | yes | yes | Maximum degree-day factor for impervious areas [m degC-1] |
| [max_degree_day_factor_pervious](#max_degree_day_factor_pervious) | type(parameter_t) | yes | yes | Maximum degree-day factor for pervious areas [m degC-1] |

## Field details

### snow_threshold_temperature

Threshold for rain and snow partitioning [degC] `snow_threshold_temperature`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: -2.0, upper_bound: 2.0}`

Components:
- `snow_threshold_temperature%value`: `real(dp)`; declared required yes; input required yes
- `snow_threshold_temperature%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `snow_threshold_temperature%lower_bound`: `real(dp)`; declared required yes; input required no; default `-2.0` (object default)
- `snow_threshold_temperature%upper_bound`: `real(dp)`; declared required yes; input required no; default `2.0` (object default)

### degree_day_factor_forest

Degree-day factor for forest [m degC-1] `degree_day_factor_forest`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.0001, upper_bound: 4.0}`

Components:
- `degree_day_factor_forest%value`: `real(dp)`; declared required yes; input required yes
- `degree_day_factor_forest%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `degree_day_factor_forest%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.0001` (object default)
- `degree_day_factor_forest%upper_bound`: `real(dp)`; declared required yes; input required no; default `4.0` (object default)

### degree_day_factor_impervious

Degree-day factor for impervious areas [m degC-1] `degree_day_factor_impervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.0, upper_bound: 1.0}`

Components:
- `degree_day_factor_impervious%value`: `real(dp)`; declared required yes; input required yes
- `degree_day_factor_impervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `degree_day_factor_impervious%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `degree_day_factor_impervious%upper_bound`: `real(dp)`; declared required yes; input required no; default `1.0` (object default)

### degree_day_factor_pervious

Degree-day factor for pervious areas [m degC-1] `degree_day_factor_pervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.0, upper_bound: 2.0}`

Components:
- `degree_day_factor_pervious%value`: `real(dp)`; declared required yes; input required yes
- `degree_day_factor_pervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `degree_day_factor_pervious%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `degree_day_factor_pervious%upper_bound`: `real(dp)`; declared required yes; input required no; default `2.0` (object default)

### degree_day_factor_precipitation

Precipitation-dependent degree-day factor increase [degC-1] `degree_day_factor_precipitation`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.1, upper_bound: 0.9}`

Components:
- `degree_day_factor_precipitation%value`: `real(dp)`; declared required yes; input required yes
- `degree_day_factor_precipitation%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `degree_day_factor_precipitation%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.1` (object default)
- `degree_day_factor_precipitation%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.9` (object default)

### max_degree_day_factor_forest

Maximum degree-day factor for forest [m degC-1] `max_degree_day_factor_forest`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.0, upper_bound: 8.0}`

Components:
- `max_degree_day_factor_forest%value`: `real(dp)`; declared required yes; input required yes
- `max_degree_day_factor_forest%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `max_degree_day_factor_forest%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `max_degree_day_factor_forest%upper_bound`: `real(dp)`; declared required yes; input required no; default `8.0` (object default)

### max_degree_day_factor_impervious

Maximum degree-day factor for impervious areas [m degC-1] `max_degree_day_factor_impervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.0, upper_bound: 8.0}`

Components:
- `max_degree_day_factor_impervious%value`: `real(dp)`; declared required yes; input required yes
- `max_degree_day_factor_impervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `max_degree_day_factor_impervious%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `max_degree_day_factor_impervious%upper_bound`: `real(dp)`; declared required yes; input required no; default `8.0` (object default)

### max_degree_day_factor_pervious

Maximum degree-day factor for pervious areas [m degC-1] `max_degree_day_factor_pervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.0, upper_bound: 8.0}`

Components:
- `max_degree_day_factor_pervious%value`: `real(dp)`; declared required yes; input required yes
- `max_degree_day_factor_pervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `max_degree_day_factor_pervious%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `max_degree_day_factor_pervious%upper_bound`: `real(dp)`; declared required yes; input required no; default `8.0` (object default)

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
&snow_1
  snow_threshold_temperature%value = 1.0
  snow_threshold_temperature%optimize = .true.
  snow_threshold_temperature%lower_bound = -2.0
  snow_threshold_temperature%upper_bound = 2.0
  degree_day_factor_forest%value = 1.5
  degree_day_factor_forest%optimize = .true.
  degree_day_factor_forest%lower_bound = 0.0001
  degree_day_factor_forest%upper_bound = 4.0
  degree_day_factor_impervious%value = 0.5
  degree_day_factor_impervious%optimize = .true.
  degree_day_factor_impervious%lower_bound = 0.0
  degree_day_factor_impervious%upper_bound = 1.0
  degree_day_factor_pervious%value = 0.5
  degree_day_factor_pervious%optimize = .true.
  degree_day_factor_pervious%lower_bound = 0.0
  degree_day_factor_pervious%upper_bound = 2.0
  degree_day_factor_precipitation%value = 0.5
  degree_day_factor_precipitation%optimize = .true.
  degree_day_factor_precipitation%lower_bound = 0.1
  degree_day_factor_precipitation%upper_bound = 0.9
  max_degree_day_factor_forest%value = 3.0
  max_degree_day_factor_forest%optimize = .true.
  max_degree_day_factor_forest%lower_bound = 0.0
  max_degree_day_factor_forest%upper_bound = 8.0
  max_degree_day_factor_impervious%value = 3.5
  max_degree_day_factor_impervious%optimize = .true.
  max_degree_day_factor_impervious%lower_bound = 0.0
  max_degree_day_factor_impervious%upper_bound = 8.0
  max_degree_day_factor_pervious%value = 4.0
  max_degree_day_factor_pervious%optimize = .true.
  max_degree_day_factor_pervious%lower_bound = 0.0
  max_degree_day_factor_pervious%upper_bound = 8.0
/
```

