# Snow - Case 1 {#snow_1}

[TOC]

Parameters for the degree-day snow module.

**Namelist**: `snow_1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [snow_threshold_temperature](#snow_threshold_temperature) | type(parameter_t) | yes | yes | Threshold for rain and snow partitioning [degC] |
| [degree_day_factor_forest](#degree_day_factor_forest) | type(parameter_t) | yes | yes | Degree-day factor for forest [mm d-1 degC-1] |
| [degree_day_factor_impervious](#degree_day_factor_impervious) | type(parameter_t) | yes | yes | Degree-day factor for impervious areas [mm d-1 degC-1] |
| [degree_day_factor_pervious](#degree_day_factor_pervious) | type(parameter_t) | yes | yes | Degree-day factor for pervious areas [mm d-1 degC-1] |
| [degree_day_factor_precipitation](#degree_day_factor_precipitation) | type(parameter_t) | yes | yes | Precipitation-dependent degree-day factor increase [degC-1] |
| [max_degree_day_factor_forest](#max_degree_day_factor_forest) | type(parameter_t) | yes | yes | Maximum degree-day factor for forest [mm d-1 degC-1] |
| [max_degree_day_factor_impervious](#max_degree_day_factor_impervious) | type(parameter_t) | yes | yes | Maximum degree-day factor for impervious areas [mm d-1 degC-1] |
| [max_degree_day_factor_pervious](#max_degree_day_factor_pervious) | type(parameter_t) | yes | yes | Maximum degree-day factor for pervious areas [mm d-1 degC-1] |

## Field details

### snow_threshold_temperature

Threshold for rain and snow partitioning [degC] `snow_threshold_temperature`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: -2.0, max: 2.0}`

Components:
- `snow_threshold_temperature%value`: `real(dp)`; declared required yes; input required yes
- `snow_threshold_temperature%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `snow_threshold_temperature%min`: `real(dp)`; declared required yes; input required no; default `-2.0` (object default)
- `snow_threshold_temperature%max`: `real(dp)`; declared required yes; input required no; default `2.0` (object default)

### degree_day_factor_forest

Degree-day factor for forest [mm d-1 degC-1] `degree_day_factor_forest`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0001, max: 4.0}`

Components:
- `degree_day_factor_forest%value`: `real(dp)`; declared required yes; input required yes
- `degree_day_factor_forest%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `degree_day_factor_forest%min`: `real(dp)`; declared required yes; input required no; default `0.0001` (object default)
- `degree_day_factor_forest%max`: `real(dp)`; declared required yes; input required no; default `4.0` (object default)

### degree_day_factor_impervious

Degree-day factor for impervious areas [mm d-1 degC-1] `degree_day_factor_impervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 1.0}`

Components:
- `degree_day_factor_impervious%value`: `real(dp)`; declared required yes; input required yes
- `degree_day_factor_impervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `degree_day_factor_impervious%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `degree_day_factor_impervious%max`: `real(dp)`; declared required yes; input required no; default `1.0` (object default)

### degree_day_factor_pervious

Degree-day factor for pervious areas [mm d-1 degC-1] `degree_day_factor_pervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 2.0}`

Components:
- `degree_day_factor_pervious%value`: `real(dp)`; declared required yes; input required yes
- `degree_day_factor_pervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `degree_day_factor_pervious%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `degree_day_factor_pervious%max`: `real(dp)`; declared required yes; input required no; default `2.0` (object default)

### degree_day_factor_precipitation

Precipitation-dependent degree-day factor increase [degC-1] `degree_day_factor_precipitation`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.1, max: 0.9}`

Components:
- `degree_day_factor_precipitation%value`: `real(dp)`; declared required yes; input required yes
- `degree_day_factor_precipitation%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `degree_day_factor_precipitation%min`: `real(dp)`; declared required yes; input required no; default `0.1` (object default)
- `degree_day_factor_precipitation%max`: `real(dp)`; declared required yes; input required no; default `0.9` (object default)

### max_degree_day_factor_forest

Maximum degree-day factor for forest [mm d-1 degC-1] `max_degree_day_factor_forest`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 8.0}`

Components:
- `max_degree_day_factor_forest%value`: `real(dp)`; declared required yes; input required yes
- `max_degree_day_factor_forest%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `max_degree_day_factor_forest%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `max_degree_day_factor_forest%max`: `real(dp)`; declared required yes; input required no; default `8.0` (object default)

### max_degree_day_factor_impervious

Maximum degree-day factor for impervious areas [mm d-1 degC-1] `max_degree_day_factor_impervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 8.0}`

Components:
- `max_degree_day_factor_impervious%value`: `real(dp)`; declared required yes; input required yes
- `max_degree_day_factor_impervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `max_degree_day_factor_impervious%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `max_degree_day_factor_impervious%max`: `real(dp)`; declared required yes; input required no; default `8.0` (object default)

### max_degree_day_factor_pervious

Maximum degree-day factor for pervious areas [mm d-1 degC-1] `max_degree_day_factor_pervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 8.0}`

Components:
- `max_degree_day_factor_pervious%value`: `real(dp)`; declared required yes; input required yes
- `max_degree_day_factor_pervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `max_degree_day_factor_pervious%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `max_degree_day_factor_pervious%max`: `real(dp)`; declared required yes; input required no; default `8.0` (object default)

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
&snow_1
  snow_threshold_temperature%value = 1.0
  snow_threshold_temperature%optimize = .true.
  snow_threshold_temperature%min = -2.0
  snow_threshold_temperature%max = 2.0
  degree_day_factor_forest%value = 1.5
  degree_day_factor_forest%optimize = .true.
  degree_day_factor_forest%min = 0.0001
  degree_day_factor_forest%max = 4.0
  degree_day_factor_impervious%value = 0.5
  degree_day_factor_impervious%optimize = .true.
  degree_day_factor_impervious%min = 0.0
  degree_day_factor_impervious%max = 1.0
  degree_day_factor_pervious%value = 0.5
  degree_day_factor_pervious%optimize = .true.
  degree_day_factor_pervious%min = 0.0
  degree_day_factor_pervious%max = 2.0
  degree_day_factor_precipitation%value = 0.5
  degree_day_factor_precipitation%optimize = .true.
  degree_day_factor_precipitation%min = 0.1
  degree_day_factor_precipitation%max = 0.9
  max_degree_day_factor_forest%value = 3.0
  max_degree_day_factor_forest%optimize = .true.
  max_degree_day_factor_forest%min = 0.0
  max_degree_day_factor_forest%max = 8.0
  max_degree_day_factor_impervious%value = 3.5
  max_degree_day_factor_impervious%optimize = .true.
  max_degree_day_factor_impervious%min = 0.0
  max_degree_day_factor_impervious%max = 8.0
  max_degree_day_factor_pervious%value = 4.0
  max_degree_day_factor_pervious%optimize = .true.
  max_degree_day_factor_pervious%min = 0.0
  max_degree_day_factor_pervious%max = 8.0
/
```

