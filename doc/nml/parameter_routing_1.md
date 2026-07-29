# Routing - Case 1 {#routing_1}

[TOC]

Parameters for Muskingum routing.

**Namelist**: `routing_1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [travel_time_constant](#travel_time_constant) | type(parameter_t) | yes | yes | Muskingum travel-time constant |
| [travel_time_river_length](#travel_time_river_length) | type(parameter_t) | yes | yes | River-length contribution to Muskingum travel time |
| [travel_time_river_slope](#travel_time_river_slope) | type(parameter_t) | yes | yes | River-slope contribution to Muskingum travel time |
| [travel_time_impervious](#travel_time_impervious) | type(parameter_t) | yes | yes | Impervious-area contribution to Muskingum travel time |
| [attenuation_river_slope](#attenuation_river_slope) | type(parameter_t) | yes | yes | River-slope contribution to Muskingum attenuation |

## Field details

### travel_time_constant

Muskingum travel-time constant `travel_time_constant`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.31, upper_bound: 0.35}`

Components:
- `travel_time_constant%value`: `real(dp)`; declared required yes; input required yes
- `travel_time_constant%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `travel_time_constant%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.31` (object default)
- `travel_time_constant%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.35` (object default)

### travel_time_river_length

River-length contribution to Muskingum travel time `travel_time_river_length`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.07, upper_bound: 0.08}`

Components:
- `travel_time_river_length%value`: `real(dp)`; declared required yes; input required yes
- `travel_time_river_length%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `travel_time_river_length%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.07` (object default)
- `travel_time_river_length%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.08` (object default)

### travel_time_river_slope

River-slope contribution to Muskingum travel time `travel_time_river_slope`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 1.95, upper_bound: 2.1}`

Components:
- `travel_time_river_slope%value`: `real(dp)`; declared required yes; input required yes
- `travel_time_river_slope%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `travel_time_river_slope%lower_bound`: `real(dp)`; declared required yes; input required no; default `1.95` (object default)
- `travel_time_river_slope%upper_bound`: `real(dp)`; declared required yes; input required no; default `2.1` (object default)

### travel_time_impervious

Impervious-area contribution to Muskingum travel time `travel_time_impervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.09, upper_bound: 0.11}`

Components:
- `travel_time_impervious%value`: `real(dp)`; declared required yes; input required yes
- `travel_time_impervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `travel_time_impervious%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.09` (object default)
- `travel_time_impervious%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.11` (object default)

### attenuation_river_slope

River-slope contribution to Muskingum attenuation `attenuation_river_slope`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.01, upper_bound: 0.5}`

Components:
- `attenuation_river_slope%value`: `real(dp)`; declared required yes; input required yes
- `attenuation_river_slope%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `attenuation_river_slope%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.01` (object default)
- `attenuation_river_slope%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.5` (object default)

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
&routing_1
  travel_time_constant%value = 0.325
  travel_time_constant%optimize = .true.
  travel_time_constant%lower_bound = 0.31
  travel_time_constant%upper_bound = 0.35
  travel_time_river_length%value = 0.075
  travel_time_river_length%optimize = .true.
  travel_time_river_length%lower_bound = 0.07
  travel_time_river_length%upper_bound = 0.08
  travel_time_river_slope%value = 2.0
  travel_time_river_slope%optimize = .true.
  travel_time_river_slope%lower_bound = 1.95
  travel_time_river_slope%upper_bound = 2.1
  travel_time_impervious%value = 0.1
  travel_time_impervious%optimize = .true.
  travel_time_impervious%lower_bound = 0.09
  travel_time_impervious%upper_bound = 0.11
  attenuation_river_slope%value = 0.3
  attenuation_river_slope%optimize = .true.
  attenuation_river_slope%lower_bound = 0.01
  attenuation_river_slope%upper_bound = 0.5
/
```

