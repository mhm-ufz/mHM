# River temperature - Case 1 {#river_temperature_1}

[TOC]

Parameters for river-temperature routing.

**Namelist**: `river_temperature_1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [albedo_water](#albedo_water) | type(parameter_t) | yes | yes | Albedo of open water [-] |
| [pt_a_water](#pt_a_water) | type(parameter_t) | yes | yes | Priestley-Taylor coefficient for open water [-] |
| [emissivity_water](#emissivity_water) | type(parameter_t) | yes | yes | Emissivity of open water [-] |
| [turbulent_heat_exchange_coefficient](#turbulent_heat_exchange_coefficient) | type(parameter_t) | yes | yes | Turbulent heat-exchange coefficient [W m-2 K-1] |

## Field details

### albedo_water

Albedo of open water [-] `albedo_water`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.03, upper_bound: 0.2}`

Components:
- `albedo_water%value`: `real(dp)`; declared required yes; input required yes
- `albedo_water%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `albedo_water%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.03` (object default)
- `albedo_water%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.2` (object default)

### pt_a_water

Priestley-Taylor coefficient for open water [-] `pt_a_water`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 1.2, upper_bound: 2.0}`

Components:
- `pt_a_water%value`: `real(dp)`; declared required yes; input required yes
- `pt_a_water%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `pt_a_water%lower_bound`: `real(dp)`; declared required yes; input required no; default `1.2` (object default)
- `pt_a_water%upper_bound`: `real(dp)`; declared required yes; input required no; default `2.0` (object default)

### emissivity_water

Emissivity of open water [-] `emissivity_water`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.95, upper_bound: 0.99}`

Components:
- `emissivity_water%value`: `real(dp)`; declared required yes; input required yes
- `emissivity_water%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `emissivity_water%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.95` (object default)
- `emissivity_water%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.99` (object default)

### turbulent_heat_exchange_coefficient

Turbulent heat-exchange coefficient [W m-2 K-1] `turbulent_heat_exchange_coefficient`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 10.0, upper_bound: 50.0}`

Components:
- `turbulent_heat_exchange_coefficient%value`: `real(dp)`; declared required yes; input required yes
- `turbulent_heat_exchange_coefficient%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `turbulent_heat_exchange_coefficient%lower_bound`: `real(dp)`; declared required yes; input required no; default `10.0` (object default)
- `turbulent_heat_exchange_coefficient%upper_bound`: `real(dp)`; declared required yes; input required no; default `50.0` (object default)

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
&river_temperature_1
  albedo_water%value = 0.15
  albedo_water%optimize = .false.
  albedo_water%lower_bound = 0.03
  albedo_water%upper_bound = 0.2
  pt_a_water%value = 1.26
  pt_a_water%optimize = .false.
  pt_a_water%lower_bound = 1.2
  pt_a_water%upper_bound = 2.0
  emissivity_water%value = 0.96
  emissivity_water%optimize = .false.
  emissivity_water%lower_bound = 0.95
  emissivity_water%upper_bound = 0.99
  turbulent_heat_exchange_coefficient%value = 20.0
  turbulent_heat_exchange_coefficient%optimize = .false.
  turbulent_heat_exchange_coefficient%lower_bound = 10.0
  turbulent_heat_exchange_coefficient%upper_bound = 50.0
/
```

