# Interflow - Case 1 {#interflow_1}

[TOC]

Parameters for parallel interflow reservoirs.

**Namelist**: `interflow_1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [storage_capacity_factor](#storage_capacity_factor) | type(parameter_t) | yes | yes | Interflow storage-capacity factor |
| [recession_slope](#recession_slope) | type(parameter_t) | yes | yes | Slope multiplier for interflow recession |
| [fast_recession_forest](#fast_recession_forest) | type(parameter_t) | yes | yes | Forest multiplier for fast interflow recession |
| [slow_recession_ks](#slow_recession_ks) | type(parameter_t) | yes | yes | Saturated-conductivity multiplier for slow interflow recession |
| [slow_recession_exponent](#slow_recession_exponent) | type(parameter_t) | yes | yes | Slow interflow exponent |

## Field details

### storage_capacity_factor

Interflow storage-capacity factor `storage_capacity_factor`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 75.0, upper_bound: 200.0}`

Components:
- `storage_capacity_factor%value`: `real(dp)`; declared required yes; input required yes
- `storage_capacity_factor%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `storage_capacity_factor%lower_bound`: `real(dp)`; declared required yes; input required no; default `75.0` (object default)
- `storage_capacity_factor%upper_bound`: `real(dp)`; declared required yes; input required no; default `200.0` (object default)

### recession_slope

Slope multiplier for interflow recession `recession_slope`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.0, upper_bound: 10.0}`

Components:
- `recession_slope%value`: `real(dp)`; declared required yes; input required yes
- `recession_slope%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `recession_slope%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `recession_slope%upper_bound`: `real(dp)`; declared required yes; input required no; default `10.0` (object default)

### fast_recession_forest

Forest multiplier for fast interflow recession `fast_recession_forest`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 1.0, upper_bound: 3.0}`

Components:
- `fast_recession_forest%value`: `real(dp)`; declared required yes; input required yes
- `fast_recession_forest%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `fast_recession_forest%lower_bound`: `real(dp)`; declared required yes; input required no; default `1.0` (object default)
- `fast_recession_forest%upper_bound`: `real(dp)`; declared required yes; input required no; default `3.0` (object default)

### slow_recession_ks

Saturated-conductivity multiplier for slow interflow recession `slow_recession_ks`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 1.0, upper_bound: 30.0}`

Components:
- `slow_recession_ks%value`: `real(dp)`; declared required yes; input required yes
- `slow_recession_ks%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `slow_recession_ks%lower_bound`: `real(dp)`; declared required yes; input required no; default `1.0` (object default)
- `slow_recession_ks%upper_bound`: `real(dp)`; declared required yes; input required no; default `30.0` (object default)

### slow_recession_exponent

Slow interflow exponent `slow_recession_exponent`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.05, upper_bound: 0.3}`

Components:
- `slow_recession_exponent%value`: `real(dp)`; declared required yes; input required yes
- `slow_recession_exponent%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `slow_recession_exponent%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.05` (object default)
- `slow_recession_exponent%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.3` (object default)

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
&interflow_1
  storage_capacity_factor%value = 85.0
  storage_capacity_factor%optimize = .true.
  storage_capacity_factor%lower_bound = 75.0
  storage_capacity_factor%upper_bound = 200.0
  recession_slope%value = 7.0
  recession_slope%optimize = .true.
  recession_slope%lower_bound = 0.0
  recession_slope%upper_bound = 10.0
  fast_recession_forest%value = 1.5
  fast_recession_forest%optimize = .true.
  fast_recession_forest%lower_bound = 1.0
  fast_recession_forest%upper_bound = 3.0
  slow_recession_ks%value = 15.0
  slow_recession_ks%optimize = .true.
  slow_recession_ks%lower_bound = 1.0
  slow_recession_ks%upper_bound = 30.0
  slow_recession_exponent%value = 0.125
  slow_recession_exponent%optimize = .true.
  slow_recession_exponent%lower_bound = 0.05
  slow_recession_exponent%upper_bound = 0.3
/
```

