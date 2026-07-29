# Direct runoff - Case 1 {#direct_runoff_1}

[TOC]

Parameters for sealed-area runoff.

**Namelist**: `direct_runoff_1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [impervious_storage_capacity](#impervious_storage_capacity) | type(parameter_t) | yes | yes | Capacity of impervious storage [mm] |

## Field details

### impervious_storage_capacity

Capacity of impervious storage [mm] `impervious_storage_capacity`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.0, upper_bound: 5.0}`

Components:
- `impervious_storage_capacity%value`: `real(dp)`; declared required yes; input required yes
- `impervious_storage_capacity%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `impervious_storage_capacity%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `impervious_storage_capacity%upper_bound`: `real(dp)`; declared required yes; input required no; default `5.0` (object default)

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
&direct_runoff_1
  impervious_storage_capacity%value = 0.5
  impervious_storage_capacity%optimize = .true.
  impervious_storage_capacity%lower_bound = 0.0
  impervious_storage_capacity%upper_bound = 5.0
/
```

