# Soil moisture - Case 1 {#soil_moisture_1}

[TOC]

Feddes ET reduction with multi-layer infiltration and Brooks-Corey parameterization.

**Namelist**: `soil_moisture_1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [organic_matter_forest](#organic_matter_forest) | type(parameter_t) | yes | yes | Organic matter content for forest |
| [organic_matter_impervious](#organic_matter_impervious) | type(parameter_t) | yes | yes | Organic matter content for impervious areas |
| [organic_matter_pervious](#organic_matter_pervious) | type(parameter_t) | yes | yes | Organic matter content for pervious areas |
| [ptf_lower_66_5_constant](#ptf_lower_66_5_constant) | type(parameter_t) | yes | yes | PTF constant below 66.5 percent sand |
| [ptf_lower_66_5_clay](#ptf_lower_66_5_clay) | type(parameter_t) | yes | yes | PTF clay multiplier below 66.5 percent sand |
| [ptf_lower_66_5_bulk_density](#ptf_lower_66_5_bulk_density) | type(parameter_t) | yes | yes | PTF bulk-density multiplier below 66.5 percent sand |
| [ptf_upper_66_5_constant](#ptf_upper_66_5_constant) | type(parameter_t) | yes | yes | PTF constant above 66.5 percent sand |
| [ptf_upper_66_5_clay](#ptf_upper_66_5_clay) | type(parameter_t) | yes | yes | PTF clay multiplier above 66.5 percent sand |
| [ptf_upper_66_5_bulk_density](#ptf_upper_66_5_bulk_density) | type(parameter_t) | yes | yes | PTF bulk-density multiplier above 66.5 percent sand |
| [ptf_ks_constant](#ptf_ks_constant) | type(parameter_t) | yes | yes | PTF constant for saturated hydraulic conductivity |
| [ptf_ks_sand](#ptf_ks_sand) | type(parameter_t) | yes | yes | PTF sand multiplier for saturated hydraulic conductivity |
| [ptf_ks_clay](#ptf_ks_clay) | type(parameter_t) | yes | yes | PTF clay multiplier for saturated hydraulic conductivity |
| [root_fraction_forest](#root_fraction_forest) | type(parameter_t) | yes | yes | Root-fraction coefficient for forest |
| [root_fraction_impervious](#root_fraction_impervious) | type(parameter_t) | yes | yes | Root-fraction coefficient for impervious areas |
| [root_fraction_pervious](#root_fraction_pervious) | type(parameter_t) | yes | yes | Root-fraction coefficient for pervious areas |
| [infiltration_shape_factor](#infiltration_shape_factor) | type(parameter_t) | yes | yes | Shape factor partitioning effective precipitation into runoff and infiltration |

## Field details

### organic_matter_forest

Organic matter content for forest `organic_matter_forest`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 5.0, upper_bound: 10.0}`

Components:
- `organic_matter_forest%value`: `real(dp)`; declared required yes; input required yes
- `organic_matter_forest%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `organic_matter_forest%lower_bound`: `real(dp)`; declared required yes; input required no; default `5.0` (object default)
- `organic_matter_forest%upper_bound`: `real(dp)`; declared required yes; input required no; default `10.0` (object default)

### organic_matter_impervious

Organic matter content for impervious areas `organic_matter_impervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.0, upper_bound: 1.0}`

Components:
- `organic_matter_impervious%value`: `real(dp)`; declared required yes; input required yes
- `organic_matter_impervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `organic_matter_impervious%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `organic_matter_impervious%upper_bound`: `real(dp)`; declared required yes; input required no; default `1.0` (object default)

### organic_matter_pervious

Organic matter content for pervious areas `organic_matter_pervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 1.0, upper_bound: 5.0}`

Components:
- `organic_matter_pervious%value`: `real(dp)`; declared required yes; input required yes
- `organic_matter_pervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `organic_matter_pervious%lower_bound`: `real(dp)`; declared required yes; input required no; default `1.0` (object default)
- `organic_matter_pervious%upper_bound`: `real(dp)`; declared required yes; input required no; default `5.0` (object default)

### ptf_lower_66_5_constant

PTF constant below 66.5 percent sand `ptf_lower_66_5_constant`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.75, upper_bound: 0.8}`

Components:
- `ptf_lower_66_5_constant%value`: `real(dp)`; declared required yes; input required yes
- `ptf_lower_66_5_constant%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_lower_66_5_constant%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.75` (object default)
- `ptf_lower_66_5_constant%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.8` (object default)

### ptf_lower_66_5_clay

PTF clay multiplier below 66.5 percent sand `ptf_lower_66_5_clay`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.0008, upper_bound: 0.0012}`

Components:
- `ptf_lower_66_5_clay%value`: `real(dp)`; declared required yes; input required yes
- `ptf_lower_66_5_clay%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_lower_66_5_clay%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.0008` (object default)
- `ptf_lower_66_5_clay%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.0012` (object default)

### ptf_lower_66_5_bulk_density

PTF bulk-density multiplier below 66.5 percent sand `ptf_lower_66_5_bulk_density`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: -0.27, upper_bound: -0.25}`

Components:
- `ptf_lower_66_5_bulk_density%value`: `real(dp)`; declared required yes; input required yes
- `ptf_lower_66_5_bulk_density%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_lower_66_5_bulk_density%lower_bound`: `real(dp)`; declared required yes; input required no; default `-0.27` (object default)
- `ptf_lower_66_5_bulk_density%upper_bound`: `real(dp)`; declared required yes; input required no; default `-0.25` (object default)

### ptf_upper_66_5_constant

PTF constant above 66.5 percent sand `ptf_upper_66_5_constant`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.8, upper_bound: 0.9}`

Components:
- `ptf_upper_66_5_constant%value`: `real(dp)`; declared required yes; input required yes
- `ptf_upper_66_5_constant%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_upper_66_5_constant%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.8` (object default)
- `ptf_upper_66_5_constant%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.9` (object default)

### ptf_upper_66_5_clay

PTF clay multiplier above 66.5 percent sand `ptf_upper_66_5_clay`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: -0.0012, upper_bound: -0.0008}`

Components:
- `ptf_upper_66_5_clay%value`: `real(dp)`; declared required yes; input required yes
- `ptf_upper_66_5_clay%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_upper_66_5_clay%lower_bound`: `real(dp)`; declared required yes; input required no; default `-0.0012` (object default)
- `ptf_upper_66_5_clay%upper_bound`: `real(dp)`; declared required yes; input required no; default `-0.0008` (object default)

### ptf_upper_66_5_bulk_density

PTF bulk-density multiplier above 66.5 percent sand `ptf_upper_66_5_bulk_density`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: -0.35, upper_bound: -0.3}`

Components:
- `ptf_upper_66_5_bulk_density%value`: `real(dp)`; declared required yes; input required yes
- `ptf_upper_66_5_bulk_density%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_upper_66_5_bulk_density%lower_bound`: `real(dp)`; declared required yes; input required no; default `-0.35` (object default)
- `ptf_upper_66_5_bulk_density%upper_bound`: `real(dp)`; declared required yes; input required no; default `-0.3` (object default)

### ptf_ks_constant

PTF constant for saturated hydraulic conductivity `ptf_ks_constant`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: -1.2, upper_bound: -0.285}`

Components:
- `ptf_ks_constant%value`: `real(dp)`; declared required yes; input required yes
- `ptf_ks_constant%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_ks_constant%lower_bound`: `real(dp)`; declared required yes; input required no; default `-1.2` (object default)
- `ptf_ks_constant%upper_bound`: `real(dp)`; declared required yes; input required no; default `-0.285` (object default)

### ptf_ks_sand

PTF sand multiplier for saturated hydraulic conductivity `ptf_ks_sand`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.006, upper_bound: 0.026}`

Components:
- `ptf_ks_sand%value`: `real(dp)`; declared required yes; input required yes
- `ptf_ks_sand%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_ks_sand%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.006` (object default)
- `ptf_ks_sand%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.026` (object default)

### ptf_ks_clay

PTF clay multiplier for saturated hydraulic conductivity `ptf_ks_clay`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.003, upper_bound: 0.013}`

Components:
- `ptf_ks_clay%value`: `real(dp)`; declared required yes; input required yes
- `ptf_ks_clay%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_ks_clay%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.003` (object default)
- `ptf_ks_clay%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.013` (object default)

### root_fraction_forest

Root-fraction coefficient for forest `root_fraction_forest`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.9, upper_bound: 0.999}`

Components:
- `root_fraction_forest%value`: `real(dp)`; declared required yes; input required yes
- `root_fraction_forest%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `root_fraction_forest%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.9` (object default)
- `root_fraction_forest%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.999` (object default)

### root_fraction_impervious

Root-fraction coefficient for impervious areas `root_fraction_impervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.9, upper_bound: 0.95}`

Components:
- `root_fraction_impervious%value`: `real(dp)`; declared required yes; input required yes
- `root_fraction_impervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `root_fraction_impervious%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.9` (object default)
- `root_fraction_impervious%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.95` (object default)

### root_fraction_pervious

Root-fraction coefficient for pervious areas `root_fraction_pervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 0.001, upper_bound: 0.09}`

Components:
- `root_fraction_pervious%value`: `real(dp)`; declared required yes; input required yes
- `root_fraction_pervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `root_fraction_pervious%lower_bound`: `real(dp)`; declared required yes; input required no; default `0.001` (object default)
- `root_fraction_pervious%upper_bound`: `real(dp)`; declared required yes; input required no; default `0.09` (object default)

### infiltration_shape_factor

Shape factor partitioning effective precipitation into runoff and infiltration `infiltration_shape_factor`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{lower_bound: 1.0, upper_bound: 4.0}`

Components:
- `infiltration_shape_factor%value`: `real(dp)`; declared required yes; input required yes
- `infiltration_shape_factor%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `infiltration_shape_factor%lower_bound`: `real(dp)`; declared required yes; input required no; default `1.0` (object default)
- `infiltration_shape_factor%upper_bound`: `real(dp)`; declared required yes; input required no; default `4.0` (object default)

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
&soil_moisture_1
  organic_matter_forest%value = 7.0
  organic_matter_forest%optimize = .true.
  organic_matter_forest%lower_bound = 5.0
  organic_matter_forest%upper_bound = 10.0
  organic_matter_impervious%value = 0.5
  organic_matter_impervious%optimize = .true.
  organic_matter_impervious%lower_bound = 0.0
  organic_matter_impervious%upper_bound = 1.0
  organic_matter_pervious%value = 2.5
  organic_matter_pervious%optimize = .true.
  organic_matter_pervious%lower_bound = 1.0
  organic_matter_pervious%upper_bound = 5.0
  ptf_lower_66_5_constant%value = 0.788
  ptf_lower_66_5_constant%optimize = .true.
  ptf_lower_66_5_constant%lower_bound = 0.75
  ptf_lower_66_5_constant%upper_bound = 0.8
  ptf_lower_66_5_clay%value = 0.001
  ptf_lower_66_5_clay%optimize = .true.
  ptf_lower_66_5_clay%lower_bound = 0.0008
  ptf_lower_66_5_clay%upper_bound = 0.0012
  ptf_lower_66_5_bulk_density%value = -0.263
  ptf_lower_66_5_bulk_density%optimize = .true.
  ptf_lower_66_5_bulk_density%lower_bound = -0.27
  ptf_lower_66_5_bulk_density%upper_bound = -0.25
  ptf_upper_66_5_constant%value = 0.8907
  ptf_upper_66_5_constant%optimize = .true.
  ptf_upper_66_5_constant%lower_bound = 0.8
  ptf_upper_66_5_constant%upper_bound = 0.9
  ptf_upper_66_5_clay%value = -0.001
  ptf_upper_66_5_clay%optimize = .true.
  ptf_upper_66_5_clay%lower_bound = -0.0012
  ptf_upper_66_5_clay%upper_bound = -0.0008
  ptf_upper_66_5_bulk_density%value = -0.322
  ptf_upper_66_5_bulk_density%optimize = .true.
  ptf_upper_66_5_bulk_density%lower_bound = -0.35
  ptf_upper_66_5_bulk_density%upper_bound = -0.3
  ptf_ks_constant%value = -0.585
  ptf_ks_constant%optimize = .true.
  ptf_ks_constant%lower_bound = -1.2
  ptf_ks_constant%upper_bound = -0.285
  ptf_ks_sand%value = 0.0125
  ptf_ks_sand%optimize = .true.
  ptf_ks_sand%lower_bound = 0.006
  ptf_ks_sand%upper_bound = 0.026
  ptf_ks_clay%value = 0.0063
  ptf_ks_clay%optimize = .true.
  ptf_ks_clay%lower_bound = 0.003
  ptf_ks_clay%upper_bound = 0.013
  root_fraction_forest%value = 0.97
  root_fraction_forest%optimize = .true.
  root_fraction_forest%lower_bound = 0.9
  root_fraction_forest%upper_bound = 0.999
  root_fraction_impervious%value = 0.93
  root_fraction_impervious%optimize = .true.
  root_fraction_impervious%lower_bound = 0.9
  root_fraction_impervious%upper_bound = 0.95
  root_fraction_pervious%value = 0.02
  root_fraction_pervious%optimize = .true.
  root_fraction_pervious%lower_bound = 0.001
  root_fraction_pervious%upper_bound = 0.09
  infiltration_shape_factor%value = 1.75
  infiltration_shape_factor%optimize = .true.
  infiltration_shape_factor%lower_bound = 1.0
  infiltration_shape_factor%upper_bound = 4.0
/
```

