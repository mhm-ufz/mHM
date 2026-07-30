# Soil moisture - Case 3 {#soil_moisture_3}

[TOC]

Jarvis ET reduction with global field-capacity dependency on root fraction.

**Namelist**: `soil_moisture_3`

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
| [root_fraction_sand](#root_fraction_sand) | type(parameter_t) | yes | yes | Sand contribution to the root-fraction coefficient |
| [root_fraction_clay](#root_fraction_clay) | type(parameter_t) | yes | yes | Clay contribution to the root-fraction coefficient |
| [field_capacity_min](#field_capacity_min) | type(parameter_t) | yes | yes | Global minimum field capacity |
| [field_capacity_delta](#field_capacity_delta) | type(parameter_t) | yes | yes | Global field-capacity range |
| [jarvis_sm_threshold](#jarvis_sm_threshold) | type(parameter_t) | yes | yes | Jarvis soil-moisture threshold |

## Field details

### organic_matter_forest

Organic matter content for forest `organic_matter_forest`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 20.0}`

Components:
- `organic_matter_forest%value`: `real(dp)`; declared required yes; input required yes
- `organic_matter_forest%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `organic_matter_forest%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `organic_matter_forest%max`: `real(dp)`; declared required yes; input required no; default `20.0` (object default)

### organic_matter_impervious

Organic matter content for impervious areas `organic_matter_impervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 1.0}`

Components:
- `organic_matter_impervious%value`: `real(dp)`; declared required yes; input required yes
- `organic_matter_impervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `organic_matter_impervious%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `organic_matter_impervious%max`: `real(dp)`; declared required yes; input required no; default `1.0` (object default)

### organic_matter_pervious

Organic matter content for pervious areas `organic_matter_pervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 4.0}`

Components:
- `organic_matter_pervious%value`: `real(dp)`; declared required yes; input required yes
- `organic_matter_pervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `organic_matter_pervious%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `organic_matter_pervious%max`: `real(dp)`; declared required yes; input required no; default `4.0` (object default)

### ptf_lower_66_5_constant

PTF constant below 66.5 percent sand `ptf_lower_66_5_constant`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.6462, max: 0.9506}`

Components:
- `ptf_lower_66_5_constant%value`: `real(dp)`; declared required yes; input required yes
- `ptf_lower_66_5_constant%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_lower_66_5_constant%min`: `real(dp)`; declared required yes; input required no; default `0.6462` (object default)
- `ptf_lower_66_5_constant%max`: `real(dp)`; declared required yes; input required no; default `0.9506` (object default)

### ptf_lower_66_5_clay

PTF clay multiplier below 66.5 percent sand `ptf_lower_66_5_clay`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0001, max: 0.0029}`

Components:
- `ptf_lower_66_5_clay%value`: `real(dp)`; declared required yes; input required yes
- `ptf_lower_66_5_clay%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_lower_66_5_clay%min`: `real(dp)`; declared required yes; input required no; default `0.0001` (object default)
- `ptf_lower_66_5_clay%max`: `real(dp)`; declared required yes; input required no; default `0.0029` (object default)

### ptf_lower_66_5_bulk_density

PTF bulk-density multiplier below 66.5 percent sand `ptf_lower_66_5_bulk_density`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: -0.3727, max: -0.1871}`

Components:
- `ptf_lower_66_5_bulk_density%value`: `real(dp)`; declared required yes; input required yes
- `ptf_lower_66_5_bulk_density%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_lower_66_5_bulk_density%min`: `real(dp)`; declared required yes; input required no; default `-0.3727` (object default)
- `ptf_lower_66_5_bulk_density%max`: `real(dp)`; declared required yes; input required no; default `-0.1871` (object default)

### ptf_upper_66_5_constant

PTF constant above 66.5 percent sand `ptf_upper_66_5_constant`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.5358, max: 1.1232}`

Components:
- `ptf_upper_66_5_constant%value`: `real(dp)`; declared required yes; input required yes
- `ptf_upper_66_5_constant%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_upper_66_5_constant%min`: `real(dp)`; declared required yes; input required no; default `0.5358` (object default)
- `ptf_upper_66_5_constant%max`: `real(dp)`; declared required yes; input required no; default `1.1232` (object default)

### ptf_upper_66_5_clay

PTF clay multiplier above 66.5 percent sand `ptf_upper_66_5_clay`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: -0.0055, max: 0.0049}`

Components:
- `ptf_upper_66_5_clay%value`: `real(dp)`; declared required yes; input required yes
- `ptf_upper_66_5_clay%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_upper_66_5_clay%min`: `real(dp)`; declared required yes; input required no; default `-0.0055` (object default)
- `ptf_upper_66_5_clay%max`: `real(dp)`; declared required yes; input required no; default `0.0049` (object default)

### ptf_upper_66_5_bulk_density

PTF bulk-density multiplier above 66.5 percent sand `ptf_upper_66_5_bulk_density`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: -0.5513, max: -0.0913}`

Components:
- `ptf_upper_66_5_bulk_density%value`: `real(dp)`; declared required yes; input required yes
- `ptf_upper_66_5_bulk_density%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_upper_66_5_bulk_density%min`: `real(dp)`; declared required yes; input required no; default `-0.5513` (object default)
- `ptf_upper_66_5_bulk_density%max`: `real(dp)`; declared required yes; input required no; default `-0.0913` (object default)

### ptf_ks_constant

PTF constant for saturated hydraulic conductivity `ptf_ks_constant`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: -1.2, max: -0.285}`

Components:
- `ptf_ks_constant%value`: `real(dp)`; declared required yes; input required yes
- `ptf_ks_constant%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_ks_constant%min`: `real(dp)`; declared required yes; input required no; default `-1.2` (object default)
- `ptf_ks_constant%max`: `real(dp)`; declared required yes; input required no; default `-0.285` (object default)

### ptf_ks_sand

PTF sand multiplier for saturated hydraulic conductivity `ptf_ks_sand`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.006, max: 0.026}`

Components:
- `ptf_ks_sand%value`: `real(dp)`; declared required yes; input required yes
- `ptf_ks_sand%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_ks_sand%min`: `real(dp)`; declared required yes; input required no; default `0.006` (object default)
- `ptf_ks_sand%max`: `real(dp)`; declared required yes; input required no; default `0.026` (object default)

### ptf_ks_clay

PTF clay multiplier for saturated hydraulic conductivity `ptf_ks_clay`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.003, max: 0.013}`

Components:
- `ptf_ks_clay%value`: `real(dp)`; declared required yes; input required yes
- `ptf_ks_clay%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ptf_ks_clay%min`: `real(dp)`; declared required yes; input required no; default `0.003` (object default)
- `ptf_ks_clay%max`: `real(dp)`; declared required yes; input required no; default `0.013` (object default)

### root_fraction_forest

Root-fraction coefficient for forest `root_fraction_forest`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.97, max: 0.985}`

Components:
- `root_fraction_forest%value`: `real(dp)`; declared required yes; input required yes
- `root_fraction_forest%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `root_fraction_forest%min`: `real(dp)`; declared required yes; input required no; default `0.97` (object default)
- `root_fraction_forest%max`: `real(dp)`; declared required yes; input required no; default `0.985` (object default)

### root_fraction_impervious

Root-fraction coefficient for impervious areas `root_fraction_impervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.97, max: 0.985}`

Components:
- `root_fraction_impervious%value`: `real(dp)`; declared required yes; input required yes
- `root_fraction_impervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `root_fraction_impervious%min`: `real(dp)`; declared required yes; input required no; default `0.97` (object default)
- `root_fraction_impervious%max`: `real(dp)`; declared required yes; input required no; default `0.985` (object default)

### root_fraction_pervious

Root-fraction coefficient for pervious areas `root_fraction_pervious`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.97, max: 0.985}`

Components:
- `root_fraction_pervious%value`: `real(dp)`; declared required yes; input required yes
- `root_fraction_pervious%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `root_fraction_pervious%min`: `real(dp)`; declared required yes; input required no; default `0.97` (object default)
- `root_fraction_pervious%max`: `real(dp)`; declared required yes; input required no; default `0.985` (object default)

### infiltration_shape_factor

Shape factor partitioning effective precipitation into runoff and infiltration `infiltration_shape_factor`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 1.0, max: 4.0}`

Components:
- `infiltration_shape_factor%value`: `real(dp)`; declared required yes; input required yes
- `infiltration_shape_factor%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `infiltration_shape_factor%min`: `real(dp)`; declared required yes; input required no; default `1.0` (object default)
- `infiltration_shape_factor%max`: `real(dp)`; declared required yes; input required no; default `4.0` (object default)

### root_fraction_sand

Sand contribution to the root-fraction coefficient `root_fraction_sand`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.001, max: 0.09}`

Components:
- `root_fraction_sand%value`: `real(dp)`; declared required yes; input required yes
- `root_fraction_sand%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `root_fraction_sand%min`: `real(dp)`; declared required yes; input required no; default `0.001` (object default)
- `root_fraction_sand%max`: `real(dp)`; declared required yes; input required no; default `0.09` (object default)

### root_fraction_clay

Clay contribution to the root-fraction coefficient `root_fraction_clay`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.9, max: 0.999}`

Components:
- `root_fraction_clay%value`: `real(dp)`; declared required yes; input required yes
- `root_fraction_clay%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `root_fraction_clay%min`: `real(dp)`; declared required yes; input required no; default `0.9` (object default)
- `root_fraction_clay%max`: `real(dp)`; declared required yes; input required no; default `0.999` (object default)

### field_capacity_min

Global minimum field capacity `field_capacity_min`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.1, max: 0.2}`

Components:
- `field_capacity_min%value`: `real(dp)`; declared required yes; input required yes
- `field_capacity_min%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `field_capacity_min%min`: `real(dp)`; declared required yes; input required no; default `0.1` (object default)
- `field_capacity_min%max`: `real(dp)`; declared required yes; input required no; default `0.2` (object default)

### field_capacity_delta

Global field-capacity range `field_capacity_delta`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.1, max: 0.4}`

Components:
- `field_capacity_delta%value`: `real(dp)`; declared required yes; input required yes
- `field_capacity_delta%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `field_capacity_delta%min`: `real(dp)`; declared required yes; input required no; default `0.1` (object default)
- `field_capacity_delta%max`: `real(dp)`; declared required yes; input required no; default `0.4` (object default)

### jarvis_sm_threshold

Jarvis soil-moisture threshold `jarvis_sm_threshold`

A model parameter with optional calibration metadata.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0, max: 1.0}`

Components:
- `jarvis_sm_threshold%value`: `real(dp)`; declared required yes; input required yes
- `jarvis_sm_threshold%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `jarvis_sm_threshold%min`: `real(dp)`; declared required yes; input required no; default `0.0` (object default)
- `jarvis_sm_threshold%max`: `real(dp)`; declared required yes; input required no; default `1.0` (object default)

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
&soil_moisture_3
  organic_matter_forest%value = 3.4
  organic_matter_forest%optimize = .true.
  organic_matter_forest%min = 0.0
  organic_matter_forest%max = 20.0
  organic_matter_impervious%value = 0.1
  organic_matter_impervious%optimize = .true.
  organic_matter_impervious%min = 0.0
  organic_matter_impervious%max = 1.0
  organic_matter_pervious%value = 0.6
  organic_matter_pervious%optimize = .true.
  organic_matter_pervious%min = 0.0
  organic_matter_pervious%max = 4.0
  ptf_lower_66_5_constant%value = 0.76
  ptf_lower_66_5_constant%optimize = .true.
  ptf_lower_66_5_constant%min = 0.6462
  ptf_lower_66_5_constant%max = 0.9506
  ptf_lower_66_5_clay%value = 0.0009
  ptf_lower_66_5_clay%optimize = .true.
  ptf_lower_66_5_clay%min = 0.0001
  ptf_lower_66_5_clay%max = 0.0029
  ptf_lower_66_5_bulk_density%value = -0.264
  ptf_lower_66_5_bulk_density%optimize = .true.
  ptf_lower_66_5_bulk_density%min = -0.3727
  ptf_lower_66_5_bulk_density%max = -0.1871
  ptf_upper_66_5_constant%value = 0.89
  ptf_upper_66_5_constant%optimize = .true.
  ptf_upper_66_5_constant%min = 0.5358
  ptf_upper_66_5_constant%max = 1.1232
  ptf_upper_66_5_clay%value = -0.001
  ptf_upper_66_5_clay%optimize = .true.
  ptf_upper_66_5_clay%min = -0.0055
  ptf_upper_66_5_clay%max = 0.0049
  ptf_upper_66_5_bulk_density%value = -0.324
  ptf_upper_66_5_bulk_density%optimize = .true.
  ptf_upper_66_5_bulk_density%min = -0.5513
  ptf_upper_66_5_bulk_density%max = -0.0913
  ptf_ks_constant%value = -0.585
  ptf_ks_constant%optimize = .true.
  ptf_ks_constant%min = -1.2
  ptf_ks_constant%max = -0.285
  ptf_ks_sand%value = 0.0125
  ptf_ks_sand%optimize = .true.
  ptf_ks_sand%min = 0.006
  ptf_ks_sand%max = 0.026
  ptf_ks_clay%value = 0.0063
  ptf_ks_clay%optimize = .true.
  ptf_ks_clay%min = 0.003
  ptf_ks_clay%max = 0.013
  root_fraction_forest%value = 0.975
  root_fraction_forest%optimize = .true.
  root_fraction_forest%min = 0.97
  root_fraction_forest%max = 0.985
  root_fraction_impervious%value = 0.975
  root_fraction_impervious%optimize = .true.
  root_fraction_impervious%min = 0.97
  root_fraction_impervious%max = 0.985
  root_fraction_pervious%value = 0.975
  root_fraction_pervious%optimize = .true.
  root_fraction_pervious%min = 0.97
  root_fraction_pervious%max = 0.985
  infiltration_shape_factor%value = 1.75
  infiltration_shape_factor%optimize = .true.
  infiltration_shape_factor%min = 1.0
  infiltration_shape_factor%max = 4.0
  root_fraction_sand%value = 0.09
  root_fraction_sand%optimize = .true.
  root_fraction_sand%min = 0.001
  root_fraction_sand%max = 0.09
  root_fraction_clay%value = 0.98
  root_fraction_clay%optimize = .true.
  root_fraction_clay%min = 0.9
  root_fraction_clay%max = 0.999
  field_capacity_min%value = 0.15
  field_capacity_min%optimize = .false.
  field_capacity_min%min = 0.1
  field_capacity_min%max = 0.2
  field_capacity_delta%value = 0.25
  field_capacity_delta%optimize = .false.
  field_capacity_delta%min = 0.1
  field_capacity_delta%max = 0.4
  jarvis_sm_threshold%value = 0.5
  jarvis_sm_threshold%optimize = .true.
  jarvis_sm_threshold%min = 0.0
  jarvis_sm_threshold%max = 1.0
/
```

