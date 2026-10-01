# Lakes - Case 1 {#lakes_1}

[TOC]

Parameters for the mHM lake module (mLM) controlling spill over the dam crest and the
seepage and baseflow from the lake groundwater store. Each parameter is given per lake
per domain as the four components value, optimize-flag, min, max (min/max are the
calibration bounds). Reservoir release-regulation parameters are not included at this
stage (minimal natural-lake scope).

**Namelist**: `lakes_1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [spill_coefficient_of_discharge](#spill_coefficient_of_discharge) | type(parameter_t) array | yes | yes | Spill coefficient of discharge |
| [spill_exponent](#spill_exponent) | type(parameter_t) array | yes | yes | Spill exponent |
| [longterm_baseflow](#longterm_baseflow) | type(parameter_t) array | yes | yes | Long-term baseflow |
| [exponent_for_percolation](#exponent_for_percolation) | type(parameter_t) array | yes | yes | Exponent for percolation |

## Field details

### spill_coefficient_of_discharge

Spill coefficient of discharge `spill_coefficient_of_discharge`

Coefficient of discharge C_d of the spillway in the generalized spill equation
Q_spill = C_d * L_s * deltaH^spill_exponent, where deltaH is the head of the lake
level above the dam crest level (DCL) and L_s is the effective spillway length.
Spill is produced only when the lake level rises above the crest. Based on the ogee
free-overflow weir form (USBR, 1998). Dimensionless.

Summary:
- Type: `type(parameter_t), dimension(n_domains, n_lakes)`
- Declared required: yes
- Input required: yes
- Default: `{min: 1.4, max: 3.0}` (broadcast item default)

Components:
- `spill_coefficient_of_discharge%value`: `real(dp)`; declared required yes; input required yes
- `spill_coefficient_of_discharge%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `spill_coefficient_of_discharge%min`: `real(dp)`; declared required yes; input required no; default `1.4` (item default)
- `spill_coefficient_of_discharge%max`: `real(dp)`; declared required yes; input required no; default `3.0` (item default)

### spill_exponent

Spill exponent `spill_exponent`

Exponent on the head term deltaH in the generalized spill equation
Q_spill = C_d * L_s * deltaH^spill_exponent. The textbook ogee free-overflow value
is 3/2; it is exposed as a calibration parameter here. Dimensionless.

Summary:
- Type: `type(parameter_t), dimension(n_domains, n_lakes)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.5, max: 2.5}` (broadcast item default)

Components:
- `spill_exponent%value`: `real(dp)`; declared required yes; input required yes
- `spill_exponent%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `spill_exponent%min`: `real(dp)`; declared required yes; input required no; default `0.5` (item default)
- `spill_exponent%max`: `real(dp)`; declared required yes; input required no; default `2.5` (item default)

### longterm_baseflow

Long-term baseflow `longterm_baseflow`

Long-term baseflow feed [m3 s-1] representing the sustained groundwater contribution
to lake outflow. Used as the baseline baseflow term alongside the recession-driven
baseflow from the lake groundwater storage.

Summary:
- Type: `type(parameter_t), dimension(n_domains, n_lakes)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.0001, max: 5.0}` (broadcast item default)

Components:
- `longterm_baseflow%value`: `real(dp)`; declared required yes; input required yes
- `longterm_baseflow%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `longterm_baseflow%min`: `real(dp)`; declared required yes; input required no; default `0.0001` (item default)
- `longterm_baseflow%max`: `real(dp)`; declared required yes; input required no; default `5.0` (item default)

### exponent_for_percolation

Exponent for percolation `exponent_for_percolation`

Exponential constant chi controlling the fraction kfrac of the head above the lake
bed that is available for percolation (seepage): kfrac = (h - bed_level) / exp(chi).
Percolation from each lake-bed slice is p = sum(kp * kfrac * H * As/A). A larger head
makes more water available to percolate; a larger chi reduces that fraction.
Dimensionless.

Summary:
- Type: `type(parameter_t), dimension(n_domains, n_lakes)`
- Declared required: yes
- Input required: yes
- Default: `{min: 8.0, max: 15.0}` (broadcast item default)

Components:
- `exponent_for_percolation%value`: `real(dp)`; declared required yes; input required yes
- `exponent_for_percolation%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `exponent_for_percolation%min`: `real(dp)`; declared required yes; input required no; default `8.0` (item default)
- `exponent_for_percolation%max`: `real(dp)`; declared required yes; input required no; default `15.0` (item default)

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
&lakes_1
  spill_coefficient_of_discharge(:,:)%value = 2.0
  spill_coefficient_of_discharge(:,:)%optimize = .true.
  spill_coefficient_of_discharge(:,:)%min = 1.4
  spill_coefficient_of_discharge(:,:)%max = 3.0
  spill_exponent(:,:)%value = 1.5
  spill_exponent(:,:)%optimize = .true.
  spill_exponent(:,:)%min = 0.5
  spill_exponent(:,:)%max = 2.5
  longterm_baseflow(:,:)%value = 4.5
  longterm_baseflow(:,:)%optimize = .true.
  longterm_baseflow(:,:)%min = 0.0001
  longterm_baseflow(:,:)%max = 5.0
  exponent_for_percolation(:,:)%value = 11.0
  exponent_for_percolation(:,:)%optimize = .true.
  exponent_for_percolation(:,:)%min = 8.0
  exponent_for_percolation(:,:)%max = 15.0
/
```

