# Lakes - Case 1 {#lakes_1}

[TOC]

Parameters for the mHM lake module (mLM) controlling spill over the dam crest and the
seepage and baseflow from the lake groundwater store. The spill parameters are given per
lake per domain as the four components value, optimize-flag, min, max (min/max are the
calibration bounds). The seepage leakance and the baseflow recession are global (pan-domain,
shared across all lakes, like the routing parameters), so they are scalars. 

**Namelist**: `lakes_1`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [spill_coefficient_of_discharge](#spill_coefficient_of_discharge) | type(parameter_t) array | yes | yes | Spill coefficient of discharge |
| [spill_exponent](#spill_exponent) | type(parameter_t) array | yes | yes | Spill exponent |
| [b](#b) | type(parameter_t) | yes | yes | Effective lakebed thickness |
| [ks_lakebed](#ks_lakebed) | type(parameter_t) | yes | yes | Lakebed saturated hydraulic conductivity |
| [baseflow_recession_coefficient](#baseflow_recession_coefficient) | type(parameter_t) | yes | yes | Baseflow recession coefficient |

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

### b

Effective lakebed thickness `b`

Effective lakebed (colmation / clogging-layer) thickness b [m] in the Darcy lakebed
leakance Ks/b that drives lake seepage: the slice percolation is
L_i = (Ks/b) * (lake_level - bed_level) * (slice_area / lake_area), linear in head.
b is a single global (pan-domain) calibrated constant - there is no global lakebed-thickness
dataset, so only the leakance ratio Ks/b is identifiable; Ks carries the spatial variability.
Lower bound strictly positive.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 0.1, max: 10.0}`

Components:
- `b%value`: `real(dp)`; declared required yes; input required yes
- `b%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `b%min`: `real(dp)`; declared required yes; input required no; default `0.1` (object default)
- `b%max`: `real(dp)`; declared required yes; input required no; default `10.0` (object default)

### ks_lakebed

Lakebed saturated hydraulic conductivity `ks_lakebed`

Saturated hydraulic conductivity Ks of the lakebed [cm d-1] in the Darcy lakebed leakance
Ks/b. Global (pan-domain) stand-in value used until the detached lake-MPR path regionalises
Ks per lake from the Cosby pedotransfer function over the lake footprint; at that point this
scalar is replaced by the per-lake field. Converted to model units (m per time step) and
divided by b to form the per-step leakance.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 1.1, max: 1000.0}`

Components:
- `ks_lakebed%value`: `real(dp)`; declared required yes; input required yes
- `ks_lakebed%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `ks_lakebed%min`: `real(dp)`; declared required yes; input required no; default `1.1` (object default)
- `ks_lakebed%max`: `real(dp)`; declared required yes; input required no; default `1000.0` (object default)

### baseflow_recession_coefficient

Baseflow recession coefficient `baseflow_recession_coefficient`

Lake groundwater-storage recession time constant kbf [d]. Each hydrological time step
the per-step recession coefficient is formed as (tsHydroSecs / 86400) / kbf (the SI
form of the v5 c2TSTu / kBaseFlow conversion); the baseflow is that coefficient times
the lake groundwater storage, so a larger value gives a slower, more sustained baseflow.
Global (pan-domain), shared across all lakes.

Summary:
- Type: `type(parameter_t)`
- Declared required: yes
- Input required: yes
- Default: `{min: 1.0, max: 1000.0}`

Components:
- `baseflow_recession_coefficient%value`: `real(dp)`; declared required yes; input required yes
- `baseflow_recession_coefficient%optimize`: `logical`; declared required no; input required no; default `.false.` (component default)
- `baseflow_recession_coefficient%min`: `real(dp)`; declared required yes; input required no; default `1.0` (object default)
- `baseflow_recession_coefficient%max`: `real(dp)`; declared required yes; input required no; default `1000.0` (object default)

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
  b%value = 1.0
  b%optimize = .true.
  b%min = 0.1
  b%max = 10.0
  ks_lakebed%value = 50.0
  ks_lakebed%optimize = .true.
  ks_lakebed%min = 1.1
  ks_lakebed%max = 1000.0
  baseflow_recession_coefficient%value = 100.0
  baseflow_recession_coefficient%optimize = .true.
  baseflow_recession_coefficient%min = 1.0
  baseflow_recession_coefficient%max = 1000.0
/
```

