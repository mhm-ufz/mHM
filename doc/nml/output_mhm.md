# mHM output configuration {#output_mhm}

[TOC]

Output configuration for mHM.

**Namelist**: `output_mhm`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [output_deflate_level](#output_deflate_level) | integer | no | no | Output deflate level |
| [output_double_precision](#output_double_precision) | logical | no | no | Output double precision |
| [output_time_reference](#output_time_reference) | integer | no | no | Output time reference |
| [output_frequency](#output_frequency) | integer | no | no | Output time step |
| [out_interception](#out_interception) | logical | no | no | Interception |
| [out_snowpack](#out_snowpack) | logical | no | no | Snowpack |
| [out_SWC](#out_swc) | logical | no | no | Layered Soil Water Content |
| [out_SM](#out_sm) | logical | no | no | Layered Volumetric Soil Moisture |
| [out_SM_all](#out_sm_all) | logical | no | no | Mean Volumetric Soil Moisture |
| [out_sealedSTW](#out_sealedstw) | logical | no | no | Reservoir of Sealed areas |
| [out_unsatSTW](#out_unsatstw) | logical | no | no | Reservoir of Unsaturated areas |
| [out_satSTW](#out_satstw) | logical | no | no | Reservoir of Saturated areas |
| [out_PET](#out_pet) | logical | no | no | Potential Evapotranspiration |
| [out_aET_all](#out_aet_all) | logical | no | no | Mean actual Evapotranspiration |
| [out_Q](#out_q) | logical | no | no | Total Discharge |
| [out_QD](#out_qd) | logical | no | no | Direct Runoff |
| [out_QIf](#out_qif) | logical | no | no | Fast Interflow |
| [out_QIs](#out_qis) | logical | no | no | Slow Interflow |
| [out_QB](#out_qb) | logical | no | no | Baseflow |
| [out_recharge](#out_recharge) | logical | no | no | Groundwater Recharge |
| [out_soil_infil](#out_soil_infil) | logical | no | no | Soil Infiltration |
| [out_neutrons](#out_neutrons) | logical | no | no | Neutrons |
| [out_aET_layer](#out_aet_layer) | logical | no | no | Actual Evapotranspiration from Soil Layers |
| [out_preEffect](#out_preeffect) | logical | no | no | Effective Precipitation |
| [out_Qsm](#out_qsm) | logical | no | no | Snow Melt |

## Field details

### output_deflate_level

Output deflate level `output_deflate_level`

Compression level for netCDF4 outputs (0: no compression, 9: maximum compression).

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `6`
- Minimum: `>= 0`
- Maximum: `<= 9`

### output_double_precision

Output double precision `output_double_precision`

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### output_time_reference

Output time reference `output_time_reference`

Location of reference time point in outputs:
- if 0 : start of the time interval (i.e. 1990-01-01 00:00 for daily output on 1990-01-01)
- if 1 : center of the time interval (i.e. 1990-01-01 12:00 for daily output on 1990-01-01)
- if 2 : end of the time interval (i.e. 1990-01-02 00:00 for daily output on 1990-01-01)

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `2`
- Allowed values: `0`, `1`, `2`

### output_frequency

Output time step `output_frequency`

switch to control write out frequency of the gridded model outputs below
- if >0 : fixed output interval in hours; it must be a whole multiple of the global model time step
- if 0 : only at end of run
- if -1 : daily
- if -2 : monthly
- if -3 : yearly
Hydrological amounts are written with units `mm`. For accumulated outputs,
the represented interval is described by the output time bounds rather than
encoded in the units.

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `-2`
- Minimum: `>= -3`

### out_interception

Interception `out_interception`

Canopy interception storage (L1_inter) [mm]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_snowpack

Snowpack `out_snowpack`

Height of snowpack (L1_snowpack) [mm]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_SWC

Layered Soil Water Content `out_SWC`

soil water content in the single layers (L1_soilMoist)

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_SM

Layered Volumetric Soil Moisture `out_SM`

volumetric soil moisture in the single layers (L1_soilMoist / L1_soilMoistSat) [mm/mm]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_SM_all

Mean Volumetric Soil Moisture `out_SM_all`

mean volumetric soil moisture averaged over all soil layers (L1_soilMoist / L1_soilMoistSat) [mm/mm]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_sealedSTW

Reservoir of Sealed areas `out_sealedSTW`

waterdepth in reservoir of sealed areas (L1_sealSTW) [mm]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_unsatSTW

Reservoir of Unsaturated areas `out_unsatSTW`

waterdepth in reservoir of unsat. soil zone (L1_unsatSTW) [mm]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_satSTW

Reservoir of Saturated areas `out_satSTW`

Water depth in reservoir of sat. soil zone (L1_satSTW) [mm]
--> level of GW reservoir

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_PET

Potential Evapotranspiration `out_PET`

Potential evapotranspiration PET [mm/T]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_aET_all

Mean actual Evapotranspiration `out_aET_all`

Actual evapotranspiration aET [mm/T]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_Q

Total Discharge `out_Q`

Total discharge generated per cell (L1_total_runoff) [mm/T]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_QD

Direct Runoff `out_QD`

Direct runoff generated per cell (L1_runoffSeal) [mm/T]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_QIf

Fast Interflow `out_QIf`

Fast interflow generated per cell (L1_fastRunoff) [mm/T]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_QIs

Slow Interflow `out_QIs`

Slow interflow generated per cell (L1_slowRunoff) [mm/T]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_QB

Baseflow `out_QB`

Baseflow generated per cell (L1_baseflow) [mm/T]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_recharge

Groundwater Recharge `out_recharge`

Groundwater Recharge (L1_percol) [mm/T]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_soil_infil

Soil Infiltration `out_soil_infil`

Infiltration (L1_infilSoil) [mm/T]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_neutrons

Neutrons `out_neutrons`

Ground albedo Neutrons related to soil moisture
THIS IS WORK IN PROGRESS, DO NOT USE FOR RESEARCH

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_aET_layer

Actual Evapotranspiration from Soil Layers `out_aET_layer`

Actual evapotranspiration from the soil layers [mm/T]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_preEffect

Effective Precipitation `out_preEffect`

Effective Precipitation (L1_preEffect) [mm/T]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_Qsm

Snow Melt `out_Qsm`

Snow melt (L1_melt) [mm/T]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

## Example

```fortran
&output_mhm
  output_deflate_level = 6
  output_double_precision = .false.
  output_time_reference = 2
  output_frequency = -2
  out_interception = .false.
  out_snowpack = .false.
  out_SWC = .false.
  out_SM = .false.
  out_SM_all = .false.
  out_sealedSTW = .false.
  out_unsatSTW = .false.
  out_satSTW = .false.
  out_PET = .false.
  out_aET_all = .false.
  out_Q = .false.
  out_QD = .false.
  out_QIf = .false.
  out_QIs = .false.
  out_QB = .false.
  out_recharge = .false.
  out_soil_infil = .false.
  out_neutrons = .false.
  out_aET_layer = .false.
  out_preEffect = .false.
  out_Qsm = .false.
/
```

