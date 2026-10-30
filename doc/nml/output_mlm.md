# mLM output configuration {#output_mlm}

[TOC]

Output configuration for mLM lake points.

**Namelist**: `output_mlm`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [output_deflate_level](#output_deflate_level) | integer | no | no | Output deflate level |
| [output_double_precision](#output_double_precision) | logical | no | no | Output double precision |
| [output_time_reference](#output_time_reference) | integer | no | no | Output time reference |
| [output_frequency](#output_frequency) | integer | no | no | Output time step |
| [output_timeseries](#output_timeseries) | logical | no | no | Use station time-series layout for lake-point output |
| [out_lake_level](#out_lake_level) | logical | no | no | Lake level |
| [out_lake_surface_area](#out_lake_surface_area) | logical | no | no | Lake surface area |
| [out_lake_volume](#out_lake_volume) | logical | no | no | Lake volume |
| [out_lake_gw_storage](#out_lake_gw_storage) | logical | no | no | Lake groundwater storage |
| [out_lake_inflow](#out_lake_inflow) | logical | no | no | Lake inflow |
| [out_lake_precipitation](#out_lake_precipitation) | logical | no | no | Lake precipitation |
| [out_lake_evaporation](#out_lake_evaporation) | logical | no | no | Lake evaporation |
| [out_lake_seepage](#out_lake_seepage) | logical | no | no | Lake seepage |
| [out_lake_baseflow](#out_lake_baseflow) | logical | no | no | Lake baseflow |
| [out_lake_spill](#out_lake_spill) | logical | no | no | Lake spill |
| [out_lake_outflow](#out_lake_outflow) | logical | no | no | Lake outflow |

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
- if 0 : start of the time interval
- if 1 : center of the time interval
- if 2 : end of the time interval

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `2`
- Allowed values: `0`, `1`, `2`

### output_frequency

Output time step `output_frequency`

Switch controlling the lake-point output frequency:
- if >0 : fixed output interval in hours; it must be a whole multiple of the one-hour mLM step
- if 0 : only at end of run
- if -1 : daily
- if -2 : monthly
- if -3 : yearly

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `-1`
- Minimum: `>= -3`

### output_timeseries

Use station time-series layout for lake-point output `output_timeseries`

Use the mo_points_io time-series layout for temporal lake-point variables.
- if .false. : (time, lake) dimension order (classic layout; read directly by CDO).
- if .true.  : (lake, time) CF timeSeries featureType layout.

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_lake_level

Lake level `out_lake_level`

Lake surface level [m a.s.l.].

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_lake_surface_area

Lake surface area `out_lake_surface_area`

Lake surface area [m2].

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_lake_volume

Lake volume `out_lake_volume`

Lake storage volume [m3].

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_lake_gw_storage

Lake groundwater storage `out_lake_gw_storage`

Lake groundwater storage [m3].

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_lake_inflow

Lake inflow `out_lake_inflow`

Lake inflow from mRM [m3 s-1].

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_lake_precipitation

Lake precipitation `out_lake_precipitation`

Precipitation onto the lake surface [mm].

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_lake_evaporation

Lake evaporation `out_lake_evaporation`

Open-water evaporation from the lake surface [mm].

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_lake_seepage

Lake seepage `out_lake_seepage`

Seepage (percolation) from the lake to groundwater storage [mm].

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_lake_baseflow

Lake baseflow `out_lake_baseflow`

Baseflow released from the lake groundwater storage [m3 s-1].

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_lake_spill

Lake spill `out_lake_spill`

Crest spill from the lake [m3 s-1]. Lake outflow = spill + baseflow.

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### out_lake_outflow

Lake outflow `out_lake_outflow`

Hourly mLM release supplied to mRM [m3 s-1].

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

## Example

```fortran
&output_mlm
  output_deflate_level = 6
  output_double_precision = .false.
  output_time_reference = 2
  output_frequency = -1
  output_timeseries = .false.
  out_lake_level = .false.
  out_lake_surface_area = .false.
  out_lake_volume = .false.
  out_lake_gw_storage = .false.
  out_lake_inflow = .false.
  out_lake_precipitation = .false.
  out_lake_evaporation = .false.
  out_lake_seepage = .false.
  out_lake_baseflow = .false.
  out_lake_spill = .false.
  out_lake_outflow = .false.
/
```

