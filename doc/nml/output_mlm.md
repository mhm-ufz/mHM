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
  out_lake_outflow = .false.
/
```

