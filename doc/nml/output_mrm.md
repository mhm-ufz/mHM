# mRM output configuration {#output_mrm}

[TOC]

Output configuration for mRM.

**Namelist**: `output_mrm`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [output_deflate_level](#output_deflate_level) | integer | no | no | Output deflate level |
| [output_double_precision](#output_double_precision) | logical | no | no | Output double precision |
| [output_time_reference](#output_time_reference) | integer | no | no | Output time reference |
| [output_frequency](#output_frequency) | integer | no | no | Output time step |
| [out_Qrouted](#out_qrouted) | logical | no | no | Routed Streamflow |
| [out_RivTemp](#out_rivtemp) | logical | no | no | Routed Temperature |

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
- if >0 : after every 'n' time steps
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

### out_Qrouted

Routed Streamflow `out_Qrouted`

Routed streamflow (Qrouted in output flux) (L11_qMod) [m3 s-1]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`
- Examples: `.true.`

### out_RivTemp

Routed Temperature `out_RivTemp`

Routed temperature (only if River Temperature is enabled) [K]

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`
- Examples: `.true.`

## Example

```fortran
&output_mrm
  output_deflate_level = 6
  output_double_precision = .false.
  output_time_reference = 2
  output_frequency = -1
  out_Qrouted = .true.
  out_RivTemp = .true.
/
```

