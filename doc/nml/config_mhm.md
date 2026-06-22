# mHM model configuration {#config_mhm}

[TOC]

Configuration for the mHM model setup.

**Namelist**: `config_mhm`

## Fields

| Name | Type | Required | Info |
| --- | --- | --- | --- |
| [output_path](#output_path) | string array | no | Output path |
| [read_restart](#read_restart) | logical array | no | Read restart |
| [restart_input_path](#restart_input_path) | string array | no | Restart input path |
| [write_restart](#write_restart) | logical array | no | Write restart |
| [restart_output_path](#restart_output_path) | string array | no | Restart output path |
| [evap_coeff](#evap_coeff) | real array | no | Evaporation coefficients |
| [share_evap_coeff](#share_evap_coeff) | logical | no | Share evaporation coefficients between domains |

## Field details

### output_path

Output path `output_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Examples: `["mhm_output.nc"]`

### read_restart

Read restart `read_restart`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### restart_input_path

Restart input path `restart_input_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Examples: `["mhm_restart_in.nc"]`

### write_restart

Write restart `write_restart`

Summary:
- Type: `logical, dimension(n_domains)`
- Required: no
- Default: `.false.`

### restart_output_path

Restart output path `restart_output_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Examples: `["mhm_restart_out.nc"]`

### evap_coeff

Evaporation coefficients `evap_coeff`

Monthly free pan evaporation coefficients for free-water surfaces.
(array dimension 1: month, dimension 2: domain)

Summary:
- Type: `real(dp), dimension(12, n_domains)`
- Required: no
- Examples: `[1.3, 1.2, 0.72, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.5]`

### share_evap_coeff

Share evaporation coefficients between domains `share_evap_coeff`

Share the same evaporation coefficient values between all domains taking the values from the first domain.
If set to false, different evaporation coefficient values can be set for each domain.

Summary:
- Type: `logical`
- Required: no
- Default: `.true.`

## Example

```fortran
&config_mhm
  output_path(:) = "mhm_output.nc"
  read_restart(:) = .false.
  restart_input_path(:) = "mhm_restart_in.nc"
  write_restart(:) = .false.
  restart_output_path(:) = "mhm_restart_out.nc"
  evap_coeff(:,1) = 1.3, 1.2, 0.72, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.5
  share_evap_coeff = .true.
/
```

