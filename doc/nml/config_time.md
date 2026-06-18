# Time configuration {#config_time}

[TOC]

Configuration for simulation and evaluation time periods in mHM.

**Namelist**: `config_time`

## Fields

| Name | Type | Required | Info |
| --- | --- | --- | --- |
| [sim_start](#sim_start) | string array | no | Simulation start |
| [eval_start](#eval_start) | string array | no | Evaluation start |
| [sim_end](#sim_end) | string array | no | Simulation end |
| [share_time_period](#share_time_period) | logical | no | Share time period between domains |
| [time_step](#time_step) | integer array | no | Time step of the simulation |
| [share_time_step](#share_time_step) | logical | no | Share time step between domains |

## Field details

### sim_start

Simulation start `sim_start`

Start date of the simulation. Only hours are supported; minutes and seconds are set to zero.
Format: YYYY-MM-DD [hh[:mm[:ss]]]

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Examples: `"2020-06-01 00:00"`, `"2025-01-01"`

### eval_start

Evaluation start `eval_start`

Start date of the evaluation period. Only hours are supported; minutes and seconds are set to zero.
Format: YYYY-MM-DD [hh[:mm[:ss]]]

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Examples: `"2021-01-01 00:00"`, `"2025-02-15"`

### sim_end

Simulation end `sim_end`

End date of the simulation. Only hours are supported; minutes and seconds are set to zero.
Format: YYYY-MM-DD [hh[:mm[:ss]]]

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Required: no
- Examples: `"2022-01-01 00:00"`, `"2026-01-01"`

### share_time_period

Share time period between domains `share_time_period`

Share the same simulation and evaluation time periods between all domains taking the values from the first domain.
If set to false, different time periods can be set for each domain.

Summary:
- Type: `logical`
- Required: no
- Default: `.false.`

### time_step

Time step of the simulation `time_step`

Time step of the simulation (1 or 24, in hours).

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Required: no
- Default: `1`

### share_time_step

Share time step between domains `share_time_step`

Share the same time step value between all domains taking the value from the first domain.
If set to false, different time step values can be set for each domain.

Summary:
- Type: `logical`
- Required: no
- Default: `.true.`

## Example

```fortran
&config_time
  sim_start(:) = "2020-06-01 00:00"
  eval_start(:) = "2021-01-01 00:00"
  sim_end(:) = "2022-01-01 00:00"
  share_time_period = .false.
  time_step(:) = 1
  share_time_step = .true.
/
```

