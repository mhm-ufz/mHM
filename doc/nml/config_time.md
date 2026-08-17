# Time configuration {#config_time}

[TOC]

Configuration for simulation and evaluation time periods in mHM.

**Namelist**: `config_time`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [sim_start](#sim_start) | string array | no | no | Simulation start |
| [eval_start](#eval_start) | string array | no | no | Evaluation start |
| [sim_end](#sim_end) | string array | no | no | Simulation end |
| [share_time_period](#share_time_period) | logical | no | no | Share time period between domains |
| [time_step](#time_step) | integer array | no | no | Global model time step |
| [share_time_step](#share_time_step) | logical | no | no | Share time step between domains |

## Field details

### sim_start

Simulation start `sim_start`

Start date of the simulation. Only hours are supported; minutes and seconds are set to zero.
Format: YYYY-MM-DD [hh[:mm[:ss]]]

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `date-time`
- Declared required: no
- Input required: no
- Examples: `"2020-06-01 00:00"`, `"2025-01-01"`

### eval_start

Evaluation start `eval_start`

Start date of the evaluation period. Only hours are supported; minutes and seconds are set to zero.
Format: YYYY-MM-DD [hh[:mm[:ss]]]

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `date-time`
- Declared required: no
- Input required: no
- Examples: `"2021-01-01 00:00"`, `"2025-02-15"`

### sim_end

Simulation end `sim_end`

End date of the simulation. Only hours are supported; minutes and seconds are set to zero.
Format: YYYY-MM-DD [hh[:mm[:ss]]]

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `date-time`
- Declared required: no
- Input required: no
- Examples: `"2022-01-01 00:00"`, `"2026-01-01"`

### share_time_period

Share time period between domains `share_time_period`

Share the same simulation and evaluation time periods between all domains taking the values from the first domain.
If set to false, different time periods can be set for each domain.

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
- Default: `.false.`

### time_step

Global model time step `time_step`

Global model update cadence in hours. The supported values are whole divisors of one day.
This controls how often process containers are updated and is independent of the support
interval of input data. Individual process components may support only a subset of these
globally valid cadences.

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `1`
- Allowed values: `1`, `2`, `3`, `4`, `6`, `8`, `12`, `24`

### share_time_step

Share time step between domains `share_time_step`

Share the same time step value between all domains taking the value from the first domain.
If set to false, different time step values can be set for each domain.

Summary:
- Type: `logical`
- Declared required: no
- Input required: no
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

