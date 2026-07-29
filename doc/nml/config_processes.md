# Processes configuration {#config_processes}

[TOC]

Configuration for process case selection in mHM.

**Namelist**: `config_processes`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [interception](#interception) | integer | no | no | Interception process case |
| [snow](#snow) | integer | no | no | Snow process case |
| [soil_moisture](#soil_moisture) | integer | no | no | Soil moisture process case |
| [direct_runoff](#direct_runoff) | integer | no | no | Direct runoff process case |
| [pet](#pet) | integer | no | no | Potential evapotranspiration (PET) process case |
| [interflow](#interflow) | integer | no | no | Interflow process case |
| [percolation](#percolation) | integer | no | no | Percolation process case |
| [baseflow](#baseflow) | integer | no | no | Baseflow process case |
| [neutrons](#neutrons) | integer | no | no | Ground albedo of cosmic-ray neutrons process case |
| [routing](#routing) | integer | no | no | Routing process case |
| [temperature_routing](#temperature_routing) | integer | no | no | River temperature routing process case |

## Field details

### interception

Interception process case `interception`

Process case for interception.
- -1: pass through precipitation as throughfall
- 0: deactivated
- 1: maximum Interception

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `0`
- Allowed values: `-1`, `0`, `1`
- Examples: `1`

### snow

Snow process case `snow`

Process case for snow.
- -1: pass through throughfall as effective precipitation
- 0: deactivated
- 1: degree-day approach

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `0`
- Allowed values: `-1`, `0`, `1`
- Examples: `1`

### soil_moisture

Soil moisture process case `soil_moisture`

Process case for soil moisture.
- 0: deactivated
- 1: Feddes equation for ET reduction, multi-layer infiltration capacity approach, Brooks-Corey like
- 2: Jarvis equation for ET reduction, multi-layer infiltration capacity approach, Brooks-Corey like
- 3: Jarvis equation for ET reduction and global FC dependency on root fraction coefficient
- 4: Feddes equation for ET reduction and global FC dependency on root fraction coefficient

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `0`
- Allowed values: `0`, `1`, `2`, `3`, `4`
- Examples: `1`

### direct_runoff

Direct runoff process case `direct_runoff`

Process case for direct runoff.
- 0: deactivated
- 1: linear reservoir exceedance approach

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `0`
- Allowed values: `0`, `1`
- Examples: `1`

### pet

Potential evapotranspiration (PET) process case `pet`

Process case for potential evapotranspiration (PET).
- -2: PET is input, aspect driven correction
- -1: PET is input, LAI driven correction
- 0: deactivated
- 1: Hargreaves-Samani method
- 2: Priestley-Taylor method
- 3: Penman-Monteith method

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `0`
- Allowed values: `-2`, `-1`, `0`, `1`, `2`, `3`
- Examples: `-2`

### interflow

Interflow process case `interflow`

Process case for interflow.
- 0: deactivated
- 1: storage reservoir with one outflow threshold and nonlinear response

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `0`
- Allowed values: `0`, `1`
- Examples: `1`

### percolation

Percolation process case `percolation`

Process case for percolation.
- 0: deactivated
- 1: GW assumed as linear reservoir

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `0`
- Allowed values: `0`, `1`
- Examples: `1`

### baseflow

Baseflow process case `baseflow`

Process case for baseflow.
- 0: deactivated
- 1: recession parameters (not regionalized yet)

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `0`
- Allowed values: `0`, `1`
- Examples: `1`

### neutrons

Ground albedo of cosmic-ray neutrons process case `neutrons`

Ground albedo of cosmic-ray neutrons process case.
Work in progress, do not use for research.
- 0: deactivated
- 1: inverse N0 based on Desilets et al. 2010
- 2: COSMIC forward operator by Shuttleworth et al. 2013

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `0`
- Allowed values: `0`, `1`, `2`
- Examples: `0`

### routing

Routing process case `routing`

Process case for routing.
- 0: deactivated
- 1: Muskingum approach
- 2: adaptive timestep, constant celerity
- 3: adaptive timestep, spatially varying celerity

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `0`
- Allowed values: `0`, `1`, `2`, `3`
- Examples: `3`

### temperature_routing

River temperature routing process case `temperature_routing`

River temperature routing process case (needs routing).
- 0: deactivated
- 1: following: Beek et al., 2012

Summary:
- Type: `integer(i4)`
- Declared required: no
- Input required: no
- Default: `0`
- Allowed values: `0`, `1`
- Examples: `0`

## Example

```fortran
&config_processes
  interception = 1
  snow = 1
  soil_moisture = 1
  direct_runoff = 1
  pet = -2
  interflow = 1
  percolation = 1
  baseflow = 1
  neutrons = 0
  routing = 3
  temperature_routing = 0
/
```

