# Project configuration {#config_project}

[TOC]

Configuration for the overall project setup in mHM.

**Namelist**: `config_project`

## Fields

| Name | Type | Required | Info |
| --- | --- | --- | --- |
| [project_details](#project_details) | string | no | Project name |
| [setup_description](#setup_description) | string | no | Description of the setup |
| [simulation_type](#simulation_type) | string | no | Type of simulation |
| [Conventions](#conventions) | string | no | Convention used for dataset |
| [contact](#contact) | string | no | Contact details, incl. PI name, modellers |
| [mHM_details](#mhm_details) | string | no | Developing institution |
| [history](#history) | string | no | Some details on data/model run version. |
| [n_domains](#n_domains) | integer | no | Number of domains |
| [n_geo_units](#n_geo_units) | integer | no | Number of geological units |
| [max_layers](#max_layers) | integer | no | Maximum number of soil layers |
| [read_domains_from_dirs](#read_domains_from_dirs) | logical | no | Flag for separate domains |

## Field details

### project_details

Project name `project_details`

Summary:
- Type: `character(len=buf)`
- Required: no
- Default: `"mHM project"`

### setup_description

Description of the setup `setup_description`

Summary:
- Type: `character(len=buf)`
- Required: no
- Default: `"Model run"`

### simulation_type

Type of simulation `simulation_type`

e.g. hindcast simulation, seasonal forecast, climate projection

Summary:
- Type: `character(len=buf)`
- Required: no
- Default: `"Simulation"`

### Conventions

Convention used for dataset `Conventions`

Summary:
- Type: `character(len=buf)`
- Required: no
- Default: `"None"`

### contact

Contact details, incl. PI name, modellers `contact`

Summary:
- Type: `character(len=buf)`
- Required: no
- Default: `"Developer"`

### mHM_details

Developing institution `mHM_details`

Summary:
- Type: `character(len=buf)`
- Required: no
- Default: `"Research unit"`

### history

Some details on data/model run version. `history`

Summary:
- Type: `character(len=buf)`
- Required: no
- Default: `"Model run version 1"`

### n_domains

Number of domains `n_domains`

Number of domains to be simulated.

Summary:
- Type: `integer(i4)`
- Required: no
- Default: `1`

### n_geo_units

Number of geological units `n_geo_units`

Number of geological units in the global parameter set.

Summary:
- Type: `integer(i4)`
- Required: no
- Default: `25`

### max_layers

Maximum number of soil layers `max_layers`

Maximum number of soil-layer entries stored per domain in this namelist file.

Summary:
- Type: `integer(i4)`
- Required: no
- Default: `10`
- Minimum: `>= 1`

### read_domains_from_dirs

Flag for separate domains `read_domains_from_dirs`

Read domain configurations from namelists in separate directories.
Only "config_project", "config_domain", "config_processes" and "config_optimize" are read from this file then.

Summary:
- Type: `logical`
- Required: no
- Default: `.false.`

## Example

```fortran
&config_project
  project_details = "mHM project"
  setup_description = "Model run"
  simulation_type = "Simulation"
  Conventions = "None"
  contact = "Developer"
  mHM_details = "Research unit"
  history = "Model run version 1"
  n_domains = 1
  n_geo_units = 25
  max_layers = 10
  read_domains_from_dirs = .false.
/
```

