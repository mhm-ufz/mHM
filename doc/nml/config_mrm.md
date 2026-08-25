# mRM configuration {#config_mrm}

[TOC]

Configuration for the multi-scale routing model (mRM) in mHM.

**Namelist**: `config_mrm`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [river_net_order_root_based](#river_net_order_root_based) | logical array | no | no | Flag for root based river network ordering. |
| [river_net_omp_level_min](#river_net_omp_level_min) | integer array | no | no | Minimum level size for OpenMP parallelization. |
| [max_route_step](#max_route_step) | integer array | no | no | Maximum numerical routing substep in seconds. |
| [upscale_mode](#upscale_mode) | integer array | no | no | River upscaling mode. |
| [length_percentile](#length_percentile) | real array | no | no | Percentile for the minimum upscaled link length. |
| [scc_gauges_path](#scc_gauges_path) | string array | no | no | Path for SCC gauges NetCDF file. |
| [scc_gauges_as_poi](#scc_gauges_as_poi) | logical array | no | no | Use SCC gauges as POIs. |
| [poi_gauges_path](#poi_gauges_path) | string array | no | no | Path for POI gauges NetCDF file. |
| [output_path](#output_path) | string array | no | no | Path for output file. |
| [output_node_path](#output_node_path) | string array | no | no | Path for node based output file. |
| [output_poi_path](#output_poi_path) | string array | no | no | Path for POI based output file. |
| [poi_timeseries](#poi_timeseries) | logical array | no | no | Use station time-series layout for POI output. |
| [node_timeseries](#node_timeseries) | logical array | no | no | Use point time-series layout for node output. |
| [read_restart](#read_restart) | logical array | no | no | Read restart |
| [read_restart_fluxes](#read_restart_fluxes) | logical array | no | no | Read restart fluxes |
| [restart_input_path](#restart_input_path) | string array | no | no | Restart input path |
| [write_restart](#write_restart) | logical array | no | no | Write restart |
| [restart_output_path](#restart_output_path) | string array | no | no | Restart output path |
| [diagnostics_path](#diagnostics_path) | string array | no | no | Diagnostics output path |

## Field details

### river_net_order_root_based

Flag for root based river network ordering. `river_net_order_root_based`

Flag to indicate if the river network is ordered in root based order.
If false, the ordering is leaf based.
Root based ordering results in more equal distributed level sizes for parallelization.
Leaf based ordering has huge levels of nodes at the headwaters, which can lead to load balancing issues.

Summary:
- Type: `logical, dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `.false.`
- Examples: `[.true.]`

### river_net_omp_level_min

Minimum level size for OpenMP parallelization. `river_net_omp_level_min`

Minimum level size in the river network where OpenMP parallelization starts.
Levels smaller than this size are always computed serially.
This can be used to avoid parallelization overhead for small levels of the river network,
especially when ordering is leaf based.

Special values are:
- -1 : lets the model choose a default based on the number of threads.
-  1 : forces parallelization on all levels.
-  0 : disables parallelization.
By default: threads * 8 (indicated by -1)

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `-1`
- Minimum: `>= -1`
- Examples: `[100]`

### max_route_step

Maximum numerical routing substep in seconds. `max_route_step`

Maximum allowed numerical routing substep in seconds.
This parameter can be used to limit the substep in case of very large CFL time steps
due to very low flow velocities.
This is useful in coupling scenarios to match the time step of other models.
If the routing substep is smaller than the model time step, multiple routing iterations are performed per model update.
If it is larger, the completed routing step is enlarged to a compatible whole-hour interval.
Valid values range from 1 minute (60s) to 1 day (86400s).
Value needs to be a divisor of 3600 or a multiple of 3600 and a divisor of 86400.

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `86400`
- Allowed values: `60`, `120`, `180`, `240`, `300`, `360`, `600`, `720`, `900`, `1200`, `1800`, `3600`, `7200`, `10800`, `14400`, `21600`, `28800`, `43200`, `86400`
- Examples: `[3600]`

### upscale_mode

River upscaling mode. `upscale_mode`

Method used to construct the coarse river graph from the fine D8 network.
- 0 : legacy topology connecting to the first entered neighboring coarse cell.
- 1 : FLOW-like topology connecting directly to the next downstream link endpoint.

Summary:
- Type: `integer(i4), dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `1`
- Allowed values: `0`, `1`
- Examples: `[1]`

### length_percentile

Percentile for the minimum upscaled link length. `length_percentile`

Percentile of the non-sink link-length distribution used as a lower bound
for routing stability, in percent.

Summary:
- Type: `real(dp), dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `40.0`
- Minimum: `>= 0.0`
- Maximum: `<= 100.0`
- Examples: `[40.0]`

### scc_gauges_path

Path for SCC gauges NetCDF file. `scc_gauges_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no
- Examples: `["scc_gauges.nc"]`

### scc_gauges_as_poi

Use SCC gauges as POIs. `scc_gauges_as_poi`

Use the exact coarse river nodes retained for SCC gauges as point-of-interest output locations.

Summary:
- Type: `logical, dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `.false.`

### poi_gauges_path

Path for POI gauges NetCDF file. `poi_gauges_path`

Point locations used to select the nearest mRM river nodes.

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no
- Examples: `["poi_gauges.nc"]`

### output_path

Path for output file. `output_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no
- Examples: `["mrm_output.nc"]`

### output_node_path

Path for node based output file. `output_node_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no
- Examples: `["mrm_node_output.nc"]`

### output_poi_path

Path for POI based output file. `output_poi_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no
- Examples: `["mrm_poi_output.nc"]`

### poi_timeseries

Use station time-series layout for POI output. `poi_timeseries`

Use the mo_points_io time-series layout for temporal POI variables.

Summary:
- Type: `logical, dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `.true.`

### node_timeseries

Use point time-series layout for node output. `node_timeseries`

Use the mo_points_io time-series layout for temporal node variables.

Summary:
- Type: `logical, dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `.false.`

### read_restart

Read restart `read_restart`

Summary:
- Type: `logical, dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `.false.`

### read_restart_fluxes

Read restart fluxes `read_restart_fluxes`

Summary:
- Type: `logical, dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `.true.`

### restart_input_path

Restart input path `restart_input_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no
- Examples: `["mrm_restart_in.nc"]`

### write_restart

Write restart `write_restart`

Summary:
- Type: `logical, dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `.false.`

### restart_output_path

Restart output path `restart_output_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no
- Examples: `["mrm_restart_out.nc"]`

### diagnostics_path

Diagnostics output path `diagnostics_path`

Path for diagnostics output file containing information about the river upscaling and SCC.

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no
- Examples: `["mrm_diagnostics.nc"]`

## Example

```fortran
&config_mrm
  river_net_order_root_based(:) = .true.
  river_net_omp_level_min(:) = 100
  max_route_step(:) = 3600
  upscale_mode(:) = 1
  length_percentile(:) = 40.0
  scc_gauges_path(:) = "scc_gauges.nc"
  scc_gauges_as_poi(:) = .false.
  poi_gauges_path(:) = "poi_gauges.nc"
  output_path(:) = "mrm_output.nc"
  output_node_path(:) = "mrm_node_output.nc"
  output_poi_path(:) = "mrm_poi_output.nc"
  poi_timeseries(:) = .true.
  node_timeseries(:) = .false.
  read_restart(:) = .false.
  read_restart_fluxes(:) = .true.
  restart_input_path(:) = "mrm_restart_in.nc"
  write_restart(:) = .false.
  restart_output_path(:) = "mrm_restart_out.nc"
  diagnostics_path(:) = "mrm_diagnostics.nc"
/
```

