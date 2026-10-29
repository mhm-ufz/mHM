# Developing mHM v6 {#v6_development}

The v6 structure makes components usable across independent domains, repeated
parameter runs, and restart workflows. Clear ownership and explicit inputs let
you extend one process while keeping those workflows predictable. This guide
follows data through configuration, connection, stepping, and restart so you can
add a process or container. `AGENTS.md` contains detailed coding rules and local
validation commands. For build and run instructions, see
[Getting started](@ref v6_getting_started).

[TOC]

## Follow the domain

Start with `app/driver.F90` and `src/domain/mo_domain.F90`. Each `domain_t` owns
its exchange and component state. The domain orchestrates components in this
order: **Input → MPR → meteo → mHM → mLM → mRM**.

Orchestration and exchange contracts live in `src/domain/`; file-input handling
lives in `src/input/`; process code lives in its component folder. Keep shared,
process-neutral utilities in `src/common/`. Scientific kernels receive their
inputs explicitly, with mutable state owned by a container. Group related state
and caches in small derived types so their purpose and lifetime remain clear.

The exchange holds the shared clock, configuration, parameter registry, spatial
support, and field contracts. It makes data available between components; each
scientific calculation belongs to its process component.

| Component | What it owns and does |
| --- | --- |
| [Input](@ref f_input) | Opens configured datasets, establishes native grids and point/river support, publishes static fields, and selects or reads forcing chunks each step. |
| [MPR](@ref f_mpr) | Owns LUTs, the soil database, scaling and temporal caches. `connect` prepares static inputs; `initialize` regionalizes current parameters; `update` selects active land-cover/LAI slices. |
| [Meteo](@ref f_meteo) | Remaps raw forcing with count weights, calculates/corrects PET, and disaggregates to the model interval. Land forcing and lake-point forcing use their respective support. Model steps are 1 or 24 hours. |
| [mHM](@ref f_mhm) | Owns selected land hydrology states and fluxes and runs their kernels in process order. Disabled subprocesses may receive externally provided fields through the exchange. |
| [mLM](@ref f_mLM) | Owns lake definition, natural-lake state and delayed outflow. Cases implement pass-through, precipitation/PET balance, or bathymetry, spill, seepage, groundwater and open-water evaporation. |
| [mRM](@ref f_mrm) | Owns routing topology, hydraulic/thermal state, routing accumulators, discharge, and lake inflow. It completes routing intervals within the shared domain clock. Routing-only runs also support intermediate model cadences. |

For a compact component example, read [mo_mlm_container](@ref mo_mlm_container).
Use the relevant component's kernels as the scientific reference when adding a
process case.

## Give each lifecycle phase one job

| Phase | Responsibility |
| --- | --- |
| `configure` | Read and validate configuration, select active processes, declare provided fields, and register parameters. |
| `prepare_restart` | Restore model definition needed by consumers, such as grids and horizon bounds. |
| `connect` | Resolve ownership and dependencies; prepare parameter-independent inputs and caches. |
| `initialize` | Fetch current runtime parameters, rebuild derived fields and numerical scratch, and initialize or restore mutable state. |
| `update` | Advance the current model interval using prepared storage. |
| `finalize` | Finish I/O. |
| `destroy` | Release owned storage after the domain clears exchange publications. |

The domain finishes each phase for all selected components before starting the
next phase. `create` precedes these phases and wires component exchange pointers
and namelist dimensions. Consumers may therefore see a declared provider before
its storage or values are available.

Process selections come from `exchange%config%processes`. Each component sets
its `active` flag in `configure`; the domain gates subsequent calls. Keep
parameter input configuration separate from the runtime registry
`exchange%parameters`, which is sealed after configuration. Repeated
`initialize` calls must apply changed parameters without rereading static inputs
or invalidating published pointers. Geometry, LUTs, and connection maps live for
the connected domain; parameter-dependent caches and state are rebuilt/reset for
each initialization; forcing slices and fluxes represent the current interval.
Changing model structure requires a new setup through the full lifecycle.

Keep the domain instance a valid Fortran `target` throughout its lifetime, as in
`mo_mhm_runner`. A connected container contains pointers and cannot be copied as
an independent domain. Finalization leaves publications valid for I/O. Destruction
first clears the exchange, then releases components in reverse order:
mRM, mLM, mHM, meteo, MPR, Input.

The domain advances one shared clock. Updates cover
`[exchange%time_step_start, exchange%time]`: use interval start for selections
such as land cover and LAI, and endpoint time for input lookup and output
timestamps. Prepare numerical work arrays before stepping and keep thread-called
kernels reentrant.

## Separate geometry from packed values

A `grid_t` defines geometry, coordinates, orientation, and an active-cell mask.
A gridded exchange field stores values only for those active cells. Its cell
axis is **flat in Fortran array order**: traverse `(x, y)` with the first index
varying fastest, skipping inactive cells. `grid%ncells` is the packed length;
`grid%cell_ij(k, :)` identifies the original grid position of packed cell `k`.
For a mask keeping `(1,1)`, `(3,1)`, and `(2,2)`, that is also their packed order.

Use the existing grid/I/O helpers to preserve this mapping:

```fortran
! field_xy has shape (grid%nx, grid%ny); packed has length grid%ncells.
call grid%pack_into(field_xy, packed)
call grid%unpack_into(packed, field_xy)
```

Allocate those buffers before use. Unpacking fills inactive cells with nodata.
Respect the grid's stored axis direction; do not infer a geographic row order
or flatten a masked raster with `reshape`. FORCES input readers already produce
this packed representation when given a packed destination.

| Field support | In-memory shape and ordering |
| --- | --- |
| Grid | `(ncells)` in the selected grid's mask order. |
| Grid plus horizons | `(ncells, n_horizons)`, with the packed cell axis first. Each column is one horizon. |
| Points | `(n_points)` in the published point-set order. Stable IDs identify entities across files. |
| River | `(n_nodes)` in the published river's node order, which is distinct from raster packing. |

Full, land, and lake grids can have different masks and packed indices. Hydrology
uses `level1_land`; use named selectors such as `l1_land` and `model_horizons` in
field metadata. Reuse or validate existing geometry; the component that first
establishes it retains ownership. Changing support requires an explicit mapping
or scaler. A row `data(k, :)` of a layered field is strided; pass it to an
assumed-shape kernel without requiring contiguity or copying it per thread.

## Establish an exchange contract

The following excerpts follow the existing mHM runoff handoff to mRM. They
belong to the declaration and lifecycle scopes named below.

### Define metadata and the provider

Add the typed field to `exchange_t`, initialize its canonical metadata in
`exchange_create`, and clear it in `exchange_destroy`:

```fortran
! exchange_t declaration
type(var_dp) :: runoff_total
! exchange_create
self%runoff_total = var_dp(grid=l1_land, name="Q", units="mm", &
  long_name="total runoff", standard_name="runoff_amount")
! exchange_destroy
call self%runoff_total%destroy()
```

Select exactly one primary support: grid, points, or river. A layered field uses
a `var2d_*` type and an explicit layer selector. Declare units, meaning, and
whether it is static. Share only fields needed by other components; keep local
scratch and private states in the owning container.

When mHM owns total runoff, its `configure` declares:

```fortran
call self%exchange%runoff_total%provide("mHM")
```

`provide` records ownership; it does not allocate data. A second declaration is
an error, even from the same provider. Resolve alternatives from selected
processes during configuration: an internal producer and a file/coupled producer
must not both claim a field. Do not redeclare ownership during initialization.

### Bind stable storage

For container-owned data, allocate the packed target during `connect`, prepare
its timing metadata, and publish a pointer:

```fortran
! After allocating self%runoff%total_runoff on level1_land:
call self%exchange%runoff_total%prepare_data("mHM", self%exchange%step_hours)
self%exchange%runoff_total%data => self%runoff%total_runoff
```

The container remains the owner. Updates and repeated initialization assign into
the existing allocation; consumers access it through the exchange. Never publish
a routine-local temporary or free/reallocate a target while it remains published.
Type-bound procedures use `self` with `target`, backed by a persistent caller.

Static file inputs can instead use exchange-owned storage. For example, Input
reads aspect directly into its final packed destination:

```fortran
call self%exchange%aspect%prepare_storage("Input", no_time, [self%aspect%grid%ncells])
call self%aspect%read_static(self%exchange%aspect%data)
```

`prepare_storage` takes `i8` extents and preserves an existing matching allocation;
`set_storage` copies values into prepared storage. An exchange alias uses
`associate_alias` only after the source has data and retains canonical metadata.

### Check requirements when they are needed

A consumer checks that a producer was configured before assuming values exist:

```fortran
! mRM connect: runoff may come from dynamic Input.
call self%exchange%runoff_total%check_provided("mRM")
! Before consuming its bound values:
call self%exchange%check_data(self%exchange%runoff_total, "mRM")
```

`check_provided` checks the declaration. `check_data` also checks association,
temporal metadata, support, and shape derived from that support. It does not
guarantee that a field meets a consumer's particular cadence; validate that
requirement explicitly. Missing required fields fail with `log_fatal` and
`error stop`. Static inputs needed to build maps are checked in `connect`;
dynamic inputs can become available during `initialize` or the first `update`.

This distinction also handles later providers. mLM publishes outflow in
`connect`, then mRM connects and publishes inflow. mLM validates that inflow in
`initialize`. During stepping, mLM consumes the preceding completed hourly inflow;
mRM then uses current outflow and produces the next inflow. Preserve this delay
and the established component order when adding coupling.

Set `stepping` to the interval the field represents: `no_time` (0) for static,
positive hours for fixed intervals, or `daily`, `monthly`, `yearly`, or `varying`
from `mo_grid_io`. Input retains native timing; processed fields publish their
model interval. Temporal caches can switch the published pointer, so read through
`exchange%field%data` when consuming a field; do not keep a stale consumer alias.

## Keep restart responsibilities with the owner

A component's restart stores the definition and model memory needed to continue
its work. Persist physical state, required topology/IDs, timing, and necessary
parameter-derived fields. Reconstruct scratch, scalers, file handles, and pointers
from that definition. Each component writes its own restart.

Restart paths and enable flags are configured explicitly. `prepare_restart`
restores definitions before any consumer connects; `connect` validates existing
support and binds storage; `initialize` restores mutable state into that storage.
MPR also loads its regionalized parameter fields with its restart definition.
mRM makes restart lake inflow available while connecting so earlier mLM
initialization can check it, and restores it again when resetting routing state.

| Component | Restart scope |
| --- | --- |
| Input | Reopens configured forcing sources; owns no restart. |
| Meteo | Optional versioned lake support: stable IDs, ragged sparse L2 indices/weights, full L2 geometry/mask/packing and active lake disaggregation settings. Land processing rebuilds its caches. |
| MPR | Regionalized fields, grid, horizon definitions, and LAI/land-cover timing. A restart-backed run bypasses original LUTs/maps and rejects overrides of the parameters behind those fields. |
| mHM | Active hydrological stores, grid/horizon definition, process cases, and timestamp. Scratch and derived diagnostics are rebuilt. |
| mLM | Lake IDs/points, required static definition (including footprint areas), natural-lake state and delayed outflow, generating parameters and timestamp. Optional L0 topology. Restore per-lake fields by stable ID. |
| mRM | Routing grid/network, hydraulic and thermal memory, public discharge, lake coupling memory, selected output points, and timestamp. |

For lake cases 1 and -2, `config_mlm%restart_save_topology(:)` defaults to false.
It affects restart writing only. Reading detects topology by its presence. mLM
requires its current version-1 format and a complete definition for the saved process
case, including footprint areas and effective generating parameters where applicable.
Unversioned, unsupported or incomplete lake restarts must be regenerated; no legacy
definition reconstruction or parameter fallback is supported.
Explicitly requesting topology without available L0 geometry fails during connection.
`prepare_restart` publishes restart-owned lake IDs and points before meteo connects;
Input retains ownership of any supplied definition.

`config_meteo` has per-domain `read_restart`/`write_restart` flags (default false)
and requires explicit `restart_input_path`/`restart_output_path` when enabled.
Initially these options require lake forcing (case 1 or -2). The version-1 payload
contains no forcing values, file handles, pointers or scratch. Its creation time is
informational: matching support can be reused across simulation periods. Validate
L2 origin, orientation, mask and packing as well as resolution. Reorder by unique
stable lake ID. Validate finite, nonnegative, normalized sparse and temporal weights.
If supplied L0 geometry produces a different mapping on the same L2 grid, reject it.
If L2 changes, rebuild from available L0 geometry; weighted rebuilding needs source
weight files for the current grid. Without compatible support or L0 geometry, fail.
Configure active disaggregation settings consistently with the restart: the effective
`read_meteo_weights` value (including its false default) must match the saved mode,
and active day/night fractions must be configured and match saved values. No
assignment-presence detection or automatic inheritance is needed. Cached temporal
weights remain authoritative while the support is reused, so compatible cached
support does not require reopening source weight files.

Lake regionalization caches static maps, LUT mappings and catchment areas during
connection. Each initialization fetches the registry's conditional lake vector and
borrowed soil/baseflow parameters. Restored Ks, recession and spillway length freeze
their generating values, recorded in every current-format natural-lake restart.
State and forcing targets remain stable during repeated initialization; diagnostics
and output resources are reset.
Lake evaporation kernels and meteorological coefficients live in meteo, while mLM
calls them with explicit current state in water-balance order. Regionalization
coefficients live in MPR constants. Reservoir regulation remains outside this scope.

Validate timestamps, process cases, geometry/masks, layer interfaces, and stable
IDs as applicable. A restart-owned definition must be usable without the original
static inputs; any supplied definition must agree. Restore only states owned by
the selected process, keeping externally supplied fields under their provider.

`finalize` writes requested restarts and closes output before destruction.
mHM writes only at its output boundaries. mRM requires a completed routing
boundary and, when output is active, an output boundary; partial routing
accumulators are not persisted. A new process must identify all memory needed
for an uninterrupted-versus-restarted comparison, including coupling delays.

Use existing grid/point/river writers and canonical exchange metadata for I/O.
Gridded files expand packed arrays back onto their geometry. Layered output opts
into `as_output_var(layered=.true.)`, including one-layer fields, and supplies
layer tops and bounds. Runtime `(ncells, n_horizons)` becomes Fortran `(x, y, z)`
at the I/O boundary; the NetCDF output CDL order is `(time, z, y, x)` (or
`lat, lon` for the horizontal axes). Preserve state averaging and flux
accumulation. Append NetCDF variable definitions with `add_var` from
`mo_netcdf_utils`. See [horizon/output conventions](@ref v6_getting_started) for
the variable-bottom marker and output names.

## Size namelists before reading them

Namelists are generated from `nml-schemas/*.yml` and `nml-config.toml`. Their
runtime dimensions must be known before dependent arrays can be read:

1. `exchange%create` reads `config_project` and establishes project dimensions.
2. `domain%create` calls `set_dims` on the exchange and components.
3. Each component passes the required dimensions to its generated config type.
4. `configure` reads and validates the dependent namelists.

For example, `config_mpr` declares `soil_depth` with schema shape
`[max_layers, n_domains]`. The project supplies `max_layers`, and `mpr_set_dims`
passes `exchange%max_layers` and `exchange%nml_n_domains` to the generated
`set_dims` method. Each domain's `n_layers` then selects how many entries it uses.
The [project reference](@ref config_project) and
[MPR reference](@ref config_mpr) describe these options.

When adding a configurable namelist dimension, declare it in `nml-config.toml`,
expose its value in `config_project`, propagate it through the exchange, and wire
it into the owning component's `set_dims`. Generator defaults alone do not
express the project's required capacity. Counts obtained from datasets, such as
the actual number of lake points, remain data-derived unless needed to size a
namelist array.

Use `exchange%nml_n_domains` and `exchange%nml_domain_id` for domain-indexed
component namelists. With separate domain directories, these are both one;
the top-level project still describes all domains and supplies sufficient
allocation capacity for them.

For a new parameter block, register its schema/file profile in `nml-config.toml`
and its generated config type in `parameter_namelists_t`. The owning component
preserves an API-configured instance or reads its selected block from the
parameter file, validates it, and calls `exchange%parameters%add_process` with
deterministically ordered values and names. Ignore unselected blocks. Fetch
current values with `get_process` during `initialize`; keep parameter-dependent
calculations there so repeated runs use the new values.

## Make file locations explicit

Configure complete paths to component input, output, and restart files. Do not
invent filenames, assume standard subfolders, or search directories for likely
inputs. For example, precipitation can be configured as:

```fortran
&config_input
  pre_path(1) = '../forcing/rainfall_2020.nc'
/
```

Resolve configured paths with `exchange%get_path`. Relative paths use the
domain's configured working directory; `root=.true.` selects the immutable run
root. Absolute paths are supported. Moving a namelist into a subfolder does not
automatically change the base for its data paths. In separate-directory runs,
`config_domain` explicitly selects each domain directory and namelist.

Check required paths when the selected feature needs them, and report missing
active inputs with `log_fatal` and `error stop`. An omitted optional path must
have a documented meaning. The entrypoint's documented defaults for `mhm.nml`,
`parameters.nml`, and `outputs.nml` are separate from component file-path rules.

## Add a feature in a small, complete change

1. Identify the component, state owner, exchanged inputs/outputs, and lifecycle
   phases involved. A process case belongs in its existing component; a new
   container needs its own cohesive state and responsibilities. State any
   intended scientific behavior change. Extend `config_processes` for new
   selections and derive component activity and field ownership from them.
2. Extend the schemas and dimension wiring, then run
   `nml-tools generate --config nml-config.toml`. Include generated Fortran,
   templates, and references together.
3. Implement declarations, dependency checks, initialization, updates, and
   restart handling in their respective phases. Preserve units, calculation order,
   and stable bindings. Declare any new exchange fields, metadata, and destruction.
4. For a new container, add it to `domain_t`, wire its exchange pointer and
   `set_dims`, and dispatch every lifecycle method. Determine its dependency order
   and activity gating explicitly; update the architecture contract if it changes.
5. Register new sources and tests explicitly in `src/CMakeLists.txt` and
   `src/tests/CMakeLists.txt`. Add API comments and the component's Doxygen group
   from `doc/directories.dox`.
6. Verify the behavior with numerical references and useful failure cases. For
   ownership changes, cover repeated initialization and two independent domains;
   for file handling, exercise explicitly configured nonstandard paths. Compare
   uninterrupted and restarted runs, and cover packed masks and layered data where
   used. Use `AGENTS.md` to select the required local and OpenMP checks.
