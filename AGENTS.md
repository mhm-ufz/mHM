# mHM v6 development guide

This is the repository-wide guide for developers and AI agents. Apply it to new or
touched code; keep unrelated cleanup separate. Commands below run from the repository root.

## Architecture and lifecycle

- Before substantial changes, map affected ownership, lifecycle, dependencies, and compatibility
  requirements; state which behavior changes intentionally.
- Use [app/driver.F90](app/driver.F90) as the execution entrypoint. Domain orchestration
  belongs in `src/domain/`; process implementations belong in their component folders.
- Each `domain_t` owns its `exchange_t` and component state; each component references that exchange.
  Pass container-owned state explicitly to scientific kernels; introduce no mutable module-global state.
- Group cohesive process/cache state in nested derived types instead of widening top-level containers.
- Every process container exposes `configure -> prepare_restart -> connect -> initialize -> update -> finalize -> destroy`.
  Preserve this lifecycle and component order: `input -> MPR -> meteo -> mHM -> mLM -> mRM`.
- `exchange%config%processes` is the sole process-selection source. Components set their own `active`
  flag in `configure`; use `any([selections] /= 0_i4)` for aggregate activation. The domain gates later calls.
- `prepare_restart` restores model definition before binding; `initialize` restores mutable restart state afterward.
- Finalizers perform I/O only. Domain destruction calls `exchange%destroy()` first to clear publications,
  release exchange-owned storage, and reset declarations; then destroys mRM, mLM, mHM, meteo, MPR, Input.
- Preserve scientific units, equations, and calculation order during structural refactors.
- When intentionally changing a documented architectural contract, update this guide and its
  regression coverage in the same change.

## Exchange, grids, and ownership

- Each selected provider calls `variable%provide("Component")` exactly once in `configure`;
  duplicate declarations are errors, including declarations from the same component.
- Centralize ownership/source resolution in `connect`. Use the exchange directly for couplable
  subprocess handoffs; avoid a second pointer layer, consumer registries, lifecycle states, or global usage tracking.
- Consumers use `check_provided` for dynamic inputs and `exchange%check_data` when bound data is needed.
  Check dependencies in the phase that needs them; fail clearly with `log_fatal` and `error stop`.
- Metadata declares exactly one primary support (grid, points, or river), optional layers, and canonical
  timing/output metadata. Derive dimensions in `check_data`; do not pass expected shapes to helpers.
- Borrowed producers call `prepare_data(component, stepping)` then bind `data` to a stable target.
  Static file inputs use `prepare_storage` and read directly into final support; DEM/slope work arrays may use `set_storage`.
- Never reallocate published targets while pointers refer to them; reset in place during repeated initialization.
  Use `associate_alias` only after its source has data, refresh aliases after rebinding, and retain alias canonical metadata.
- File inputs own their detected native stepping; processed/model fields publish the interval they represent.
  Stepping is 0 static, -1/-2/-3 daily/monthly/yearly, positive hours fixed, or `varying` irregular.
- Grid publications are non-owning. The first component establishing a grid retains ownership;
  later providers validate it rather than replace it. Use named selector constants, never arithmetic on selectors.
- Full/land/lake selector ranges are 0--99/100--199/200--299. Hydrology uses `level1_land`;
  `level1` aliases it without lakes and remains optional with lakes until explicitly provided.
- `0.0_dp` means unset resolution. Components require positive resolution only when deriving grids;
  validate positive configured resolutions against supplied/restart grids without requiring redundant resolution inputs.
- Externally owned mHM subprocess fields are allowed only when the owning subprocess is disabled.
  `InitState*` constants are initialize/reset-only defaults; use `reset_fields` for both states and flux/diagnostic arrays.

## Paths, time, and parameters

- Resolve configurable paths with `exchange%get_path`. The model interval is
  `[exchange%time_step_start, exchange%time]`; `exchange%step` controls global cadence, with no divergent private time loops.
- Use interval start for LAI, land cover, day/night, and monthly parameters; use endpoint time for
  input lookup, output/restart timestamps, aggregation flushing, diagnostics, and simulation-end checks.
- Positive output frequencies are fixed hours: compatible with mHM's model step and whole multiples
  in either direction of mRM's completed routing step. Only mRM backward-fills finer records with completed routing means.
- Write mHM restarts at its output boundaries; write mRM restarts at completed routing boundaries
  and its output boundaries when output is active. Never persist partial routing accumulators.
- Keep data-only `parameter_t` application-owned in `mo_parameter_types`; schemas reference it with
  `x-fortran-module`. Preserve field/schema order `value, optimize, min, max`; do not duplicate the type in `nml_helper`.
- Keep parameter input instances in `exchange%config%parameters`, separate from the flattened runtime
  registry `exchange%parameters`. For selected blocks: preserve interface-configured instances, otherwise read
  optional `parameter_file`, validate, and register named values deterministically. Ignore unselected blocks; seal after all configuration.
- MPR `connect` caches parameter-independent grids, LUTs, inputs, indices, and temporal source slices;
  every `initialize` fetches current runtime parameters and rebuilds derived fields without rereading static inputs.
- MPR `update` switches active slices only. Keep caches in `mpr_t`, publish active 1D/2D slices, and
  use interval start and half-open periods `[start, end)`; validate timing, dimensions, and bounds before switching pointers.
- MPR initializes its scaler with count weights and cached bounds. Pass that scaler, selected cases,
  and process-specific parameter arrays directly to kernels; do not reconstruct process matrices, padded parameters, or grid-bound copies.

## Component scaling, workspaces, and soil horizons

- Meteo accepts 1-hour and 24-hour model steps; routing-only intermediate cadences remain supported.
  Use FORCES scalers directly with explicit count weights and cached bounds for land forcing and weight maps.
  Keep lake support integration and the forcing/PET/disaggregation order intact.
- PET kernels and coefficients belong in meteo; hydrology/neutron constants in mHM, regionalization
  coefficients in MPR, neutral run defaults/file-format identifiers in common, and exchange contracts in domain.
- Prepare numerical scratch during initialization. Pass strided soil-horizon slices directly to
  assumed-shape kernel arguments; do not require contiguity or copy them into per-thread vectors.
  Neutron kernels share one container-owned, read-only geometry cache and use local scalar accumulators;
  keep timestep calls allocation-free and independent of team size. Rebuild parameter-dependent caches
  on each initialize.
- MPR owns model horizon bounds; exchange publications are non-owning. Derive layer tops from the boundary
  array. Mode 0 retains soil-specific local depths for regionalization and repeats the final layer top as
  its terminal bound, without an L1 depth field. This mHM-specific zero-thickness marker means the bottom
  depth is unspecified, not that the physical layer has zero thickness; never use it to calculate physical
  bottom-layer thickness or volume. Mode 1 uses configured fixed interfaces.
- Derive the fixed/variable bottom flag from `soil_db_mode`; do not maintain terminal-boundary provenance.
  Restore bounds before connecting consumers; validate counts and strictly increasing internal interfaces,
  comparing the terminal boundary only for fixed-depth mode. Restore older v6 axes through bounds and
  normalize variable nominal endpoints to the repeated final top. Explain the mHM convention on output axes.
- Reuse general axis metadata supplied by NetCDF writers. Restart coordinates use
  `write_horizon_metadata` for units, layer-top descriptions, and the variable-bottom explanation.
  Use the coordinate variable's `comment` attribute for the explanation.
- Both neutron formulations exclude the final soil layer. COSMIC L3 must be finite and positive;
  mode-1 class properties map by horizon before harmonic upscaling.

## Fortran conventions and pitfalls

- Follow surrounding formatting; use explicit interfaces, `implicit none`, `intent`, `mo_kind` kinds,
  and typed literals. Prefer `use ..., only:` and private-by-default modules with explicit public exports.
- Prefix type-bound implementations by their owning type (`mpr_*` for `mpr_t`). Declare `self` with
  `target`; keep dummy declarations together before local declarations. Prefer compact declaration
  sections, especially in short routines and interface specifications; do not add blank lines solely
  to separate imports, arguments, or locals. Keep useful block separation in longer routines.
  Declare each dummy on its own line; document each except `self`/`this` on its declaration with `!<`.
  Remove superseded prepended comments and parameter descriptions instead of duplicating them.
- Import `mo_logging` without `only`; its logging macros depend on module names not visible at the call site.
- Extend NetCDF variable metadata arrays with `mo_netcdf_utils:add_var`; avoid in-place array constructors
  involving derived types with allocatable components.
- An unallocated allocatable or disassociated pointer actual is absent for a compatible `optional`
  dummy that is neither allocatable nor a pointer. Check `present()` before use; this does not apply to allocatable/pointer dummies.
- Do not guard invalid access with `.and.`, `.or.`, or `merge()`; evaluation need not short-circuit.
  Use separate presence/allocation/association guards before dereferencing or querying array size/bounds.
- Guard `log_*` macros with block `if ... then ... end if`, never a one-line `if`.
  They expand to their own conditional `write`; nesting one-line logical `if` statements is invalid Fortran.
- A dummy's `target` alone does not ensure a persistent caller-owned target. Keep domain instances
  and persistent publication targets valid; detach aliases before freeing storage. Dangling pointers are not disassociated pointers.
- Mutable local declaration initialization implies `SAVE`; initialize scratch in executable statements.
  Derived-type component defaults remain valid. Keep thread-called routines reentrant (use `recursive` where needed by supported compilers).
- `intent(out)` can deallocate allocatables/components on entry; intrinsic assignment can reallocate
  arrays. Reset published arrays in place. Pointer `intent(in)` protects association, not target values.
- `=>` associates pointers; `=` assigns values. Derived-type assignment copies pointer associations,
  not their targets; do not assume copying a connected container creates independent state.
- Use constructs supported by the CI compiler matrix; do not rely on one compiler accepting questionable interfaces or syntax.

## OpenMP

- Use OpenMP for substantial independent cell/node work, including large array passes, when thread
  overhead can be amortized. Keep small vectors and I/O serial; preserve subprocess, timestep, and routing-level dependencies.
- Prefer `schedule(static)` for regular loops, disjoint writes, private scratch, and explicit data sharing.
  Reuse outer parallel regions with worksharing loops; avoid nested teams and per-cell team creation.
- All team threads must encounter worksharing/barriers consistently. Preserve required synchronization
  and summation order where practical; do not add `nowait` across dependent work.
- Follow the existing kernel OpenMP policy for meteo PET and output preparation; do not add cell-count
  cutoffs without workload evidence. Preserve existing routing eligibility guards. Keep serial
  builds functional, including conditional `omp_lib` imports/calls. Validate one and multiple threads,
  exercise the parallel path, and provide representative timings for claimed performance improvements.

## Schemas and documentation

- `nml-schemas/*.yml` and `nml-config.toml` are authoritative. Run
  `nml-tools generate --config nml-config.toml`; do not hand-edit generated `src/nml/*.f90`, templates, or namelist references.
- Include generated modules, root templates, and `doc/nml/` references together with schema changes.
- Follow the file/module header in [mo_river_flow.f90](src/mRM/mo_river_flow.f90): adapt filenames,
  module names, factual authorship, dates, version/changelog, copyright, and license. Preserve donor provenance; never invent attribution.
- Every production module declares its appropriate `\ingroup` from [doc/directories.dox](doc/directories.dox);
  preserve exact group spelling (including `f_mLM`). Centralize directory/group definitions there,
  qualify `\dir` paths (e.g. `src/nml`), and use `\model_sources` for the root source directory.
- Document types with `!> \class type_name` and procedures with `\brief`. Explain non-obvious flow
  with short block comments; avoid line-by-line commentary on trivial logic.
- Build Doxygen from the repository root with `PROJECT_NUMBER="$(cat version.txt)" doxygen doc/doxygen.config`.
  See [v6 getting started](doc/v6-getting-started.md) for HTML/PDF commands and dependencies.

## Validation

| Change | Mandatory local checks |
| --- | --- |
| Guidance/prose | Review diff, links, and commands; no model runs solely for prose edits. |
| Fortran behavior | Debug build, relevant CTest tests, and affected component scenarios. |
| Schema | Regenerate outputs, build, affected configuration/parameter tests and scenarios. |
| Exchange/lifecycle/ownership/shared timing | Full serial CTest and curated v6 smoke scenarios; focused contract regression coverage. |
| Parallel behavior | Add OpenMP checks at one and multiple threads; exercise parallel work and compare numerical results with justified tolerances. |
| Build/packaging | Verify affected driver, library-only, or Python targets. |
| Doxygen structure/configuration | Generate docs; inspect affected warnings/pages. Check PDF generation when its infrastructure changes. |

Serial Debug setup and full checks (use `ctest ... -R '<affected-test-pattern>'` for focused tests):

```bash
cmake -S . -B debug -DCMAKE_BUILD_TYPE=Debug -DMHM_BUILD_DRIVER=ON \
  -DMHM_BUILD_TESTING=ON -DBUILD_TESTING=ON -DMHM_WITH_OpenMP=OFF
cmake --build debug --parallel 4
ctest --test-dir debug --output-on-failure
bash CI-scripts/run_v6_namelists.sh --log-dir logs_v6 ./debug/app/mhm
```

For parallel changes, configure `openmp_debug` with the same options and `-DMHM_WITH_OpenMP=ON`;
build it, run affected CTest tests with `OMP_NUM_THREADS=1` and `OMP_NUM_THREADS=4`, and run the smoke
runner with `--threads 1` and `--threads 4` against `./openmp_debug/app/mhm`.

- CTest and the smoke runner stage isolated run roots. Keep restart write/read sequences together;
  retain the runner's output configurations, including `mhm-output-template.nml` for required `output_mrm`.
- Every pFUnit `@assert...` stays on one physical source line without continuations; compute long expressions in locals first.
- Add meaningful regressions for changed behavior; do not add tests merely mirroring implementation.
  Missing active inputs must fail with `log_fatal` + `error stop`; failure tests check exit status and the expected diagnostic.
- For scientific refactors, compare affected values against pre-change behavior or an independent
  reference, covering relevant boundary and configuration variants. Successful smoke execution alone
  does not establish numerical equivalence.
- Ownership changes cover domain isolation, repeated initialization, stable bindings, and destruction.
  Report unavailable pFUnit coverage explicitly; standalone tests do not replace skipped `.pf` tests.
- Keep compiler/Release matrices, coverage, and memory tooling in CI; do not repeat the entire CI matrix for every local edit.

## Hygiene

- Create commits only when explicitly requested by the user. Preserve rename recognition by keeping
  pure moves/renames separate from content edits, and group commits by component or topic.
- Register new/renamed sources and tests explicitly in `src/CMakeLists.txt` and `src/tests/CMakeLists.txt`.
- Keep diffs focused and fixtures immutable. Commit schema-generated outputs with their source changes;
  exclude runtime outputs/restarts (`.nc`), logs, build products, and generated HTML/PDF artifacts.
- Keep `cmake/`, `hpc-module-loads/`, and `doc/doxygen-awesome-css/` unchanged; they are external projects.
- Preserve functional documentation infrastructure, headers, figures, scientific overview, licenses,
  references, and release history. Keep current workflow guidance aligned with v6.
