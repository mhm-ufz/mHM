# Example domains

From the repository root, run:

```bash
./debug/app/mhm examples/domain_01
./debug/app/mhm examples/domain_02
./debug/app/mhm examples/
./debug/app/mhm -n combine.nml examples/
```

The top-level `mhm.nml` configures both domains in one file. `combine.nml`
uses `read_domains_from_dirs = .true.` and reads each domain's local
configuration. Inputs are immutable; outputs and restarts belong to each run.

The nml-tools 0.6 project profiles describe these file sets: `standard` uses
`main` + `parameter` + `output`, while `combined` uses `meta` + `parameter` +
`output`. The combined profile suggests `combine.nml`; the model loads the
domain-local files separately at runtime. The `member` profile contains only
`main`, for a domain within a combined run whose parameter and output settings
are supplied by that combined run. These profiles also support GUI project
selection without changing the model's file-loading rules.

CLI validation infers array dimensions from explicit `config_project` values.
For a separate parameter file, provide the main or meta file with `--dim-file`:

```bash
nml-tools validate --config nml-config.toml --profile main examples/mhm.nml
nml-tools validate --config nml-config.toml --profile meta examples/combine.nml
nml-tools validate --config nml-config.toml --profile main examples/domain_01/mhm.nml
nml-tools validate --config nml-config.toml --profile parameter \
  --dim-file examples/mhm.nml examples/parameters.nml
nml-tools validate --config nml-config.toml --profile output examples/outputs.nml
```

An existing output file may omit `output_mhm`, `output_mlm`, or `output_mrm`.
Each omitted group uses its empty-group defaults; configured output paths still
control stream activation, and output variable switches retain their defaults.

`max_layers` is the per-file allocation cap for `config_mpr%soil_depth`.
The active layer count is `config_mpr%n_layers(:)`. Domain-local allocation
caps must fit within the combined configuration's cap.

`fixtures/` contains additional inputs for the v6 validation scenarios.
