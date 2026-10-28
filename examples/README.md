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

`max_layers` is the per-file allocation cap for `config_mpr%soil_depth`.
The active layer count is `config_mpr%n_layers(:)`. Domain-local allocation
caps must fit within the combined configuration's cap.

`fixtures/` contains additional inputs for the v6 validation scenarios.
