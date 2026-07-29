# Example Domains

Run the commands from the repository root:

```bash
mhm_driver examples/
mhm_driver -n combine.nml examples/
mhm_driver examples/domain_01
mhm_driver examples/domain_02
```

`examples/` is a repository-root symlink to the tracked `example/` directory.

`example/combine.nml` uses `read_domains_from_dirs = .true.` and reads
`config_domain` to find the per-domain namelists in `example/domain_01` and
`example/domain_02`.

`max_layers` in `config_project` is only the per-file allocation cap for
`config_mpr%soil_depth`; the active layer count is `config_mpr%n_layers(:)`.
For combined runs, each domain-local `max_layers` must fit inside the top-level
cap, and the top-level cap is used when reading each domain file.
