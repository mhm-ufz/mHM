# Example Domains

Run the commands from the repository root:

```bash
mhm_driver examples/
mhm_driver examples/domain_01
mhm_driver examples/domain_02
```

`examples/` is a repository-root symlink to the tracked `example/` directory.

`example/combine.nml` stages the intended `read_domains_from_dirs = .true.`
layout for `example/domain_01` and `example/domain_02`, but it is not expected
to run yet because `domain_dirs` are not wired through the current driver path.
