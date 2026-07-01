# Domain configuration {#config_domain}

[TOC]

Domain-indexed configuration for combined mHM runs with domain files.

**Namelist**: `config_domain`

## Fields

| Name | Type | Required | Info |
| --- | --- | --- | --- |
| [domain_dirs](#domain_dirs) | string array | no | Domain directories |
| [domain_nmls](#domain_nmls) | string array | no | Domain namelists |

## Field details

### domain_dirs

Domain directories `domain_dirs`

Paths to domain directories for combined runs with "read_domains_from_dirs" enabled.

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `directory-path`
- Required: no

### domain_nmls

Domain namelists `domain_nmls`

Paths to domain-specific namelists when "read_domains_from_dirs" is true.
The path will be interpreted as relative to the given domain directory.

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Required: no
- Default: `"mhm.nml"`

## Example

```fortran
&config_domain
  domain_dirs(:) = ""
  domain_nmls(:) = "mhm.nml"
/
```

