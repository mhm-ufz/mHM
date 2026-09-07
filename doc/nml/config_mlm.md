# mLM configuration {#config_mlm}

[TOC]

Configuration for the v6 lake component (mLM).

**Namelist**: `config_mlm`

## Fields

| Name | Type | Declared required | Input required | Info |
| --- | --- | --- | --- | --- |
| [output_path](#output_path) | string array | no | no | mLM lake-point output path |
| [read_restart](#read_restart) | logical array | no | no | Read mLM restart |
| [restart_input_path](#restart_input_path) | string array | no | no | mLM restart input path |
| [write_restart](#write_restart) | logical array | no | no | Write mLM restart |
| [restart_output_path](#restart_output_path) | string array | no | no | mLM restart output path |

## Field details

### output_path

mLM lake-point output path `output_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no

### read_restart

Read mLM restart `read_restart`

Summary:
- Type: `logical, dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `.false.`

### restart_input_path

mLM restart input path `restart_input_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no

### write_restart

Write mLM restart `write_restart`

Summary:
- Type: `logical, dimension(n_domains)`
- Declared required: no
- Input required: no
- Default: `.false.`

### restart_output_path

mLM restart output path `restart_output_path`

Summary:
- Type: `character(len=buf), dimension(n_domains)`
- Item format: `file-path`
- Declared required: no
- Input required: no

## Example

```fortran
&config_mlm
  output_path(:) = ""
  read_restart(:) = .false.
  restart_input_path(:) = ""
  write_restart(:) = .false.
  restart_output_path(:) = ""
/
```

