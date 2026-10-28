# mHM Python package

The package runs mHM v6 forward simulations through the same domain runner as
the native executable. It also exposes version metadata and example downloads:

```python
import mhm

print(mhm.__version__)  # Python distribution version
print(mhm.f_version)    # Compiled Fortran model version
mhm.download_examples(domain=1, path="my-domain")
mhm.run_forward(cwd="my-domain")
```

```bash
mhm my-domain
python -m mhm my-domain
python -m mhm -n custom.nml -p parameters.nml -o outputs.nml -q my-domain

mhm-download --domain 2 --path my-domain
python -m mhm.download --branch v6alpha --domain 1
python -m mhm.download --source git.ufz --domain 2
```

The run defaults are `mhm.nml`, `parameters.nml`, `outputs.nml`, and the current
directory, shared with the native driver. Model paths are relative to the supplied
run directory; running the model does not change Python's working directory.
`mhm.run_forward` accepts optional `main_file`, `para_file`, `out_file`, and `cwd`
arguments as strings or pathlike objects. Logging options include repeated
`-q`/`-v`, `--log-file`, `--log-date`, `--log-time`, `--log-line`, and `--log-scope`;
the Python function accepts the corresponding keyword arguments. Log files are
relative to the process working directory. Use `mhm --help` for all options,
`mhm --version` for the compiled version, and `mhm --pyversion` for the distribution
version. Optimization and stepwise model-control bindings will follow later.

Development distributions download from `v6alpha`; releases use the matching
version tag. Downloads contain the selected example domain's configuration
and inputs, ready for a forward run. The downloader tries GitHub first and falls
back to UFZ GitLab if the transfer fails or the archive is unusable. Select
`--source github` or `--source git.ufz` to use one host exclusively. The Python
`download_examples` function accepts the same `source` argument, defaulting to
`"auto"`.

To build and test locally with a Fortran/C compiler and NetCDF available:

```bash
python -m pip install '.[test]'
python -m pytest pybind/tests/
```
