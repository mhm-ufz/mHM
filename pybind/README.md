# mHM Python package

This package retains the build and distribution infrastructure for future v6
bindings. It currently exposes version metadata and example-domain downloads:

```python
import mhm

print(mhm.__version__)  # Python distribution version
print(mhm.f_version)    # Compiled Fortran model version
mhm.download_examples(domain=1, path="my-domain")
```

```bash
mhm-download --domain 2 --path my-domain
mhm-download --branch v6alpha --domain 1
```

Development distributions download from `v6alpha`; releases use the matching
version tag. Downloads contain the selected example domain's configuration
and inputs, ready for the standalone Fortran `mhm` executable.

Python model execution and field-access APIs are intentionally absent while the
new wrapper is developed. The package does not install a competing `mhm` CLI.

To build and test locally with a Fortran/C compiler and NetCDF available:

```bash
python -m pip install '.[test]'
python -m pytest pybind/tests/
```
