# mHM - Python bindings

Python bindings to control mHH.

The wrapper (`src/mhm/wrapper.f90`) is just a small layer on top of these
interfaces to be compatible with [f2py](https://numpy.org/doc/stable/f2py/index.html).

To compile everything locally
([editable install](https://pip.pypa.io/en/stable/cli/pip_install/#install-editable)),
you can use pip:

```bash
pip install -e .
```

## Example

Assuming to have the [mhm repository](https://git.ufz.de/mhm/mhm), with the default `develop` branch checked out,
cloned next to this parents folder.

Then you can do the following to simply run mhm without optimization:
```python
import mhm

# assuming the mhm repo to be in the parent dir
mhm.model.init(
    namelist_mhm="mhm.nml",
    namelist_mhm_param="mhm_parameter.nml",
    namelist_mhm_output="mhm_outputs.nml",
    namelist_mrm_output="mrm_outputs.nml",
    cwd="../mhm",
)
mhm.model.run()
mhm.model.finalize()
```

Or you can do the following to control each timestep:
```python
import mhm

# assuming the mhm repo to be in the parent dir
mhm.model.init(
    namelist_mhm="mhm.nml",
    namelist_mhm_param="mhm_parameter.nml",
    namelist_mhm_output="mhm_outputs.nml",
    namelist_mrm_output="mrm_outputs.nml",
    cwd="../mhm",
)

mhm.run.prepare() # global_parameters(:, 3) by default
ndomians = mhm.run.get_ndomains()
for i in range(1, ndomians + 1):
    mhm.run.prepare_domain(domain=i) # 0 by default
    while not mhm.run.finished():
        mhm.run.do_time_step()
        mhm.run.write_output()
    mhm.run.finalize_domain()
mhm.run.finalize()

mhm.model.finalize()
```
