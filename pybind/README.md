# mHM - Python bindings

[TOC]

Python bindings to control mHM.

The wrapper (`mhm/wrapper.f90`) is just a small layer on top of the
interfaces provided by mHM to be compatible with [f2py](https://numpy.org/doc/stable/f2py/index.html).

## Installation

To compile everything after cloning/downloading, you can use pip:

```bash
pip install -v .
```

To install it directly from the git repository you can type:

```bash
pip install -v git+https://git.ufz.de/mhm/mhm.git
```

There will be a PyPI package in the future to install the latest release with:

```bash
pip install mhm
```

Installing the mHM Python package will provide the `mhm` command to execute mHM the traditional way.

## Documentation

See `mhm.tools` and `wrapper.f90` for further information on the provided routines.

## Examples

If you have cloned the repository, you can do the following to simply run mhm without optimization:

```python
import mhm

# assuming the mhm repo to be in the parent dir
mhm.model.init(
    namelist_mhm="mhm.nml",
    namelist_mhm_param="mhm_parameter.nml",
    namelist_mhm_output="mhm_outputs.nml",
    namelist_mrm_output="mrm_outputs.nml",
    cwd=".",
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
    cwd=".",
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

See also the `examples` folder.
