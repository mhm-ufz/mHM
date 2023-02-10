"""!
Python bindings of mHM.

@copyright Copyright 2005-@today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
    mHM is released under the LGPLv3+ license @license_note
"""

##
# @dir pybind
# @brief @copybrief mhm
# @details @copydetails mhm
#
# @dir mhm
# @brief mhm Python package
# @details @copydetails mhm
#
# @defgroup   mhm mhm - Python bindings
# @brief      Python wrapper to control mHM.
# @details    The mhm python package provides a wrapper module to control mHM from Python.
#             This includes:
#             - initialization, running and finalizing a model run
#             - control of time stepping
#             - access to internal variables

from . import cli, download
from .download import download_test
from .tools import get_parameter, get_runoff, get_runoff_eval, get_variable
from .wrapper import get, model, run, set

try:
    from ._version import __version__
except ModuleNotFoundError:  # pragma: no cover
    # package is not installed
    __version__ = "0.0.0.dev0"


def __getattr__(name):
    """Magic method to provide 'f_version' in Python."""
    if name == "f_version":
        return model.version().decode("utf-8").strip()
    raise AttributeError(f"module {__name__} has no attribute {name}")


__all__ = [
    "cli",
    "model",
    "get",
    "set",
    "run",
    "download",
    "download_test",
    "f_version",
    "get_runoff",
    "get_variable",
    "get_runoff_eval",
    "get_parameter",
    "__version__",
]
