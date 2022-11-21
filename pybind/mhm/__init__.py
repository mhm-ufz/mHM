"""Python bindings of mHM."""
from . import cli
from .tools import get_runoff, get_variable
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
    "f_version",
    "get_runoff",
    "get_variable",
    "__version__",
]
