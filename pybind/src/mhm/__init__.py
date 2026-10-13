r"""mHM version metadata and example downloads; model bindings will be added later.

\ingroup py_mhm
"""

from .download import download_examples
from .wrapper import metadata as _metadata

try:
    from ._version import __version__
except ModuleNotFoundError:  # source trees without generated packaging metadata
    __version__ = "0.0.0.dev0"


def __getattr__(name):
    """Expose the compiled Fortran version on demand."""
    if name == "f_version":
        return _metadata.version().decode("utf-8").strip()
    raise AttributeError(f"module {__name__} has no attribute {name}")


__all__ = ["__version__", "f_version", "download_examples"]
