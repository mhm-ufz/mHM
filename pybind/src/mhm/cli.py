"""!
Module to provide a script to execute mHM.

@copyright Copyright 2005-@today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
    mHM is released under the LGPLv3+ license @license_note
"""

import argparse

from .wrapper import model

try:
    from ._version import __version__
except ModuleNotFoundError:  # pragma: no cover
    # package is not installed
    __version__ = "0.0.0.dev0"


def f_version():
    return model.version().decode("utf-8").strip()


def mhm():  # pragma: no cover
    """Execute mhm as a command line program."""
    parser = argparse.ArgumentParser(
        prog="mhm",
        description="The mesoscale hydrological model - mHM with Python bindings",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f_version(),
        help="show mHM version and exit",
    )
    parser.add_argument(
        "--pyversion",
        action="version",
        version=__version__,
        help="show mHM Python bindings version and exit",
    )
    parser.add_argument(
        "-n", "--nml", default="mhm.nml", help="The mHM configuration namelist."
    )
    parser.add_argument(
        "-p",
        "--parameter",
        default="mhm_parameter.nml",
        help="The mHM parameter namelist.",
    )
    parser.add_argument(
        "-o", "--mhm_output", default="mhm_outputs.nml", help="The mHM output namelist."
    )
    parser.add_argument(
        "-r", "--mrm_output", default="mrm_outputs.nml", help="The mRM output namelist."
    )
    parser.add_argument(
        "-q", "--quiet", action="count", default=0, help="Decrease verbosity level."
    )
    parser.add_argument(
        "cwd", nargs="?", default=".", help="The desired working directory."
    )

    args = parser.parse_args()
    # set verbosity
    model.set_verbosity(level=3 - args.quiet)
    # init model
    model.init(
        namelist_mhm=args.nml,
        namelist_mhm_param=args.parameter,
        namelist_mhm_output=args.mhm_output,
        namelist_mrm_output=args.mrm_output,
        cwd=args.cwd,
    )
    # simple run or optimization
    model.run_or_optimize()
    # finalize will deallocate all variables.
    model.finalize()
