r"""Command line forward execution for mHM v6.

\ingroup py_mhm
"""

import argparse

from . import __version__, f_version
from ._run import _defaults, run_forward


def mhm(argv=None):
    """Execute the forward runner using the native driver's command line options."""
    main_file, para_file, out_file, cwd = _defaults()
    parser = argparse.ArgumentParser(
        prog="mhm",
        description="The mesoscale hydrological model - mHM v6 (forward run)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-V", "--version", action="version", version=f_version)
    parser.add_argument("--pyversion", action="version", version=__version__)
    parser.add_argument(
        "-n", "--nml", default=main_file, help="The mHM configuration namelist."
    )
    parser.add_argument(
        "-p", "--parameters", default=para_file, help="The mHM parameter namelist."
    )
    parser.add_argument(
        "-o", "--output", default=out_file, help="The mHM output namelist."
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase logging verbosity level.",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="count",
        default=0,
        help="Decrease logging verbosity level.",
    )
    parser.add_argument(
        "--log-date", action="store_true", help="Output date while logging."
    )
    parser.add_argument(
        "--log-time", action="store_true", help="Output time while logging."
    )
    parser.add_argument(
        "--log-line",
        action="store_true",
        help="Output file and line information while logging.",
    )
    parser.add_argument(
        "--log-scope",
        default="#",
        help="Comma-separated logging scopes; '#' selects root.",
    )
    parser.add_argument("--log-file", help="Write to a log file.")
    parser.add_argument(
        "cwd", nargs="?", default=cwd, help="The desired working directory."
    )
    args = parser.parse_args(argv)
    run_forward(
        main_file=args.nml,
        para_file=args.parameters,
        out_file=args.output,
        cwd=args.cwd,
        verbose=args.verbose,
        quiet=args.quiet,
        log_date=args.log_date,
        log_time=args.log_time,
        log_line=args.log_line,
        log_scope=args.log_scope,
        log_file=args.log_file,
    )
