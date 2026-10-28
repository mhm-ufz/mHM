r"""Forward execution through the v6 domain runner.

\ingroup py_mhm
"""

import os

from .wrapper import metadata, runner


def _defaults():
    return tuple(value.decode("utf-8").strip() for value in metadata.run_defaults())


def run_forward(
    main_file=None,
    para_file=None,
    out_file=None,
    cwd=None,
    *,
    verbose=0,
    quiet=0,
    log_date=False,
    log_time=False,
    log_line=False,
    log_scope="#",
    log_file=None,
):
    """Execute a forward run, resolving model paths against ``cwd``.

    Omitted paths use the native driver's defaults: mhm.nml, parameters.nml,
    outputs.nml, and the current directory. Paths accept strings or pathlike
    objects. Logging options match the CLI; log_file is relative to the process
    working directory. Model optimization and stepwise control are not provided.
    """
    paths = tuple(
        default if value is None else os.fsdecode(value)
        for value, default in zip((main_file, para_file, out_file, cwd), _defaults())
    )
    runner.run_forward(
        *paths,
        verbose=verbose,
        quiet=quiet,
        log_date=log_date,
        log_time=log_time,
        log_line=log_line,
        log_scope=log_scope,
        log_file="" if log_file is None else os.fsdecode(log_file),
    )
