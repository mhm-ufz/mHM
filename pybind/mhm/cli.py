import os
import subprocess
import sys


def _run(name):
    executable = os.path.join(os.path.dirname(__file__), name)
    return subprocess.call([executable] + sys.argv[1:])


def mhm():
    """Execute mhm as a command line program."""
    raise SystemExit(_run("mhm"))
