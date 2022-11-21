"""Module to provide a script to execute mHM."""
import os
import subprocess
import sys


def mhm():
    """Execute mhm as a command line program."""
    exe = os.path.join(os.path.dirname(__file__), "mhm")
    raise SystemExit(subprocess.call([exe] + sys.argv[1:]))
