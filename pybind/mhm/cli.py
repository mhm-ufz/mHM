"""!
Module to provide a script to execute mHM.

@copyright Copyright 2005-@today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
    mHM is released under the LGPLv3+ license @license_note
"""
import os
import subprocess
import sys


def mhm():
    """Execute mhm as a command line program."""
    exe = os.path.join(os.path.dirname(__file__), "mhm")
    if not os.path.exists(exe):
        raise RuntimeError("mhm: python bindings were installed without driver.")
    raise SystemExit(subprocess.call([exe] + sys.argv[1:]))
