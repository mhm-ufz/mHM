# -*- coding: utf-8 -*-
"""
Clean up the mHM output in the check case folders.

Author
------
    Sebastian Mueller

Created
-------
    Feb 2020
"""
import glob
import os
import shutil

# output folder and reference folder
OUT = "output_b1"


if __name__ == "__main__":
    # checking path
    cases_path = os.path.dirname(os.path.realpath(__file__))
    # get all cases folders (in the cases_path)
    raw_paths = glob.glob(os.path.join(cases_path, "case*"))
    # sort out the mrm check cases
    cases = [path for path in raw_paths if "mrm" not in os.path.basename(path)]
    # sort the cases by name
    cases.sort()
    for case in cases:
        base = os.path.basename(case)
        print(base)
        # get the output directory
        out_dir = os.path.join(case, OUT)
        # remove it
        shutil.rmtree(out_dir, ignore_errors=True)
        logs = glob.glob(os.path.join(case, "*_" + base + "_log.txt"))
        for log in logs:
            os.remove(log)
    # final result
    print(" ..cleaned.")
