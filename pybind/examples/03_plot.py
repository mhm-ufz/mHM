from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np

import mhm

plt.ion()
fig = plt.figure()
# assuming the mhm repo to be in the parent dir
mhm.model.init(
    namelist_mhm="mhm.nml",
    namelist_mhm_param="mhm_parameter.nml",
    namelist_mhm_output="mhm_outputs.nml",
    namelist_mrm_output="mrm_outputs.nml",
    cwd="../mhm",
)
mhm.run.prepare()  # global_parameters(:, 3) by default
ndomians = mhm.run.get_ndomains()
for i in range(1, ndomians + 1):
    first_time = True
    mhm.run.prepare_domain(i)  # 1 by default
    while not mhm.run.finished():
        mhm.run.do_time_step()
        mhm.run.write_output()
        run_off = mhm.get_variable("L1_total_runoff")
        time = datetime(*mhm.run.current_time())
        if first_time:
            fig.clf()
            ax = fig.add_subplot(1, 1, 1)
            run_im = ax.imshow(run_off, vmin=0, vmax=0.5, origin="lower")
            fig.colorbar(run_im, ax=ax, label="runoff")
            first_time = False
        else:
            run_im.set_data(run_off.T)
        if time.weekday() == 0 and time.hour == 0:
            fig.suptitle(f"{time}")
            fig.canvas.flush_events()
    mhm.run.finalize_domain()
mhm.run.finalize()
mhm.model.finalize()
