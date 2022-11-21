import matplotlib.pyplot as plt
import numpy as np

import mhm

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
    mhm.run.prepare_domain(i)  # 1 by default
    while not mhm.run.finished():
        mhm.run.do_time_step()
        mhm.run.write_output()
    mhm.run.finalize_domain()
mhm.run.finalize()

runoff = mhm.get_runoff()

# this will deallocate all internal variables
mhm.model.finalize()

plt.plot(runoff[:, 0])
plt.show()
