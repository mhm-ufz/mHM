import matplotlib.pyplot as plt

import mhm

# assuming the mhm repo to be in the parent dir
mhm.model.init(
    namelist_mhm="mhm.nml",
    namelist_mhm_param="mhm_parameter.nml",
    namelist_mhm_output="mhm_outputs.nml",
    namelist_mrm_output="mrm_outputs.nml",
    cwd="../mhm",
)
mhm.model.run()
runoff = mhm.get_runoff()

# finalize will deallocate all variables.
mhm.model.finalize()

plt.plot(runoff[:, 0])
plt.show()
