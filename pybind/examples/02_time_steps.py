from pathlib import Path

import matplotlib.pyplot as plt

import mhm

here = Path(__file__).parent

# assuming to be in the mhm repo
mhm.model.init(cwd=here / ".." / "..")
mhm.run.prepare()
ndomians = mhm.run.get_ndomains()
for i in range(1, ndomians + 1):
    mhm.run.prepare_domain(i)  # 1 by default
    while not mhm.run.finished():
        mhm.run.do_time_step()
        mhm.run.write_output()
    mhm.run.finalize_domain()
mhm.run.finalize()
# get runoff before deallocation
runoff = mhm.get_runoff()
# this will deallocate all internal variables
mhm.model.finalize()
# plot runoff
plt.plot(runoff[:, 0])
plt.show()
