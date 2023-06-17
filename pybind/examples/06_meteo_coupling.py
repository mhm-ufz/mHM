from datetime import datetime as dt
from datetime import timedelta as td
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import mhm

here = Path(__file__).parent
test_domain = here / ".." / ".." / "test_domain"
pre = xr.open_dataset(test_domain / "input" / "meteo" / "pre" / "pre.nc")

# determine meteo timestep from first two time stamps
times = [dt.combine(d, t) for d, t in zip(pre.time.dt.date.data, pre.time.dt.time.data)]
meteo_ts = (times[1] - times[0]) // td(hours=1)  # 1 or 24

# only show errors
mhm.model.set_verbosity(level=1)
# configure coupling and init model
mhm.model.config_coupling(meteo_expect_pre=True, meteo_timestep=meteo_ts)
mhm.model.init(cwd=test_domain)
# coupled run only supports one domain
mhm.run.prepare()
mhm.run.prepare_domain()  # domain=1 by default
# mask in mHM is the opposite as in numpy and in xy order (transpose it)
shp = mhm.get.l1_domain_shape()
sel = np.array(mhm.get.l1_domain_mask(m=shp[0], n=shp[1]), dtype=bool).T
# run model
while not mhm.run.finished():
    time = dt(*mhm.run.current_time())
    if time.hour % meteo_ts == 0:
        data = pre["pre"].sel(time=time).data[sel]
        mhm.set.meteo(data, "PRE", time.year, time.month, time.day, time.hour)
    mhm.run.do_time_step()
    mhm.run.write_output()
mhm.run.finalize_domain()
mhm.run.finalize()
# get runoff before deallocation
runoff_sim_obs = mhm.get_runoff_eval(gauge_id=1)
# this will deallocate all internal variables
mhm.model.finalize()
pre.close()

# plotting
plt.plot(runoff_sim_obs[:, 0])
plt.plot(runoff_sim_obs[:, 1], linestyle="--")
plt.show()
