from datetime import datetime as dt
from datetime import timedelta as td
from pathlib import Path

import matplotlib.pyplot as plt
import xarray as xr

import mhm

here = Path(__file__).parent
test_domain = here / ".." / ".." / "test_domain"
pre = xr.open_dataset(test_domain / "input" / "meteo" / "pre" / "pre.nc")
temp = xr.open_dataset(test_domain / "input" / "meteo" / "tavg" / "tavg.nc")
pet = xr.open_dataset(test_domain / "input" / "meteo" / "pet" / "pet.nc")

# determine meteo timestep from first two time stamps (in pre)
times = [dt.combine(d, t) for d, t in zip(pre.time.dt.date.data, pre.time.dt.time.data)]
meteo_ts = (times[1] - times[0]) // td(hours=1)  # 1 or 24

# configure coupling and init model
mhm.model.config_coupling(
    couple_case=1,
    meteo_expect_pre=True,
    meteo_expect_temp=True,
    meteo_expect_pet=True,
    meteo_timestep=meteo_ts,
)
mhm.model.init(cwd=test_domain)
# coupled run only supports one domain
mhm.run.prepare()
mhm.run.prepare_domain()  # domain=1 by default
# run model
while not mhm.run.finished():
    time = dt(*mhm.run.current_time())
    if time.hour % meteo_ts == 0:
        mhm.set_meteo(time=time, pre=pre["pre"].sel(time=time).data)
        mhm.set_meteo(time=time, temp=temp["tavg"].sel(time=time).data)
        mhm.set_meteo(time=time, pet=pet["pet"].sel(time=time).data)
    mhm.run.do_time_step()
    mhm.run.write_output()
mhm.run.finalize_domain()
mhm.run.finalize()
# get runoff before deallocation
runoff_sim_obs = mhm.get_runoff_eval(gauge_id=1)
# this will deallocate all internal variables
mhm.model.finalize()
pre.close()
temp.close()
pet.close()

# plotting
plt.plot(runoff_sim_obs[:, 0])
plt.plot(runoff_sim_obs[:, 1], linestyle="--")
plt.show()
