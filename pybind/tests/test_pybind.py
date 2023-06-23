import os
import shutil
import unittest
from datetime import datetime as dt
from datetime import timedelta as td
from pathlib import Path

import numpy as np
import xarray as xr
from numpy.testing import assert_allclose

import mhm

HERE = Path(__file__).parent


class TestPybind(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        os.chdir(HERE)
        self.test_domain = HERE / "test_domain"
        mhm.download_test(path=self.test_domain)

    def tearDown(self):
        # move out of test-domain directory after each test
        os.chdir(HERE)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.test_domain)

    def test_coupling(self):
        ref_runoff = np.loadtxt(HERE / "test_files" / "ref_runoff.txt")
        pre = xr.open_dataset(self.test_domain / "input" / "meteo" / "pre" / "pre.nc")
        temp = xr.open_dataset(
            self.test_domain / "input" / "meteo" / "tavg" / "tavg.nc"
        )
        pet = xr.open_dataset(self.test_domain / "input" / "meteo" / "pet" / "pet.nc")

        # determine meteo timestep from first two time stamps (in pre)
        times = [
            dt.combine(d, t)
            for d, t in zip(pre.time.dt.date.data, pre.time.dt.time.data)
        ]
        meteo_ts = (times[1] - times[0]) // td(hours=1)  # 1 or 24

        # configure coupling and init model
        mhm.model.config_coupling(
            couple_case=1,
            meteo_expect_pre=True,
            meteo_expect_temp=True,
            meteo_expect_pet=True,
            meteo_timestep=meteo_ts,
        )
        mhm.model.set_verbosity(level=1)
        mhm.model.init(cwd=self.test_domain)
        mhm.model.disable_output()
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
        assert_allclose(runoff_sim_obs, ref_runoff)

    def test_get_variable(self):
        ref_stream = np.loadtxt(HERE / "test_files" / "ref_stream_mean.txt")
        stream = 0.0
        n = 0
        mhm.model.set_verbosity(level=1)
        mhm.model.init(cwd=self.test_domain)
        mhm.model.disable_output()
        mhm.run.prepare()
        mhm.run.prepare_domain()
        while not mhm.run.finished():
            mhm.run.do_time_step()
            stream = (stream * n + mhm.get_variable("L11_QMOD")) / (n + 1)
            n += 1  # increase counter
        mhm.run.finalize_domain()
        mhm.run.finalize()
        mhm.model.finalize()
        assert_allclose(stream.filled(fill_value=np.nan), ref_stream, equal_nan=True)

    def test_version(self):
        fv = mhm.f_version
        pv = mhm.__version__
        self.assertTrue(fv.split(".")[:2] == pv.split(".")[:2])


if __name__ == "__main__":
    unittest.main()
