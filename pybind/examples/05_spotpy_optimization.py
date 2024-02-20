"""
Copyright 2023 by Sebastian Müller and Tobias Houska
"""

from pathlib import Path

import matplotlib.pyplot as plt
import spotpy
from spotpy.parameter import Constant, Uniform

import mhm


class spot_setup:
    """
    Spotpy setup for optimizing mHM with Spotpy.

    An existing mHM setup is optimized for one gauge with KGE as
    objective function.

    The `optimize` variable needs to be set to `.false.`, to prevent
    internal optimization triggers in mHM.

    All runtime output will be disabled.

    Parameters
    ----------
    namelist_mhm : str, optional
        Main namelist file path for mHM, by default "mhm.nml"
    namelist_mhm_param : str, optional
        Parameter namelist file path for mHM, by default "mhm_parameter.nml"
    namelist_mhm_output : str, optional
        Output namelist file path for mHM, by default "mhm_outputs.nml"
    namelist_mrm_output : str, optional
        Output namelist file path for mRM, by default "mrm_outputs.nml"
    cwd : str, optional
        Working directory to execute mHM, by default "."
    gauge_id : int, optional
        Gauge ID for optimization (starting at 1), by default 1
    minimize : bool, optional
        Flag to change the sign of the objective function
        for minimizing algorithms, by default False
    """

    def __init__(
        self,
        namelist_mhm="mhm.nml",
        namelist_mhm_param="mhm_parameter.nml",
        namelist_mhm_output="mhm_outputs.nml",
        namelist_mrm_output="mrm_outputs.nml",
        cwd=".",
        gauge_id=1,
        minimize=False,
    ):
        # only show errors
        mhm.model.set_verbosity(level=1)
        # initialize model
        mhm.model.init(
            namelist_mhm=namelist_mhm,
            namelist_mhm_param=namelist_mhm_param,
            namelist_mhm_output=namelist_mhm_output,
            namelist_mrm_output=namelist_mrm_output,
            cwd=cwd,
        )
        # disable output during optimization
        mhm.model.disable_output()
        # set the gauge ID to optimize
        self.gauge_id = gauge_id
        # factor for target function when minimizing
        self.factor = -1 if minimize else 1

        # get parameter configuration of mHM
        para_names, para_config = mhm.get_parameter()
        self.params = []
        for name, config in zip(para_names, para_config):
            # names shouldn't have "," (like "GeoParam(1,:)")
            name = name.replace(",:", "")
            use_para = config[3] > 0.5  # Flag: 0 or 1
            if use_para:
                para = Uniform(name, low=config[0], high=config[1], optguess=config[2])
            else:
                para = Constant(name, config[2])
            self.params.append(para)

        # get observed runoff read by mHM
        self.runoff_obs = mhm.get_runoff_eval(gauge_id=self.gauge_id)[:, 1]

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, parameter):
        mhm.model.run_with_parameter(parameter)
        # runoff = [sim, obs]
        return mhm.get_runoff_eval(gauge_id=self.gauge_id)[:, 0]

    def evaluation(self):
        return self.runoff_obs

    def objectivefunction(self, simulation, evaluation):
        mask = evaluation > 0.0
        like = spotpy.objectivefunctions.kge(evaluation[mask], simulation[mask])
        return self.factor * like


if __name__ == "__main__":
    here = Path(__file__).parent
    # only use test domain 1
    test_domain = here / ".." / ".." / "test_domain"
    spot_setup = spot_setup(gauge_id=1, minimize=False, cwd=test_domain)

    # Start a spotpy analysis
    sampler = spotpy.algorithms.lhs(spot_setup, dbname="LHS_mHM", dbformat="csv")
    sampler.sample(repetitions=200)
    # finalize mHM
    mhm.model.finalize()

    # Load the results gained with the latin hypercube sampler
    results = spotpy.analyser.load_csv_results("LHS_mHM")
    plt.plot(results["like1"])
    plt.show()
