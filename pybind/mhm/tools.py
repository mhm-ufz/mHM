"""!
Tools to interact with mHM.

@copyright Copyright 2005-@today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
    mHM is released under the LGPLv3+ license @license_note
@ingroup mhm
"""
import numpy as np

from .wrapper import get


def get_parameter():
    """
    Get parameter names and configuration.

    @retval names (List[Str]): Names of all used parameters in mHM.
    @retval config (numpy.ndarray): Configuration for all used parameters (min, max, value, flag, scale).
    """
    para_n = get.parameter_length()
    names = [
        get.parameter_name(i).decode("utf-8").strip() for i in range(1, para_n + 1)
    ]
    config = get.parameter_config(para_n)
    return names, config


def get_runoff():
    """
    Get 2D array of runoff time-series for all gauges.

    @retval runoff (numpy.ndarray): The runoff for all gauges with dims (time, gauge).
    """
    shp = get.runoff_shape()
    return get.runoff(*shp)


def get_runoff_eval(gauge_id):
    """
    Get 2D array of simulated and observed runoff time-series for selected gauges.

    @retval runoff (numpy.ndarray(TS, 2)): The runoff for selected gauges with dims (time-steps, 2).
    """
    length = get.runoff_eval_length(gauge_id)
    return get.runoff_eval(gauge_id, length)


def get_variable(name, index=1, indexing="ij"):
    """
    Get a specific variable from mHM in the current time-step.

    @param name (str): Name of the variable
    @param index (int, optional): If the variable has an additional dimension,
        one needs to specify an index, by default 1
    @param indexing (str, optional): Indexing for the 2D variable,
        either "xy" or "ij", by default "ij"
    @retval variable (numpy.ndarray): Numpy array holding the desired variable.
    @throws ValueError: If the variable name doesn't start with "L0", "L1", "L11" or "L2".
    """
    name = name.upper()  # convention
    grid = name.split("_")[0].lower()
    if grid not in ["l0", "l1", "l11", "l2"]:
        raise ValueError(f"Unknown variable: {name}")
    n = getattr(get, grid + "_domain_size")()
    shp = getattr(get, grid + "_domain_shape")()
    # mask in mHM is the opposite in numpy
    mask = np.array(getattr(get, grid + "_domain_mask")(m=shp[0], n=shp[1]), dtype=bool)
    # ncols, nrows, ncells, xll, yll, cell_size, no_data
    grid_info = getattr(get, grid + "_domain_info")()
    var = getattr(get, grid + "_variable")(name=name, n=n, idx=index)
    # reshaping
    mask = mask.ravel(order="F")
    output = np.ma.empty_like(mask, dtype=float)
    output.fill_value = grid_info[-1]
    output.mask = ~mask
    output[mask] = var
    output = output.reshape((grid_info[0], grid_info[1]), order="C")
    return output.T if indexing == "xy" else output
