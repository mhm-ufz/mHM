"""!
Tools to interact with mHM.

@copyright Copyright 2005-@today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
    mHM is released under the LGPLv3+ license @license_note
@ingroup mhm
"""

import numpy as np

from . import wrapper as wr


def get_parameter():
    """
    Get parameter names and configuration.

    @retval names (List[Str]): Names of all used parameters in mHM.
    @retval config (numpy.ndarray): Configuration for all used parameters (min, max, value, flag, scale).
    """
    para_n = wr.get.parameter_length()
    names = [
        wr.get.parameter_name(i).decode("utf-8").strip() for i in range(1, para_n + 1)
    ]
    config = wr.get.parameter_config(para_n)
    return names, config


def get_runoff():
    """
    Get 2D array of runoff time-series for all gauges.

    @retval runoff (numpy.ndarray): The runoff for all gauges with dims (time, gauge).
    """
    shp = wr.get.runoff_shape()
    return wr.get.runoff(*shp)


def get_runoff_eval(gauge_id):
    """
    Get 2D array of simulated and observed runoff time-series for selected gauges.

    @retval runoff (numpy.ndarray(TS, 2)): The runoff for selected gauges with dims (time-steps, 2).
    """
    length = wr.get.runoff_eval_length(gauge_id)
    return wr.get.runoff_eval(gauge_id, length)


def get_mask(level, indexing="ij", selection=False):
    """
    Get mask for a certain mHM level.

    @param level (str): Name level ("L0", "L1", "L11", "L2")
    @param indexing (str, optional): Indexing for the 2D mask,
        either "xy" or "ij" (yx order), by default "ij"
    @param selection (bool): A masked value in mHM indicates cells inside the domain.
        In numpy, masked values are outside the domain. That means, by default False
    @retval mask (numpy.ndarray): Boolean numpy array holding the mask.
    @throws ValueError: If the level is not in ["L0", "L1", "L11" or "L2"].
    """
    level = level.lower()
    shp = getattr(wr.get, level + "_domain_shape")()
    # mask in mHM is the opposite in numpy
    mask = np.array(
        getattr(wr.get, level + "_domain_mask")(m=shp[0], n=shp[1]), dtype=bool
    )
    if not selection:
        mask = ~mask
    return mask.T if indexing == "ij" else mask


def get_variable(name, index=1, indexing="ij", compressed=False):
    """
    Get a specific variable from mHM in the current time-step.

    @param name (str): Name of the variable
    @param index (int, optional): If the variable has an additional dimension
        (e.g. the horizon), one needs to specify an index, by default 1
    @param indexing (str, optional): Indexing for the 2D variable,
        either "xy" or "ij", by default "ij"
    @param compressed (bool): Whether the data should be flattened and only contain
        values for each unmasked domain cell. By default False
    @retval variable (numpy.ndarray): Numpy array holding the desired variable.
    @throws ValueError: If the variable name doesn't start with "L0", "L1", "L11" or "L2".
    """
    name = name.upper()  # convention
    grid = name.split("_")[0].lower()
    if grid not in ["l0", "l1", "l11", "l2"]:
        raise ValueError(f"Unknown variable: {name}")
    n = getattr(wr.get, grid + "_domain_size")()
    var = getattr(wr.get, grid + "_variable")(name=name, n=n, idx=index)
    if compressed:
        return np.asarray(var, dtype=float)
    # ncols, nrows, ncells, xll, yll, cell_size, no_data
    grid_info = getattr(wr.get, grid + "_domain_info")()
    # mask in mHM is the opposite in numpy (selection)
    sel = get_mask(grid, indexing="xy", selection=True)
    sel = sel.ravel(order="F")
    output = np.ma.empty_like(sel, dtype=float)
    output.fill_value = grid_info[-1]
    output.mask = ~sel
    output[sel] = var
    output = output.reshape((grid_info[0], grid_info[1]), order="C")
    return output.T if indexing == "xy" else output


def set_meteo(
    time,
    pre=None,
    temp=None,
    pet=None,
    tmin=None,
    tmax=None,
    netrad=None,
    absvappress=None,
    windspeed=None,
    ssrd=None,
    strd=None,
    tann=None,
    compressed=False,
    indexing="ij",
):
    """
    Set meteo data with a time stamp in mHM.

    @param time (datetime.datetime): Timestamp
    @param pre (numpy.ndarray, optional):         [mm]      Precipitation
    @param temp (numpy.ndarray, optional):        [degC]    Air temperature
    @param pet (numpy.ndarray, optional):         [mm TS-1] Potential evapotranspiration
    @param tmin (numpy.ndarray, optional):        [degC]    minimum daily air temperature
    @param tmax (numpy.ndarray, optional):        [degC]    maximum daily air temperature
    @param netrad (numpy.ndarray, optional):      [W m2]    net radiation
    @param absvappress (numpy.ndarray, optional): [Pa]      absolute vapour pressure
    @param windspeed (numpy.ndarray, optional):   [m s-1]   windspeed
    @param ssrd (numpy.ndarray, optional):        [W m2]    short wave radiation
    @param strd (numpy.ndarray, optional):        [W m2]    long wave radiation
    @param tann (numpy.ndarray, optional):        [degC]    annual mean air temperature
    @param compressed (bool): Whether the data is flattened and only contains
        values for each unmasked domain cell. By default False
    @param indexing (str, optional): Indexing for 2D arrays if data is not compressed,
        either "xy" or "ij" (yx order), by default "ij"
    """
    sel = slice(None) if compressed else get_mask("L1", indexing, selection=True)
    if pre is not None:
        wr.set.meteo(pre[sel], "PRE", time.year, time.month, time.day, time.hour)
    if temp is not None:
        wr.set.meteo(temp[sel], "TEMP", time.year, time.month, time.day, time.hour)
    if pet is not None:
        wr.set.meteo(pet[sel], "PET", time.year, time.month, time.day, time.hour)
    if tmin is not None:
        wr.set.meteo(tmin[sel], "TMIN", time.year, time.month, time.day, time.hour)
    if tmax is not None:
        wr.set.meteo(tmax[sel], "TMAX", time.year, time.month, time.day, time.hour)
    if netrad is not None:
        wr.set.meteo(netrad[sel], "NETRAD", time.year, time.month, time.day, time.hour)
    if absvappress is not None:
        wr.set.meteo(
            absvappress[sel], "ABSVAPPRESS", time.year, time.month, time.day, time.hour
        )
    if windspeed is not None:
        wr.set.meteo(
            windspeed[sel], "WINDSPEED", time.year, time.month, time.day, time.hour
        )
    if ssrd is not None:
        wr.set.meteo(ssrd[sel], "SSRD", time.year, time.month, time.day, time.hour)
    if strd is not None:
        wr.set.meteo(strd[sel], "STRD", time.year, time.month, time.day, time.hour)
    if tann is not None:
        wr.set.meteo(tann[sel], "TANN", time.year, time.month, time.day, time.hour)
