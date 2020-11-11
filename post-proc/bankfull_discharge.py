#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author:
    Lennart Schueler

Purpose:
    Calculate the river discharge at bankfull conditions.
"""

from __future__ import division, absolute_import, print_function

import argparse
import pandas as pd
import xarray as xr
import numpy as np
import netcdf4 as nc
from netcdf4 import NcDataset


def find_nearest_idx(array, value):
    return (np.abs(array-value)).argmin()
def find_nearest(array, value):
    return array[find_nearest_idx(array, value)]


def get_cmdline_args(default_return_period):
    parser = argparse.ArgumentParser(description=
                'Calculate the river discharge at bankfull conditions.')
    parser.add_argument('ncin_path',
            help='The path of the mRM NetCDF file with the discharge data')
    parser.add_argument('ncout_path',
            help='The path of the output NetCDF file')
    parser.add_argument('-p', '--return_period', type=float,
            help='The return period of the flood, default: {:1} years'.
                format(default_return_period))
    parser.add_argument(
        "-w",
        "--wetted_perimeter",
        action="store_true",
        default=False,
        dest="peri_bkfl",
        help="Also estimate the wetted perimeter. (default: no)",
    )
    args = parser.parse_args()
    ncin_path = args.ncin_path
    ncout_path = args.ncout_path
    return_period = args.return_period
    peri_bkfl = args.peri_bkfl
    if return_period is None:
        return_period = default_return_period
    return ncin_path, ncout_path, return_period, peri_bkfl


def read_discharge(filename, var_name='Qrouted'):
    """Reads in the discharge from a previous mHM run.

    Assumes that the time variable is named 'time'.
    Converts the time into an array of datetime objects.

    Args:
        filename (str): name of the mHM output file
        var_name (str): name of the variable to be read in
    Returns:
        t (1d ndarray): array of datetime objects
        Q (3d ndarray): the mHM data
    """
    rootgrp = NcDataset(filename, 'r')
    t_nc = rootgrp['time']
    t = t_nc[:]
    t_units = t_nc.units
    try:
        t_cal = t_nc.calendar
    except AttributeError:
        t_cal = 'standard'

    t = nc.num2date(t, units = t_units, calendar = t_cal)
    Q = rootgrp[var_name][:]
    return t, Q

def write_Q_bkfl(Q_bkfl, discharge_filename, ncout_filename, peri_bkfl=False):
    """Copies dims and attrs from given file and writes the bankfull discharge

    Args:
        Q_bkfl (2d ndarray): the bankfull discharge
        discharge_filename (str): the filename of the discharge data
        ncout_filename (str): the output filename
        peri_bkfl (bool): whether to calculate the wetted perimeter derived
            from bankful discharge
    """
    ncin = NcDataset(discharge_filename, 'r')
    ncout = NcDataset(ncout_filename, 'w')

    dims = nc.getDimensions(ncin)
    variables = nc.getVariables(ncin)

    nc.copyDimensions(ncout, dims)
    nc.copyVariables(ncout, variables, skip='Qrouted')

    Q_bkfl_nc = ncout.createVariable('Q_bkfl', 'f8', (dims['northing'].name, dims['easting'].name))
    set_nc_attrs(Q_bkfl_nc, 'Discharge at bankfull conditions')
    Q_bkfl_nc[:] = Q_bkfl
    if peri_bkfl:
        P_bkfl_nc = ncout.createVariable('P_bkfl', 'f8', (dims['northing'].name, dims['easting'].name))
        set_nc_attrs(P_bkfl_nc, 'Perimeter at bankfull conditions', units="m")
        P_bkfl_nc[:] = 4.8 * np.sqrt(Q_bkfl)

def set_nc_attrs(nc_var, long_name, units='m3 s-1'):
    nc_var.setncattr('FillValue', -9999.)
    nc_var.setncattr('long_name', long_name)
    nc_var.setncattr('units', units)
    nc_var.setncattr('scale_factor', 1.)
    nc_var.setncattr('missing_value', -9999.)
    nc_var.setncattr('coordinates', 'lat lon')


def calc_monthly_means(t, Q):
    """Calculates the monthly mean of Q

    Args:
        t (1d ndarray): the time
        Q (3d ndarray): the discharge on the mRM grid
    """
    dates = pd.Series(t, name='time')
    ds = xr.Dataset({'Q': (['time', 'y', 'x'], Q)},
                    coords={'time': t})
    # ds_mon = ds.resample('M', dim='time')
    ds_mon = ds.resample(time='1M').mean()
    return ds_mon['Q']

def calc_Q_bkfl( Q, return_period):
    """Calculates the discharge at bankfull conditions for a single time series

    Args:
        t (1d ndarray): the time
        Q (1d ndarray): the discharge
        return_period (float, opt.): the return period of bankfull conditions
    """
    # exceedance probability
    ex_prob = np.linspace(0, 1, len(Q), endpoint=False)
    # empirical CDF
    Q_sort = np.sort(Q)[::-1]
    # plotting Q_sort against ex_prob gives the exceedance probability
    # pt.plot(Q_sort, ex_prob)
    # X-year flood is defined as a flood which has a
    # 1/x% chance to occur during a year
    idx_bkfl = find_nearest_idx(ex_prob, 1/return_period)
    return Q_sort[idx_bkfl]

def process_grid(Q, return_period):
    """Calculates the discharge at bankfull conditions for a complete grid

    Args:
        t (1d ndarray): the time
        Q (3d ndarray): the discharge on the mRM grid
        return_period (float, opt.): the return period in years
    """
    Q_bkfl = np.ma.empty(Q.shape[1:], Q.dtype)
    for i in range(Q.shape[1]):
        for j in range(Q.shape[2]):
            if not np.all(Q[:,i,j]):
                Q_bkfl[i,j] = calc_Q_bkfl(Q[:,i,j], return_period)
    return Q_bkfl


if __name__ == '__main__':
    ncin_filename, ncout_filename, return_period, peri_bkfl = get_cmdline_args(1.5)
    t, Q = read_discharge(ncin_filename)
    Q_mon = calc_monthly_means(t, Q)

    Q_bkfl = process_grid(Q_mon, return_period)
    write_Q_bkfl(Q_bkfl, ncin_filename, ncout_filename, peri_bkfl)
