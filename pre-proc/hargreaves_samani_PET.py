#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""

Author
------
David Schaefer

History
-------
- Rohini Kumar, Mar 2016 - Added comments

Purpose
-------
HARGREAVES and SAMANI (1985) based PET estimates

Arguments
---------
- tavg     (3D) [degree Celsius]
- tmax     (3D) [degree Celsius]
- tmin     (3D) [degree Celsius]
- latitude (2D) [degree decimals]

all stored in one or several NetCDF file(s). The input data files are
expected to hold a time variable following the CF Conventions
(http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#time-coordinate)

Output
------
- pet (3D) [mm/timestep]

written into a NetCDF file.

Requirements
------------
numpy
netCDF4
argparse
netcdf4

Usage
-----
1) Satisfy all dependencies
2) Adapt your PYTHONPATH to include the module netcdf4 from PYTHON_chs_lib.
   Skip this step if you received this program as a part of mHM.
3) Run the program using one of the following possibilities
3.1) python hargreaves_samani_PET.py -g tavg.nc -n tmin.nc -x tmax.nc -l lat.nc pet.nc
3.2) python hargreaves_samani_PET.py --tavg tavg.nc --tmin tmin.nc --tmax tmax.nc --lat lat.nc pet.nc
3.3) python hargreaves_samani_PET.py --tavg tavg.nc,tavg --tmin tmin.nc,tmin --tmax tmax.nc,tmax --lat lat.nc,lat pet.nc

Note
----
When running the program as in 3.1) and 3.2) the NetCDF variables holding the data are expected to be named:
- tavg: average temperature
- tmin: minimum temperature
- tmax: maximum temperature
- lat : latitude

If the input has a differing naming scheme, 3.3) gives the possibility to pass the variables names
explicitly, after the filename and separated by comma (no additional whitespaces allowed)

"""

import numpy as np
from netcdf4 import NcDataset
from argparse import ArgumentParser

VARNAMES = {
    "tavg" : "tavg",
    "tmin" : "tmin",
    "tmax" : "tmax",
    "lat"  : "lat"
}

parser = ArgumentParser(
    description="Calculate PET according to HARGREAVES",
    usage="calc_pet.py [-h] -g TAVGFILE[,varname] -n TMINFILE[,varname] -x TMAXFILE[,varname] -l LATFILE[,varname] outfile"
)
parser.add_argument(
    "-g","--tavg", dest="tavg", required=True,
    help="NetCDF file holding average temperature data. The variable name should either be 'tavg' or given explicitly."
)
parser.add_argument(
    "-n","--tmin", dest="tmin", required=True,
    help="NetCDF file holding minimum temperature data. The variable name should either be 'tmin' or given explicitly."
)
parser.add_argument(
    "-x","--tmax", dest="tmax", required=True,
    help="NetCDF file holding maximum temperature data. The variable name should either be 'tmax' or given explicitly."
)
parser.add_argument(
    "-l","--lat", dest="lat", required=True,
    help="NetCDF file holding latitude data. The variable name should either be 'lat' or given explicitly."
)
parser.add_argument('outfile', help="Output file")


def epotHargreaves(tavg, tmin, tmax, lat, julian_dates):
    """
    Arguments
    ---------
    tavg         : np.ndarray(3D) -> average tmperature
    tmin         : np.ndarray(3D) -> minimum temperature
    tmax         : np.ndarray(3D) -> maximum tmperature
    lat          : np.ndarray(2D) -> latitude 
    julian_dates : np.ndarray(1D) -> day of the year for every timestep in tas/tasmin/tasmax

    Returns
    -------
    ndarray(3D)

    Purpose
    -------
    Calculates potential evapotranspiration in [mm/timestep] according to HARGREAVES.           

    Note
    ----
    For large argument arrays it is well possible to run into MemoryErrors.
    Split your problem along the timeaxis of the input arguments, in order to trade
    speed for results.
    """
    
    lat = np.radians(np.abs(lat[None,:]))
    julian_dates = julian_dates[:,None,None]
    # relative distance between sun and earth
    dr = 1 + 0.033 * (np.cos(((2*np.pi*julian_dates) / 365.25)))
    # solar inclination
    delta = 0.4093 * np.sin( 2*np.pi*julian_dates/365.25 - 1.39)
    # sunset hour angle [degree]
    ws = np.arccos(np.clip(-np.tan(lat) * np.tan(delta), -1, 1))
    # extraterrestrial radiation (DUFFIE and BECKMAN, 1980)
    ra = 15.3351 * dr * (ws * np.sin(lat) * np.sin(delta) + np.cos(lat) * np.cos(delta) * np.sin(ws))
    pet = 0.0023 * ra * (tavg + 17.8) * np.sqrt(np.maximum(tmax-tmin,0))
    return np.maximum(pet, 0)

def getData(files, varmap):
    out = {} 
    for fname in files:
        with NcDataset(fname,"r") as nc:
            for k, v in varmap.items():
                if v in nc.variables:
                    out[k] = np.squeeze(nc.variables[v][:])
    return out

def getJulianDates(fname):
    with NcDataset(fname, "r") as nc:
        dates = nc.getDates()
        try:
            out = np.array([d.timetuple().tm_yday for d in dates])
        except AttributeError:
            out = np.array([d.timetuple()[-2] for d in dates])
    return out

def parseArgument(arg, sep=","):
    return ([e.strip() for e in arg.split(sep)] + [None, ])[:2]

if __name__== "__main__":

    args = parser.parse_args()

    tavgfile, tavg = parseArgument(args.tavg)
    tminfile, tmin = parseArgument(args.tmin)
    tmaxfile, tmax = parseArgument(args.tmax)
    latfile, lat   = parseArgument(args.lat)
    
    varmap = {
        "tavg" : tavg if tavg else VARNAMES["tavg"],
        "tmin" : tmin if tmin else VARNAMES["tmin"],
        "tmax" : tmax if tmax else VARNAMES["tmax"],
        "lat"  : lat  if lat  else VARNAMES["lat"],
    }

    pet = epotHargreaves(
        julian_dates=getJulianDates(tavgfile),
        **getData([tavgfile, tminfile, tmaxfile, latfile], varmap)
    )

    with NcDataset(tavgfile, "r") as ncin:
        with NcDataset(args.outfile, "w") as ncout:
            basevar = varmap["tavg"]
            ncout.copyDataset(ncin, skipvars=(basevar, "time_bnds"))
            var = ncout.createVariable(
                "pet",
                pet.dtype,
                ncin.filterDimensions(pet.shape).keys(),
                fill_value=ncin.variables[basevar].fill_value,
            )
            var[:] = pet
            var.createAttributes({
                "units"       : "mm/timestep",
                "long_name"   : "daily potential evapotranspiration",
            })

