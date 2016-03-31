#! /usr/bin/env python
# -*- coding: utf-8 -*-
#===========================================================================================
# Created by   David Schafer (david.schaefer@ufz.de) on 23/03/2016
# Commented by Rohini Kumar (rohini.kumar@ufz.de)    on 23/03/2016
#------------------------------------------------------------------------
#
# HARGREAVES AND SAMANI (1985) BASED DAILY PET ESTIMATES
# Requires: 3D tavg, tmax, and tmin [all in degree celcius] and
#           2D latitude grid [degree decimals] 
# Output will be written in highest precision among all input datasets
#
# To execute this script, you need to activate 
# 1) our CHS python library scripts
#    by e.g., export PYTHONPATH=../PYTHON_chs_lib/:$PYTHONPATH
# 2) virtual C-binding virtual environment
#    by e.g., module load /global/apps/chs-virtualenv/chspython/2.7.6
#
# After loding these environment run command
#  by default the variable in the 'avg. temp. file' should be named as 'tavg'
#                          in the 'max. temp. file' should be named as 'tmax'
#                          in the 'min. temp. file' should be named as 'tmin'
#                          in the 'lat        file' should be named as 'lat'
# python hargreaves_samani_PET.py -g tavg.nc -n tmin.nc -x tmax.nc -l lat.nc pet.nc
#
# another way is to explicitly spell out the variable names
# python hargreaves_samani_PET.py -g tavg.nc,tavg -n tmin.nc,tmin -x tmax.nc,tmin -l lat.nc,latlon.nc pet.nc
#
#===========================================================================================
import numpy as np
from ufz.netcdf4 import NcDataset
from argparse import ArgumentParser

VARNAMES = {
    "tavg" : "tavg",
    "tmin" : "tmin",
    "tmax" : "tmax",
    "lat"  : "lat"
}

parser = ArgumentParser(description="Calculate PET according to HARGREAVES",
                        usage="calc_pet.py [-h] -g TAVGFILE[,varname] -n TMINFILE[,varname] -x TMAXFILE[,varname] -l LATFILE[,varname] outfile")
parser.add_argument("-g","--tavg", dest="tavg", required=True,
                    help="NetCDF file holding average temperature data. The variable name should either be 'tavg' or given explictly.")
parser.add_argument("-n","--tmin", dest="tmin", required=True,
                    help="NetCDF file holding minimum temperature data. The variable name should either be 'tmin' or given explictly.")
parser.add_argument("-x","--tmax", dest="tmax", required=True,
                    help="NetCDF file holding maximum temperature data. The variable name should either be 'tmax' or given explictly.")
parser.add_argument("-l","--lat", dest="lat", required=True,
                    help="NetCDF file holding latitude data. The variable name should either be 'lat' or given explictly.")
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
    For large argument arrays it is however well possible to run into MemoryErrors.
    Split your problem along the timeaxis of tas and co. in order to trade speed
    for results.
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

