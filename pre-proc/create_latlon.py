#!/usr/bin/env python
#############################################
# calculate latitude and longitude coordinates
# given header file specifications
#
# created by Matthias Zink 04.04.2012
# modified Stephan Thober, Nov 2013 - mHM confirmative
#
#############################################
#
# COORDINATE SYSTEM
#   code 31463 defines GK4 (DHDN3 Zone 3)
#   [website: http://www.spatialreference.org/ref/epsg/31463/]
#   equal +proj=tmerc +ellps=bessel +lon_0=12 +x_0=3,500,000 +y_0=0
coord_sys = 'epsg:31463'

# HEADER FILE
#   specifies the grid properties
#   for example, use a copy of the input/meteo/pre/header.txt
#   and adapt cellsize, ncols, nrows to your hydrologic resolution
headerfile = 'header.txt'

# OUTPUT FILE
#   path to the output file, latlon.nc is hard-coded in mHM
latlonfile = 'latlon.nc'

#############################################

import numpy as np                       # array manipulation
import netCDF4 as nc                     # netCDF interphase
import ufz
import time, os, sys                     # call current time for timestamp
from pyproj      import Proj
from writenetcdf import *

# check input files

if not os.path.isfile(headerfile):
  sys.exit("No "+headerfile+" file found here, are you in the right directory?")

# read header information
header_info = np.loadtxt( headerfile, dtype='|S20')
ncols       = np.int(header_info[0,1])
nrows       = np.int(header_info[1,1])
xllcorner   = np.float(header_info[2,1])
yllcorner   = np.float(header_info[3,1])
cs          = np.int(header_info[4,1])
missVal     = header_info[5,1]

# create x and y grid
xx          = np.arange( xllcorner + cs/2, xllcorner + cs/2 + ncols*cs, cs)
yy          = np.arange( yllcorner + cs/2, yllcorner + cs/2 + nrows*cs, cs)
xx, yy      = np.meshgrid(xx,yy)

#
# determine latitude and longitude of the Aimgrid
projAim    = Proj(init=coord_sys)
lons, lats = projAim(xx, yy, inverse=True)
#
# write netCDF
fName   = latlonfile
fhandle = nc.Dataset(fName, 'w', format='NETCDF4')
#
FiAtt  = ([['description', 'lat lon file'],
           ['history','Created ' + time.ctime(time.time()) ]])
writenetcdf(fhandle, fileattributes=FiAtt)
#
varAtt  = ([['axis'     , 'X']])
varName = 'xc'
dims    = ncols
var     = xx[1,:] #np.arange(ncols)+1
writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, isdim=True)
#
varAtt  = ([['axis'     , 'Y']])
varName = 'yc'
dims    = nrows
var     = yy[:,1] #np.arange(nrows)+1
writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, isdim=True)
# Variables
# lon
varAtt  = ([['units'        , 'degrees_east'],
           ['long_name'     , 'longitude'],
           ['missing_value' , missVal],
           ['fill_value'    , missVal]])
#
varName = 'lon'
var     = lons
dims    = ['yc','xc']
writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, vartype='f8', comp=True)
# lat
varAtt  = ([['units'        , 'degrees_north'],
           ['long_name'     , 'latitude'],
           ['missing_value' , missVal],
           ['fill_value'    , missVal]])
varName = 'lat'
var     = lats
dims    = ['yc','xc']
writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, vartype='f8', comp=True)
#
fhandle.close()
print latlonfile + " created."
