#!/usr/bin/env python
#
# purpose: conversion of binary files to netcdf files
#
# written: Juliane Mai & Matthias Zink Feb 2014
#
from __future__ import print_function

import numpy   as np
import netCDF4 as nc  
import struct
from writenetcdf import writenetcdf # from ufz
from readnetcdf  import readnetcdf  # from ufz
from fread       import fread       # from ufz
from autostring  import astr        # from ufz

# -------------------------------------------------------------------------
# Command line arguments
# -------------------------------------------------------------------------

headerfile         = '../test_basin/input/meteo/pre/header.txt'
outfile            = '../test_basin/input/meteo/pre/out.nc'
indir              = '../test_basin/input/meteo/pre/'
years              = '1989,1993'
#
latlonfile         = '../test_basin/input/latlon/latlon.nc'
variable           = 'pre'
variable_long_name = 'daily sum of precipitation'
variable_unit      = 'mm d-1'

import optparse
parser = optparse.OptionParser(usage='%prog [options]',
                               description="Converts binary input files for mHM within <years> stored under <indir> with header information from <headerfile> into NetCDF file <outfile>.")

# usage example with command line arguments
# -------------------------------------------------------------------------
# Precipitation:
#
# python bin2nc.py -f ../test_basin/input/meteo/pre/header.txt -o ../test_basin/input/meteo/pre/pre.nc -i ../test_basin/input/meteo/pre/ -y 1989,1993 -c ../test_basin/input/latlon/latlon.nc -v pre -l 'daily sum of precipitation' -u 'mm d-1'
#
# Temperature
#
# python bin2nc.py -f '../test_basin/input/meteo/tavg/header.txt' -o '../test_basin/input/meteo/tavg/out.nc' -i '../test_basin/input/meteo/tavg/' -y '1989,1993' -c '../test_basin/input/latlon/latlon.nc' -v 'tavg' -l 'daily mean air temperature' -u 'degC'
#
# PET
#
# python bin2nc.py -f '../test_basin/input/meteo/pet/header.txt' -o '../test_basin/input/meteo/pet/out.nc' -i '../test_basin/input/meteo/pet/' -y '1989,1993' -c '../test_basin/input/latlon/latlon.nc' -v 'pet' -l 'daily sum of potential evapotranspiration' -u 'mm d-1'
# -------------------------------------------------------------------------

parser.add_option('-f', '--header', action='store', dest='headerfile', type='string',
                  default=headerfile, metavar='Header file.',
                  help='Header file containing information about e.g. number of rows and columns. (default: header.txt).')
parser.add_option('-o', '--outfile', action='store', dest='outfile', type='string',
                  default=outfile, metavar='NetCDF file.',
                  help='Name of NetCDF file. (default: outfile=../test_basin/input/meteo/out.nc).')
parser.add_option('-i', '--indir', action='store', dest='indir', type='string',
                  default=indir, metavar='Directory of *.bin files.',
                  help='Directory containing annual files [YYYY].bin (default: ../test_basin/input/meteo/pre/).')
parser.add_option('-y', '--years', action='store', dest='years', type='string',
                  default=years, metavar='Years of *.bin files.',
                  help='Yearly files to be converted via -y <from year>,<to year> (default: "1989,1993").')
parser.add_option('-c', '--coordinates', action='store', dest='latlonfile', type='string',
                  default=latlonfile, metavar='Latlon file.',
                  help='NetCDF containing latitude and longitide grid. You may use create_latlon.py. ' +
                        '(default: ../test_basin/input/latlon/latlon.nc)')
parser.add_option('-v', '--varname', action='store', dest='variable', type='string',
                  default=variable, metavar='NetCDF Variable.',
                  help='Variable name used in NetCDF and has to be exactly one of {pre, tavg, pet} (default: pre).')
parser.add_option('-l', '--varnamelong', action='store', dest='variable_long_name', type='string',
                  default=variable_long_name, metavar='Long_name in NetCDF.',
                  help='Attribute long_name of variable in NetCDF (default: "daily sum of precipitation").')
parser.add_option('-u', '--variable_unit', action='store', dest='variable_unit', type='string',
                  default=variable_unit, metavar='Unit in NetCDF.',
                  help='The attribute unit of variable in NetCDF (default: "mm d-1").')

(opts, args) = parser.parse_args()

headerfile         = opts.headerfile
outfile            = opts.outfile   
indir              = opts.indir 
years              = opts.years
#
latlonfile         = opts.latlonfile
variable           = opts.variable
variable_long_name = opts.variable_long_name
variable_unit      = opts.variable_unit

del parser, opts, args

# -------------------------------------------------------------------------
# Read header file
# -------------------------------------------------------------------------

# info = (ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)
info         = fread(headerfile,cskip=1)

ncols        = int(info[0][0])
nrows        = int(info[1][0])
xllcorner    = info[2][0]
yllcorner    = info[3][0]
cellsize     = int(info[4][0])
NODATA_value = info[5][0]

# init some things
startyear    = int(years.split(',')[0])
endyear      = int(years.split(',')[1])
times        = []
# -------------------------------------------------------------------------
# Open and set NetCDF file
# -------------------------------------------------------------------------
#
# Dimensions
fhandle      = nc.Dataset(outfile, 'w', format='NETCDF4')
startTime    = ('days since '  + str(startyear) + '-' + str(1).zfill(2) + '-' + str(1).zfill(2) + ' ' 
                  + str(0).zfill(2) + ':' + str(0).zfill(2) + ':' + str(0).zfill(2))   
varName = 'time'
dims    = None
varAtt  = ([['units'    , startTime],
            ['calendar' , 'standard']])
thand   = writenetcdf(fhandle, name=varName, dims=dims, attributes=varAtt, isdim=True, vartype='i4')
#
varAtt  = ([['axis'     , 'X']])
varName = 'xc'
dims    = ncols
var     = np.arange(xllcorner, xllcorner + ncols * cellsize, cellsize)
writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, isdim=True)
#
varAtt  = ([['axis'     , 'Y']])
varName = 'yc'
dims    = nrows
var     = np.arange(yllcorner + nrows * cellsize, yllcorner, -cellsize)
writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, isdim=True)
# Variables
# lon
varAtt  = ([['units'         , 'degrees_east'],
            ['long_name'     , 'longitude'   ]])
#
varName = 'lon'
var     = readnetcdf(latlonfile,var=varName)
dims    = ['yc','xc']
writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, comp=True)
# lat
varAtt  = ([['units'        , 'degrees_north'],
            ['long_name'     , 'latitude'    ]])
varName = 'lat'
var     = readnetcdf(latlonfile,var=varName)
dims    = ['yc','xc']
writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, comp=True)
#
# create variable and write attributes of it to NetCDF file
varAtt  = ([['long_name'     , variable_long_name           ],
            ['units'         , variable_unit                ],
            ['_FillValue'    , float(info[5][0])            ],
            ['coordinates'   , 'lon lat'                    ]])
varName = variable
#var     = mask_upscale
dims    = ['time','yc','xc']
vhandle = writenetcdf(fhandle, name=varName, dims=dims, attributes=varAtt, vartype='f8',comp=True)
# -------------------------------------------------------------------------
# Read binary files
# -------------------------------------------------------------------------

# Format        C Type          Python type     Standard size   
#   i           int             integer         4               
#   I           unsigned int    integer         4               
#   l           long            integer         4               
#   L           unsigned long   integer         4               
#   f           float           float           4               
#   d           double          float           8               
#   s           char[]          string          
from_day   = 0
to_day     = 0
#
for iyear in range(startyear, endyear+1):
  leap       = (((iyear % 4) == 0) & ((iyear % 100) != 0)) | ((iyear % 400) == 0)
  days_year  = 365 + leap
  bindata    = open(indir + str(iyear) + '.bin', "rb").read()
  values     = np.array(struct.unpack(astr(ncols*nrows*days_year)+'f',bindata[0:4*ncols*nrows*days_year]))
  # create daily fields and roll axis, because for netcdf writer first dimension has to be time
  var        = np.rollaxis(np.reshape(values, (nrows, ncols, days_year), order='Fortran'),2,0)
  from_day   = to_day
  to_day     = to_day + days_year
  # write data to nc file
  writenetcdf(fhandle, vhandle, var=var, time=np.arange(from_day, to_day))
  # write time steps to nc
  times      = np.arange(from_day, to_day)
  writenetcdf(fhandle, thand  , time=np.arange(times[0], times[-1]+1), var=times)
# close netcdf
fhandle.close()

