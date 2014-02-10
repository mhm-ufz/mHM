#!/usr/bin/env python
#
# purpose: conversion of binary files to netcdf files
#
# written: Juliane Mai & Matthias Zink Feb 2014
#
from __future__ import print_function

import ufz
import numpy   as np
import netCDF4 as nc  
import struct  as struct

# -------------------------------------------------------------------------
# Command line arguments
# -------------------------------------------------------------------------

headerfile         = '../test/input/meteo/pre/header.txt'
outfile            = '../test/input/meteo/pre/out.nc'
indir              = '../test/input/meteo/pre/'
years              = '1980,1981'
#
latlonfile         = '../test/input/latlon/latlon.nc'
variable           = 'pre'
variable_long_name = 'daily sum of precipitation'
variable_unit      = 'mm d-1'

import optparse
parser = optparse.OptionParser(usage='%prog [options]',
                               description="Converts binary input files for mHM within <years> stored under <indir> with header information from <headerfile> into NetCDF file <outfile>.")

# usage example with command line arguments
# -------------------------------------------------------------------------
# python bin2nc.py -f '../test/input/meteo/tavg/header.txt' -o '../test/input/meteo/tavg/out.nc' -i '../test/input/meteo/tavg/' \
#                  -y '1980,1981' -c '../test/input/latlon/latlon.nc' -v 'tavg' -l 'daily mean air temperature' -u 'degC'
# -------------------------------------------------------------------------

parser.add_option('-f', '--header', action='store', dest='headerfile', type='string',
                  default=headerfile, metavar='Header file.',
                  help='Header file containing information about e.g. number of rows and columns. (default: header.txt).')
parser.add_option('-o', '--outfile', action='store', dest='outfile', type='string',
                  default=outfile, metavar='NetCDF file.',
                  help='Name of NetCDF file. (default: outfile=../test/input/meteo/out.nc).')
parser.add_option('-i', '--indir', action='store', dest='indir', type='string',
                  default=indir, metavar='Directory containing *.bin files',
                  help='Directory containing annual files [YYYY].bin (default: ../test/input/meteo/pre/).')
parser.add_option('-y', '--years', action='store', dest='years', type='string',
                  default=years, metavar='Yearly files to be converted.',
                  help='Yearly files to be converted via -y <from year>,<to year> (e.g. defualt= -y 1980,1981).')
parser.add_option('-c', '--coordinates', action='store', dest='latlonfile', type='string',
                  default=latlonfile, metavar='NetCDF containing latitude and longitide grid.',
                  help='NetCDF containing latitude and longitide grid. <You may use create_latlon.py>. ' +
                        'defualt= ../test/input/latlon/latlon.nc')
parser.add_option('-v', '--varname', action='store', dest='variable', type='string',
                  default=variable, metavar='Variable name in NetCDF (e.g. pre, tavg, pet)',
                  help='Variable name used in NetCDF (e.g. pre, tavg, pet).')
parser.add_option('-l', '--varnamelong', action='store', dest='variable_long_name', type='string',
                  default=variable_long_name, metavar='Attribute long_name of variable in NetCDF (e.g. <daily sum of precipitation>).',
                  help='Attribute long_name of variable in NetCDF (e.g. <daily sum of precipitation>).')
parser.add_option('-u', '--variable_unit', action='store', dest='variable_unit', type='string',
                  default=variable_unit, metavar='The attribute unit of variable in NetCDF (e.g. <mm d-1>).',
                  help='The attribute unit of variable in NetCDF (e.g. <mm d-1>).')

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
info         = ufz.fread(headerfile,cskip=1)

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
thand   = ufz.writenetcdf(fhandle, name=varName, dims=dims, attributes=varAtt, isdim=True, vartype='i4')
#
varAtt  = ([['axis'     , 'X']])
varName = 'xc'
dims    = ncols
var     = np.arange(xllcorner, xllcorner + ncols * cellsize, cellsize)
ufz.writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, isdim=True)
#
varAtt  = ([['axis'     , 'Y']])
varName = 'yc'
dims    = nrows
var     = np.arange(yllcorner + nrows * cellsize, yllcorner, -cellsize)
ufz.writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, isdim=True)
# Variables
# lon
varAtt  = ([['units'         , 'degrees_east'],
            ['long_name'     , 'longitude'   ]])
#
varName = 'lon'
var     = ufz.readnetcdf(latlonfile,var=varName)
dims    = ['yc','xc']
ufz.writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, comp=True)
# lat
varAtt  = ([['units'        , 'degrees_north'],
            ['long_name'     , 'latitude'    ]])
varName = 'lat'
var     = ufz.readnetcdf(latlonfile,var=varName)
dims    = ['yc','xc']
ufz.writenetcdf(fhandle, name=varName, dims=dims, var=var, attributes=varAtt, comp=True)
#
# create variable and write attributes of it to NetCDF file
varAtt  = ([['long_name'     , variable_long_name           ],
            ['units'         , variable_unit                ],
            ['missing_value' , float(info[5][0])            ],
            ['fill_value'    , float(info[5][0])            ],
            ['coordinates'   , 'lon lat'                    ]])
varName = variable
#var     = mask_upscale
dims    = ['time','yc','xc']
vhandle = ufz.writenetcdf(fhandle, name=varName, dims=dims, attributes=varAtt, vartype='f8',comp=True)
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
  values     = np.array(struct.unpack(ufz.astr(ncols*nrows*days_year)+'f',bindata[0:4*ncols*nrows*days_year]))
  # create daily fields and roll axis, because for netcdf writer first dimension has to be time
  var        = np.rollaxis(np.reshape(values, (nrows, ncols, days_year), order='Fortran'),2,0)
  from_day   = to_day
  to_day     = to_day + days_year
  # write data to nc file
  ufz.writenetcdf(fhandle, vhandle, var=var, time=np.arange(from_day, to_day))
  # write time steps to nc
  times      = np.arange(from_day, to_day)
  ufz.writenetcdf(fhandle, thand  , time=np.arange(times[0], times[-1]+1), var=times)
# close netcdf
fhandle.close()

