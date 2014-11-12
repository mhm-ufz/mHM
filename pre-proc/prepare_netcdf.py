#!/usr/bin/env python
from __future__ import print_function
"""


History
-------
Written  Matthias Cuntz & Juliane Mai Nov 2014 - write a netcdf file in L0
"""

# -------------------------------------------------------------------------
# Command line arguments
#

import argparse

addargs = []
parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description=('''!!!.'''))
args    = parser.parse_args()
del parser, args

# import packages needed after help so that help with command line -h is fast
import numpy as np
import netCDF4 as nc

# -------------------------------------------------------------------------
# Prepare data
#

# test_basin
ncol = 288
nrow = 432

# # test_basin_2
# ncol = 240
# nrow = 240

ntime = 1826
timeunit = "days since 1989-01-01 00:00:00"
varname = 'lai'
longname= 'LAI from Modis'

dat = np.arange(nrow*ncol*ntime).reshape((ntime, nrow, ncol))
dat = np.sin(dat)*3 + 3.

# -------------------------------------------------------------------------
# Write netcdf file
#

# NETCDF - out
filename = 'lai.nc'
print('Create netcdf file ', filename)
f = nc.Dataset(filename, 'w', format='NETCDF4')

# Structure
x      = f.createDimension('x', ncol)
y      = f.createDimension('y', nrow)
time   = f.createDimension('time', None)

lai                = f.createVariable(varname, 'f8', ('time','y','x',), fill_value=-9999., zlib=True, least_significant_digit=3)
lai.long_name      = longname

time               = f.createVariable('time', 'i4', ('time',))
time.units         = timeunit
time.calendar      = "standard"

f.Production = '2014/11/12'

# Data
lai[:,:,:]      = dat
time[:]         = np.arange(ntime)

f.close()

# -------------------------------------------------------------------------
# Finish
