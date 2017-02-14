#!/usr/bin/env python
#############################################
# create restart file from mHM output
#
# created by Stephan Thober 17.04.2014
#
#############################################
#

# INIT FILE
#   contains state variables as they are
#   obtained from mHM output
states_file = "mHM_states.nc"

# DATE
#   date at which states shall be extracted from states_file
date = ''

# RESTART FILE
#   reference restart file that has 
#   been obtained for this basin
restart_file = "ref_states.nc"

# OUTPUT FILE
#   path to the output file, latlon.nc is hard-coded in mHM
out_file = 'test.nc'

import optparse
parser = optparse.OptionParser(usage='%prog [options]',
                               description="Creates a restart file <outfile> given the states of an mHM output file <States File>, a reference restart file  <Restart File>")

# usage example with command line arguments
# -------------------------------------------------------------------------
#
# python create_restart_from_output.py -s mHM_Fluxes_States.nc -r 001_states.nc -o states.nc
#
# -------------------------------------------------------------------------

parser.add_option('-s', '--states', action='store', dest='states_file', type='string',
                  default=states_file, metavar='States File',
                  help='States as obtained from mHM output file')
parser.add_option('-d', '--date', action='store', dest='date', type='string',
                  default=date, metavar='Date',
                  help='Date as YYYY-MM-DD HH at which states shall be extracted from States File.')
parser.add_option('-r', '--restart', action='store', dest='restart_file', type='string',
                  default=restart_file, metavar='Restart File',
                  help='Restart File obtained for exactly the same domain as the States File.')
parser.add_option('-o', '--outfile', action='store', dest='out_file', type='string',
                  default=out_file, metavar='Output File',
                  help='Name of NetCDF file. (default: test.nc).')

(opts, args) = parser.parse_args()

states_file  = opts.states_file
out_file     = opts.out_file   
restart_file = opts.restart_file
date         = opts.date

#############################################

import numpy as np                     # array manipulation
import os                              # call current time for timestamp
from date2dec import date2dec          # from ufz
from readnetcdf import readnetcdf      # from ufz

#############################################

# translate date to time slice index
if date != '':
    # check whether ref dates is in hours
    if (readnetcdf( states_file, 'time', attributes = True)['units'][:5] == 'hours') and (np.diff(readnetcdf( states_file, 'time'))[0] != 1):
        print( '***ERROR: time axis in states file is not in hours' )
        exit
    # get reference time from states_file
    ref_date = readnetcdf( states_file, 'time', attributes = True)['units'].split(' ')[2]
    ref_year = np.int( ref_date.split('-')[0])
    ref_mo   = np.int( ref_date.split('-')[1])
    ref_day  = np.int( ref_date.split('-')[2])
    ref_julian = date2dec( yr = ref_year, mo = ref_mo, dy = ref_day)
    # get the date of the first time step in states_file
    sta_julian = readnetcdf( states_file, 'time')[0]
    if readnetcdf( states_file, 'time', attributes = True)['units'].split(' ')[0] == 'hours':
        # transform to days
        sta_julian = (sta_julian + 1.) / 24.
    # add julian days up to reference date
    ref_julian = ref_julian + sta_julian - 1

    # get julian date of actual time step
    ref_year = np.int( date.split(' ')[0].split('-')[0])
    ref_mo   = np.int( date.split(' ')[0].split('-')[1])
    ref_day  = np.int( date.split(' ')[0].split('-')[2])
    ref_hour = np.int( date.split(' ')[1])
    act_julian = date2dec( yr = ref_year, mo = ref_mo, dy = ref_day)
    dd = np.int(act_julian - ref_julian)
    tt = dd * 24 + ref_hour
else:
    tt = 0

print( 'selecting timestep: '+ str(tt) )
#############################################

# name of states in states_file
in_states = [ 'interception', 'snowpack', 'SWC', 'sealedSTW', 'unsatSTW', 'satSTW' ]

# name of states in restart file
out_states = [ 'L1_Inter', 'L1_snowPack', 'L1_soilMoist', 'L1_sealSTW', 'L1_unsatSTW', 'L1_satSTW' ]

#############################################
# COPY THE STATES ###########################
#############################################

# copy restart file to out_file
os.system( 'cp ' + restart_file + ' ' + out_file )

# overwrite states in outfile
for ii in np.arange( len(in_states) ):

    print 'processing: ' + in_states[ii]

    if in_states[ii] == 'SWC':

        # read out_states first to obtain number of soil layers
        fo, r_arr = readnetcdf( out_file, out_states[ii], overwrite = True )
        
        # get number of soil layers
        n_soil = r_arr.shape[0]

        for nn in np.arange( n_soil ):
            s_name = in_states[ii]+'_L'+str(nn+1).zfill(2)
            fi, s_arr  = readnetcdf( states_file, s_name, pointer = True )
            r_arr[nn,:,:] = s_arr[tt,:,:]
            fi.close()

    else:
        # read states from mHM
        fi, s_arr = readnetcdf( states_file, in_states[ii], pointer = True )
    
        # copy states to out_file
        fo, r_arr = readnetcdf( out_file, out_states[ii], overwrite = True )

        r_arr[:] = s_arr[tt,:,:]

        # close input file
        fi.close()
        
    # close output file
    fo.close()
