#!/usr/bin/env python
#############################################
# cut mhm input from gauge location in in river network
# given header file specifications
#
# created by Stephan Thober 17.02.2016
#
#############################################

from __future__ import print_function
import optparse
parser = optparse.OptionParser(usage='%prog [options]',
                               description="Cretaes latitude and longitude grids for <coord_sys> with the domain defined in  <headerfile> into NetCDF file <outfile>.")

# usage example with command line arguments
# -------------------------------------------------------------------------
#
# py create_latlon.py -c 'epsg:31463' -f header.txt -o latlon.nc
#
# -------------------------------------------------------------------------

def check_dir(dir):
    if not isdir(dir):
        print(' create directory: ' + dir)
        makedirs(dir)
    

# default variables
indir = '../test_basin/input/'
outdir = './test_cut_mhm_input/'
gauge_id = -9999.

parser.add_option('-i', '--indir', action='store', dest='indir', type='string',
                  default=indir, metavar='Property',
                  help='input directory containing mhm gauge, land use, and morphological information')
parser.add_option('-o', '--outdir', action='store', dest='outdir', type='string',
                  default=outdir, metavar='Property',
                  help='directory where output should be written')
parser.add_option('-g', '--gauge_id', action='store', dest='gauge_id', type='float',
                  default=gauge_id, metavar='Property',
                  help='gauge_id from where subbasin should be cutted')

(opts, args) = parser.parse_args()

indir = opts.indir
outdir = opts.outdir
gauge_id = opts.gauge_id

if indir[-1] != '/':
    indir += '/'
if outdir[-1] != '/':
    outdir += '/'

del opts, args
#############################################
## START PROCESSING
#############################################
## mHM files to process
files = ['gauge/gauge_loc.asc',
         'morph/aspect.asc',
         'morph/idgauges.asc',
         'morph/dem.asc',
         'morph/slope.asc',
         'morph/facc.asc',
         'morph/fdir.asc',
         'luse/lc_2001.asc']
#############################################


import numpy as np # array manipulation
from os.path import isdir
from os import makedirs
from ufz import fread, river_network, fwrite # from ufz

# read gauge_id file
idgauges = fread(indir + '/gauge/gauge_loc.asc', skip=6)
if gauge_id == -9999.:
    gauges_loc = np.where(idgauges > 0.)
    gauges_id = idgauges[gauges_loc]
else:
    gauges_loc = np.where(idgauges == gauge_id)
    gauges_id = gauge_id
    if gauges_loc[0].shape[0] == 0:
        raise ValueError('***ERROR: gauge with ID {:2.0f} does not exist'.format(gauge_id))    

# read header
header = fread(indir + '/morph/fdir.asc', nc=2, skip=6, header=True)    
# read flow direction of entire basin
fdir = fread(indir + '/morph/fdir.asc', skip=6)

for gg in np.arange(len(gauges_id)):
    print(' process gauge with gauge id ' + str(gauges_id[gg]))
    # recalculate flow accumulation from gauge location
    sn = river_network(fdir=fdir, do_fa=True, sinks=(np.array([gauges_loc[0][gg]]), np.array([gauges_loc[1][gg]])))
    # set mask
    mask = sn.fa.mask

    # check whether output directory exist
    outpath = outdir + str(gauges_id[gg]) + '/'
    check_dir(outpath)

    for ff in files:
        print('   process ' + ff)
        # read
        arr = fread(indir + ff, skip=6)
        # set mask
        arr[mask] = -9999
        # check dir
        check_dir(outpath + '/'.join(ff.split('/')[:-1]))
        # write
        fwrite(outpath + ff, arr, header=header)

    # copy gauge file
    gfile = '{:05.0f}.txt'.format(gauges_id[gg])
    print('   you have to copy gauge file ' + gfile + ' manually using')
    print('   cp ' + indir + 'gauge/' + gfile + ' ' + outpath + 'gauge/' + gfile)

# done
print('Done')

