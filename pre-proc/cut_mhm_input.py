#!/usr/bin/env python
#############################################
# cut mhm input from gauge location in in river network
# given header file specifications
#
# created by Stephan Thober 17.02.2016
#
#############################################

from __future__ import print_function
import argparse
import textwrap
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
          description:
            This is the python pre-processing script for cutting a subbasin within a given original basin. The rectangle overlying the catchment remains unchanged. This allows to use the same meteorological forcing for the subbasin and the given original basin.

          Example:
            python cut_mhm_basin.py -i ../test_basin/input/ -o cutted_subbasin -g 398

          Note:
            All morphological input data of the given original basin have to be processed and are assumed to be inside the input directory. The files are: gauge/gauge_loc.asc, morph/aspect.asc, morph/idgauges.asc, morph/dem.asc, morph/slope.asc, morph/facc.asc, morph/fdir.asc, luse/lc_2001.asc
          
          Author:
            Stephan Thober
          Created:
            31 May 2016
          '''))

def check_dir(dir):
    if not isdir(dir):
        print(' create directory: ' + dir)
        makedirs(dir)
    

# default variables
indir = '../test_basin/input/'
outdir = './test_cut_mhm_input/'
lc_file = ['lc_1991.asc','lc_2000.asc','lc_2006.asc']
gauge_id = -9999.

parser.add_argument('-g', '--gauge_id', action='store', dest='gauge_id', default=gauge_id, 
                  help='gauge_id from where subbasin should be cutted')
parser.add_argument('-i', '--indir', action='store', dest='indir', default=indir, 
                  help='input directory containing mhm gauge, land use, and morphological information')
parser.add_argument('-l', '--lc_file', action='store', dest='lc_file', default=lc_file, 
                  help='land cover file to process')
parser.add_argument('-o', '--outdir', action='store', dest='outdir', default=outdir, 
                  help='directory where output should be written')

args = parser.parse_args()
gauge_id = int(args.gauge_id)
indir = args.indir
lc_file = args.lc_file
outdir = args.outdir

if indir[-1] != '/':
    indir += '/'
if outdir[-1] != '/':
    outdir += '/'

del args
#############################################
## START PROCESSING
#############################################
## mHM files to process

files = ['morph/aspect.asc',
         'morph/idgauges.asc',
         'morph/dem.asc',
         'morph/slope.asc',
         'morph/facc.asc',
         'morph/fdir.asc',
         'morph/soil_class.asc',
         'morph/geology_class.asc',
         'morph/LAI_class.asc',]

## append landuse files    
for lcf in lc_file:    
    files.append('luse/' + lcf)
#############################################


import numpy as np # array manipulation
from os.path import isdir
from os import makedirs
from fread import fread
from ufz import river_network, fwrite
#from river_network import river_network
#from fwrite import fwrite # from ufz
from subprocess import call # for system call

# read gauge_id file
idgauges = fread(indir + '/morph/idgauges.asc', skip=6)
if gauge_id == -9999.:
    gauges_loc = np.where(idgauges > 0.)
    gauges_id = idgauges[gauges_loc]
else:
    gauges_loc = np.where(idgauges == gauge_id)
    gauges_id = [gauge_id]
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
    print('   This script copies the gauge file using system call')
    print('   cp ' + indir + 'gauge/' + gfile + ' ' + outpath + 'gauge/' + gfile)
    check_dir(outpath + 'gauge/')
    call(['cp', indir + 'gauge/' + gfile, outpath + 'gauge/' + gfile])

# done
print('Done')

