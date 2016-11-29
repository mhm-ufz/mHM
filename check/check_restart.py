from ufz import readnc
import numpy as np

new_file = 'output_b1/mRM_restart_001.nc'
save_file = 'output_save/mRM_restart_001.nc'

import argparse
import textwrap
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
          Description:
            Check whether variables in two netcdf files are identical.

          Author:
            Stephan Thober

          Created:
            Nov 2016 

          Example:
            python check_restart.py -n <new_file> -s <save_file>

          Note:
          '''))
parser.add_argument('-n', '--new_file', action='store',
                    default=new_file, dest='new_file', metavar='new_file',
                    help='new file (default: output_b1/mRM_restart_001.nc)')
parser.add_argument('-s', '--save_file', action='store',
                    default=save_file, dest='save_file', metavar='save_file',
                    help='save file that new file is compared against, variables are taken from this file (default: output_save/mRM_restart_001.nc)')

args      = parser.parse_args()
new_file  = args.new_file
save_file = args.save_file
del parser, args

# get variables
v_names = readnc(save_file, variables=True)

for var in v_names:
    new_arr = readnc(new_file, var)
    save_arr = readnc(save_file, var)

    diff = np.ma.abs(new_arr - save_arr)

    if np.ma.amax(diff) > 0.:
        print('processing: ' + var)
        print('  max diff: ', np.ma.amax(diff))
        print('')

print('Done!')
