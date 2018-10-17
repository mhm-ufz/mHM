#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : check_restart_file
Project Name: scripts
Description : insert your description here, if applicable
Author      : ottor
Created     : 07.05.18 15:10
"""

# IMPORTS
import xarray as xr
from numpy import allclose
import argparse
import textwrap

# GLOBAL VARIABLES
MATCH_VARS = {'L1_basin_mask': 'L1_basin_Mask', 'L1_basin_cellarea': 'L1_areaCell',
              'L11_basin_mask': 'L11_basin_Mask', 'L11_basin_cellarea': 'L11_areaCell',
              'L11_nLinkFracFPimp': 'L11_FracFPimp'}
IGNORE_VARS = ['L1_basin_lat', 'L11_basin_lat',
               'L1_basin_lon', 'L11_basin_lon',
               'LC_year_start', 'LC_year_end', ]
NEW_FILE_DEFAULT = './case_00/output_b1/b1_mRM_restart_001.nc'
REF_FILE_DEFAULT = './case_00/output_save/backup/b1_mRM_restart_001.nc'
RTOL = 1e-03
ATOL = 1e-03

# FUNCTIONS
def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
              Description:
                Check whether variables in two netcdf files are identical.

              Author:
                Stephan Thober, Robert Schweppe

              Created:
                Nov 2016 

              Example:
                python check_files.py -n <new_file> -s <save_file>

              Note:
              '''))
    parser.add_argument('-n', '--new_file', action='store',
                        default=NEW_FILE_DEFAULT, dest='new_file', metavar='new_file',
                        help='new file (default: {})'.format(NEW_FILE_DEFAULT))
    parser.add_argument('-s', '--save_file', action='store',
                        default=REF_FILE_DEFAULT, dest='save_file', metavar='save_file',
                        help='reference file that new file is compared against, '
                             'variables are taken from this file (default: {})'.format(REF_FILE_DEFAULT))
    args = parser.parse_args()
    if args.new_file.endswith('.nc'):
        return args.new_file, args.save_file, 'nc'
    if args.new_file.endswith('discharge.out'):
        return args.new_file, args.save_file, 'discharge'
    if args.new_file.endswith('FinalParam.out'):
        return args.new_file, args.save_file, 'param'

def compare_arrays(new_array, ref_array, rtol=RTOL, atol=ATOL):
    if not allclose(ref_array, new_array, equal_nan=True, rtol=rtol, atol=atol):
        return 1
    return 0

def read_nc_file(path):
    with xr.open_dataset(path) as ds:
        ds.load()
        return ds

def message(n_diff, n_big_diff, n_total):
    print(' {} of {} records differ'.format(n_diff, n_total))
    print(' {} of {} records differ more than {}'.format(n_big_diff, n_total, ATOL))

def compare_restarts(new_file, ref_file):
    elim_dim = ['LCoverScenes', 'LAI_timesteps']
    with xr.open_dataset(new_file) as ds_new:
        if elim_dim[0] in ds_new.dims:
            ds_new = next(iter(ds_new.groupby(elim_dim[0])))[1]
        if elim_dim[1] in ds_new.dims:
            it = iter(ds_new.groupby(elim_dim[1], squeeze=False))
            try:
                while True:
                    ds_new = next(it)
            except StopIteration:
                ds_new = ds_new[1].squeeze(elim_dim[1])
        ds_new.load()
    with xr.open_dataset(ref_file) as ds_ref:
        ds_ref.load()
    diff_cnt = 0
    big_diff_cnt = 0
    for var_name in ds_new.data_vars:
        if var_name in ds_ref.data_vars or var_name in MATCH_VARS:
            if not ds_ref[MATCH_VARS.get(var_name, var_name)].equals(ds_new[var_name]):
                diff_cnt += 1
                # print(var_name, 'differs')
                # pb = big_diff_cnt
                big_diff_cnt += compare_arrays(ds_new[var_name].values,
                                               ds_ref[MATCH_VARS.get(var_name, var_name)].values)
                # if pb != big_diff_cnt:
                    # print(var_name, 'differs much')
        elif var_name not in IGNORE_VARS:
            print(var_name, 'is not contained in reference')
            diff_cnt += 1
            big_diff_cnt += 1
    message(diff_cnt, big_diff_cnt, len(ds_ref.data_vars))


# CLASSES

# SCRIPT
if __name__ == '__main__':
    new_file, ref_file, ftype = parse_args()
    if ftype == 'nc':
        compare_restarts(new_file, ref_file)
