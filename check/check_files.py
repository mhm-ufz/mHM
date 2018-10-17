#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : check_files
Project Name: mHM check case repo
Description : insert your description here, if applicable
Author      : Robert Schweppe, Stephan Thober
Created     : 07.05.18 15:10
"""

# IMPORTS
import xarray as xr
from numpy import allclose, loadtxt, array_equal
import argparse
import textwrap

# GLOBAL VARIABLES
# NEW_FILE_DEFAULT = './case_00/output_b1/b1_mRM_restart_001.nc'
# REF_FILE_DEFAULT = './case_00/output_save/backup/b1_mRM_restart_001.nc'
NEW_FILE_DEFAULT = './case_10/output_b1/b1_mRM_Fluxes_States.nc'
REF_FILE_DEFAULT = './case_10/output_save/b1_mRM_Fluxes_States.nc'
RTOL = 1e-03
ATOL = 1e-04

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

def compare_nc_files(new_file, ref_file):
    ds_new = read_nc_file(new_file)
    ds_ref = read_nc_file(ref_file)
    diff_cnt = 0
    big_diff_cnt = 0
    for var_name in ds_ref.data_vars:
        if var_name in ds_new.data_vars:
            if not ds_new[var_name].equals(ds_ref[var_name]):
                diff_cnt += 1
                big_diff_cnt += compare_arrays(ds_new[var_name].values, ds_ref[var_name].values)
        else:
            print(var_name, 'is not contained in reference')
            diff_cnt += 1
            big_diff_cnt += 1
    message(diff_cnt, big_diff_cnt, len(ds_ref.data_vars))

def compare_csv_files(new_file, ref_file):
    header_new, ds_new = read_csv_file(new_file)
    header_ref, ds_ref = read_csv_file(ref_file)
    diff_cnt = 0
    big_diff_cnt = 0
    if not array_equal(ds_ref, ds_new):
        diff_cnt += 1
        big_diff_cnt += compare_arrays(ds_new, ds_ref)
    message(diff_cnt, big_diff_cnt, len(header_ref))


def compare_arrays(new_array, ref_array, rtol=RTOL, atol=ATOL):
    if not allclose(ref_array, new_array, equal_nan=True, rtol=rtol, atol=atol):
        return 1
    return 0

def read_nc_file(path):
    with xr.open_dataset(path) as ds:
        ds.load()
        return ds

def read_csv_file(path):
    with open(path) as f_in:
        # read the first lines with meta information
        header = f_in.readline().strip().split()
        # read the rest of the thing into a numpy array
        return header, loadtxt(f_in, dtype=float)

def message(n_diff, n_big_diff, n_total):
    print(' {} of {} records differ'.format(n_diff, n_total))
    print(' {} of {} records differ more than {}'.format(n_big_diff, n_total, ATOL))

# CLASSES

# SCRIPT
if __name__ == '__main__':
    new_file, ref_file, ftype = parse_args()
    if ftype == 'nc':
        compare_nc_files(new_file, ref_file)
    if ftype in ['param', 'discharge']:
        compare_csv_files(new_file, ref_file)
