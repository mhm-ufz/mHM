#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : asc2nc
Project Name: mhm
Description : converts all ascii files of a typical (test_basin) mhm <v.6.0 environment and converts them to
              netcdf files ready to be used in MPR from version 6.0 onwards
Author      : ottor
Created     : 28.08.17 14:53
"""

# IMPORTS
import os
import shutil
import pathlib
from itertools import product

import numpy
import pandas as pd
import xarray as xr
import argparse


# GLOBAL VARIABLES
# default input directory
IN_DIR = '../test_basin/input'
# default output directory
OUT_DIR = '../../MPR/reference/test_basin/input/temp'
# all file types scanned in input directory
POSSIBLE_SUFFIXES = ['.nc', '.asc']
# all folders scanned in input directory
FOLDER_LIST = ['lai', 'luse', 'morph']

# FUNCTIONS
def parse_args():
    """
    parses all command line arguments

    Returns
    -------
    Namespace object with args
    """

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='''This is a Python script to convert the input of mHM < v6.0 into 
                                     input for mHM >= v6.0.

    author: Robert Schweppe
    created: Oct 2018''')
    parser.add_argument('-i', '--in_dir', action='store',
                        default=IN_DIR, dest='input_dir', metavar='input_dir',
                        help='input directories, default: {}'.format(IN_DIR))
    parser.add_argument('-o', '--out_dir', action='store',
                        default=OUT_DIR, dest='output_dir', metavar='output_dir',
                        help='Output directory, where netcdf files for mHM >6.0 are stored'
                             ', default: {}'.format(OUT_DIR))

    return parser.parse_args()


def get_all_subfiles(path, relation=None):
    """
    returns generator, yielding all files found in any subfolder of a given path considering POSSIBLE_SUFFIXES and FOLDER_LIST
    relative to the given path

    Parameters
    ----------
    path: pathlib.Path
        the path to look for files
    relation:
        a parent folder to be prefixed to the resulting file path

    Returns
    -------
    generator
    """
    if path.is_file():
        if path.suffix in POSSIBLE_SUFFIXES:
            yield path.relative_to(relation or path)
    else:
        for sub_path in path.iterdir():
            if sub_path.name not in FOLDER_LIST and sub_path.is_dir():
                continue
            yield from get_all_subfiles(sub_path, relation or path)


# CLASSES
class MyAsciiToNetcdfConverter(object):
    POSSIBLE_LOOKUP_MODES = ['dims_as_col', 'dims_in_col']
    def __init__(self, input_file, output_file, lookup=None, sel=None, name=None, attrs=None, values_dtype=int,
                 iterate={}, csv_lut=False, lookup_mode='dims_in_col', lookup_nrows=None):
        """
        initializes an object that is capable of handling the voncersion from ascii to netcdf

        Parameters
        ----------
        input_file: str
            path to input file (ascii file with header information)
        output_file: str
            path to nc output file (can optionally include placeholders for formatting values)
            e.g. './h_{horizon:0>2d}.nc' with arg iterate={'horizon': x}
        lookup: str
            path to lookup table mapping values to the indices in the ascii file
            it is (whitespace-seperated) and may look like this:
            '''
            ID HORIZON UD[mm] LD[mm] CLAY[%] SAND[%] Bd[gcm-3] MUID
            1  1   0   50   33.   33.   1.5   DUMMY
            1   2   50   100   33. 33.   1.5   DUMMY
            1   3   100   200   33.  33.   1.5   DUMMY
            '''
        sel: str or list of str
            only applies if lookup is given and selects columns in the lookup file by name
        name: str
            if only one variable is created, this is the name set for this
        attrs: dict
            if only one variable is created, these are the attributes set for this
        values_dtype: type
            dtype set for values in array, to support np.nan, set to float
        iterate: dict
            dict containing additional dimensions to create in netcdf file
            they are grabbed from lookup file, created in dataset and used for creating output
        csv_lut: bool
            flag specify format of lookup table, either comma- or whitespace-separated
        lookup_mode: bool
            the way the lookup table is parsed, either:
                dims_as_col - all columns are used as a dimension, requires iterate (e.g. for LAI_class)
                              iterate must be the new dimension, e.g. {'month_of_year': range(1, 13)}
                dims_in_col - all unique values in column are used as a dimension,
                              requires iterate (e.g. for soil_class)
                              iterate must be the column name and its index, e.g. {'horizon': 1}
        lookup_nrows: int
            number of rows to read in lookup file
        """
        self.lookup_nrows = lookup_nrows
        self.lookup_mode = lookup_mode
        if self.lookup_mode not in self.POSSIBLE_LOOKUP_MODES:
            raise Exception('The provided lookup_mode "' + self.lookup_mode + '" is not valid. I must be one of: '
                            ', '.join(self.POSSIBLE_LOOKUP_MODES))
        self.csv_lut = csv_lut
        self.iterate = iterate
        self.values_dtype = values_dtype
        self.attrs = attrs
        if self.attrs is None:
            self.attrs = {}
        self.name = name
        self.sel = sel
        self.lookup = lookup
        if self.lookup is not None:
            self.lookup = self._format_path(self.lookup)
        self.input_file = self._format_path(input_file)
        self.output_file = self._format_path(output_file)

        self.raw_data = None
        self.data = None
        self.iterators = {}

    def read(self):
        """
        reads the ascii file and optionally the lookup table and creates a xr.Dataset at self.data
        """

        data_meta = {}

        with open(self.input_file) as f_in:
            # read the first lines with meta information
            for line in range(6):
                key, value = f_in.readline().split()
                data_meta[key] = float(value)
            # read the rest of the thing into a numpy array
            self.raw_data = numpy.loadtxt(f_in, dtype=self.values_dtype)

        # apply nan
        self.raw_data[self.raw_data == data_meta['NODATA_value']] = numpy.nan

        # span the coord arrays
        lats = self._make_coord(data_meta, 'lat')[::-1]
        lat_attrs = {'standard_name': 'latitude',
                     'long_name': 'latitude',
                     'units': 'degrees_north',
                     'axis': 'Y'}
        lons = self._make_coord(data_meta, 'lon')
        lon_attrs = {'standard_name': 'longitude',
                     'long_name': 'longitude',
                     'units': 'degrees_east',
                     'axis': 'X'}
        if self.lookup is not None:
            # read the lookup table with the specified columns, set the ID as index, skip initial row
            if self.csv_lut:
                delim_whitespace = False
                skiprows = 0
            else:
                delim_whitespace = True
                skiprows = 1
            index_col = [0]
            if self.lookup_mode == self.POSSIBLE_LOOKUP_MODES[1]:
                index_col.extend([_ for _ in list(self.iterate.values()) if _ is not None])
            lookup_data = pd.read_csv(self.lookup, skiprows=skiprows, usecols=self.sel, nrows=self.lookup_nrows,
                                      delim_whitespace=delim_whitespace,
                                      index_col=index_col)
            # iteratively build the nc Dataset
            coords = {'lat': lats, 'lon': lons}
            if self.lookup_mode == self.POSSIBLE_LOOKUP_MODES[1]:
                # update coord by unique values occuring in additional index col
                self.iterators = {key: lookup_data.index.get_level_values(value).unique().values for key, value in
                                  self.iterate.items()}
                # TODO: the horizon coord needs to get the value of the lower boundary of each layer
                coords.update(self.iterators)
                self.data = xr.Dataset({col_name.split('[')[0]: xr.DataArray(data=self._convert_raw(lookup_data[col_name]),
                                                                             coords=coords,
                                                                             dims=list(coords.keys()),
                                                                             name=col_name.split('[')[0],
                                                                             attrs={'units': col_name.split('[')[-1].rstrip(
                                                                                 ']')}) for col_name in
                                        lookup_data.columns},
                                       attrs=self.attrs)
            elif self.lookup_mode == self.POSSIBLE_LOOKUP_MODES[0]:
                # update coord by additional columns
                self.iterators = self.iterate
                coords.update(self.iterators)
                self.data = xr.DataArray(data=self._convert_raw(lookup_data),
                                         coords=coords,
                                         dims=list(coords.keys()),
                                         name=self.name,
                                         attrs=self.attrs).sortby(['lon', 'lat'])

        else:
            self.data = xr.DataArray(data=self.raw_data,
                                     coords={'lon': lons, 'lat': lats},
                                     dims=['lat', 'lon'],
                                     name=self.name,
                                     attrs=self.attrs,
                                     ).sortby(['lon', 'lat'])
        self.data.lon.attrs = lon_attrs
        self.data.lat.attrs = lat_attrs

    def _convert_raw(self, series):
        """
        applies the lookup table values to array

        Parameters
        ----------
        series: pd.Series
            Series with index as array indices and replacement values

        Returns
        -------
        np.array
        """
        converted_clone = self.raw_data.copy()
        # here all new dimensions are added to temporary field converted_clone
        for new_axes in self.iterators.values():
            # dynamically append a dimension to data to accommodate new values
            # build slice object for array
            old_slice = tuple(
                [slice(None) for x in converted_clone.shape] + [numpy.newaxis for x in self.iterators.keys()])
            # build slice object for new dimension
            new_slice = tuple([numpy.newaxis for x in converted_clone.shape] + [slice(None)])
            # use numpys broadcasting to combine
            converted_clone = converted_clone[old_slice] * numpy.ones_like(new_axes)[new_slice]
        if self.iterators:
            if self.lookup_mode == self.POSSIBLE_LOOKUP_MODES[1]:
                # loop over product of new dims
                for axes in product(*list(self.iterators.values())):
                    # https://stackoverflow.com/questions/5036816/numpy-lookup-map-or-point
                    # select the slice in the array, convert axes values to their corresponding index values
                    array_slice = tuple(
                        [slice(None) for x in self.raw_data.shape] + [list(list(self.iterators.values())[i_a]).index(a) for
                                                                      i_a, a in enumerate(axes)])
                    # select the slice in the lookup data
                    lookup_slice = tuple([slice(None)] + list(axes))
                    # select the lookup data
                    lookup_sel = series[lookup_slice]
                    # get lookup values and index them by index array
                    # append a nan value because searchsort returns len(lookup_sel.index) if nan is found
                    converted_clone[array_slice] = lookup_sel.append(pd.Series(numpy.nan)).values[
                        lookup_sel.index.searchsorted(converted_clone[array_slice])]
            elif self.lookup_mode == self.POSSIBLE_LOOKUP_MODES[0]:
                for index in numpy.unique(self.raw_data[~numpy.isnan(self.raw_data)]):
                    converted_clone[self.raw_data == index] = series.loc[index, :]
        else:
            #converted_clone = series.values[series.index.searchsorted(converted_clone)]
            for idx, val in series.iteritems():
                converted_clone[converted_clone == idx] = val
        return converted_clone

    def _make_coord(self, metas, which='lat'):
        """
        creates array with coordinates based on meta information from asc file header

        Parameters
        ----------
        metas: dict
            must contain entries 'xllcorner', 'yllcorner', 'nrows', 'ncols', 'cellsize'
        which: str
            one of 'lat' or 'lon'

        Returns
        -------
        np.array with corresponding latitude or longitude values

        """
        coord_dict = {'lat': ('yllcorner', 'nrows'), 'lon': ('xllcorner', 'ncols')}
        ll = metas[coord_dict[which][0]]
        n = metas[coord_dict[which][1]]
        w = metas['cellsize']
        return numpy.linspace(start=ll + 0.5 * w,
                              stop=ll + (n + 0.5) * w,
                              num=int(n),
                              endpoint=False)

    def _format_path(self, path):
        """
        formats a path

        Parameters
        ----------
        path: str
            raw path

        Returns
        -------
        str
        """
        return pathlib.Path(os.path.abspath(os.path.expandvars(path)))

    def write(self):
        """
        dumps the data to netcdf file
        """
        if hasattr(self.data, 'data_vars'):
            # if multiple data arrays are contained in dataset, write each in seperate file
            for data_var in self.data.data_vars:
                output_file = pathlib.Path(self.output_file.parent, data_var + self.output_file.suffix)
                self.data[data_var].to_netcdf(output_file)
        else:
            # if only dataarray is contained, then set name properly
            self.data.name = self.output_file.stem
            self.data.to_netcdf(self.output_file)

# SCRIPT
if __name__ == '__main__':
    args = parse_args()

    input_dir = pathlib.Path(args.input_dir)
    if input_dir.is_file():
        raise Exception("Input directory must be a directory, it is a file")
    path_list = get_all_subfiles(input_dir)

    for path in path_list:
        print('working on file: {}'.format(path))
        kwargs = {}
        if path.stem.endswith('_class'):
            kwargs['lookup'] = pathlib.Path(input_dir, path.parent, path.stem + 'definition.txt')
            if path.stem == 'LAI_class':
                #kwargs['csv_lut'] = False
                kwargs['sel'] = ['ID', 'Jan.', 'Feb.', 'Mar.', 'Apr.', 'May', 'Jun.', 'Jul.', 'Aug.', 'Sep.',
                                 'Oct.', 'Nov.', 'Dec.']
                kwargs['lookup_mode'] = 'dims_as_col'
                kwargs['attrs'] = {'units': '-'}
                kwargs['iterate'] = {'month_of_year': list(range(1,13))}
            elif path.stem == 'soil_class':
                kwargs['iterate'] = {'horizon': 1}
            elif path.stem == 'geology_class':
                kwargs['sel'] = ['GeoParam(i)', 'ClassUnit', 'Karstic']
                kwargs['lookup_nrows'] = 12
        elif '_class_horizon_' in path.stem:
            kwargs['lookup'] = pathlib.Path(input_dir, path.parent, 'soil_classdefinition_iFlag_soilDB_1.txt')
        if path.suffix == '.nc':
            shutil.copy(pathlib.Path(input_dir, path), pathlib.Path(args.output_dir, path.name))
        else:
            my_conv = MyAsciiToNetcdfConverter(
                input_file=pathlib.Path(input_dir, path),
                output_file=pathlib.Path(args.output_dir, path.stem + '.nc'),
                values_dtype=float,
                **kwargs
            )
            my_conv.read()
            my_conv.write()
