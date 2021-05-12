#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : asc2nc
Project Name: mhm
Description : converts all ascii files of a typical (test_domain) mhm <v.6.0 environment and converts them to
              netcdf files ready to be used in MPR from version 6.0 onwards
Author      : ottor
Created     : 28.08.17 14:53
"""

# IMPORTS
import os
import pathlib
import shutil
from itertools import product
import re
from packaging import version

import numpy
import pandas as pd
import xarray as xr
import argparse
import sys

# we rely on dict preserving key order
MIN_PYTHON = (3, 7)
if sys.version_info < MIN_PYTHON:
    sys.exit("Python %s.%s or later is required.\n" % MIN_PYTHON)

# GLOBAL VARIABLES
# default input directory
HORIZON_COORD_NAME = 'horizon'
IN_DIR = '../test_domain/input'
IN_DIR = '/Users/ottor/nc/Home/local_libs/fortran/mhm_original/test_domain/input'
# default output directory
OUT_DIR = '../../MPR/reference/test_domain/input/temp'
OUT_DIR = '/Users/ottor/temp/test_domain'
# all file types scanned in input directory
POSSIBLE_SUFFIXES = ['.nc', '.asc']
# all folders scanned in input directory
FOLDER_LIST = ['lai', 'luse', 'morph', 'latlon', 'optional_data', 'meteo', 'pet', 'pre', 'tavg']

LAT_ATTRS = {'standard_name': 'latitude',
             'long_name': 'latitude',
             #'units': 'degrees_north',
             'axis': 'Y'}
LON_ATTRS = {'standard_name': 'longitude',
             'long_name': 'longitude',
             #'units': 'degrees_east',
             'axis': 'X'}
LC_ATTRS = {'standard_name': 'land cover periods',
             'long_name': 'periods of common land cover',
             'units': 'years',
             'axis': 'T'}

SOIL_ATTRS = {'standard_name': 'horizon',
             'long_name': 'soil horizons',
             'units': 'm',
             'axis': 'Z',
             'positive': 'down',
             'bounds': 'horizon_bnds'}

PROPERTIES_MAPPING = {
    'aspect': ('mpr', 'int32', 1.0, -9999),
    'bd': ('mpr', 'float64', 0.01, -9999.0),
    'classunit': ('mpr', 'int32', 1.0, -9999),
    'clay': ('mpr', 'float64', 0.01, -9999.0),
    'dem': ('routing', 'int32', 1.0, -9999),
    'facc': ('routing', 'int32', 1.0, -9999),
    'fdir': ('routing', 'int32', 1.0, -9999),
    'idgauges': ('routing', 'int32', 1.0, -9999),
    'karstic': ('mpr', 'int32', 1.0, -9999),
    'land_cover': ('mpr', 'int32', 1.0, -9999),
    'lai_class': ('mpr', 'float64', 0.01, -9999.0),
    'ld': ('mpr', 'int32', 1.0, -9999),
    'sand': ('mpr', 'float64', 0.01, -9999.0),
    'slope': ('mpr', 'int32', 1.0, -9999),
    'ud': ('mpr', 'int32', 1.0, -9999),
    'lat_l0': ('mpr', 'float64', 1.0, -9999.0),
    'et': ('optional_data', 'float64', 1.0, -9999.0),
    'neutrons': ('optional_data', 'float64', 1.0, -9999.0),
    'Q_bkfl': ('optional_data', 'float64', 1.0, -9999.0),
    'sm': ('optional_data', 'float64', 1.0, -9999.0),
    'twsa': ('optional_data', 'float64', 1.0, -9999.0),
    'eabs': ('meteo', 'float64', 1.0, -9999.0),
    'mask': ('meteo', 'float64', 1.0, -9999.0),
    'net_rad': ('meteo', 'float64', 1.0, -9999.0),
    'ssrd': ('meteo', 'float64', 1.0, -9999.0),
    'strd': ('meteo', 'float64', 1.0, -9999.0),
    'tmax': ('meteo', 'float64', 1.0, -9999.0),
    'tmin': ('meteo', 'float64', 1.0, -9999.0),
    'windspeed': ('meteo', 'float64', 1.0, -9999.0),
    'tann': ('meteo', 'float64', 1.0, -9999.0),
    'tavg': ('meteo', 'float64', 1.0, -9999.0),
    'pre': ('meteo', 'float64', 1.0, -9999.0),
    'pet': ('meteo', 'float64', 1.0, -9999.0),
    SOIL_ATTRS['bounds']: ('mpr', 'float64', 1.0, -9999),
}

COMPRESSION_DICT = {
    'zlib': True, 'shuffle': True, 'complevel': 4,
}

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
                 iterate={}, csv_lut=False, lookup_mode='dims_in_col', infer_nrows=False, index_as_col=False):
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
        infer_nrows: bool
            number of rows to read in lookup file is inferred by parsing the first line first and then read the
            file again
        index_as_col: bool
            whether to set the first column as index also (ds.index = ds.iloc[:, 0])
        """
        self.index_as_col = index_as_col
        self.infer_nrows = infer_nrows
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
        lons = self._make_coord(data_meta, 'lon')
        coords = {'lat': lats, 'lon': lons}
        if self.lookup is not None:
            kwargs = {}
            # read the lookup table with the specified columns, set the ID as index, skip initial row
            if self.csv_lut:
                kwargs['delim_whitespace'] = False
                kwargs['skiprows'] = 0
            kwargs['delim_whitespace'] = True
            kwargs['skiprows'] = 1
            kwargs['index_col'] = [0]
            if self.lookup_mode == self.POSSIBLE_LOOKUP_MODES[1]:
                kwargs['index_col'].extend([_ for _ in list(self.iterate.values()) if _ is not None])
            if self.infer_nrows:
                with open(self.lookup) as config_file:
                    first_line = config_file.readline()
                try:
                    kwargs['nrows'] = int(first_line.split()[-1])
                except ValueError:
                    raise
            lookup_data = pd.read_csv(self.lookup, usecols=self.sel, **kwargs)
            if self.index_as_col:
                lookup_data.index = lookup_data.iloc[:, 0]
            # iteratively build the nc Dataset
            if self.lookup_mode == self.POSSIBLE_LOOKUP_MODES[1]:
                # narrow down the lookup values to ids actually occurring in the data
                use_cols = lookup_data.columns
                new_z = None
                add_info = None
                if isinstance(lookup_data.index, pd.MultiIndex):
                    unique_ids = numpy.unique(self.raw_data[~numpy.isnan(self.raw_data)])
                    lookup_data = lookup_data[lookup_data.index.isin(unique_ids, 0)]
                    # set nan to non-existing id-soil layer combinations
                    lookup_data.index = lookup_data.index.remove_unused_levels()
                    new_index = pd.MultiIndex.from_product(lookup_data.index.levels,
                                                           names=lookup_data.index.names)
                    lookup_data = lookup_data.reindex(new_index)
                    sel_cols = [col for col in lookup_data.columns if col.startswith(('UD', 'LD'))]
                    add_info = lookup_data[sel_cols]
                    new_z = numpy.unique(lookup_data[sel_cols].dropna().values)
                    self.iterators.update({HORIZON_COORD_NAME: new_z[1:]})

                    use_cols = [col for col in lookup_data.columns if col not in sel_cols]
                else:
                    # update coord by unique values occurring in additional index col
                    self.iterators.update({key: lookup_data.index.get_level_values(value).unique().values for key, value in
                                      self.iterate.items()})
                coords.update(self.iterators)
                self.data = xr.Dataset({col_name.split('[')[0]: xr.DataArray(data=self._convert_raw(lookup_data[col_name], add_info),
                                                                             coords=coords,
                                                                             dims=list(coords.keys()),
                                                                             name=col_name.split('[')[0],
                                                                             attrs={'units': col_name.split('[')[-1].rstrip(
                                                                                 ']')}) for col_name in use_cols},
                                       attrs=self.attrs)
                if new_z is not None:
                    # add the bound data to the data dict for xr.Dataset initialization
                    self.data[SOIL_ATTRS['bounds']] = ([HORIZON_COORD_NAME, 'bnds'],
                                                  numpy.stack([new_z[:-1], new_z[1:]], 1))

            elif self.lookup_mode == self.POSSIBLE_LOOKUP_MODES[0]:
                # update coord by additional columns
                self.iterators = self.iterate
                coords.update(self.iterators)
                self.data = xr.DataArray(data=self._convert_raw(lookup_data),
                                         coords=coords,
                                         dims=list(coords.keys()),
                                         name=self.name,
                                         attrs=self.attrs)

        else:
            self.data = xr.DataArray(data=self.raw_data,
                                     coords=coords,
                                     dims=['lat', 'lon'],
                                     name=self.name,
                                     attrs=self.attrs,
                                     )
        # sort the data so they are always in ascending order along every dimension
        self.data = self.data.sortby(['lon', 'lat'])
        coords['lon'] = numpy.sort(coords['lon'])
        coords['lat'] = numpy.sort(coords['lat'])

        # set coordinate attributes
        self.data.lon.attrs = LON_ATTRS
        self.data.lat.attrs = LAT_ATTRS
        if HORIZON_COORD_NAME in self.data:
            # add the nice new attributes
            self.data[HORIZON_COORD_NAME] = self.data[HORIZON_COORD_NAME] / 1000
            self.data[SOIL_ATTRS['bounds']] = self.data[SOIL_ATTRS['bounds']] / 1000
            self.data[HORIZON_COORD_NAME].attrs = SOIL_ATTRS

        # HACKS: sort them in ascending order
        if self.input_file.stem == 'fdir':
            self._rotate_fdir()

    def _convert_raw(self, series, add_info=None):
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
        def select_layer(df, lower_bound):
            # check for soil layer that matches lower_bound
            mask = (df['UD[mm]'] < lower_bound) & (df['LD[mm]'] >= lower_bound)
            try:
                return df.index[mask].values[0]
            except IndexError:
                # return that layer else the last one available (assumed to be nan)
                return df.index[-1]

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
                    if add_info is not None:
                        # now select the correct values taking into account the layer depth
                        lookup_slice = add_info.groupby(level=0).apply(select_layer, axes)
                        # select the lookup data
                        lookup_sel = series[lookup_slice]
                        lookup_sel.index = lookup_sel.index.droplevel(1)
                    else:
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
            # if multiple data arrays are contained in dataset, write each in separate file
            for data_var in self.data.data_vars:
                if data_var.endswith('_bnds'):
                    continue
                self.data = self.data.rename({data_var: data_var.lower()})
                export_vars = [data_var.lower()] + [var for var in self.data.data_vars if var.endswith('_bnds')]
                self._write_data_var(self.data[export_vars], data_var.lower())
        else:
            # if only dataarray is contained, then set name properly
            self.data.name = self.output_file.stem.lower()
            self._write_data_var(self.data, self.data.name)

    def _write_data_var(self, data_var, stem):
        output_file = pathlib.Path(self.output_file.parent, PROPERTIES_MAPPING.get(stem, ('mpr',))[0],
                                   stem + self.output_file.suffix)
        if not output_file.parent.exists():
            output_file.parent.mkdir(parents=True)
        print('writing variable {} to file: {}'.format(stem, output_file))
        data_var.to_netcdf(output_file, encoding={stem: dict(
            dtype=PROPERTIES_MAPPING.get(stem, ['int32']*2)[1],
            # scale_factor=PROPERTIES_MAPPING.get(stem, [1.0]*3)[2],
            _FillValue=PROPERTIES_MAPPING.get(stem, [-9999]*4)[3],
            **COMPRESSION_DICT
        )})

    def _rotate_fdir(self):
        masks = [self.data.values == 2**(i) for i in range(8)]
        replacements = [4, 8, 16, 32, 64, 128, 1, 2]
        for mask, replacement in zip(masks, replacements):
            self.data.values[mask] = replacement


def combine_lc_files(output_dir):
    """
    combines single land_cover_files into one file
    """
    # combine all files
    path_list = sorted(output_dir.glob('lc_*.nc'))
    # get the years as future dimension values
    years = [int(re.search('(\d){4}', path.stem).group()) for path in path_list]
    # open all files and concat them
    open_kwargs = {'parallel': False}
    if version.parse(xr.__version__) >= version.parse('0.14'):
        open_kwargs['combine'] = 'by_coords'
    ds = xr.open_mfdataset(paths=path_list, **open_kwargs)
    ds = xr.concat([ds[_] for _ in ds.data_vars], pd.Index(years, name='land_cover_period'))
    # set some attributes
    var_name = 'land_cover'
    ds.name = var_name
    bound_attr = {'bounds': 'land_cover_period_bnds'}
    # and convert to Dataset so we can add more dimensions
    ds = ds.to_dataset()
    ds['land_cover_period_bnds'] = xr.DataArray(
        numpy.array([years, years[1:] + [years[-1] + int((years[-1] - years[0]) / len(years))]]).T,
        dims=['land_cover_period', 'bnds'])
    # add all dimension attributes
    ds.lon.attrs = LON_ATTRS
    ds.land_cover_period.attrs = dict(**LC_ATTRS, **bound_attr)
    ds.lat.attrs = LAT_ATTRS

    # write new file
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    output_file = pathlib.Path(output_dir, f'{var_name}.nc')
    print('writing variable {} to file: {}'.format(var_name, output_file))
    ds.to_netcdf(output_file, encoding={var_name: dict(
            dtype=PROPERTIES_MAPPING.get(var_name, ['int32']*2)[1],
            # scale_factor=PROPERTIES_MAPPING.get(stem, [1.0]*3)[2],
            _FillValue=PROPERTIES_MAPPING.get(var_name, [-9999]*4)[3],
            **COMPRESSION_DICT)})
    # delete the old files
    for path in path_list:
        path.unlink()
    ds.close()

def sort_y_dim(filename_in, filename_out):
    """
    read in some file, check for dimensions starting with 'y' and sort by this dimension
    """
    ds = xr.open_dataset(filename_in)
    sortby_dims = [dim for dim in ds.dims if dim.startswith('y')]
    if sortby_dims:
        ds = ds.sortby(sortby_dims, ascending=True)
    for var in ds.data_vars:
        missing_value = ds[var].attrs.get('missing_value')
        if missing_value is not None:
            if 'i' in ds['bd'].dtype.str:
                ds[var].encoding['_FillValue'] = int(missing_value)
            else:
                ds[var].encoding['_FillValue'] = float(missing_value)
    if not filename_out.parent.exists():
        filename_out.parent.mkdir(parents=True)
    print('writing variable {} to file: {}'.format(filename_out.stem, filename_out))
    ds.to_netcdf(filename_out, encoding={data_var: COMPRESSION_DICT for data_var in ds.data_vars})
    ds.close()


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
                kwargs['infer_nrows'] = True
            elif path.stem == 'soil_class':
                kwargs['iterate'] = {HORIZON_COORD_NAME: 1}
            elif path.stem == 'geology_class':
                kwargs['sel'] = ['GeoParam(i)', 'ClassUnit', 'Karstic']
                # put the number of GeoParam(i) rows here, usually 10
                kwargs['infer_nrows'] = True
                kwargs['index_as_col'] = True
        elif '_class_horizon_' in path.stem:
            kwargs['lookup'] = pathlib.Path(input_dir, path.parent, 'soil_classdefinition_iFlag_soilDB_1.txt')
        if path.suffix == '.nc':
            output_file = pathlib.Path(args.output_dir, PROPERTIES_MAPPING.get(path.stem, ('mpr',))[0],
                          path.stem + path.suffix)
            sort_y_dim(pathlib.Path(input_dir, path), output_file)
        else:
            my_conv = MyAsciiToNetcdfConverter(
                input_file=pathlib.Path(input_dir, path),
                output_file=pathlib.Path(args.output_dir, path.stem + '.nc'),
                values_dtype=float,
                **kwargs
            )
            my_conv.read()
            my_conv.write()

    combine_lc_files(pathlib.Path(args.output_dir, PROPERTIES_MAPPING.get('land_cover', ('mpr',))[0]))
