#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : reformat_test_files
Project Name: mhm_fork
Description : insert your description here, if applicable
Author      : ottor
Created     : 17.05.21 10:02
"""

# IMPORTS
import pathlib
from asc2nc import get_all_subfiles, sort_y_dim
import shutil
import f90nml
from os import sep
from datetime import date
import re


# GLOBAL VARIABLES
USE_SYMLINKS = False
ROOT_TARGET_FOLDER = '../../test_domain/input_new/'
# ROOT_TARGET_FOLDER = './output'
STATIC_TARGET_FOLDER = 'mpr'
FORCING_TARGET_FOLDER = 'meteo'
CONFIG_TARGET_FOLDER = 'config'
OUTPUT_TARGET_FOLDER = 'output'
TARGET_MPR_NML_NAME = 'mpr.nml'
# input paths (environment from Zink et al., 2017)
ROOT_SOURCE_FOLDER = '/data/stohyd/data/processed/Germany/basins/'
FORCING_SOURCE_FOLDER = 'meteo'
FORCING_SOURCE_SEL = ['pre', 'tavg', 'pet']
# local mHM repository
ROOT_REPO_FOLDER = '/home/ottor/lib/mhm_mpr/'
# ROOT_REPO_FOLDER = '/Users/ottor/nc/Home/local_libs/fortran/mhm_mpr/'
CONFIG_SOURCE_SEL = [
    'mhm.nml', 'mhm_outputs.nml', 'mhm_parameter.nml', 'mrm_outputs.nml',
    # 'mrm.nml',  'mrm_parameter.nml',
    'mpr_fso.nml'
]


# commented lines denote keywords for upcoming mHM version
MHM_NML_REPLACE_DICT = {
    ('project_description', 'contact') : 'mHM developers (email:mhm-developers@ufz.de)',
    ('project_description', 'conventions') : 'XXX',
    ('project_description', 'history') : 'model run version 1',
    ('project_description', 'mhm_details') : 'Helmholtz Center for Environmental Research - UFZ, Department Computational Hydrosystems, Stochastic Hydrology Group',
    ('project_description', 'project_details') : 'mHM test domain project',
    ('project_description', 'setup_description') : 'model run for the Mosel domain, forced with the E-OBS meteorologic data',
    ('project_description', 'simulation_type') : 'historical simulation',
    ('mainconfig', 'iflag_cordinate_sys'): 0,
    ('mainconfig', 'ndomains'): 1,
    #('mainconfig', 'nbasins'): 1,
    ('mainconfig', 'resolution_hydrology'): [4000],
    ('mainconfig', 'l0domain'): [1],
    #('mainconfig', 'l0basin'): [1],
    ('mainconfig', 'write_restart'): True,
    ('mainconfig', 'read_opt_domain_data'): [0],
    ('mainconfig_mhm_mrm', 'dir_restartin'): [''],
    ('mainconfig_mhm_mrm', 'resolution_routing'): [4000],
    ('mainconfig_mhm_mrm', 'timestep'): 1,
    ('mainconfig_mhm_mrm', 'read_restart'): False,
    ('mainconfig_mhm_mrm', 'optimize'): False,
    ('mainconfig_mhm_mrm', 'optimize_restart'): False,
    ('mainconfig_mhm_mrm', 'opti_method'): 1,
    ('mainconfig_mhm_mrm', 'opti_function'): 10,
    ('mainconfig_mrm', 'alma_convention'): True,
    ('mainconfig_mrm', 'varnametotalrunoff'): 'total_runoff',
    ('mainconfig_mrm', 'filenametotalrunoff'): 'total_runoff',
    ('mainconfig_mrm', 'gw_coupling'): False,
    ('directories_general', 'dirconfigout'): OUTPUT_TARGET_FOLDER + sep,
    ('directories_general', 'dir_restartout'): [OUTPUT_TARGET_FOLDER + sep],
    ('directories_general', 'dir_out'): [OUTPUT_TARGET_FOLDER + sep],
    ('directories_mhm', 'path_mpr_nml'): [TARGET_MPR_NML_NAME],
    ('directories_mhm', 'dir_precipitation'): [str(pathlib.Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER)) + sep],
    ('directories_mhm', 'dir_temperature'): [str(pathlib.Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER)) + sep],
    ('directories_mhm', 'dir_referenceet'): [str(pathlib.Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER)) + sep],
    ('directories_mhm', 'dir_mintemperature'): [str(pathlib.Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER)) + sep],
    ('directories_mhm', 'dir_maxtemperature'): [str(pathlib.Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER)) + sep],
    ('directories_mhm', 'dir_netradiation'): [str(pathlib.Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER)) + sep],
    ('directories_mhm', 'dir_absvappressure'): [str(pathlib.Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER)) + sep],
    ('directories_mhm', 'dir_windspeed'): [str(pathlib.Path(ROOT_TARGET_FOLDER, FORCING_TARGET_FOLDER)) + sep],
    ('directories_mhm', 'time_step_model_inputs'): [0],
    ('directories_mhm', 'timestep_lai_input'): [0],
    ('directories_mrm', 'dir_gauges'): [str(pathlib.Path(ROOT_TARGET_FOLDER, 'routing')) + sep],
    ('directories_mrm', 'dir_total_runoff'): [str(pathlib.Path(ROOT_TARGET_FOLDER, 'routing')) + sep],
    ('directories_mrm', 'dir_bankfull_runoff'): [str(pathlib.Path(ROOT_TARGET_FOLDER, 'routing')) + sep],
    ('optional_data', 'dir_soil_moisture'): [str(pathlib.Path(ROOT_TARGET_FOLDER, 'optional_data')) + sep],
    ('optional_data', 'nsoilhorizons_sm_input'): 1,
    ('optional_data', 'timestep_sm_input'): -2,
    ('optional_data', 'dir_neutrons'): [str(pathlib.Path(ROOT_TARGET_FOLDER, 'optional_data')) + sep],
    ('optional_data', 'dir_evapotranspiration'): [str(pathlib.Path(ROOT_TARGET_FOLDER, 'optional_data')) + sep],
    ('optional_data', 'timestep_et_input'): -2,
    ('optional_data', 'dir_tws'): [str(pathlib.Path(ROOT_TARGET_FOLDER, 'optional_data')) + sep],
    #('optional_data', 'file_tws'): [str(pathlib.Path(ROOT_TARGET_FOLDER, 'optional_data', 'tws_basin_1.txt'))],
    ('optional_data', 'timestep_tws_input'): -2,
    ('processselection', 'processcase'): [1, 1, 1, 1, 0, 1, 1, 3, 1, 0],
    #('processselection', 'processcase'): [1, 1, 1, 1, 0, 1, 1, 2, 1, 1],
    ('time_periods', 'warming_days'): [0],
    ('time_periods', 'eval_per'): [
        {
            'ystart': 1990,
            'mstart': 1,
            'dstart': 1,
            'yend': 1991,
            'mend': 12,
            'dend': 31,
        }
    ],
    ('lcover', 'nlandcoverperiods'): 3,
    ('evaluation_gauges', 'ngaugestotal'): 1,
    ('evaluation_gauges', 'nogauges_domain'): [1],
    #('evaluation_gauges', 'nogauges_basin'): [1],
    ('evaluation_gauges', 'gauge_id'): [[1234]],
    ('evaluation_gauges', 'gauge_filename'): [['1234.txt']],
    ('inflow_gauges', 'ninflowgaugestotal'): 0,
    ('inflow_gauges', 'noinflowgauges_domain'): [0],
    #('inflow_gauges', 'noinflowgauges_basin'): [0],
    ('inflow_gauges', 'inflowgauge_id'): [[-9]],
    ('inflow_gauges', 'inflowgauge_filename'): [['']],
    ('inflow_gauges', 'inflowgauge_headwater'): [[False]],
    ('panevapo', 'evap_coeff'): [1.3, 1.2, 0.72, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.5],
    ('nightdayratio', 'read_meteo_weights'): False,
    ('nightdayratio', 'fnight_prec'): [0.46, 0.5, 0.52, 0.51, 0.48, 0.5, 0.49, 0.48, 0.52, 0.56, 0.5, 0.47],
    ('nightdayratio', 'fnight_pet'): [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
    ('nightdayratio', 'fnight_temp'): [-0.76, -1.3, -1.88, -2.38, -2.72, -2.75, -2.74, -3.04, -2.44, -1.6, -0.94,
                                       -0.53],
    ('optimization', 'niterations'): 7,
    ('optimization', 'seed'): 1235876,
    ('optimization', 'dds_r'): 0.2,
    ('optimization', 'sa_temp'): -9.0,
    ('optimization', 'sce_ngs'): 2,
    ('optimization', 'sce_npg'): -9,
    ('optimization', 'sce_nps'): -9,
    ('optimization', 'mcmc_opti'): False,
    ('optimization', 'mcmc_error_params'): [0.01, 0.6],
}


# FUNCTIONS
def _read_namelist(path):
    parser = f90nml.Parser()
    parser.global_start_index = 1
    return parser.read(path)


def prepare_config(file_name, out_file_name):
    if file_name.name == 'mhm.nml':
        _configure_mhm_nml(file_name, out_file_name,
                           output_path=output_target_root,
                           forcing_path=pathlib.Path('../../test_domain/input_new/meteo'),
                           routing_path=pathlib.Path('../../test_domain/input_new/routing'),
                           # this is a custom hack, make more flexible if needed
                           gauge_id=int(re.search(r'(\d)+', basin).group()),
                           start_date=date(2000, 1, 1),
                           end_date=date(2004, 12, 31),
                           spinup=1825,
                           )
    elif file_name.name == 'mpr_fso.nml':
        data_arrays_props = {var: {
            'from_file': str(pathlib.Path(ROOT_TARGET_FOLDER, STATIC_TARGET_FOLDER, basin, 'mpr', '{}.nc'.format(var)))}
                 for var in MPR_READ_VARS}
        data_arrays_props['lat_l0'] = {
            'from_file': str(pathlib.Path(ROOT_TARGET_FOLDER, STATIC_TARGET_FOLDER, basin, 'mpr', 'latlon.nc'))}
        _configure_mpr_nml(config_path, pathlib.Path(config_target_root, TARGET_MPR_NML_NAME),
                           coords={
                               'land_cover_period_out': {
                                   'coord_from_values': [2000, 2006, 2011],
                                   'coord_from_values_bound': 1950,
                                   'coord_stagger': 'end',
                               },
                               'lat_out': {
                                   'coord_from_range_step': 4000
                               },
                               'lon_out': {
                                   'coord_from_range_step': 4000
                               },
                               'horizon_out': {
                                   'coord_from_values': [0.05, 0.25, 2.0],
                                   'coord_from_values_bound': 0.0,
                                   'coord_stagger': 'end',
                               },
                               'horizon_till': {
                                   'coord_from_values': [0.05, 0.25],
                                   'coord_from_values_bound': 0.0,
                                   'coord_stagger': 'end',
                               },
                               'horizon_notill': {
                                   'coord_from_values': [2.0],
                                   'coord_from_values_bound': 0.25,
                                   'coord_stagger': 'end',
                               },
                           },
                           data_arrays=data_arrays_props,
                           out_filename=str(pathlib.Path(output_target_root, 'mHM_parameters.nc'))
                           )
    else:
        # create new file
        print('Creating {} file: {}'.format('config', pathlib.Path(config_target_root, file_name)))
        _copy_or_link_file(config_path, pathlib.Path(config_target_root, file_name))


def _copy_or_link_file(source, target):
    if USE_SYMLINKS:
        if target.is_symlink():
            target.unlink()
        target.symlink_to(source)
    else:
        if target.exists():
            target.unlink()
        shutil.copy(source, target)


def _configure_mhm_nml(in_path, out_path, output_path=None, forcing_path=None, routing_path=None, gauge_id=None,
                       start_date=None, end_date=None, spinup=None):
    local_mhm_nml_replace_dict = {}
    if output_path is not None:
        local_mhm_nml_replace_dict[('directories_general', 'dirconfigout')] = str(output_path) + sep
        local_mhm_nml_replace_dict[('directories_general', 'dir_restartout')] = [str(output_path) + sep]
        local_mhm_nml_replace_dict[('directories_general', 'dir_out')] = [str(output_path) + sep]
    if forcing_path is not None:
        local_mhm_nml_replace_dict[('directories_mhm', 'dir_precipitation')] = [str(forcing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mhm', 'dir_temperature')] = [str(forcing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mhm', 'dir_referenceet')] = [str(forcing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mhm', 'dir_mintemperature')] = [str(forcing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mhm', 'dir_maxtemperature')] = [str(forcing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mhm', 'dir_netradiation')] = [str(forcing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mhm', 'dir_absvappressure')] = [str(forcing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mhm', 'dir_windspeed')] = [str(forcing_path) + sep]
    if routing_path is not None:
        local_mhm_nml_replace_dict[('directories_mrm', 'dir_gauges')] = [str(routing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mrm', 'dir_total_runoff')] = [str(routing_path) + sep]
        local_mhm_nml_replace_dict[('directories_mrm', 'dir_bankfull_runoff')] = [str(routing_path) + sep]
    if gauge_id is not None:
        local_mhm_nml_replace_dict[('evaluation_gauges', 'gauge_id')] = [[gauge_id]]
        local_mhm_nml_replace_dict[('evaluation_gauges', 'gauge_filename')] = [['{}.txt'.format(gauge_id)]]
    if start_date is not None and end_date is not None:
        local_mhm_nml_replace_dict[('time_periods', 'warming_days')] = [spinup]
        local_mhm_nml_replace_dict[('time_periods', 'eval_per')] = [
            {
                'ystart': start_date.year,
                'mstart': start_date.month,
                'dstart': start_date.day,
                'yend': end_date.year,
                'mend': end_date.month,
                'dend': end_date.day,
            }
        ]

    nml = _read_namelist(in_path)
    # set global settings
    for key, value in MHM_NML_REPLACE_DICT.items():
        nml[key[0]][key[1]] = value
    # set "local" settings
    for key, value in local_mhm_nml_replace_dict.items():
        nml[key[0]][key[1]] = value
    print('Creating {} file: {}'.format('config', out_path))
    nml.write(out_path, force=True)


def _get_item_index(nml, item_type, name, item):
    names = nml[item_type][name]
    # it fails if item not in names
    return names.index(item)


def _configure_mpr_nml(in_path, out_path, coords=None, data_arrays=None, out_filename=''):
    if data_arrays is None:
        data_arrays = {}
    if coords is None:
        coords = {}
    nml = _read_namelist(in_path)
    for item, item_type, item_label in zip(
            [coords, data_arrays],
            ['coordinates', 'data_arrays'],
            ['coord_name', 'name']):
        for item_name, item_props in item.items():
            item_index = _get_item_index(nml, item_type=item_type, name=item_label, item=item_name)
            for prop, value in item_props.items():
                nml[item_type][prop][item_index] = value
    if out_filename:
        nml['main']['out_filename'] = out_filename
    print('Creating {} file: {}'.format('config', out_path))
    nml.write(out_path, force=True)



# CLASSES

# SCRIPT
if __name__ == '__main__':
    input_dir = pathlib.Path('../../mhm_fork2/check')
    output_dir = pathlib.Path('../check_new')
    if input_dir.is_file():
        raise Exception("Input directory must be a directory, it is a file")
    nc_files = get_all_subfiles(input_dir, suffixes=['.nc'])
    # nc_files = []
    for nc_file in nc_files:
        print('working on file: {}'.format(nc_file))
        output_file = pathlib.Path(output_dir, nc_file)
        sort_y_dim(pathlib.Path(input_dir, nc_file), output_file, y_dims=('ncols', 'northing'))
    nml_files = get_all_subfiles(input_dir, suffixes=['.nml', '.out'])
    for nml_file in nml_files:
        print('working on file: {}'.format(nml_file))
        output_file = pathlib.Path(output_dir, nml_file)
        shutil.copy(src=pathlib.Path(input_dir, nml_file), dst=output_file)



