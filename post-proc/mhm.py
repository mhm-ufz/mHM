#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from readnetcdf import readnetcdf

class MHM(object):
    
    """
        MHM python class


        Example use
        --------
        # import class
        >>> from mhm import MHM

        # init class
        >>> model_path = path_to_my_mhm_model
        >>> mhm = MHM(model_path)

        # import data
        >>> mhm.import_data('lat_lon')
        >>> mhm.import_data('dem')
        >>> mhm.import_data('states_and_fluxes', 'all')
        
        # process data (getting the average of the actual evapotranspiration)
        >>> import numpy as np
        >>> mhm.combine_variables(['aET'])
        >>> mean_ET = np.mean(mhm.aET, axis=0)

        # plot data (digital elevation model)
        >>> import matplotlib as plt
        >>> plt.imshow(mhm.dem)
        >>> plt.show()


        License
        -------
        This file is part of the UFZ Python package.

        The UFZ Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The UFZ Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2015 Falk Hesse


        History
        -------
        Written,  FH, Nov 2015
    """
    
##-- init function ------------------------------------------------------------

    def __init__(self, model_path, file_name = 'mhm.nml', 
                 rel_path_name = False):
        
        # to do: - option for restart file
        #        - both implicite and explicite
        #        - adapt for use with several basins
        #        - onyl simple upscaling
        #        - documentation
        #        - doctests
        
        if model_path[-1] == '/':
            self.model_path = model_path
        else:
            self.model_path = model_path + '/'
        
        self.nml = {'nml_path' : self.model_path + file_name}
        self.output_nml = {'nml_path' : self.model_path + 'mhm_outputs.nml'}
        self.parameter_nml = {'nml_path' : self.model_path + 'mhm_parameter.nml'}
        
        self.cellsize = {}
        self.cellsize_ratio = {}
        self.ncols = {}
        self.nrows = {}
        self.mask = {}
        self.lat = {}
        self.lon = {}
        
        # reading data from the 'mhm.nml' file for general information
         
        f_id = open(self.model_path + file_name)
        for line_i in f_id:
            line = line_i.strip()
            if not (line.startswith("!") or line.startswith("&") ):
                if ('=' in line):
                    line = line.replace('"', '')
                    line = line.split()
                    if len(line) == 3:
                        self.nml[line[0]] = line[2]
        f_id.close()
        
        self.nBasins = int(self.nml['nBasins'])

        
        if   self.nml['timeStep_sm_input'] == '-1':
            self.t_stepsize = 'daily'
        elif self.nml['timeStep_sm_input'] == '-2':
            self.t_stepsize = 'monthly'
        elif self.nml['timeStep_sm_input'] == '-3':
            self.t_stepsize = 'yearly'
        
        if rel_path_name:
            for k,v in self.nml.items():
                if ('dir_' in k) or ('file_' in k):
                    self.nml[k] = model_path + self.nml[k]
        
        # reading data from the 'dem.asc' file for domain information   
             
        f_id = open(self.nml['dir_Morpho(1)'] + 'dem.asc')
        
        for i, line in enumerate(f_id):
            if i == 4:
                self.cellsize['L0'] = int(line[9:])
        self.cellsize['L1'] = int(self.nml['resolution_Hydrology(1)'])
        self.cellsize['L2'] = 4000
        f_id.close()
        
        self.cellsize_ratio['L0'] = 1
        self.cellsize_ratio['L1'] = float(self.nml['resolution_Hydrology(1)'])/self.cellsize['L0']
        self.cellsize_ratio['L2'] = self.cellsize['L2']/self.cellsize['L0']
        
        f_id = open(self.nml['dir_Morpho(1)'] + 'dem.asc')
        
        for line_i in range(6):
            line = f_id.next().strip()
            line = line.split()
            if   line[0] == 'ncols':
                self.ncols['L0'] = int(line[1])
                self.ncols['L1'] = int(float(line[1])/self.cellsize_ratio['L1'])
                self.ncols['L2'] = int(float(line[1])/self.cellsize_ratio['L2'])          
            elif line[0] == 'nrows':
                self.nrows['L0'] = int(line[1])
                self.nrows['L1'] = int(float(line[1])/self.cellsize_ratio['L1'])
                self.nrows['L2'] = int(float(line[1])/self.cellsize_ratio['L2'])
            elif line[0] == 'xllcorner':
                self.xllcorner = float(line[1])
            elif line[0] == 'yllcorner':
                self.yllcorner = float(line[1])
            elif line[0] == 'NODATA_value':
                self.NODATA_value = float(line[1])    
        
        self.mask['L0'] = np.zeros((self.nrows['L0'], self.ncols['L0']))
        self.mask['L1'] = np.zeros((self.nrows['L1'], self.ncols['L1']))
        self.mask['L2'] = np.zeros((self.nrows['L2'], self.ncols['L2']))
        
        for line_i in range(6, 6 + self.nrows['L0']):
            line = f_id.next().strip()
            line = line.split()
            for col_i in range(0, self.ncols['L0']):
                if float(line[col_i]) == float(self.NODATA_value):
                    self.mask['L0'][line_i - 6][col_i] = 1
        for row_i in range(0, self.nrows['L1']):
            for col_i in range(0, self.ncols['L1']):
                row_a = int(row_i*self.cellsize_ratio['L1'])
                row_e = int((row_i + 1)*self.cellsize_ratio['L1'])
                col_a = int(col_i*self.cellsize_ratio['L1'])
                col_e = int((col_i + 1)*self.cellsize_ratio['L1'])
                if np.mean(self.mask['L0'][row_a:row_e,col_a:col_e]) == 1: 
                    self.mask['L1'][row_i][col_i] = 1
        for row_i in range(0, self.nrows['L2']):
            for col_i in range(0, self.ncols['L2']):
                row_a = int(row_i*self.cellsize_ratio['L2'])
                row_e = int((row_i + 1)*self.cellsize_ratio['L2'])
                col_a = int(col_i*self.cellsize_ratio['L2'])
                col_e = int((col_i + 1)*self.cellsize_ratio['L2'])          
                if np.mean(self.mask['L0'][row_a:row_e,col_a:col_e]) == 1: 
                    self.mask['L2'][row_i][col_i] = 1            
        f_id.close()
        
##-- importing functions ------------------------------------------------------
     
    def import_data(self, data_type, *kwargs ):
        if data_type == 'states_and_fluxes':            
            self.import_states_and_fluxes( *kwargs )
        elif data_type == 'lat_lon_L0':
            data_path = str(kwargs[0])
            self.lon_L0 = self.import_lat_lon(data_path, 'lon')
            self.lat_L0 = self.import_lat_lon(data_path, 'lat')
        elif data_type == 'lat_lon':
            data_path = self.nml['file_LatLon(1)']
            self.lon_L1 = self.import_lat_lon(data_path, 'lon')
            self.lat_L1 = self.import_lat_lon(data_path, 'lat')
        elif data_type == 'lat_lon_L2':
            data_path = str(kwargs[0])
            self.lon_L2 = self.import_lat_lon(data_path, 'lon')
            self.lat_L2 = self.import_lat_lon(data_path, 'lat') 
        elif data_type == 'landcover':
            data_path = self.nml['dir_LCover(1)'] + str(kwargs[0])
            setattr(self, data_type, self.import_L0_data(data_path))
        elif data_type == 'lai_class':
            data_path = self.nml['dir_Morpho(1)'] + 'LAI_class.asc'
            self.lai_class = self.import_L0_data( data_path )
        elif data_type == 'slope':
            data_path = self.nml['dir_Morpho(1)'] + 'slope.asc'
            setattr(self, data_type, self.import_L0_data(data_path))
        elif data_type == 'soil_class':
            data_path = self.nml['dir_Morpho(1)'] + 'soil_class.asc'
            setattr(self, data_type, self.import_L0_data(data_path))
        elif data_type == 'pre':
            data_path = 'meteo/pre/pre.nc'
            self.pre = self.import_precipitation( data_path )
        elif data_type == 'dem':
            data_path = self.nml['dir_Morpho(1)'] + 'dem.asc'
            setattr(self, data_type, self.import_L0_data(data_path)) 
        elif data_type == 'facc':
            data_path = self.nml['dir_Morpho(1)'] + 'facc.asc'
            setattr(self, data_type, self.import_L0_data(data_path)) 
        
    def import_lat_lon(self, f_path, var):
        return readnetcdf(f_path, var=var) 
        
    def import_precipitation(self, f_path):
        path = self.nml['dir_Precipitation(1)'] + 'pre.nc'
        return readnetcdf(path, var='pre')
    
    def import_states_and_fluxes(self, *kwargs):
        
        f_path = self.nml['dir_Out(1)'] + 'mHM_Fluxes_States.nc'
#        print(f_path)
    
        if kwargs[0] == 'all':
            var_list = [str(i) for i in readnetcdf(f_path, variables=True)]
        elif kwargs[0] == 'default':
            var_list = ['time', 'northing', 'easting', 'lon', 'lat', 
                        'SWC_L01', 'SWC_L02', 'SWC_L03',
                        'SM_L01', 'SM_L02', 'SM_L03', 
                        'unsatSTW', 'satSTW', 'aET', 'QIf', 'QIs', 'QB',
                        'recharge', 'aET_L01', 'aET_L02', 'aET_L03',
                        'preEffect']
        else:
            var_list = kwargs[0]
#        print(var_list)
        for var_i in range(0, len(var_list)):
            var = var_list[var_i]
            setattr(self, var, readnetcdf(f_path, var=var))
        
        self.t_no = len( getattr(self, 'time') )
        self.x_no = self.lon.shape[1]
        self.y_no = self.lon.shape[0]       
             
        self.t = np.linspace(1, self.t_no, self.t_no)

    def import_L0_data(self, f_path):
               
        f_id = open(f_path)
        f_lines = f_id.readlines()
        f_id.close()
        
        tmp = np.zeros( ( self.nrows['L0'], self.ncols['L0'] ) )
        
        for line_i in range(0, self.nrows['L0']):
            tmp[line_i] = f_lines[line_i + 6].split()    
        
        return np.ma.array(tmp, mask = self.mask['L0'])
        
##-- importing restart file ---------------------------------------------------

    def import_restart_file(self, *kwargs):
        
        f_path = self.nml['dir_Out(1)'] + 'mHM_restart_001.nc'
        var_list = [str(i) for i in readnetcdf(f_path, variables=True)]
#        print(var_list)
        for var_i in range(0, len(var_list)):
            var = var_list[var_i]
            setattr(self, var, readnetcdf(f_path, var=var))

##-- upscaling and downscaling functions --------------------------------------

    def upscale_data(self, data_type, data_level = 'L1', flag = 'any'):      
        setattr(self, data_type + '_L1',
                self.upscale_L1_data(getattr(self, data_type), flag))
        
    def upscale_L1_data(self, L0_array, flag):
        
        L1_array = np.ma.array(np.zeros((self.nrows['L1'], 
                                         self.ncols['L1'] )), 
                                         mask = self.mask['L1'])
#        print(self.nrows['L1']*self.ncols['L1'] - np.sum(L1_array.mask))
        num = 0
        for row_i in range(0, self.nrows['L1']):
            for col_i in range(0, self.ncols['L1']):
                if self.mask['L1'][row_i, col_i]:
                    continue
                row_s = int(row_i*self.cellsize_ratio['L1'])
                row_e = int((row_i + 1)*self.cellsize_ratio['L1'])
                row_c = int((row_s + row_e)/2)
                col_s = int(col_i*self.cellsize_ratio['L1'])
                col_e = int((col_i + 1)*self.cellsize_ratio['L1'])
                col_c = int((col_s + col_e)/2)
                num += 1
                if flag == 'any':
                    if np.any(L0_array[row_s:row_e,col_s:col_e] == 2):
                        L1_array[row_i][col_i] = 2
                    else:
                        L1_array[row_i][col_i] = L0_array[row_c, col_c]
                elif flag == 'mean':
#                    L1_array[row_i][col_i] = L0_array[row_c, col_c]
                    L1_array[row_i][col_i] = np.mean(L0_array[row_s:row_e,col_s:col_e])
                elif flag == 'most':
                    L1_array[row_i][col_i] = np.mean(L0_array[row_s:row_e,col_s:col_e])
#        print(num)
#        print(self.nrows['L1']*self.ncols['L1'] - np.sum(L1_array.mask))
        return L1_array
        
    def downscale_data(self, data_type):
        
        pre_L1 = np.zeros((self.nrows['L1'], self.ncols['L1']))
        pre = np.mean(self.pre, axis=0)
        L2 = self.cellsize_ratio['L2']
#        L2 = 10
        ratio = int(L2/self.cellsize_ratio['L1'])
        
        for row_i in range(0, self.nrows['L1']):
            for col_i in range(0, self.ncols['L1']):
                pre_L1[row_i, col_i] = pre[row_i/ratio, col_i/ratio]
        
        self.pre_L1 = np.ma.array(pre_L1, mask = self.mask['L1'])
#        self.pre_L1 = np.ma.array( np.mean(self.preEffect, axis = 0), mask = self.mask['L1'])
        
        
##-- flow direction functions -------------------------------------------------

        
    def flow_dir(self, col, row):
        
        self.flow_path = {}
        self.flow_path_backward = {}
        
        self.dem = self.import_L0_data(self.nml['dir_Morpho(1)'] + 'dem.asc')
        self.facc = self.import_L0_data(self.nml['dir_Morpho(1)'] + 'facc.asc')
        self.fdir = self.import_L0_data(self.nml['dir_Morpho(1)'] + 'fdir.asc')
        
        self.mask['basin_L0'] = np.ones((self.nrows['L0'], self.ncols['L0']))
        
        outlet = np.where(self.facc == np.max(self.facc))        
        
        for row_i in range(0, self.nrows['L0']):
            for col_i in range(0, self.ncols['L0']):
                if self.mask['L0'][row_i][col_i]:
                    continue
                pos_i = (row_i, col_i)
                self.flow_path[pos_i] = self.get_flow_path(pos_i, outlet)
        
        basin_col = []
        basin_row = []
        for row_i in range(0, self.nrows['L0']):
            for col_i in range(0, self.ncols['L0']):
                if self.mask['L0'][row_i][col_i]:
                    continue
                path = self.flow_path[ (row_i, col_i) ]
                rows = np.where( path[0,:] == row )
                cols = np.where( path[1,:] == col )
                tmp = np.intersect1d( rows, cols)
                if tmp:
                    basin_col.append(col_i)
                    basin_row.append(row_i)
                    self.mask['basin_L0'][row_i][col_i] = 0
            
        
    def get_flow_path(self, start, end):
        
        flow_path_x = []
        flow_path_y = []
        y_i = start[0]
        x_i = start[1]
        
#        print(start)
        while (y_i != end[0]) and (x_i != end[1]):
            if self.fdir[y_i, x_i] == 1:
                y_i += 0
                x_i += 1
            elif self.fdir[y_i, x_i] == 2:
                y_i += 1
                x_i += 1
            elif self.fdir[y_i, x_i] == 4:
                y_i += 1
                x_i += 0
            elif self.fdir[y_i, x_i] == 8:
                y_i += 1
                x_i += -1
            elif self.fdir[y_i, x_i] == 16:
                y_i += 0
                x_i += -1
            elif self.fdir[y_i, x_i] == 32:
                y_i += -1
                x_i += -1
            elif self.fdir[y_i, x_i] == 64:
                y_i += -1
                x_i += 0
            elif self.fdir[y_i, x_i] == 128:
                y_i += -1
                x_i += 1
            flow_path_x.append(x_i)
            flow_path_y.append(y_i)
        flow_path_x.append(end[1])
        flow_path_y.append(end[0])
        
        f_path = np.zeros((2, len(flow_path_x)))
        f_path[0,:] = flow_path_y
        f_path[1,:] = flow_path_x
        
        return f_path

##-- misc ---------------------------------------------------------------------

    def combine_variables(self, my_list):
        
        f_path = self.nml['dir_Out(1)'] + 'mHM_Fluxes_States.nc'
        var_list = [str(i) for i in readnetcdf(f_path, variables=True)]
#        print(var_list)
        var_dict = {}
        for elem in my_list:
            var_dict[elem] = [s for s in var_list if elem + '_L0' in s]
            tmp = 0
            for sub_elem in var_dict[elem]:
                tmp += getattr(self, sub_elem )
            setattr(self, elem, tmp )
    