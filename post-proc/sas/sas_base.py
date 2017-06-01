
import numpy as np

from mhm import MHM
from sas.get_p import get_p_forward
from sas.get_p import get_p_backward
from sas.get_theta import get_theta

class SAS(MHM):
        
    def __init__(self, model_path, file_name = 'mhm.nml', 
                 rel_path_name = False):
        MHM.__init__(self, model_path, file_name, rel_path_name)
        
    def map_variables(self):
        
        var_dict = {}
        var_dict['x_3_array'] = ['SWC']
        var_dict['x_5_array'] = ['unsatSTW']
        var_dict['x_6_array'] = ['satSTW']
        var_dict['ET_array'] = ['aET']
        var_dict['QI'] = ['QIf', 'QIs']
        var_dict['BF_array'] = ['QB']
        var_dict['R_array'] = ['recharge']
#        var_dict['root'] = ['SWC']
#        var_dict['soil'] = ['SWC', 'unsatSTW']
#        var_dict['subsurface'] = ['SWC', 'unsatSTW', 'satSTW']
#        var_dict['Q_in_array'] = ['preEffect']

        for elem in var_dict:
            tmp = 0
            for sub_elem in var_dict[elem]:
                tmp += getattr(self, sub_elem )
            setattr(self, elem, tmp )
        
##-- determining pdf's --------------------------------------------------------

    def get_p(self, p_type = 'forward', flag = 'soil', flux_type = 'Q',
              time_series = 'all', fit_p = True):
        
        self.map_variables()
        self.set_spatial_mask('L1')
        
        if time_series == 'all':
            self.time_series_flag = np.ones((self.t_no))
        elif time_series == 'storm_events':
            self.get_storm_events()
        elif time_series == 'dry_years':
            self.get_wet_dry_years(time_series)
        elif time_series == 'wet_years':
            self.get_wet_dry_years(time_series)   
             
        if flag == 'root':
            S = self.SWC
            Q_1 = self.J
            Q_2 = self.aET
        elif flag == 'soil':
            S = self.SWC + self.unsatSTW
            Q_1 = self.QIf + self.QIs + self.recharge
            Q_2 = self.aET
        elif flag == 'subsurface':
            S = self.SWC + self.unsatSTW + self.satSTW
            Q_1 = self.QIf + self.QIs + self.OB
            Q_2 = self.aET

        if p_type == 'forward':
            if flux_type == 'Q':
                self.get_p_Q_forward(S, Q_1, Q_2)
            elif flux_type == 'ET':
                self.get_p_Q_forward(S, Q_2, Q_1)
        elif p_type == 'backward':
            if flux_type == 'Q':
                self.get_p_Q_backward(flag)
            elif flux_type == 'ET':
                self.get_p_Q_backward(flag)
        elif p_type == 'theta':
            self.get_p_theta(S, Q_1, Q_2)
            
        if fit_p == True:
            self.fit_p()

    def get_p_Q_forward(self, S, Q_1, Q_2):
        
        self.mean_p_Q = np.zeros((self.t_no, self.y_no, self.x_no))
        self.my_p_Q = np.zeros((self.t_no))                

        for x_i in range(0, self.x_no):
            for y_i in range(0, self.y_no):
                if self.mask['L1_special'][y_i, x_i]:
                    continue
                p_Q = np.zeros((self.t_no))            
                counter = 0
                for t_in in range(0, self.t_no):
                    T = self.t_no - t_in
                    if self.time_series_flag[t_in] > 0:
                        p_Q[:T] += get_p_forward(S[:, y_i, x_i], 
                                                 Q_1[:, y_i, x_i], 
                                                 Q_2[:, y_i, x_i], 
                                                 self.t, t_in)
                        counter += 1
                    if (t_in == 10) and (x_i == 40) and (y_i == 40):
#                        print(t_in)
                        self.my_p_Q = get_p_forward(S[:, y_i, x_i], 
                                                    Q_1[:, y_i, x_i], 
                                                    Q_2[:, y_i, x_i], 
                                                    self.t, t_in)
                if (x_i == 40) and (y_i == 40):
                    self.p_Q = p_Q/counter
                self.mean_p_Q[:,y_i,x_i] = p_Q/counter
               
    def get_p_Q_backward(self, flag):
        
        self.mean_p_Q = np.zeros((self.t_no, self.y_no, self.x_no))
        
        for x_i in range(0, self.x_no):
            for y_i in range(0, self.y_no):
                if self.mask['L1'][y_i, x_i]:
                    continue
                if flag == 'root':
                    S = self.x_3_array.data[:, y_i, x_i]
                    J = self.Q_in_array[:, y_i, x_i]
                    Q_1 = self.Q_out_33_array[:, y_i, x_i]
                    Q_2 = self.ET_array[:, y_i, x_i] 
                elif flag == 'soil':
                    S = self.SWC[:, y_i, x_i] + self.unsatSTW[:, y_i, x_i]
                    J = self.preEffect[:, y_i, x_i]
                    Q_1 = self.QI[:, y_i, x_i] + self.recharge[:, y_i, x_i]
                    Q_2 = self.aET[:, y_i, x_i]
                p_Q = np.zeros((self.t_no))          
#                counter = 0
                for t_ex in range(1, self.t_no):
                    p_Q[:t_ex+1] += get_p_backward(S, J, Q_1, Q_2, self.t, t_ex)
                self.mean_p_Q[:,y_i,x_i] = p_Q/(self.t_no - 1)
#                if (x_i == 20) and (y_i == 20):
#                    self.p_Q = p_Q/(self.t_no - 1)
                
    def get_p_rt(self):
        
        self.map_variables()
        
        self.p_rt_tot = np.zeros((self.t_no, self.t_no))
        self.p_a_tot = np.zeros((self.t_no, self.t_no))
        self.omega = np.zeros((self.t_no, self.t_no))
        
        soil_moisture = self.SWC + self.unsatSTW
        discharge = self.QI + self.recharge
#        soil_moisture_tot = np.mean(np.mean(soil_moisture, axis = 1), axis = 1)
#        self.soil_moisture_rel = np.zeros((self.t_no, self.y_no, self.x_no))
        self.soil_moisture_rel = np.ma.copy((self.SWC))
        self.discharge_rel = np.ma.copy((self.SWC))
        for t_i in range(0, self.t_no):
            soil_moisture_tot = np.sum(np.sum(soil_moisture[t_i,:,:]))
            self.soil_moisture_rel.data[t_i, :, :] = soil_moisture[t_i, :, :]/soil_moisture_tot
            discharge_tot = np.sum(np.sum(discharge[t_i,:,:]))
            self.discharge_rel.data[t_i, :, :] = discharge[t_i, :, :]/discharge_tot
#        print(soil_moisture.shape)
#        print(soil_moisture_tot.shape)
#        discharge = self.QIf + self.QIs + self.recharge
#        discharge_tot = np.mean(np.mean(discharge, axis = 1), axis = 1)
        
        for t_ex in range(1, self.t_no):
            p_Q = np.zeros((self.t_no))
            p_rt = np.zeros((self.t_no))
            p_a = np.zeros((self.t_no))
            for x_i in range(0, self.x_no):
                for y_i in range(0, self.y_no):
                    if self.mask['L1'][y_i, x_i]:
                        continue
                    S = self.SWC[:, y_i, x_i] + self.unsatSTW[:, y_i, x_i]
                    J = self.preEffect[:, y_i, x_i]
                    Q_1 = self.QI[:, y_i, x_i] + self.recharge[:, y_i, x_i]
                    Q_2 = self.aET[:, y_i, x_i]
                    
#                    T = self.t_no - t_ex
#                    p_Q[:T] = get_p_forward(S, Q_1, Q_2, self.t, t_ex)
                    p_Q[:t_ex+1] = get_p_backward(S, J, Q_1, Q_2, self.t, t_ex)
#                    ratio = soil_moisture[t_ex, y_i, x_i]/soil_moisture_tot[t_ex]
#                    print(ratio.shape)
#                    print(p_Q.shape)
                    p_rt += p_Q*self.soil_moisture_rel[t_ex,y_i,x_i]
                    p_a += p_Q*self.discharge_rel[t_ex,y_i,x_i]
                    
            self.p_rt_tot[:,t_ex] = p_rt
            self.p_a_tot[:,t_ex] = p_a
#            self.omega[:,t_ex] = p_a/p_rt
            
    def get_p_theta(self, S, Q_1, Q_2):

        self.theta = np.zeros((self.t_no, self.y_no, self.x_no))
        self.eta = np.zeros((self.t_no, self.y_no, self.x_no))
        
        for x_i in range(0, self.x_no):
            for y_i in range(0, self.y_no):
                if self.mask['L1_special'][y_i, x_i]:
                    continue
                for t_in in range(0, self.t_no):
                    self.theta[t_in, y_i, x_i] = get_theta(S[:, y_i, x_i], Q_1[:, y_i, x_i], Q_2[:, y_i, x_i], self.t, t_in)
                    self.eta[t_in, y_i, x_i] = get_theta(S[:, y_i, x_i], Q_2[:, y_i, x_i], Q_1[:, y_i, x_i], self.t, t_in)

        
                   
##-- fitting pdf to analytical functions --------------------------------------
                    
    def fit_p(self, func = 'exp'):
        
        from scipy.optimize import curve_fit
        from scipy.special import gamma
                  
        def exp_dist(x, a, b):
            return a*np.exp(-x/b)
        def gamma_dist(x, a, b):
            return (1/b)**a/gamma(a)*x**(a-1)*np.exp(-x/b)

        func = func[:3].lower()
        tmp = np.zeros((self.nrows['L1'], self.ncols['L1'] ))
        self.mean_tt_opt = np.ma.array(tmp, mask = self.mask['L1_special'])
   
        for x_i in range(0, self.x_no):
            for y_i in range(0, self.y_no):
                if self.mask['L1_special'][y_i][x_i]:
                    continue
                if func == 'exp':
                    popt, pcov = curve_fit(exp_dist, self.t[:100],
                                           self.mean_p_Q[:100, y_i, x_i])
                if func == 'gam':
                    popt, pcov = curve_fit(gamma_dist, self.t[:100],
                                           self.mean_p_Q[:100, y_i, x_i])                           
                self.mean_tt_opt[y_i, x_i] = popt[1]
                
##-- determining mask for evaluating the TTDs ---------------------------------
    
    def set_spatial_mask(self, mask_array, value = 0):
        
        data_level = 'L1'
        if 'L1_special' not in self.mask:
            self.mask['L1_special'] = np.copy(self.mask[data_level])
        
#        index = np.where( self.landcover_L1 == mask_id)
        if mask_array == 'non_urban':
            data_type = 'landcover'
            self.upscale_data(data_type, data_level)
            self.mask['L1_special'][np.where(self.landcover_L1 == 2)] = 1
        if mask_array == 'urban':
            data_type = 'landcover'
            self.upscale_data(data_type, data_level, flag = 'mean')
#            self.mask['L1_special'][np.where(self.landcover_L1 == 1)] = 1
#            self.mask['L1_special'][np.where(self.landcover_L1 == 3)] = 1
            self.mask['L1_special'][np.where(self.mask['L1_special'] == 0)] = 1
            self.mask['L1_special'][np.where(self.landcover_L1 == 2)] = 0
        if mask_array == 'forest':
            data_type = 'landcover'
            self.upscale_data(data_type, data_level)
            self.mask['L1_special'][np.where(self.landcover_L1 == 2)] = 1
            self.mask['L1_special'][np.where(self.landcover_L1 == 3)] = 1   
        if mask_array == 'grass':
            data_type = 'landcover'
            self.upscale_data(data_type, data_level)
            self.mask['L1_special'][np.where(self.landcover_L1 == 1)] = 1
            self.mask['L1_special'][np.where(self.landcover_L1 == 2)] = 1
        if mask_array == 'lai_class_1':
            data_type = 'lai_class'
            self.upscale_data(data_type, data_level)
            self.mask['L1_special'][:,:] = 1
            self.mask['L1_special'][np.where(self.lai_class_L1 == 1)] = 0   
        if mask_array == 'lai_class_10':
            data_type = 'lai_class'
            self.upscale_data(data_type, data_level)
            self.mask['L1_special'][:,:] = 1
            self.mask['L1_special'][np.where(self.lai_class_L1 == 10)] = 0
        if mask_array == 'soil_class':
#            self.mask['L1_special'] = np.copy(self.mask[data_level])
            data_type = 'soil_class'
            self.upscale_data(data_type, data_level)
            self.mask['L1_special'][:,:] = 1
            self.mask['L1_special'][np.where(self.soil_class_L1 == value)] = 0
            
            
#        self.mask['L1_special'] = np.copy(self.mask['L1'])
        
##-- determining weather extremes ---------------------------------------------
                             
    def get_wet_dry_years(self, time_series):
        
        time_step = 12
        self.averaged_x_3 = np.zeros((self.t_no))
#        self.dry_year = np.zeros((self.t_no))
#        self.wet_year = np.zeros((self.t_no))
        self.time_series_flag = np.zeros((self.t_no))
        
#        tmp = np.mean(self.SWC, axis=1)
        self.mean_soil_moisture = np.mean(np.mean(self.SWC, axis=1), axis=1)
        my_soil_moisture = (self.mean_soil_moisture - np.mean(self.mean_soil_moisture))/np.std(self.mean_soil_moisture)
        
        for t_i in range(0, self.t_no/time_step):
            a = t_i*time_step
            e = t_i*time_step + time_step
            tmp = np.mean(my_soil_moisture[a:e])
            tmp = np.ones(( time_step ))*tmp
            self.averaged_x_3[a:e] = tmp
            if time_series == 'wet_years':
                if tmp[0] > 0.605:
                    self.time_series_flag[a:e] = np.ones(( time_step ))
                else:
                    self.time_series_flag[a:e] = np.zeros(( time_step ))
            elif time_series == 'dry_years':
                if tmp[0] < -0.426:
                    self.time_series_flag[a:e] = np.ones(( time_step ))
                else:
                    self.time_series_flag[a:e] = np.zeros(( time_step ))

    def get_storm_events(self):

        self.time_series_flag = np.zeros((self.t_no))
        self.mean_preEffect = np.mean(np.mean(self.preEffect, axis=1), axis=1)

        for t_i in range(0, self.t_no):
            if self.mean_preEffect[t_i] > 100:
                self.time_series_flag[t_i] = 1

##-- compressing time series --------------------------------------------------

    def map_to_month(self, flux_data_list, state_data_list):
        chunk_no = 33
        self.compress_flux(flux_data_list, chunk_no)
        self.compress_state(state_data_list, chunk_no)
        self.t_no = self.t_no/chunk_no
        self.t = range(1, self.t_no + 1)


    def compress_flux(self, data_list, chunk_no):
        for data_i in data_list:
            setattr(self, data_i, self.sum_chunk(getattr(self, data_i),
                                                 chunk_no, 0))
        
    def compress_state(self, data_list, chunk_no):
        for data_i in data_list:
            setattr(self, data_i, self.mean_chunk(getattr(self, data_i),
                                                  chunk_no, 0))
        
    def sum_chunk(self, x, chunk_size, axis):
        
        shape = x.shape
        if axis < 0:
            axis += x.ndim
        shape = shape[:axis] + (-1, chunk_size) + shape[axis+1:]
        x = x.reshape(shape)
        return x.sum(axis=axis+1)
        
    def mean_chunk(self, x, chunk_size, axis):
        
        shape = x.shape
        if axis < 0:
            axis += x.ndim
        shape = shape[:axis] + (-1, chunk_size) + shape[axis+1:]
        x = x.reshape(shape)
        return x.mean(axis=axis+1)
        
    
    