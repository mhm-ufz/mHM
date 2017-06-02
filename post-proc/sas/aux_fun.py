
import numpy as np
from get_theta import *
from get_p import *

class MapPoint(object):

    def __init__(self, x_array, y_array, lon, lat, x_pos, y_pos):
        
        self.x_i = (np.abs(x_array - x_pos)).argmin()
        self.y_i = (np.abs(y_array - y_pos)).argmin()
        self.x = x_array[self.x_i]
        self.y = y_array[self.y_i]
        self.lon = lon[self.y_i, self.x_i]
        self.lat = lat[self.y_i, self.x_i]
    
    def set_soil_moisture(self, numpy_array):
    
        self.soil_moisture = numpy_array[:, self.y_i, self.x_i]
        
    def set_ET(self, numpy_array):
    
        self.ET = numpy_array[:, self.y_i, self.x_i]
    
    def set_I(self, numpy_array):
    
        self.I = numpy_array[:, self.y_i, self.x_i]    
            
    def set_groundwater(self, numpy_array):
    
        self.groundwater = numpy_array[:, self.y_i, self.x_i]
        
    def get_ttd(self, U_array, Q_out_array, ET_array, t):
            
        U = U_array.data[:, self.y_i, self.x_i]
        Q_out = Q_out_array.data[:, self.y_i, self.x_i]
        ET = ET_array.data[:, self.y_i, self.x_i]
        
        
        t_no = len(t)
        tau = np.linspace(1, t_no, t_no)
        tau_no = len(tau)
        
        self.theta = np.zeros(tau_no)
        self.eta = np.zeros(tau_no)
        self.p_tt_mean = np.zeros(tau_no)        
        self.p_tt_var = np.zeros(tau_no) 
        self.mean_p_tt = np.zeros(tau_no)

        for t_in in range(0, t_no):
            
            self.theta[t_in] = get_theta(U, Q_out, ET, t, t_in)
            self.eta[t_in] = get_theta(U, ET, Q_out, t, t_in)
            p_tt = get_pdf_tt(U, Q_out, ET, t, t_in)
            self.p_tt_mean[t_in] = np.sum(p_tt*np.linspace(1, t_no - t_in, t_no - t_in))
            self.p_tt_var[t_in] = np.sum(p_tt*np.linspace(1, t_no - t_in, t_no - t_in)**2) - self.p_tt_mean[t_in]**2
            self.mean_p_tt[0:t_no-t_in] = self.mean_p_tt[0:t_no-t_in] + p_tt
        
        self.mean_p_tt = self.mean_p_tt/t_no      
        
#        t_in = 10
#        p_tt = get_pdf_tt(U, Q_out, ET, t, t_in)
#        print( np.sum(p_tt*np.linspace(1, t_no - t_in, t_no - t_in)) )
#        plot(self.theta)
#        show()
            
class DataStruct( object ):
    
    def __init__(self, numpy_array):
        
        tmp = numpy_array.shape
        y_no = tmp[0]
        x_no = tmp[1]
        
        mask = np.zeros( ( y_no, x_no ), dtype=bool)
        
        for x_i in range(0, x_no):
            for y_i in range(0, y_no):
                if numpy_array[y_i][x_i] == -9999:
                    mask[y_i][x_i] = True
                else:
                    mask[y_i][x_i] = False
        
        self.field = np.ma.array(numpy_array, mask = mask)
#        print(mask)
#        print(tmp)

def get_nearest_value(numpy_array, value):
    
    nearest_index = (np.abs(numpy_array - value)).argmin()
    nearest_value = numpy_array[nearest_index]
    
    return nearest_value
    
def get_nearest_index(numpy_array, value):
    
    nearest_index = (np.abs(numpy_array - value)).argmin()
    
    return nearest_index
    
    
def write_data(data):
    
    eol = '\n'
    eov = '\t'
    
    f_path  = 'elevation.txt'
    shape = data.shape

    row_no = shape[1]    
    col_no = shape[0]
    
        
    f_id    = open(f_path, 'w')
    
    for row_i in range (0, row_no):
        
        output_line = ''
        for col_i in range(0, col_no):
#            if data[col_i, row_i] == False:
#                output_line = output_line + str(1) + eov
#            elif data[col_i, row_i] == True:
#                output_line = output_line + str(0) + eov
            output_line = output_line + str(data[col_i, row_i]) + eov

        f_id.write( output_line + eol)

    f_id.close()
    