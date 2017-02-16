
from pylab import *
from sas.get_p import get_p_forward
from sas.get_p import get_p_backward

class MHM_Cell(object):

    def __init__(self, mhm, x_pos, y_pos):
        
        self.x_i = (np.abs(mhm.x - x_pos)).argmin()
        self.y_i = (np.abs(mhm.y - y_pos)).argmin()
        self.x = mhm.x[self.x_i]
        self.y = mhm.y[self.y_i]
        self.lon = mhm.L1_lon[self.y_i, self.x_i]
        self.lat = mhm.L1_lat[self.y_i, self.x_i]
        
        
        self.t = mhm.t
        self.t_no = mhm.t_no        
        self.x_3 = mhm.x_3_array[:, self.y_i, self.x_i]
        self.J = mhm.Q_in_array[:, self.y_i, self.x_i]
        self.Q_out_3 = mhm.Q_out_33_array[:, self.y_i, self.x_i]
        self.ET = mhm.ET_array[:, self.y_i, self.x_i]
        
    def get_p_Q_root(self, p_type):
        if p_type == 'forward':
            self.get_p_Q_root_forward()
        elif p_type == 'backward':
            self.get_p_Q_root_backward()
            
    def get_p_Q_root_forward(self):
        
        self.p_Q_root = np.zeros( (self.t_no) )
        p_Q_root = np.zeros((self.t_no))
        
        for t_in in range(0, self.t_no):
            T = self.t_no - t_in
            p_Q_root[:T] += get_p_forward(self.x_3, self.Q_out_3, self.ET, self.t, t_in)
            
        self.p_Q_root = p_Q_root/self.t_no
        
#        t_in = 100
#        T = self.t_no - t_in
##        print(T)
#        p_Q_root[:T] = get_p_forward(self.x_3, self.Q_out_3, self.ET, self.t, t_in)
#        
#        t_in = 200
#        T = self.t_no - t_in
##        print(T)
#        p_Q_root[:T] = get_p_forward(self.x_3, self.Q_out_3, self.ET, self.t, t_in)
#                
#        t_in = 300
#        T = self.t_no - t_in
##        print(T)
#        p_Q_root[:T] = get_p_forward(self.x_3, self.Q_out_3, self.ET, self.t, t_in)
#        
#        t_in = 400
#        T = self.t_no - t_in
##        print(T)
#        p_Q_root[:T] = get_p_forward(self.x_3, self.Q_out_3, self.ET, self.t, t_in)
        
#        print(self.t)        
        
    def get_p_Q_root_backward(self):
        
        self.p_Q_root = np.zeros( (self.t_no) )
        p_Q_root = np.zeros((self.t_no))
        
        for t_ex in range(1, self.t_no):
            p_Q_root[:t_ex+1] += get_p_backward(self.x_3, self.J, self.Q_out_3, self.ET, self.t, t_ex)
                        
        self.p_Q_root = p_Q_root/(self.t_no - 1)
        
#        t_ex = 100
#        p_Q_root[:t_ex+1] = get_p_backward(self.x_3, self.J, self.Q_out_3, self.ET, self.t, t_ex)
#        
#        t_ex = 200
#        p_Q_root[:t_ex+1] = get_p_backward(self.x_3, self.J, self.Q_out_3, self.ET, self.t, t_ex)
#        
#        t_ex = 300
#        p_Q_root[:t_ex+1] = get_p_backward(self.x_3, self.J, self.Q_out_3, self.ET, self.t, t_ex)
#        
#        t_ex = 400
#        p_Q_root[:t_ex+1] = get_p_backward(self.x_3, self.J, self.Q_out_3, self.ET, self.t, t_ex)
    

def get_nearest_value(numpy_array, value):
    
    nearest_index = (np.abs(numpy_array - value)).argmin()
    nearest_value = numpy_array[nearest_index]
    
    return nearest_value
    
def get_nearest_index(numpy_array, value):
    
    nearest_index = (np.abs(numpy_array - value)).argmin()
    
    return nearest_index
    