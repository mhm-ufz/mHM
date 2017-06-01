#!/usr/bin/env python

#import numpy as np

from scipy import integrate

def get_U_num(t, Q_in, Q_out, U_0):

#    t_no = len(t)
#    t = np.linspace(1,t_no,t_no)

    return U_0 + integrate.cumtrapz( Q_in - Q_out, t, initial=0 )