#!/usr/bin/env python

import numpy as np
from scipy import integrate

from get_theta import get_theta
#from get_validity_range import get_validity_range

def get_p_forward(S, Q_out_1, Q_out_2, t, t_in):
    
    """
        Computes the forward travel-time pdf within a single reservoir with 
        one input flux and two output fluxes according to Botter et al.
        'Transport in the hydrological response: Travel time distributions, 
        soil moisture dynamics, and the old water paradox' WRR (2010), 
        Equation (35) and Botter et al. 'Catchtment residence time and travel 
        time distributions: The master equation' GRL (2011), Equation (9)


        Definition
        ----------
        def get_p_forward(S, Q_out_1, Q_out_2, t, t_in):


        Input
        -----
        S           volumetric water content of the reservoir
        Q_out_1     first output flux of reservoir (contributing to discharge)
        Q_out_2     second output flux of reservoir
        t           chronological time of reservoir dynamcis
        t_in        point in time when water is entering reservoir


        Output
        ------
        Float array of time-dependent travel-time distribution, representing 
        the PDF of the random variable T, i.e. the time from entering the
        reservoir at t_in to leaving the reservoir as discharge


        License
        -------
        This file is part of the UFZ Python package.
        
        The UFZ Python package is free software: you can redistribute it and/or 
        modify it under the terms of the GNU Lesser General Public License as 
        published by the Free Software Foundation, either version 3 of the 
        License, or (at your option) any later version.
        
        The UFZ Python package is distributed in the hope that it will be 
        useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.
        
        You should have received a copy of the GNU Lesser General Public 
        License along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.
        
        Copyright 2015 Falk He"sse


        History
        -------
        Written,  FH, Mar 2015 """
    
    theta = get_theta(S, Q_out_1, Q_out_2, t, t_in)
#    eta = get_theta(U, Q_out_2, Q_out_1, t, t_in)

    # needs to be amended later
    if theta == 0:
        theta = 0.00001

#    t_end = get_validity_range(t[t_in:], theta + eta)
    
    f = (Q_out_1[t_in:] + Q_out_2[t_in:])/S[t_in:]
    F = integrate.cumtrapz( f, t[t_in:], initial=0 )
    return Q_out_1[t_in:]/(S[t_in:]*theta)*np.exp( -F )

def get_p_backward(S, Q_in, Q_out_1, Q_out_2, t, t_ex):
    
    """
        Computes the backward travel-time pdf within a single reservoir with 
        one input flux and two output fluxes according to Botter et al.
        'Catchtment residence time and travel time distributions: The master 
        equation' GRL (2011), Equation (8)


        Definition
        ----------
        def get_p_forward(S, Q_out_1, Q_out_2, t, t_in):


        Input
        -----
        S           volumetric water content of the reservoir
        Q_in        input flux of reservoir (effective precipitation)
        Q_out_1     first output flux of reservoir (contributing to discharge)
        Q_out_2     second output flux of reservoir
        t           chronological time of reservoir dynamcis
        t_in        point in time when water is entering reservoir


        Output
        ------
        Float array of time-dependent travel-time distribution, representing 
        the PDF of the age distribution of the outflow flux at a given point in
        time.


        License
        -------
        This file is part of the UFZ Python package.
        
        The UFZ Python package is free software: you can redistribute it and/or 
        modify it under the terms of the GNU Lesser General Public License as 
        published by the Free Software Foundation, either version 3 of the 
        License, or (at your option) any later version.
        
        The UFZ Python package is distributed in the hope that it will be 
        useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.
        
        You should have received a copy of the GNU Lesser General Public 
        License along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.
        
        Copyright 2015 Falk He"sse


        History
        -------
        Written,  FH, Aug 2015 """


    f = (Q_out_1[t_ex::-1] + Q_out_2[t_ex::-1])/S[t_ex::-1]
    F = integrate.cumtrapz( f, t[t_ex::-1], initial=0 )
    
    return Q_in[t_ex::-1]/S[t_ex::-1]*np.exp(F)
