#!/usr/bin/env python

import numpy as np

from scipy import integrate
##from get_indef_integral import get_indef_integral 

def get_theta(U, Q_out_1, Q_out_2, t, t_in):
    
    """
        Computes the partitioning function theta according to Botter et al.
        'Transport in the hydrological response: Travel time distributions, 
        soil moisture dynamics, and the old water paradox' WRR (2010), 
        Equation (B2)


        Definition
        ----------
        def get_theta(U, Q_out_1, Q_out_2, t, t_in):


        Input
        -----
        U           volumetric water content of the reservoir
        Q_out_1     first output flux of reservoir (contributing to discharge)
        Q_out_2     second output flux of reservoir
        t           chronological time of reservoir dynamcis
        t_in        time when water is entering reservoir


        Output
        ------
        Float array of time-dependent partitioning function theta, representing 
        the percentage of Q_in that is entering the reservoir at t_in that is 
        leaving the reservoir as discharge eventually


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
    
#    t_num = len(t)
    tau = t[t_in:]
    
    f_1 = (Q_out_1[t_in:] + Q_out_2[t_in:])/U[t_in:]
    f_2 = Q_out_1[t_in:]/U[t_in:]
    
    aux_fun = np.exp( - integrate.cumtrapz( f_1, tau, initial=0 ))
    theta = np.trapz( f_2*aux_fun, tau )

    return theta
