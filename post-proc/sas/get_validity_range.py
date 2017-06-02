#!/usr/bin/env python

import numpy as np
#import matplotlib.pyplot as plt

def get_validity_range(t, f):
    
    """
        Computes the range of the validity of the approach of Botter et al.
        'Transport in the hydrological response: Travel time distributions, 
        soil moisture dynamics, and the old water paradox' WRR (2010), 
        Equation (35)


        Definition
        ----------
        def get_validity_range(t, f):


        Input
        -----

        t           chronological time of reservoir dynamcis
        f           function whose value should be equal to unity


        Output
        ------
        range within f is actually close enough to unity, i.e. within the 
        approach of Botter et al. is valid


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
    
    threshold = 0.8
    f = f - threshold
    f = np.sign(f)
    
    range_end = np.where(f==(-1.0))
    range_end = range_end[0]
    
#    print(range_end[0])
    
    return range_end[0]
