#!/usr/bin/env python
"""
    Python Utilities for computing StorAge-Selection (SAS) functions

    Get help on each function by typing
    >>> import sas
    >>> help(sas.function)


    License
    -------

    This file is part of the UFZ Python package.

    Not all files in the package are free software. The license is given in the
    'License' section of the docstring of each routine.

    The package is released under the GNU Lesser General Public License. The
    following applies: The SAS Python package is free software: you can 
    redistribute it and/or modify it under the terms of the GNU Lesser 
    General Public License as published by the Free Software Foundation, 
    either version 3 of the License, or (at your option) any later version.

    The SAS Python package is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
    If not, see <http://www.gnu.org/licenses/>.

    Copyright 2015 Falk He"sse


    History
    -------
    Written,  FH,   Apr 2015

"""
#from __future__ import print_function

# Routines

from sas_base import SAS
from aux_fun import *
from get_p import *
from get_theta import *
from get_U_num import *
from get_validity_range import *
#from plot_fun import *

# Information
__author__   = 'Falk Hesse'
__version__  = '0.1.0'
#__revision__ = 
__date__     = 'Date: 01.04.2015'

# Main
#if __name__ == '__main__':
#    print('\nMAD Python Package.')
#    print("Version {:s} from {:s}.".format(__version__,__date__))
#    print('\nThis is the README file. See als the license file LICENSE.\n\n')
#    f = open('README','r')
#    for line in f: print(line,end='')
#    f.close()
