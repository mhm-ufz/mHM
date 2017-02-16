#!/usr/bin/env python
"""
    Python Utilities for mHM

    Get help on each function by typing
    >>> import mhm
    >>> help(mhm.function)


    License
    -------

    This file is part of the UFZ Python package.

    The package is released under the GNU Lesser General Public License. The
    following applies: The MHM Python package is free software: you can 
    redistribute it and/or modify it under the terms of the GNU Lesser 
    General Public License as published by the Free Software Foundation, 
    either version 3 of the License, or (at your option) any later version.

    The MHM Python package is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with the UFZ makefile project (cf. gpl.txt and lgpl.txt).
    If not, see <http://www.gnu.org/licenses/>.

    Copyright 2015 Falk He"sse


    History
    -------
    Written,  FH,   Aug 2015

"""
#from __future__ import print_function

# Routines

from mhm import MHM
from readnetcdf import readnetcdf
#from mhm_cell import MHM_Cell

# Information
__author__   = 'Falk Hesse'
__version__  = '0.1.0'
#__revision__ = 
__date__     = 'Date: 10.08.2015'

# Main
#if __name__ == '__main__':
#    print('\nMAD Python Package.')
#    print("Version {:s} from {:s}.".format(__version__,__date__))
#    print('\nThis is the README file. See als the license file LICENSE.\n\n')
#    f = open('README','r')
#    for line in f: print(line,end='')
#    f.close()
