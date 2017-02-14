#!/usr/bin/env python
from __future__ import print_function
import numpy as np

def fwrite(fname, arr, header=None, precision='10.0'):
    """
        Write numbers of 2D-array to a file.

        A header can be given as optional as well as a precision.

        
        Definition
        ----------
        def fwrite(fname, data, header=None, precision='10.0'):

        
        Input
        -----
        fname        target file name
        arr          2d numpy array to write


        Optional Input Parameters
        -------------------------
        header       list of header elements: a header element is a list of two strings:
                     the first entry is the header argument, the second the value
        precision    floating point precision of array to write

        
        Examples
        --------
        >>> from ufz import fread
        >>> # Clean up doctest
        >>> filename = 'fwrite.test'
        >>> header = [['Description', 'testing'], ['author', 'ST']]
        >>> data = np.arange(10).reshape(2, 5)
        >>> fwrite(filename, data, header=header)

        >>> fread(filename, nc=2, skip=2, header=True)
        [['Description', 'testing'], ['author', 'ST']]
        >>> fread(filename, skip=2)
        array([[ 0.,  1.,  2.,  3.,  4.],
               [ 5.,  6.,  7.,  8.,  9.]])

        >>> import os
        >>> os.remove(filename)


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

        Copyright 2009-2015 Stephan Thober


        History
        -------
        Written,  ST, Feb 2016
        Modified, 
    """
    if not type(arr) is np.ndarray:
        raise ValueError('function fwrite: argument arr must be numpy.ndarray')
    fo = open(fname, 'w')
    if not header is None:
        # write header
        for ll in np.arange(len(header)):
            write_str = str(header[ll][0]) + ' ' + str(header[ll][1]) + '\n'
            fo.write(write_str)
    # write arr
    for ll in np.arange(arr.shape[0]):
        # format is mRM compatible
        write_str = ' '.join(['{:' + precision + 'f}'] * arr.shape[1]).format(*arr[ll, :]) + '\n'
        fo.write(write_str)
    fo.close()    


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
