#!/usr/bin/env python
from __future__ import print_function
import numpy as np

def autostring(num, prec=0, zero=False, set_printoptions=False, pp=False, join=False, joinall=False, sep=' '):
    """
        Format number (array) with given decimal precision.


        Definition
        ----------
        def autostring(num, prec=0, zero=False, set_printoptions=False, pp=False, join=False, joinall=False, sep=' '):
          There is a wrapper function for convenience with the short name 'astr' that calls autostring
        def astr(num, prec=0, zero=False, set_printoptions=False, pp=False, join=False, joinall=False, sep=' '):


        Input
        -----
        num                 number array


        Optional Input
        --------------
        prec                number of decimal places of formatted values
                            minimum field width for integers (default: 0)
        zero                if True, pad values with zeros rather than blanks (default: False)
        set_printoptions    if True, sets linewidth to the format times size of 1st dimension (default: False)
        pp                  shortcut for set_printoptions (default: False)
                            it will be checked for (pp | set_printoptions)
        join                if True, joins all individual strings of last (fastest) dimension into one string (default: False)
        joinall             if True, joins all individual strings into single string,
                            i.e. first flattens the array and then joins it (default: False, overwrites join)
        sep                 separator used when joining (default: space=' ')


        Output
        ------
        string (array) of formatted numbers


        Restrictions
        ------------
        None


        Examples
        --------
        >>> print(autostring(3.5967, 3))
        3.597

        >>> print(autostring(3.5967))
        4

        >>> print(autostring(3, 3))
          3

        >>> print(autostring(np.array([3.5967, 3.5964]), 3))
        ['3.597'  '3.596']

        >>> print(autostring(np.array([3.59, 1.123456e12]), 3))
        ['3.590e+00'  '1.123e+12']

        >>> print(autostring(np.array([3.59, 11.1234]), 3, zero=True))
        ['03.590'  '11.123']

        >>> print(autostring(np.array([3, 11])))
        [' 3' '11']

        >>> print(autostring(np.array([3, 11]), 3))
        ['  3' ' 11']

        >>> print(autostring(np.zeros((2,2), dtype=np.float), 1))
        [['0.0' '0.0']
         ['0.0' '0.0']]

        >>> np.set_printoptions(threshold=10)
        >>> print(autostring(np.zeros((2,10), dtype=np.float), 1))
        [['0.0' '0.0' '0.0' ..., '0.0' '0.0' '0.0']
         ['0.0' '0.0' '0.0' ..., '0.0' '0.0' '0.0']]

        >>> print(autostring(np.zeros((2,10), dtype=np.float), 1, set_printoptions=True))
        [['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']
         ['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']]

        >>> print(autostring(np.zeros((2,10), dtype=np.float), 1, pp=True))
        [['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']
         ['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']]

        >>> print(autostring(np.zeros((2,10), dtype=np.float), 1, set_printoptions=False, pp=True))
        [['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']
         ['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']]

        >>> print(autostring(np.array([3.5967, 3.5964]), 3, join=True))
        3.597 3.596

        >>> print(autostring(np.zeros((2,10), dtype=np.float), 1, join=True, sep=';'))
        ['0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0'
         '0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0']

        >>> print(autostring(np.reshape(np.arange(20,dtype=np.float),(2,10)), 1, joinall=True, sep=';'))
         0.0; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; 9.0;10.0;11.0;12.0;13.0;14.0;15.0;16.0;17.0;18.0;19.0

        >>> print(autostring(np.reshape(np.arange(20,dtype=np.float),(2,10)), 1, joinall=True, sep=';'))
         0.0; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; 9.0;10.0;11.0;12.0;13.0;14.0;15.0;16.0;17.0;18.0;19.0

        >>> print(autostring(np.array([3, 11, np.inf])))
        ['  3' ' 11' 'inf']

        >>> print(autostring(np.array([3, 11, np.nan])))
        ['  3' ' 11' 'nan']

        >>> print(autostring(np.ma.array([3, 11, np.nan], mask=[False,True,False])))
        ['  3' '-- ' 'nan']

        >>> print(autostring(np.ma.array([3, 11, np.nan], mask=[False,False,True])))
        [' 3' '11' '--']


        License
        -------
        This file is part of the UFZ Python library.

        The UFZ Python library is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The UFZ Python library is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with The UFZ Python library.  If not, see <http://www.gnu.org/licenses/>.

        Copyright 2011-2013 Matthias Cuntz


        History
        -------
        Written,  MC, Nov 2011 - from autostring.pro
        Modified, MC, May 2012 - pp
                  MC, Dec 2012 - special treatment of -0.0 on output
                  MC, Feb 2013 - nan, inf and masked arrays
                  MC, Feb 2013 - ported to Python 3
    """
    #
    # Check input
    if type(num) == type([]): num = np.array(num)
    isarr = np.ndim(num)
    if (isarr > 2):
        print("AUTOSTRING WARNING: autostring only works with scalars, 1D- and 2D arrays: return original array.")
        return num
    # Only treat int and float
    if (isarr==0):
        try:
            typ = num.dtype
        except AttributeError:
            if (type(num) == float):
                typ = np.float64
            elif (type(num) == int):
                typ = np.int32
            else:
                typ = type(num)
    else:
        typ = num.dtype
    try:
        lfloat = np.float128    # Mac/*nix
    except AttributeError:
        try:
            lfloat = np.float96 # Windows
        except AttributeError:
            lfloat = np.float64
    if np.__version__ >= "1.6":
        if (typ in [np.float16, np.float32, np.float64, lfloat]):
            isfloat = True
        elif (typ in [np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64]):
            isfloat = False
        else:
            print("AUTOSTRING WARNING: autostring cannot work with input type: return original array.")
            return num
    else:
        if (typ in [np.float32, np.float64, lfloat]):
            isfloat = True
        elif (typ in [np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64]):
            isfloat = False
        else:
            print("AUTOSTRING WARNING: autostring cannot work with input type: return original array.")
            return num
    # Scalar to array if necessary; Special treatment of -0.0
    if (isarr==0):
        if (num == 0):
            num = np.abs(num)
    else:
        num = np.where(num == 0, 0, num)
    # Zero padding
    if zero:
        nix = '0'
    else:
        nix = ''
    #
    # If we deal with an array of numbers we take the largest for the format
    # deal with inf and nan
    hasmask = False
    hasnan  = False
    if (isarr==0):
        if np.isnan(num): return 'nan'
        if np.isinf(num): return 'inf'
        abs_num = np.ma.abs(num)
        # leave room for the decimal point and the negative sign, if any
        if (num < 0.):
            num_sign_chars = 1
        else:
            num_sign_chars = 0
    else:
        if type(num) == type(np.ma.ones(1)):
            if np.sum(num.mask) > 0: hasmask = True
            if num.count() > np.ma.sum(np.isfinite(num)): hasnan = True
        else:
            if num.size > np.sum(np.isfinite(num)): hasnan = True
        inum   = np.ma.array(num, mask=~np.isfinite(num), keep_mask=True)
        abs_num = np.ma.max(np.ma.abs(inum))
        # leave room for the decimal point and the negative sign, if any
        if (np.ma.min(inum) < 0.):
            num_sign_chars = 1
        else:
            num_sign_chars = 0
    #
    # Floating point
    if isfloat: # number is a float, more or less
        if abs_num >= 1.e6:
            num_prefix_chars  = 1
            num_sci_not_chars = 4
            format_type       = 'e'
        elif ((abs_num < 1.e6) & (abs_num >= 1.)):
            nprefix = np.int_(np.log10(np.int32(abs_num)))+1
            # special treatment: the output prefix digits could
            # be one digit longer as the input prefix digits: e.g. 99.99 => 100.0
            val               = np.around(abs_num*(10.**prec))/(10.**prec)
            nprefixval        = np.int_(np.log10(val))+1
            nprefix           = np.amax(np.array([nprefix,nprefixval], dtype=np.int))
            num_prefix_chars  = nprefix
            num_sci_not_chars = 0
            format_type       = 'f'
        elif ((abs_num < 1.) & (abs_num >= 1.e-3)):
            num_prefix_chars  = 1
            num_sci_not_chars = 0
            format_type       = 'f'
        elif (abs_num == 0):
            num_prefix_chars  = 1
            num_sci_not_chars = 0
            format_type       = 'f'
        else:
            num_prefix_chars  = 1
            num_sci_not_chars = 4
            format_type       = 'e'
        #
        num_postfix_chars = prec
        num_total_chars   = num_sign_chars + num_prefix_chars + 1 + num_postfix_chars + num_sci_not_chars
        if (prec == 0): # no dot if prec=0
            num_total_chars -= 1
        if hasmask: # need space for --
            if num_total_chars < 2: num_total_chars = 2
        if hasnan: # need space for nan or inf
            if num_total_chars < 3: num_total_chars = 3
        format_string     = ("{0:s}{1:s}{2:d}{3:s}{4:d}{5:s}{6:s}".format('{0:', nix, num_total_chars,
                                                                          '.', num_postfix_chars, format_type, '}'))
    else: # number is an integer
        format_type = 'd'
        if abs_num != 0:
            num_digits = np.int_(np.log10(abs_num))+1
        else:
            num_digits = 1
        num_total_chars = np.maximum(num_digits + num_sign_chars, prec)
        if hasmask: # need space for --
            if num_total_chars < 2: num_total_chars = 2
        if hasnan: # need space for nan or inf
            if num_total_chars < 3: num_total_chars = 3
        format_string = ("{0:s}{1:s}{2:d}{3:s}{4:s}".format('{0:', nix, num_total_chars, format_type, '}'))
    #
    if (isarr == 0):
        out = format_string.format(num)
        # Special treatment of -0.0
        if np.float(out) == 0:
            out = format_string.format(0)
    else:
        fnum = num.flatten()
        nnum = fnum.size
        import sys
        if sys.hexversion > int('0x3000000',base=16):
            styp = 'U{0:d}'.format(num_total_chars)
        else:
            styp = 'S{0:d}'.format(num_total_chars)
        out = np.empty(nnum, dtype=styp)
        for i in range(nnum):
            if str(fnum[i]) == '--':
                sformat_string = ("{0:s}{1:d}s{2:s}".format('{0:', num_total_chars, '}'))
                out[i] = sformat_string.format('--')
            else:
                out[i] = format_string.format(fnum[i])
                if np.float(out[i]) == 0:
                    out[i] = format_string.format(0)
        out = np.reshape(out, num.shape)
        if (set_printoptions | pp):
            # num_total_chars+3 for '' and space, +isarr for []
            np.set_printoptions(linewidth=num.shape[-1]*(num_total_chars+3)+isarr, threshold=nnum+1)
        if (join | joinall): # There should be reduction routines in numpy
            if ((isarr == 1) | ((isarr==2) & joinall)):
                if (isarr == 2):
                    out = out.flatten()
                for i in range(out.size):
                    if (i==0):
                        outc = out[i]
                    else:
                        outc = outc+sep+out[i]
            else:
                if sys.hexversion > int('0x3000000',base=16):
                    sform = 'U{0:d}'.format((len(out[0,0])+len(sep))*out.shape[1])
                else:
                    sform = 'S{0:d}'.format((len(out[0,0])+len(sep))*out.shape[1])
                outc = np.zeros(out.shape[0], dtype=sform)
                for j in range(out.shape[0]):
                    for i in range(out.shape[1]):
                        if (i==0):
                            outc[j] = out[j,i]
                        else:
                            outc[j] = outc[j]+sep+out[j,i]
            out = outc

    # return formatted string
    return out


def astr(num, prec=0, zero=False, set_printoptions=False, pp=True, join=False, joinall=False, sep=' '):
    """
        Wrapper function for autostring with pp=True by default.
        def autostring(num, prec=0, zero=False, set_printoptions=False, pp=False, join=False, joinall=False, sep=' '):


        Examples
        --------
        >>> print(astr(3.5967, 3))
        3.597

        >>> print(astr(3.5967))
        4

        >>> print(astr(3, 3))
          3

        >>> print(astr(np.array([3.5967, 3.5964]), 3))
        ['3.597' '3.596']

        >>> print(astr(np.array([3.59, 1.123456e12]), 3))
        ['3.590e+00' '1.123e+12']

        >>> print(astr(np.array([3.59, 11.1234]), 3, zero=True))
        ['03.590' '11.123']

        >>> print(astr(np.array([3, 11])))
        [' 3' '11']

        >>> print(astr(np.array([3, 11]), 3))
        ['  3' ' 11']

        >>> print(astr(np.zeros((2,2), dtype=np.float), 1))
        [['0.0' '0.0']
         ['0.0' '0.0']]

        >>> np.set_printoptions(threshold=10)
        >>> print(astr(np.zeros((2,10), dtype=np.float), 1))
        [['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']
         ['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']]

        >>> print(astr(np.zeros((2,10), dtype=np.float), 1, set_printoptions=True))
        [['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']
         ['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']]

        >>> print(astr(np.zeros((2,10), dtype=np.float), 1, pp=True))
        [['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']
         ['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']]

        >>> print(astr(np.zeros((2,10), dtype=np.float), 1, set_printoptions=False, pp=True))
        [['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']
         ['0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0' '0.0']]

        >>> print(astr(np.array([3.5967, 3.5964]), 3, join=True))
        3.597 3.596

        >>> print(astr(np.zeros((2,10), dtype=np.float), 1, join=True, sep=';'))
        ['0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0'
         '0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0']

        >>> print(astr(np.reshape(np.arange(20,dtype=np.float),(2,10)), 1, joinall=True, sep=';'))
         0.0; 1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; 9.0;10.0;11.0;12.0;13.0;14.0;15.0;16.0;17.0;18.0;19.0

        >>> print(astr(np.array([3, 11, np.inf])))
        ['  3' ' 11' 'inf']

        >>> print(astr(np.array([3, 11, np.nan])))
        ['  3' ' 11' 'nan']

        >>> print(astr(np.ma.array([3, 11, np.nan], mask=[False,True,False])))
        ['  3' '-- ' 'nan']

        >>> print(astr(np.ma.array([3, 11, np.nan], mask=[False,False,True])))
        [' 3' '11' '--']
    """
    return autostring(num, prec=prec, zero=zero, set_printoptions=set_printoptions,
                      pp=pp, join=join, joinall=joinall, sep=sep)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # print(autostring(np.array([3, 11, np.nan])))
    # #['  3' ' 11' 'nan']

    # print(autostring(np.ma.array([3, 11, np.nan], mask=[False,True,False])))
    # #['  3' '-- ' 'nan']

    # print(autostring(np.ma.array([3, 11, np.nan], mask=[False,False,True])))
    # #['  3' ' 11' '-- ']
