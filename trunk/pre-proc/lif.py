#!/usr/bin/env python
from __future__ import print_function
import string

def lif(file, noblank=False, comment='', skip=0, maxcol=False):
    """
        Counts the numer of lines in a file.
        Blank (only whitespace) and comment lines can be excluded.


        Definition
        ----------
        def lif(file, noblank=False, comment='', skip=0):


        Input
        -----
        file         source file name


        Optional input parameters
        -------------------------
        comment      line gets excluded if first character of line is
                      in comment sequence
                     sequence can be e.g. string, list or tuple
        skip         number of lines to skip at the beginning of
                       the file (default 0)


        Options
        -------
        noblank      excludes all lines that consists only of
                       whitespace characters
        maxcol       if True, return also maximum amount of characters in one line


        Output
        ------
        if maxcol:
            number of lines in file, maximum characters in a line
        else:
            number of lines in file


        Restrictions
        ------------
        Only ascii files.


        Examples
        --------
        # Create some data
        >>> filename = 'test.dat'
        >>> file = open(filename,'w')
        >>> file.writelines('First line\\n')
        >>> file.writelines(' \\n')
        >>> file.writelines('# First comment\\n')
        >>> file.writelines('! Second comment\\n')
        >>> file.writelines('Last line\\n')
        >>> file.close()

        # Count lines
        >>> print(lif(filename))
        5
        >>> print(lif(filename,noblank=True))
        4
        >>> print(lif(filename,comment='#'))
        4
        >>> print(lif(filename,comment='#!'))
        3
        >>> print(lif(filename,comment='#S'))
        4
        >>> print(lif(filename,comment=('#','L')))
        3
        >>> print(lif(filename,comment=['#','!']))
        3
        >>> print(lif(filename,comment='#!',noblank=True))
        2
        >>> print(lif(filename,skip=2))
        3
        >>> print(lif(filename,skip=2,maxcol=True))
        (3, 16)

        # Clean up
        >>> import os
        >>> os.remove(filename)


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

        Copyright 2009-2013 Matthias Cuntz


        History
        -------
        Written,  MC, Jul 2009
        Modified, MC, Nov 2012 - maxcol
                  MC, Feb 2013 - ported to Python 3
    """
    # Open file
    try:
        f = open(file, 'r')
    except IOError:
        raise ValueError('Cannot open file '+file)
    # Count lines
    count = 0
    if skip > 0:
        iskip = 0
        while iskip < skip:
            l = f.readline()
            iskip += 1
    imax = []
    if noblank and (comment != ''):        # exclude blank, exclude comment
        for line in f:
            ll    = line
            imax += [len(ll)]
            l     = ll.strip(string.whitespace)
            if (l != ''):
                if (l[0] not in comment): count += 1
    elif noblank and (comment == ''):      # exclude blank, include comment
        for line in f:
            ll    = line
            imax += [len(ll)]
            l     = ll.strip(string.whitespace)
            if (l != ''): count += 1
    elif (not noblank) and (comment != ''):# include blank, exclude comment
        for line in f:
            ll    = line
            imax += [len(ll)]
            l     = ll.strip(string.whitespace)
            if (l == ''):
                count += 1
            else:
                if (l[0] not in comment): count += 1
    else:                                  # include blank, include comment
        for line in f:
            imax += [len(line)]
        count = len(imax)
    f.close()

    if maxcol:
        return count, max(imax)-1 # lines include \0 at the end.
    else:
        return count


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
