#!/usr/bin/env python
from __future__ import print_function
import numpy as np

def readnetcdf(file, var='', code=-1, reform=False, squeeze=False,
               variables=False, codes=False, dims=False, units=False, longnames=False,
               attributes=False, sort=False, pointer=False, overwrite=False):
    """
        Gets variables or prints information of netcdf file.


        Definition
        ----------
        def readnetcdf(file, var='', code=-1, reform=False, squeeze=False,
                       variables=False, codes=False, dims=False, units=False, longnames=False, 
                       attributes=False, sort=False, pointer=False, overwrite=False):


        Input
        -----
        file         netcdf file name


        Optional Input Parameters
        -------------------------
        var          name of variable in netcdf file
        code         code number in attribute code


        Options
        -------
        reform       if output is array then squeeze(array)
                     if codes then remove all codes==-1
                     if units or longnames then remove =''
        squeeze      same as reform
        variables    get list of variables in netcdf file
        codes        get list of codes attribute code
        dims         get list of dimensions a variable is depending on
        units        get list of units of variables from attribute units
        longnames    get list of long names of variables from
                     attribute long_name
        attributes   get dictionary of all attributes of specific variable or of file if variable is omitted
        sort         sort variable names. Codes, units and longnames will be
                     sorted accoringly so that indeces still match.
        pointer      if True, (return file pointer, variable pointer); only for reading
        overwrite    if True, (return file pointer, variable pointer); modification of file/variable possible if
                     file contains only one variable

        Output
        ------
        Either float array of variable/code or information lists
        such as list of all variables in netcdf file.


        Restrictions
        ------------
        If codes, units or longnames are reformed/squeezed,
          they do not match the variable list anymore.
        Attributes can not be sorted nor reformed/squeezed.


        Examples
        --------
        # Read varibale or code
        >>> print(readnetcdf('test_readnetcdf.nc',var='is1'))
        [[ 1.  1.  1.  1.]
         [ 1.  1.  1.  1.]]
        >>> print(readnetcdf('test_readnetcdf.nc',code=129))
        [[ 2.  2.  2.  2.]
         [ 2.  2.  2.  2.]]

        # Get variable names
        >>> print([str(i) for i in readnetcdf('test_readnetcdf.nc',variables=True)])
        ['x', 'y', 'is1', 'is2']
        >>> print([str(i) for i in readnetcdf('test_readnetcdf.nc',variables=True,sort=True)])
        ['is1', 'is2', 'x', 'y']

        # Get codes
        >>> print(readnetcdf('test_readnetcdf.nc',codes=True))
        [  -1.   -1.  128.  129.]
        >>> print(readnetcdf('test_readnetcdf.nc',codes=True,reform=True))
        [ 128.  129.]
        >>> print(readnetcdf('test_readnetcdf.nc',codes=True,sort=True))
        [128.0, 129.0, -1.0, -1.0]

        # Get special attributes units and longnames
        >>> print([str(i) for i in readnetcdf('test_readnetcdf.nc',units=True)])
        ['xx', 'yy', 'arbitrary', 'arbitrary']
        >>> print([str(i) for i in readnetcdf('test_readnetcdf.nc',units=True,sort=True)])
        ['arbitrary', 'arbitrary', 'xx', 'yy']
        >>> print([str(i) for i in readnetcdf('test_readnetcdf.nc',longnames=True)])
        ['x-axis', 'y-axis', 'all ones', 'all twos']
        >>> print([str(i) for i in readnetcdf('test_readnetcdf.nc',longnames=True,sort=True)])
        ['all ones', 'all twos', 'x-axis', 'y-axis']

        # Get dims (change from unicode to string)
        >>> print([ str(i) for i in readnetcdf('test_readnetcdf.nc',var='is1', dims=True) ])
        ['y', 'x']

        # Get attributes
        # old: {'units': 'arbitrary', 'long_name': 'all ones', 'code': 128}
        # new: {u'units': u'arbitrary', u'long_name': u'all ones', u'code': 128}
        >>> t1 = readnetcdf('test_readnetcdf.nc',var='is1',attributes=True)
        >>> print([ str(i) for i in sorted(t1)])
        ['code', 'long_name', 'units']

        # Just get file handle so that read is done later at indexing
        # useful for example to inquire remote netcdf files first
        >>> fh, var = readnetcdf('test_readnetcdf.nc',var='is1', pointer=True)
        >>> print( var.shape )
        (2, 4)
        >>> print( var[:] )
        [[ 1.  1.  1.  1.]
         [ 1.  1.  1.  1.]]
        >>> fh.close()
        
        # Change a variable in a file
        >>> print(readnetcdf('test_readnetcdf.nc',var='is1'))
        [[ 1.  1.  1.  1.]
         [ 1.  1.  1.  1.]]
        >>> fh, var = readnetcdf('test_readnetcdf.nc',var='is1', overwrite=True)
        >>> var[:] *= 2.
        >>> fh.close()
        >>> print(readnetcdf('test_readnetcdf.nc',var='is1'))
        [[ 2.  2.  2.  2.]
         [ 2.  2.  2.  2.]]
        >>> fh, var = readnetcdf('test_readnetcdf.nc',var='is1', overwrite=True)
        >>> var[:] *= 0.5
        >>> fh.close()


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
        along with the UFZ Python package (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2009-2014 Matthias Cuntz, Stephan Thober


        History
        -------
        Written,  MC, Jul 2009
        Modified, MC, Jun 2012 - removed quiet
                  MC, Feb 2013 - ported to Python 3
                  MC, Oct 2013 - netcdfread, ncread, readnc
                  ST, Apr 2014 - added overwrite flag
                  ST, May 2014 - added dims flag
                  ST, Jun 2016 - added read of file attributes
                  ST, Aug 2016 - restricted use of overwrite option
    """
    try:
        import netCDF4 as nc
    except:
        raise IOError('No NetCDF4 support available.')
    # Open netcdf file
    try:
        if overwrite:
            f = nc.Dataset(file, 'a')
            if len(f.variables) > 1:
                f.close()
                raise ValueError('ERROR: only use the overwrite option for files with one variable.')
        else:
            f = nc.Dataset(file, 'r')
    except IOError:
        raise IOError('Cannot open file for reading.'+file)
    # Variables
    vars = list(f.variables.keys())
    nvars = len(vars)
    # Sort and get sort indices
    if sort:
        svars = sorted(vars)
        ivars = list()
        for v in svars:
            ivars.append(vars.index(v))
    if variables:
        f.close()
        if sort:
            return svars
        else:
            return vars
    # Codes
    cods = np.empty(nvars)
    i=0
    for v in vars:
        attr = f.variables[v].ncattrs()
        if 'code' in attr:
            cods[i] = getattr(f.variables[v],'code')
        else:
            cods[i] = -1
        i += 1
    if codes:
        if sort:
            scods = [cods[i] for i in ivars]
        else:
            scods = cods
        if reform or squeeze:
            scods = np.compress(scods!=-1, scods)
        f.close()
        return scods
    # Get dimensions
    if dims:
        if var not in vars:
            f.close()
            raise ValueError('Variable '+var+' not in file '+file)
        dimensions = f.variables[ var ].dimensions
        f.close()
        return dimensions        
    # Get units
    if units:
        unis = list()
        for v in vars:
            attr = f.variables[v].ncattrs()
            if 'units' in attr:
                unis.append(getattr(f.variables[v],'units'))
            else:
                unis.append('')
        if sort:
            sunis = [unis[i] for i in ivars]
        else:
            sunis = unis
        if reform or squeeze:
            nn = sunis.count('')
            if nn > 0:
                for i in range(nn):
                    sunis.remove('')
        f.close()
        return sunis
    # Get long_name
    if longnames:
        longs = list()
        for v in vars:
            attr = f.variables[v].ncattrs()
            if 'long_name' in attr:
                longs.append(getattr(f.variables[v],'long_name'))
            else:
                longs.append('')
        if sort:
            slongs = [longs[i] for i in ivars]
        else:
            slongs = longs
        if reform or squeeze:
            nn = slongs.count('')
            if nn > 0:
                for i in range(nn):
                    slongs.remove('')
        f.close()
        return slongs
    # Get attributes
    if attributes:
        if var == '':
            attrs = dict()
            attr = f.ncattrs()
            for a in attr:
                attrs[a] = getattr(f, a)
            f.close()
            return attrs
        elif var not in vars:
            f.close()
            raise ValueError('Variable '+var+' not in file '+file)
        attrs = dict()
        attr = f.variables[var].ncattrs()
        for a in attr:
            attrs[a] = getattr(f.variables[var],a)
        f.close()
        return attrs
    # Get variable
    if var == '' and code==-1:
        f.close()
        raise ValueError('Variable name or code has to be given')
    if var != '':
        if var not in vars:
            f.close()
            raise ValueError('Variable '+var+' not in file '+file)
        if pointer or overwrite:
            arr = f.variables[var]
            return f, arr
        else:
            arr = f.variables[var][:]
            if reform or squeeze:
                f.close()
                return arr.squeeze()
            else:
                f.close()
                return arr
    if code != -1:
        if code not in cods:
            f.close()
            raise ValueError('Code '+str(code)+' not in file '+file)
        arr = f.variables[np.compress(cods==code, vars)[0]][:]
        if reform or squeeze:
            f.close()
            return arr.squeeze()
        else:
            f.close()
            return arr


def netcdfread(*args, **kwargs):
    """
        Wrapper for readnetcdf
        def readnetcdf(file, var='', code=-1, reform=False, squeeze=False,
                       variables=False, codes=False, dims=False, units=False, longnames=False, 
                       attributes=False, sort=False, pointer=False, overwrite=False):


        Examples
        --------
        >>> print(netcdfread('test_readnetcdf.nc',var='is1'))
        [[ 1.  1.  1.  1.]
         [ 1.  1.  1.  1.]]
        >>> print(netcdfread('test_readnetcdf.nc',code=129))
        [[ 2.  2.  2.  2.]
         [ 2.  2.  2.  2.]]
        >>> print([str(i) for i in netcdfread('test_readnetcdf.nc',variables=True)])
        ['x', 'y', 'is1', 'is2']
        >>> print([str(i) for i in netcdfread('test_readnetcdf.nc',variables=True,sort=True)])
        ['is1', 'is2', 'x', 'y']
        >>> print([str(i) for i in netcdfread('test_readnetcdf.nc',units=True)])
        ['xx', 'yy', 'arbitrary', 'arbitrary']
        >>> print([str(i) for i in netcdfread('test_readnetcdf.nc',units=True,sort=True)])
        ['arbitrary', 'arbitrary', 'xx', 'yy']
        >>> print([str(i) for i in netcdfread('test_readnetcdf.nc',longnames=True)])
        ['x-axis', 'y-axis', 'all ones', 'all twos']
        >>> print([str(i) for i in netcdfread('test_readnetcdf.nc',longnames=True,sort=True)])
        ['all ones', 'all twos', 'x-axis', 'y-axis']

        # old: {'units': 'arbitrary', 'long_name': 'all ones', 'code': 128}
        # new: {u'units': u'arbitrary', u'long_name': u'all ones', u'code': 128}
        >>> t1 = netcdfread('test_readnetcdf.nc',var='is1',attributes=True)
        >>> print([ str(i) for i in sorted(t1)])
        ['code', 'long_name', 'units']
        >>> print(netcdfread('test_readnetcdf.nc',codes=True))
        [  -1.   -1.  128.  129.]
        >>> print(netcdfread('test_readnetcdf.nc',codes=True,reform=True))
        [ 128.  129.]
        >>> print(netcdfread('test_readnetcdf.nc',codes=True,sort=True))
        [128.0, 129.0, -1.0, -1.0]
    """
    return readnetcdf(*args, **kwargs)


def ncread(*args, **kwargs):
    """
        Wrapper for readnetcdf
        def readnetcdf(file, var='', code=-1, reform=False, squeeze=False,
                       variables=False, codes=False, dims=False, units=False, longnames=False, 
                       attributes=False, sort=False, pointer=False, overwrite=False):


        Examples
        --------
        >>> print(ncread('test_readnetcdf.nc',var='is1'))
        [[ 1.  1.  1.  1.]
         [ 1.  1.  1.  1.]]
        >>> print(ncread('test_readnetcdf.nc',code=129))
        [[ 2.  2.  2.  2.]
         [ 2.  2.  2.  2.]]
        >>> print([str(i) for i in ncread('test_readnetcdf.nc',variables=True)])
        ['x', 'y', 'is1', 'is2']
        >>> print([str(i) for i in ncread('test_readnetcdf.nc',variables=True,sort=True)])
        ['is1', 'is2', 'x', 'y']
        >>> print([str(i) for i in ncread('test_readnetcdf.nc',units=True)])
        ['xx', 'yy', 'arbitrary', 'arbitrary']
        >>> print([str(i) for i in ncread('test_readnetcdf.nc',units=True,sort=True)])
        ['arbitrary', 'arbitrary', 'xx', 'yy']
        >>> print([str(i) for i in ncread('test_readnetcdf.nc',longnames=True)])
        ['x-axis', 'y-axis', 'all ones', 'all twos']
        >>> print([str(i) for i in ncread('test_readnetcdf.nc',longnames=True,sort=True)])
        ['all ones', 'all twos', 'x-axis', 'y-axis']

        # old: {'units': 'arbitrary', 'long_name': 'all ones', 'code': 128}
        # new: {u'units': u'arbitrary', u'long_name': u'all ones', u'code': 128}
        >>> t1 = ncread('test_readnetcdf.nc',var='is1',attributes=True)
        >>> print([ str(i) for i in sorted(t1)])
        ['code', 'long_name', 'units']
        >>> print(ncread('test_readnetcdf.nc',codes=True))
        [  -1.   -1.  128.  129.]
        >>> print(ncread('test_readnetcdf.nc',codes=True,reform=True))
        [ 128.  129.]
        >>> print(ncread('test_readnetcdf.nc',codes=True,sort=True))
        [128.0, 129.0, -1.0, -1.0]
    """
    return readnetcdf(*args, **kwargs)


def readnc(*args, **kwargs):
    """
        Wrapper for readnetcdf
        def readnetcdf(file, var='', code=-1, reform=False, squeeze=False,
                       variables=False, codes=False, dims=False, units=False, longnames=False, 
                       attributes=False, sort=False, pointer=False, overwrite=False):


        Examples
        --------
        >>> print(readnc('test_readnetcdf.nc',var='is1'))
        [[ 1.  1.  1.  1.]
         [ 1.  1.  1.  1.]]
        >>> print(readnc('test_readnetcdf.nc',code=129))
        [[ 2.  2.  2.  2.]
         [ 2.  2.  2.  2.]]
        >>> print([str(i) for i in readnc('test_readnetcdf.nc',variables=True)])
        ['x', 'y', 'is1', 'is2']
        >>> print([str(i) for i in readnc('test_readnetcdf.nc',variables=True,sort=True)])
        ['is1', 'is2', 'x', 'y']
        >>> print([str(i) for i in readnc('test_readnetcdf.nc',units=True)])
        ['xx', 'yy', 'arbitrary', 'arbitrary']
        >>> print([str(i) for i in readnc('test_readnetcdf.nc',units=True,sort=True)])
        ['arbitrary', 'arbitrary', 'xx', 'yy']
        >>> print([str(i) for i in readnc('test_readnetcdf.nc',longnames=True)])
        ['x-axis', 'y-axis', 'all ones', 'all twos']
        >>> print([str(i) for i in readnc('test_readnetcdf.nc',longnames=True,sort=True)])
        ['all ones', 'all twos', 'x-axis', 'y-axis']

        # old: {'units': 'arbitrary', 'long_name': 'all ones', 'code': 128}
        # new: {u'units': u'arbitrary', u'long_name': u'all ones', u'code': 128}
        >>> t1 = readnc('test_readnetcdf.nc',var='is1',attributes=True)
        >>> print([ str(i) for i in sorted(t1)])
        ['code', 'long_name', 'units']
        >>> print(readnc('test_readnetcdf.nc',codes=True))
        [  -1.   -1.  128.  129.]
        >>> print(readnc('test_readnetcdf.nc',codes=True,reform=True))
        [ 128.  129.]
        >>> print(readnc('test_readnetcdf.nc',codes=True,sort=True))
        [128.0, 129.0, -1.0, -1.0]
    """
    return readnetcdf(*args, **kwargs)


if __name__ == '__main__':
    import doctest
    try:
        import netCDF4 as nc
    except ImportError:
        raise ImportError('No NetCDF4 support available.')
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
