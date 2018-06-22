#!/usr/bin/env python
from __future__ import print_function
import numpy as np

def date2dec(calendar = 'standard', units=None,
             excelerr = True, yr=1,
             mo=1, dy=1, hr=0, mi=0, sc=0,
             ascii=None, eng=None):
    """
        Converts numpy arrays with calendar date into
        numpy arrays with decimal date. Supported calendar
        formats are standard, gregorian, julian, proleptic_gregorian,
        excel1900, excel1904, 365_day, noleap, 366_day, all_leap,
        360_day, decimal, or decimal360

        Input is year, month day, hour, minute,
        second or a combination of them. ASCII input
        is possible, too.

        Output is decimal date with day as unit.

        Requires 'netcdftime.py' from the module
        netcdftime available at:
        http://netcdf4-python.googlecode.com


        Definition
        ----------
        def date2dec(calendar = 'standard', units=None,
                     excelerr = True, yr=1,
                     mo=1, dy=1, hr=0, mi=0, sc=0,
                     ascii=None, eng=None):


        Input
        -----
        yr       -> input array with year
        mo       -> input array with month
        dy       -> input array with day
        hr       -> input array with hour
        mi       -> input array with minute
        sc       -> input array with second
        ascii    -> input array with strings of the format
                    'dd.mm.yyyy hh:mm:ss'. If hour, minutes
                    and/or seconds are missing, they will be
                    set to 00.
                    If ascii input is chosen by user,
                    other inputs will be neglected.
                    ascii and eng are mutually exclusive.
        eng      -> input array with strings of the format
                    'yyyy-mm-dd hh:mm:ss'. If hour, minutes
                    and/or seconds are missing, they will be
                    set to 00.
                    If eng input is chosen, other inputs will
                    be neglected.
                    ascii and eng are mutually exclusive.


        Parameters
        ----------
        calendar -> Input date format. Default value is
                   'standard'.

           'standard', 'gregorian'
                       =   Input date is standard format.
                           Input is in julian calendar from
                           01.01.-4712 12:00:00 (BC) until
                           05.03.1583 00:00:00 and gregorian
                           calendar from 15.03.1583 00:00:00
                           until now. Missing 10 days don't
                           exsist.
           'julian'    =   Input date is julian format.
                           Input is in julian calendar from
                           01.01.-4712 12:00:00 (BC) until now.
           'proleptic_gregorian'
                       =   Input date is gregorian format.
                           Input is in gregorian calendar from
                           01.01.0000 00:00:00 until now.
           'excel1900' =   Input date is excel 1900 format.
                           Input date is excel date with its
                           units at 01.01.1900 00:00:00 until
                           now.
           'excel1904' =   Input date is excel 1904 (lotus) format.
                           Input date is excel date with its
                           units at 01.01.1904 00:00:00 until now.
           '365_day', 'noleap'
                       =   Input date is 365 days format. Input date
                           consists of common years only (No leap years)
                           with its units at 01.01.0001 00:00:00 until now.
           '366_day', 'all_leap'
                       =   Input date is 366 days format. Input date
                           consists of leap years only (No common years)
                           with its units at 01.01.0001 00:00:00 until now.
           '360_day'   =   Input date is 360 days format.  Input
                           date consists of years with only 360 days
                           (30 days per month)with its units at
           'decimal'    =  Output is decimal year.
           'decimal360' =  Output is decimal year with a year of 360 days, i.e. 12 month with 30 days each.


        Optional Arguments
        ------------------
        units     -> Time units can be set by user. Input must be a
                     string in the format 'days since yyyy-mm-dd hh:mm:ss'.
                     Default values are set automatically.
        excelerr  -> In Excel the year 1900 is normally considered
                     as leap year, which is wrong. By default, this
                     error is taken into account (excelerr = True).
                     For excelerr = False, 1900 is considered as no
                     leap year.


        Output
        ------
        output -> Output numpy array with decimal date.


        Restrictions
        ------------
        Most versions of datetime do not support neagtive years,
        i.e. Julian days < 1721423.5 = 01.01.0001 00:00.

        There is an issue in netcdftime version < 0.9.5 in proleptic_gregorian for dates before year 301:
          ufz.dec2date(ufz.date2dec(ascii='01.01.0300 00:00:00', calendar='proleptic_gregorian'),
                       calendar='proleptic_gregorian')
            [300, 1, 2, 0, 0, 0]
          ufz.dec2date(ufz.date2dec(ascii='01.01.0301 00:00:00', calendar='proleptic_gregorian'),
                       calendar='proleptic_gregorian')
            [301, 1, 1, 0, 0, 0]

        List input is only supported up to 2 dimensions.

        Requires 'netcdftime.py' from module netcdftime available at:
        http://netcdf4-python.googlecode.com


        Examples
        --------
        #calendar = 'standard'

        # Some implementations of datetime have problems with negative years
        >>> import datetime
        >>> if datetime.MINYEAR > 0:
        ...     print('The minimum year in your datetime implementation is ', datetime.MINYEAR)
        ...     print('i.e. it does not support negative years (BC).')

        >>> if datetime.MINYEAR > 0:
        ...     year   = np.array([2000, 1810, 1630, 1510, 1271, 619, 2, 1])
        ... else:
        ...     year   = np.array([2000, 1810, 1630, 1510, 1271, 619, -1579, -4712])
        >>> month  = np.array([1, 4, 7, 9, 3, 8, 8, 1])
        >>> day    = np.array([5, 24, 15, 20, 18, 27, 23, 1])
        >>> hour   = np.array([12, 16, 10, 14, 19, 11, 20, 12])
        >>> minute = np.array([30, 15, 20, 35, 41, 8, 3, 0])
        >>> second = np.array([15, 10, 40, 50, 34, 37, 41, 0])
        >>> decimal = date2dec(calendar = 'standard', yr=year, mo=month, dy=day, hr=hour, mi=minute, sc=second)
        >>> from autostring import astr
        >>> nn = year.size
        >>> print(astr(decimal[:nn/2], 14, pp=True))
        ['2.45154902100694e+06' '2.38226217719907e+06' '2.31660093101852e+06' '2.27284810821759e+06']
        >>> print(astr(decimal[nn/2:nn-2], 14,pp=True))
        ['2.18536732053241e+06' '1.94738596431713e+06']
        >>> decimal = date2dec(calendar='standard', yr=year, mo=6, dy=15, hr=12, mi=minute, sc=second)
        >>> print(astr(decimal[:nn/2],14,pp=True))
        ['2.45171102100694e+06' '2.38231401053241e+06' '2.31657101435185e+06' '2.27275102488426e+06']
        >>> print(astr(decimal[nn/2:nn-2],14,pp=True))
        ['2.18545602886574e+06' '1.94731300598380e+06']

        # ascii input
        >>> if datetime.MINYEAR > 0:
        ...     a = np.array(['05.01.2000 12:30:15', '24.04.1810 16:15:10', '15.07.1630 10:20:40',  '20.09.1510 14:35:50',
        ...                   '18.03.1271 19:41:34', '27.08. 619 11:08:37', '23.08.0002 20:03:41', '01.01.0001 12:00:00'])
        ... else:
        ...     a = np.array(['05.01.2000 12:30:15', '24.04.1810 16:15:10', '15.07.1630 10:20:40',  '20.09.1510 14:35:50',
        ...                   '18.03.1271 19:41:34', '27.08. 619 11:08:37', '23.08.-1579 20:03:41', '01.01.-4712 12:00:00'])
        >>> decimal = date2dec(calendar='standard', ascii=a)
        >>> nn = a.size
        >>> print(astr(decimal[:nn/2],14,pp=True))
        ['2.45154902100694e+06' '2.38226217719907e+06' '2.31660093101852e+06' '2.27284810821759e+06']
        >>> print(astr(decimal[nn/2:nn-2],14,pp=True))
        ['2.18536732053241e+06' '1.94738596431713e+06']

        # calendar = 'julian'
        >>> decimal = date2dec(calendar='julian', ascii=a)
        >>> print(astr(decimal[:nn/2],14,pp=True))
        ['2.45156202100694e+06' '2.38227417719907e+06' '2.31661093101852e+06' '2.27284810821759e+06']
        >>> print(astr(decimal[nn/2:nn-2],14,pp=True))
        ['2.18536732053241e+06' '1.94738596431713e+06']

        # calendar = 'proleptic_gregorian'
        >>> decimal = date2dec(calendar='proleptic_gregorian', ascii=a)
        >>> print(astr(decimal[:nn/2], 7, pp=True))
        ['730123.5210069' '660836.6771991' '595175.4310185' '551412.6082176']
        >>> print(astr(decimal[nn/2:nn-2], 7, pp=True))
        ['463934.8205324' '225957.4643171']

        # calendar = 'excel1900' WITH excelerr=True -> 1900 considered as leap year
        >>> d = np.array(['05.01.2000 12:30:15', '27.05.1950 16:25:10', '13.08.1910 10:40:55',
        ...               '01.03.1900 00:00:00', '29.02.1900 00:00:00', '28.02.1900 00:00:00',
        ...               '01.01.1900 00:00:00'])
        >>> decimal = date2dec(calendar='excel1900', ascii=d)
        >>> nn = d.size
        >>> print(astr(decimal[:nn/2], 7, pp=True))
        ['36530.5210069' '18410.6841435' ' 3878.4450810']
        >>> print(astr(decimal[nn/2:],14,pp=True))
        ['61.00000000000000' '60.00000000000000' '59.00000000000000' ' 1.00000000000000']

        # calendar = 'excel1900' WITH excelerr = False -> 1900 is NO leap year
        >>> decimal = date2dec(calendar='excel1900', ascii=d, excelerr=False)
        >>> print(astr(decimal[:nn/2], 7, pp=True))
        ['36529.5210069' '18409.6841435' ' 3877.4450810']
        >>> print(astr(decimal[nn/2:],14,pp=True))
        ['60.00000000000000' '60.00000000000000' '59.00000000000000' ' 1.00000000000000']

        # calendar = 'excel1904'
        >>> decimal = date2dec(calendar='excel1904', ascii=d[:nn/2])
        >>> print(astr(decimal[:nn/2], 7, pp=True))
        ['35069.5210069' '16949.6841435' ' 2417.4450810']

        # calendar = '365_day'
        >>> g = np.array(['18.08.1972 12:30:15', '25.10.0986 12:30:15', '28.11.0493 22:20:40', '01.01.0001 00:00:00'])
        >>> decimal = date2dec(calendar='365_day', ascii=g)
        >>> nn = g.size
        >>> print(astr(decimal[:nn],14,pp=True))
        ['719644.52100694458932' '359822.52100694435649' '179911.93101851851679' '     0.00000000000000']

        # calendar = '366_day'
        >>> decimal = date2dec(calendar='366_day', ascii=g)
        >>> print(astr(decimal[:nn],14,pp=True))
        ['721616.52100694458932' '360808.52100694435649' '180404.93101851851679' '     0.00000000000000']

        # 360_day does not work with netcdftime.py version equal or below 0.9.2
        # calendar = '360_day'
        >>> decimal = date2dec(calendar='360_day', ascii=g)
        >>> print(astr(decimal[:nn],14,pp=True))
        ['709787.52100694458932' '354894.52100694435649' '177447.93101851851679' '     0.00000000000000']

        >>> print(astr(date2dec(yr=1992, mo=1, dy=26, hr=2, mi=0, sc=0, calendar='decimal'),14,pp=True))
       	1992.06853370763201
        >>> print(astr(date2dec(ascii='26.01.1992 02:00', calendar='decimal360'),14,pp=True))
       	1992.06967593592572
        >>> print(astr(date2dec(ascii=['26.01.1992 02:00','26.01.1992 02:00'], calendar='decimal360'),14,pp=True))
        ['1992.06967593592572' '1992.06967593592572']
        >>> print(astr(date2dec(yr=[1992,1992], mo=1, dy=26, hr=2, mi=0, sc=0, calendar='decimal360'),14,pp=True))
        ['1992.06967593592572' '1992.06967593592572']
        >>> print(astr(date2dec(yr=np.array([1992,1992]), mo=1, dy=26, hr=2, mi=0, sc=0, calendar='decimal360'),
        ...            14,pp=True))
        ['1992.06967593592572' '1992.06967593592572']
        >>> print(astr(date2dec(ascii=[['26.01.1992 02:00','26.01.1992 02:00'],
        ...                            ['26.01.1992 02:00','26.01.1992 02:00'],
        ...                            ['26.01.1992 02:00','26.01.1992 02:00']], calendar='decimal360'),14,pp=True))
        [['1992.06967593592572' '1992.06967593592572']
         ['1992.06967593592572' '1992.06967593592572']
         ['1992.06967593592572' '1992.06967593592572']]
        >>> print((date2dec(ascii='01.03.2003 00:00:00') - date2dec(ascii='01.03.2003')) == 0.)
        True


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

        Copyright 2010-2013 Arndt Piayda, Matthias Cuntz


        History
        -------
        Written  AP, Jun 2010
        Modified MC, Feb 2012 - All input can be scalar, list or array, also a mix
                              - Changed checks for easier extension
                              - decimal, decimal360
                 MC, Dec 2012 - change unit of proleptic_gregorian
                                from 'days since 0001-01-01 00:00:00'
                                to   'days since 0001-01-00 00:00:00'
                 MC, Feb 2013 - solved Excel leap year problem.
                 MC, Feb 2013 - ported to Python 3
                 MC, Jul 2013 - ascii/eng without time defaults to 00:00:00
                 MC, Oct 2013 - Excel starts at 1 not at 0
                 MC, Oct 2013 - units bugs, e.g. 01.01.0001 was substracted if Julian calendar even with units
                 MC, Nov 2013 - removed remnant of time treatment before time check in eng keyword
                 MC, Jun 2015 - adapted to new netCDF4/netcdftime (>= v1.0) and datetime (>= Python v2.7.9)
                 MC, Oct 2015 - call date2num with list instead of single netCDF4.datetime objects
                 DS, Dec 2016 - fixed import problem with newer versions of netcdftime
    """

    #
    # Constants
    calendars = ['standard', 'gregorian', 'julian', 'proleptic_gregorian',
                 'excel1900', 'excel1904', '365_day', 'noleap', '366_day',
                 'all_leap', '360_day', 'decimal', 'decimal360']
    #
    # Checks
    try:
        import netcdftime as nt
        nt.date2dec # check if method exists
        if ((nt.__version__ <= '0.9.2') & (calendar == '360_day')):
            raise ValueError("date2dec error: Your version of netcdftime.py is equal"
                             " or below 0.9.2. The 360_day calendar does not work with"
                             " arrays here. Please download a newer one.")
    except (ImportError, AttributeError):
        import netCDF4 as nt
        try:
            nt.datetime
        except AttributeError:
            # the date functions in netCDF (e.g. date2num) are only forwarded from a module
            # called cftime. At one point in time the implementation of netCDF4 changed and
            # the function datetime is not explicitly propagated anymore, so let's monkey
            # patch that here
            import cftime
            nt.datetime = cftime.datetime

    calendar = calendar.lower()
    if (calendar not in calendars):
        raise ValueError("date2dec error: Wrong calendar!"
                    " Choose: "+''.join([i+' ' for i in calendars]))
    # if ascii input is given by user, other input will be neglected
    # calculation of input size and shape
    if (ascii is not None) and (eng is not None):
        raise ValueError("date2dec error: 'ascii' and 'eng' mutually exclusive")
    if ascii is not None:
        islist = type(ascii) != type(np.array(ascii))
        isarr = np.ndim(ascii)
        if (islist & (isarr > 2)):
            raise ValueError("date2dec error: ascii input is list > 2D; Use array input")
        if isarr == 0: ascii = np.array([ascii])
        else: ascii = np.array(ascii)
        insize   = ascii.size
        outsize  = insize
        outshape = ascii.shape
        asciifl  = ascii.flatten()
        timeobj  = np.zeros(insize, dtype=object)
        # slicing of ascii strings to implement in datetime object. missing times
        # will be set to 00.
        yr = np.zeros(insize, dtype=np.int)
        mo = np.zeros(insize, dtype=np.int)
        dy = np.zeros(insize, dtype=np.int)
        hr = np.zeros(insize, dtype=np.int)
        mi = np.zeros(insize, dtype=np.int)
        sc = np.zeros(insize, dtype=np.int)
        for i in range(insize):
            aa      = asciifl[i].split('.')
            dy[i]   = int(aa[0])
            mo[i]   = int(aa[1])
            tail    = aa[2].split()
            yr[i]   = int(tail[0])
            if len(tail) > 1:
                tim     = tail[1].split(':')
                hr[i]   = int(tim[0])
                if len(tim) > 1:
                    mi[i] = int(tim[1])
                else:
                    mi[i] = 00
                if len(tim) > 2:
                    sc[i] = int(tim[2])
                else:
                    sc[i] = 00
            else:
                hr[i] = 00
                mi[i] = 00
                sc[i] = 00
            if ((yr[i]==1900) & (mo[i]==2) & (dy[i]==29)):
                timeobj[i] = nt.datetime(yr[i], 3, 1, hr[i], mi[i], sc[i])
            else:
                timeobj[i] = nt.datetime(yr[i], mo[i], dy[i], hr[i], mi[i], sc[i])
    if eng is not None:
        islist = type(eng) != type(np.array(eng))
        isarr = np.ndim(eng)
        if isarr == 0:
             eng = np.array([eng])
        else:
             eng = np.array(eng)
        if (islist & (isarr > 2)):
            raise ValueError("date2dec error: eng input is list > 2D; Use array input")
        eng = np.array(eng)
        insize   = eng.size
        outsize  = insize
        outshape = eng.shape
        engfl    = eng.flatten()
        timeobj  = np.zeros(insize, dtype=object)
        # slicing of eng strings to implement in datetime object. missing times
        # will be set to 00.
        yr = np.zeros(insize, dtype=np.int)
        mo = np.zeros(insize, dtype=np.int)
        dy = np.zeros(insize, dtype=np.int)
        hr = np.zeros(insize, dtype=np.int)
        mi = np.zeros(insize, dtype=np.int)
        sc = np.zeros(insize, dtype=np.int)
        for i in range(insize):
            ee      = engfl[i].split('-')
            yr[i]   = int(ee[0])
            mo[i]   = int(ee[1])
            tail    = ee[2].split()
            dy[i]   = int(tail[0])
            if len(tail) > 1:
                tim     = tail[1].split(':')
                hr[i]   = int(tim[0])
                if len(tim) > 1:
                    mi[i] = int(tim[1])
                else:
                    mi[i] = 00
                if len(tim) > 2:
                    sc[i] = int(tim[2])
                else:
                    sc[i] = 00
            else:
                hr[i] = 00
                mi[i] = 00
                sc[i] = 00
            if ((yr[i]==1900) & (mo[i]==2) & (dy[i]==29)):
                timeobj[i] = nt.datetime(yr[i], 3, 1, hr[i], mi[i], sc[i])
            else:
                timeobj[i] = nt.datetime(yr[i], mo[i], dy[i], hr[i], mi[i], sc[i])
    # if no ascii input, other inputs will be concidered
    # calculation of input sizes, shapes and number of axis
    if ((ascii is None) and (eng is None)):
        isnlist1 = type(yr) == type(np.array(yr))
        isarr1   = np.ndim(yr)
        if isarr1 == 0: yr = np.array([yr])
        else: yr = np.array(yr)
        isnlist2 = type(mo) == type(np.array(mo))
        isarr2   = np.ndim(mo)
        if isarr2 == 0: mo = np.array([mo])
        else: mo = np.array(mo)
        isnlist3 = type(dy) == type(np.array(dy))
        isarr3   = np.ndim(dy)
        if isarr3 == 0: dy = np.array([dy])
        else: dy = np.array(dy)
        isnlist4 = type(hr) == type(np.array(hr))
        isarr4   = np.ndim(hr)
        if isarr4 == 0: hr = np.array([hr])
        else: hr = np.array(hr)
        isnlist5 = type(mi) == type(np.array(mi))
        isarr5   = np.ndim(mi)
        if isarr5 == 0: mi = np.array([mi])
        else: mi = np.array(mi)
        isnlist6 = type(sc) == type(np.array(sc))
        isarr6   = np.ndim(sc)
        if isarr6 == 0: sc = np.array([sc])
        else: sc = np.array(sc)
        islist = not (isnlist1 | isnlist2 | isnlist3 | isnlist4 | isnlist5 | isnlist6)
        isarr  = isarr1 + isarr2 + isarr3 + isarr4 + isarr5 + isarr6
        shapes = [np.shape(yr), np.shape(mo), np.shape(dy), np.shape(hr), np.shape(mi), np.shape(sc)]
        nyr    = np.size(yr)
        nmo    = np.size(mo)
        ndy    = np.size(dy)
        nhr    = np.size(hr)
        nmi    = np.size(mi)
        nsc    = np.size(sc)
        sizes  = [nyr,nmo,ndy,nhr,nmi,nsc]
        nmax   = np.amax(sizes)
        ii     = np.argmax(sizes)
        outshape = shapes[ii]
        if (islist & (np.size(outshape) > 2)):
            raise ValueError("date2dec error: input is list > 2D; Use array input.")
        if nyr < nmax:
            if nyr == 1: yr  = np.ones(outshape,)*yr
            else: raise ValueError("date2dec error: size of yr != max input or 1.")
        if nmo < nmax:
            if nmo == 1: mo  = np.ones(outshape)*mo
            else: raise ValueError("date2dec error: size of mo != max input or 1.")
        if ndy < nmax:
            if ndy == 1: dy  = np.ones(outshape)*dy
            else: raise ValueError("date2dec error: size of dy != max input or 1.")
        if nhr < nmax:
            if nhr == 1: hr  = np.ones(outshape)*hr
            else: raise ValueError("date2dec error: size of hr != max input or 1.")
        if nmi < nmax:
            if nmi == 1: mi  = np.ones(outshape)*mi
            else: raise ValueError("date2dec error: size of mi != max input or 1.")
        if nsc < nmax:
            if nsc == 1: sc  = np.ones(outshape)*sc
            else: raise ValueError("date2dec error: size of sc != max input or 1.")
        indate  = [yr, mo, dy, hr, mi, sc]
        insize  = [np.size(i) for i in indate]
        inshape = [np.shape(i) for i in indate]
        dimnum  = [len(i) for i in inshape]
        # calculation of maximum input size and maximum number of axis for
        # reshaping the output
        indate  = [i.flatten() for i in indate]
        outsize = max(insize)
        timeobj = np.zeros(outsize, dtype=object)
        # datetime object is constructed
        for i in range(outsize):
            iyr = int(indate[0][i])
            imo = int(indate[1][i])
            idy = int(indate[2][i])
            ihr = int(indate[3][i])
            imi = int(indate[4][i])
            isc = int(indate[5][i])

            if ((iyr==1900) & (imo==2) & (idy==29)):
                timeobj[i] = nt.datetime(iyr, 3, 1, ihr, imi, isc)
            else:
                timeobj[i] = nt.datetime(iyr, imo, idy, ihr, imi, isc)
    # depending on chosen calendar and optional set of the time units
    # decimal date is calculated
    output = np.zeros(outsize)
    t0    = nt.datetime(1582, 10, 5, 00, 00, 00)
    t1    = nt.datetime(1582, 10, 15, 00, 00, 00)
    is121 = True if (min(timeobj)<t0) and (max(timeobj)>=t1) else False
    if (calendar == 'standard') or (calendar == 'gregorian'):
        if units is None:
            units = 'days since 0001-01-01 12:00:00'
            dec0 = 1721424
        else:
            dec0 = 0
        if is121 and (nt.__version__ < '1.2.2'):
            for ii, tt in enumerate(timeobj): output[ii] = nt.date2num(tt, units, calendar='gregorian')+dec0
        else:
            output = nt.date2num(timeobj, units, calendar='gregorian')+dec0
    elif calendar == 'julian':
        if units is None:
            units = 'days since 0001-01-01 12:00:00'
            dec0 = 1721424
        else:
            dec0 = 0
        if is121 and (nt.__version__ < '1.2.2'):
            for ii, tt in enumerate(timeobj): output[ii] = nt.date2num(tt, units, calendar='julian')+dec0
        else:
            output = nt.date2num(timeobj, units, calendar='julian')+dec0
    elif calendar == 'proleptic_gregorian':
        if units is None: units = 'days since 0001-01-01 00:00:00'
        if is121 and (nt.__version__ < '1.2.2'):
            for ii, tt in enumerate(timeobj): output[ii] = nt.date2num(tt, units, calendar='proleptic_gregorian')
        else:
            output = nt.date2num(timeobj, units, calendar='proleptic_gregorian')
    elif calendar == 'excel1900':
        doerr = False
        if units is None:
            units = 'days since 1899-12-31 00:00:00'
            if excelerr: doerr = True
        if is121 and (nt.__version__ < '1.2.2'):
            for ii, tt in enumerate(timeobj): output[ii] = nt.date2num(tt, units, calendar='gregorian')
        else:
            output = nt.date2num(timeobj, units, calendar='gregorian')
        if doerr:
            output = np.where(output >= 60., output+1., output)
            # date2num treats 29.02.1900 as 01.03.1990, i.e. is the same decimal number
            if np.any((output >= 61.) & (output < 62.)):
                for i in range(outsize):
                    # if (timeobj[i].year==1900) & (timeobj[i].month==2) & (timeobj[i].day==29):
                    #     output[i] -= 1.
                    if (yr[i]==1900) & (mo[i]==2) & (dy[i]==29):
                        output[i] -= 1.
    elif calendar == 'excel1904':
        if units is None: units = 'days since 1903-12-31 00:00:00'
        if is121 and (nt.__version__ < '1.2.2'):
            for ii, tt in enumerate(timeobj): output[ii] = nt.date2num(tt, units, calendar='gregorian')
        else:
            output = nt.date2num(timeobj, units, calendar='gregorian')
    elif (calendar == '365_day') or (calendar == 'noleap'):
        if units is None: units = 'days since 0001-01-01 00:00:00'
        if is121 and (nt.__version__ < '1.2.2'):
            for ii, tt in enumerate(timeobj): output[ii] = nt.date2num(tt, units, calendar='365_day')
        else:
            output = nt.date2num(timeobj, units, calendar='365_day')
    elif (calendar == '366_day') or (calendar == 'all_leap'):
        if units is None: units = 'days since 0001-01-01 00:00:00'
        if is121 and (nt.__version__ < '1.2.2'):
            for ii, tt in enumerate(timeobj): output[ii] = nt.date2num(tt, units, calendar='366_day')
        else:
            output = nt.date2num(timeobj, units, calendar='366_day')
    elif calendar == '360_day':
        if units is None: units = 'days since 0001-01-01 00:00:00'
        if is121 and (nt.__version__ < '1.2.2'):
            for ii, tt in enumerate(timeobj): output[ii] = nt.date2num(tt, units, calendar='360_day')
        else:
            output = nt.date2num(timeobj, units, calendar='360_day')
    elif calendar == 'decimal':
        ntime = np.size(yr)
        leap  = np.array((((yr%4)==0) & ((yr%100)!=0)) | ((yr%400)==0))
        tdy   = np.array(dy, dtype=np.float)
        diy   = np.array([ [-9,0, 31, 59, 90,120,151,181,212,243,273,304,334,365],
                           [-9,0, 31, 60, 91,121,152,182,213,244,274,305,335,366] ])
        for i in range(ntime):
            tdy[i] = tdy[i] + np.array(diy[leap[i],mo[i]], dtype=np.float)
        days_year = 365.
        output    = ( np.array(yr, dtype=np.float) +
                      ((tdy-1.)*24. + np.array(hr, dtype=np.float) +
                       np.array(mi, dtype=np.float)/60. +
                       np.array(sc, dtype=np.float)/3600.) /
                       ((days_year+np.array(leap, dtype=np.float))*24.) )
        # for numerical stability, i.e. back and forth transforms
        output += 1e-08 # 1/3 sec
    elif calendar == 'decimal360':
        ntime = np.size(yr)
        tdy   = np.array(dy, dtype=np.float)
        diy   = np.array([-9,  0, 30, 60, 90,120,150,180,210,240,270,300,330,360])
        for i in range(ntime):
            tdy[i] = tdy[i] + np.array(diy[mo[i]], dtype=np.float)
        days_year = 360.
        output    = ( np.array(yr, dtype=np.float) +
                      ((tdy-1.)*24. + np.array(hr, dtype=np.float) +
                       np.array(mi, dtype=np.float)/60. +
                       np.array(sc, dtype=np.float)/3600.) /
                       (days_year*24.) )
        # for numerical stability, i.e. back and forth transforms
        output += 1e-08 # 1/3 sec
    else:
        raise ValueError("date2dec error: calendar not implemented; should have been catched before.")


    # return of reshaped output
    output = np.reshape(output, outshape)
    if isarr == 0:
        output = np.float(output)
    else:
        if islist:
            ns = np.size(outshape)
            if ns == 1:
                output = [i for i in output]
            else:
                loutput = [ i for i in output[:,0]]
                for i in range(np.size(output[:,0])):
                    loutput[i] = list(np.squeeze(output[i,:]))
                output = loutput

    return output


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # year   = np.array([2000,1810,1630,1510,1271,619,1,1])#-1579,-4712])
    # month  = np.array([1,4,7,9,3,8,8,1])
    # day    = np.array([5,24,15,20,18,27,23,1])
    # hour   = np.array([12,16,10,14,19,11,20,12])
    # minute = np.array([30,15,20,35,41,8,3,0])
    # second = np.array([15,10,40,50,34,37,41,0])
    # decimal = date2dec(calendar = 'standard', yr=year, mo=month, dy=day, hr=hour, mi=minute, sc=second)
    # from autostring import astr
    # nn = year.size
    # print(astr(decimal[:nn/2],14,pp=True))
    # #    ['2.45154902100694e+06' '2.38226217719907e+06' '2.31660093101852e+06' '2.27284810821759e+06']
    # print(astr(decimal[nn/2:],14,pp=True))
    # #    ['2.18536732053241e+06' '1.94738596431713e+06' '1.14456333589120e+06' '0.00000000000000e+00']
    # decimal = date2dec(calendar='standard', yr=year, mo=6, dy=15, hr=12, mi=minute, sc=second)
    # print(astr(decimal[:nn/2],14,pp=True))
    # #    ['2.45171102100694e+06' '2.38231401053241e+06' '2.31657101435185e+06' '2.27275102488426e+06']
    # print(astr(decimal[nn/2:],14,pp=True))
    # #    ['2.18545602886574e+06' '1.94731300598380e+06' '1.14449400255787e+06' '1.66000000000000e+02']

    # from autostring import astr
    # d = np.array(['05.01.2000 12:30:15', '27.05.1950 16:25:10', '13.08.1910 10:40:55',
    #               '01.03.1900 00:00:00', '29.02.1900 00:00:00', '28.02.1900 00:00:00',
    #               '01.01.1900 00:00:00'])
    # decimal = date2dec(calendar='excel1900', ascii=d)
    # nn = d.size
    # print(astr(decimal[:nn/2], 7, pp=True))
    # #    ['36530.5210069' '18410.6841435' ' 3878.4450810']
    # print(astr(decimal[nn/2:],14,pp=True))
    # #    ['61.00000000000000' '60.00000000000000' '59.00000000000000' ' 1.00000000000000']
