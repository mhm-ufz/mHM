#!/usr/bin/env python
from __future__ import print_function
"""
usage: usgs2mrm.py [-h] [-c 1] [-l headerlines] [-u usgsidline] [-o outfile] [infile]

This is the Python script to convert tab-separated discharge values from usgs (for example,
obtained from http://waterdata.usgs.gov/nwis) to mRM discharge files. Values can also be
converted from [ft3 s-1] to [m3 s-1]. If a value is missing, nodata value -9999  will be
used. If the flag by USGS is indicating that the data is not approved, a warning for this
timestep is going to be issued, but the values are used nevertheless.

positional arguments:
  infile                Mandatory input file containing the USGS discharge data

optional arguments:
  -h, --help            show this help message and exit
  -c 1                  flag for unit conversion from ft to m (1 - yes (default), 0 - no)
  -l headerlines        number of headerlines (default 28)
  -u usgsidline         line containing the USGS ID (default 14)
  -o outfile            output file, if not given, output is written to standard out

Example:
python usgs2mrm.py -c 1 -l 28 -u 14 -o mrm.txt usgs.txt

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

Copyright 2012-2015 Stephan Thober


History
-------
Written,  Stephan Thober, Oct 2015
Modified,
         
"""
if __name__ == '__main__':
    import numpy as np
    from re import search
    from sys import stdout

    convert = 1
    infile = ''
    outfile = ''
    headlines = 28
    usgsIDline = 14
    nodata = -9999.
    # obtain arguments
    from argparse import ArgumentParser, RawDescriptionHelpFormatter
    outfile = ''
    parser  = ArgumentParser(
    formatter_class = RawDescriptionHelpFormatter,
        description = '''This is the Python script to convert tab-separated discharge values from usgs (for example,
obtained from http://waterdata.usgs.gov/nwis) to mRM discharge files. Values can also be
converted from [ft3 s-1] to [m3 s-1]. If a value is missing, nodata value -9999 will be
used. If the flag by USGS is indicating that the data is not approved, a warning for this
timestep is going to be issued, but the values are used nevertheless.
''')
    parser.add_argument('-c', '--convert',
                        action = 'store',
                        default = convert,
                        dest = 'convert',
                        metavar = 'convert',
                        help = "convert from ft to m (1 - yes (default), 0 - no)")
    parser.add_argument('-l', '--headerline',
                        action = 'store',
                        default = headlines,
                        dest = 'headlines',
                        metavar = 'headlines',
                        help = "number of header lines in input file")
    parser.add_argument('-u', '--usgsidline',
                        action = 'store',
                        default = usgsIDline,
                        dest = 'usgsIDline',
                        metavar = 'usgsIDline',
                        help = "line containing the USGS ID")
    parser.add_argument('-o', '--outfile',
                        action = 'store',
                        default = outfile,
                        dest = 'outfile',
                        metavar = 'outfile',
                        help = "output file, if not given, output is written to standard out")
    parser.add_argument('file',
                        nargs = '?',
                        default = None,
                        metavar = 'infile',
                        help = 'Mandatory input file containing the lat lon coordinates')
    args = parser.parse_args()
    convert = args.convert
    headlines = args.headlines
    usgsIDline = args.usgsIDline
    infile = args.file
    outfile = args.outfile
    del parser, args

    # factor for convert ft^3 in mm^3
    if np.int(convert) == 1:
        ft2m = 0.3048 # one ft is 0.3048 m
    else:
        ft2m = 1.

    # read file
    fi = open(infile, 'r')
    indata = fi.readlines()
    fi.close()

    # get usgs ID from indata
    usgsid = search('[0-9]{8}', indata[usgsIDline]).group()

    # create arrays
    Ndata = len(indata) - headlines
    # create list for dates 
    date = list()
    # create array for measurements
    value = np.zeros(Ndata)

    # read value and date
    for ii in np.arange(Ndata):
        readStr = indata[headlines + ii].split('\t')
        if readStr[-1] != 'A\n':
            print('***Warning: value at date ' + readStr[2] + ' is not approved!')
        date.append(readStr[2].replace('-', '  '))
        if readStr[3] == '':
            value[ii] = nodata
        else:
            value[ii] = np.float(readStr[3]) * ft2m

    # open standard out
    if outfile != '':
        fo = open(outfile, 'w')
    else:
        fo = stdout
    # write header
    fo.write(usgsid + ' Gauge 1 (daily discharge)\n')
    fo.write('nodata{:11.3f}\n'.format(nodata))
    fo.write('n       1       measurements per day [1, ' + str(Ndata) + ']\n')
    fo.write('start  ' + date[0] + ' 00 00   (YYYY MM DD HH MM)\n')
    fo.write('end    ' + date[-1] + ' 00 00   (YYYY MM DD HH MM)\n')
    # write data
    for ii in np.arange(Ndata):
        writeStr = date[ii] + '  00  00{:11.3f}\n'.format(value[ii])
        fo.write(writeStr)
    if outfile != '':
        fo.close()
    print('Done!')
