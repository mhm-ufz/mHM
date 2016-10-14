#!/usr/bin/env python
from __future__ import print_function
#
# This script reads a DEM from a netcdf file and calculates the flow direction
# following ArcGis convention:
#
# flow direction is assumed like this
#           64             y-axis
#       32      128          |
#   16      -1       1       |
#        8       2          \|/
#            4               V
#      x-axis ------>
# ORIGIN is in the upper left corner
#
# sinks are marked by -1
#
# Nomenclature:
#
# co - channel orders
# fd - flow direction
# fa - flow accumulation
# yy - row index
# xx - column index
#
# author:  Stephan Thober, David Schaefer
# created: 01.07.2015
import numpy as np

# global variables for calculating upstream cells
# yy_offset, xx_offset = np.meshgrid(np.arange(-1, 2), np.arange(-1, 2))
yy_offset, xx_offset = np.meshgrid(np.arange(1, -2, -1), np.arange(-1, 2))
yy_offset = yy_offset.ravel() # yy_offset for neighboring cells
xx_offset = xx_offset.ravel() # xx_offset for neighboring cells
# local flow direction sorted in the order of yy_offset and xx_offset
local_flow_direction = np.array([8, 16, 32, 4, -1, 64, 2, 1, 128])
# same as local_flow_direction but reverted. This means these are the
# flow directions of neighboring cells flowing into the local cell.
inflow_direction = np.array([128, 1, 2, 64, -9999, 4, 32, 16, 8])


def cal_fdir(locs, fdir, factor):
    # flow direction is assumed like this
    #           64             y-axis
    #       32      128          |
    #   16      -1       1       |
    #        8       2          \|/
    #            4               V
    #      x-axis ------>
    #
    # ORIGIN is in the upper left corner -> y axis increases to the bottom
    # loc = [nrow, ncol]
    nloc = locs[0].shape[0]
    fds = np.zeros((nloc))
    for kk in np.arange(nloc):
        # loop over all locations having maximum flow direction
        loc = [locs[0][kk], locs[1][kk]]
        fd = fdir[loc[0], loc[1]]
        if fd == 2:
            if loc[1] + 1 < factor:
                fd = 4
            elif loc[0] + 1 < factor:
                fd = 1
        elif fd == 8:
            if loc[1] > 0:
                fd = 4
            elif loc[0] + 1 < factor:
                fd = 16
        elif fd == 32:
            if loc[1] > 0:
                fd = 64
            elif loc[0] > 0:
                fd = 16
        elif fd == 128:
            if loc[1] + 1 < factor:
                fd = 64
            elif loc[0] > 0:
                fd = 1
        fds[kk] = fd
    return fds


def upscale_fdir(sn, factor, print_info=False, return_maxlocs=False, do_co=False, redo_fa=True, missing_value=-9999.):
    """
        Upscales a river network by a factor (integer > 1), that has to be a divisible of the
        resolution of the flow direction. Direction is given by the cell with the largest flow
        accumulation. If multiple of these exist, then one is chosen randomly.
    
        Definition
        ----------
        upscale_fdir(sn, factor, print_info=False, return_maxlocs=False, do_co=False, redo_fa=True)


        Input
        -----
        sn            river_network object containing flow direction, flow accumulation and sinks
        factor        integer indicating by which factor the flow direction should be upscaled

        Optional Input Parameters
        -------------------------
        print_info     flag for printing additional information
        return_maxlocs flag for return locations of cell determining flow directions at given river_network object
        do_co          flag for calculating channel order
        redo_fa        flow recalculating flow accumulation at coarser river_network

        Options
        -------

        Output
        ------
        river_network object at coarser resolution with upscaled flow direction

        Restrictions
        ------------
            Given river_network object sn has to contain flow direction, accumulation and location of
            sink following the convention below.
            
            The origin of the flow direction field is assumed to be located in the upper left corner.
            Then flow directions are following this convention:

            Flow direction is assumed like this meaning that the ORIGIN IS THE UPPER LEFT CORNER
                      64             y-axis
                  32      128          |
              16      -1       1       |
                   8       2          \|/
                       4               V
                 x-axis ------>

            Sinks are marked by -1.

        Examples
        --------
        >>> # Create some data
        >>> fd = np.ma.array([[  2,   1,   1,   2,   4,   4,   8,  8,  8], 
        ...                   [  1,   2,   1,   1,   2,   4,   4,  4,  8], 
        ...                   [128,   1, 128,   1,   1,   2,   4,  4,  4], 
        ...                   [  1, 128,  64, 128, 128,   2,   4,  4,  8], 
        ...                   [128,  64,  64,  64,   1,   2,   4,  4,  4], 
        ...                   [ 64, 128,  64,  32,   1,   1,   2,  2,  4], 
        ...                   [128,  64,  64,  64,   1,   1,   1,  1,  1], 
        ...                   [128,  64, 128,  64,  32,   1, 128, 64, 64], 
        ...                   [ 64, 128,  64,  64,  64, 128,  64, 64, 32]])
        >>> sinks = np.array([[6], [8]])
        >>> sn = river_network(fdir=fd, do_fa=True, do_co=False, sinks=sinks)
        >>> upscale_fdir(sn, 3).fdir
        masked_array(data =
         [[1.0 2.0 4.0]
         [64.0 16.0 4.0]
         [64.0 64.0 1.0]],
                     mask =
         [[False False False]
         [False False False]
         [False False False]],
               fill_value = -9999.0)
        <BLANKLINE>      
               
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
        Modified                
    """
    # consistency checks
    if not np.ma.isMaskedArray(sn.fa):
        raise ValueError('***ERROR: upscale_fdir requires flow accumulation as masked array in river_network')
    if any(np.array(sn.fdir.shape) % factor != 0):
        raise ValueError('***ERROR: factor: ' + str(factor) + ' is not a divisible of flow direction shape')
    # create upscaled arrays
    new_shape = np.array(sn.fdir.shape) / factor
    new_fd = np.ma.masked_array(np.zeros(new_shape) + missing_value,
                                mask=np.ones(new_shape), fill_value=missing_value)
    new_fa = np.ma.masked_array(np.zeros(new_shape) + missing_value,
                                mask=np.ones(new_shape), fill_value=missing_value)

    # create maxlocs list
    maxlocs = []
    for ii in np.arange(new_fd.shape[0]):
        for jj in np.arange(new_fd.shape[1]):
            # extract part of map evaluated
            tmp_fa = sn.fa[ii * factor: (ii + 1) * factor, jj * factor: (jj + 1) * factor]
            tmp_fd = sn.fdir[ii * factor: (ii + 1) * factor, jj * factor: (jj + 1) * factor]
            if np.all(tmp_fa.mask):
                if print_info:
                    print('cell is masked ', ii, jj)
                # cell is masked
                new_fd[ii, jj] = -9999.
                new_fa[ii, jj] = -9999.
            else:
                # determine locations of maximum flow accumulations
                maxloc = np.ma.where(tmp_fa == np.ma.amax(tmp_fa))
                # calculate coarse scale flow direction
                coarse_fd = cal_fdir(maxloc, tmp_fd, factor)
                if maxloc[0].shape[0] > 1:
                    # if there is more than one outflow cell, check whether flow directions are different
                    if print_info:
                        print(coarse_fd, ' flow directions of same cells')
                    # evaluate when there are more than one cell if maximum flow directions are different
                    if np.any(np.diff(coarse_fd) > 0):
                        print('***Warning: multiple cells with same flow accumulation but different flow directions found, arbitrarily choose first one')
                # store flow direction and flow accumulation
                new_fd[ii, jj] = coarse_fd[0]
                new_fa[ii, jj] = tmp_fa[maxloc[0][0], maxloc[1][0]]
                if print_info:
                    print('maxloc = ', maxloc)
                    print('tmp_fd = ', tmp_fd[maxloc], coarse_fd[0])
                    print('tmp_fa = ', tmp_fa[maxloc])
                    print('----------------')
                    print('new_fd = ', new_fd[ii, jj])
                    print('new_fa = ', new_fa[ii, jj])
                    print('================')
                # add to store maximum locations
                maxlocs.append([maxloc[0] + ii * factor, maxloc[1] + jj * factor])
    # upscale sinks
    upscale_sinks = tuple(np.array(sn.sinks)/np.int(factor))
    # return
    if return_maxlocs:
        if redo_fa:
            return maxlocs, river_network(fdir=new_fd, do_co=do_co, do_fa=True, sinks=upscale_sinks)
        else:
            return maxlocs, river_network(fdir=new_fd, do_co=do_co, fa=new_fa)
    else:
        if redo_fa:
            return river_network(fdir=new_fd, do_co=do_co, do_fa=True, sinks=upscale_sinks)
        else:
            return river_network(fdir=new_fd, do_co=do_co, fa=new_fa)

        
class river_network(object):
    def __init__(self, dem=None, fdir=None, co=None, do_co=False, fa=None, do_fa=False, print_info=False, missing_value=-9999., sinks=None):
        """
            Initializes a river_network object describing the flow path of a river through the landscape.

            Definition
            ----------
            river_network(dem=None, fdir=None, co=None, do_co=False, fa=None, do_fa=False, missing_value=-9999., sinks=None):


            Input
            -----
            dem           digital elevation model (dem), basis for calculating flow direction fdir
            fdir          flow direction following convention below.

            Optional Input Parameters
            -------------------------
            co            channel order following Strahler 1952
            do_co         flag for calculating channel order
            fa            flow accumulation
            do_fa         flag for calculating flow accumulation
            print_info    flag for printing additional information
            missing_value default: -9999.
            sinks         location of sinks as two arrays (first/second for y/x coordinate)

            Options
            -------

            Output
            ------

            Restrictions
            ------------
                Either dem or fdir has to be provided, both cannot be omitted.
                
                The origin of the flow direction field is assumed to be located in the upper left corner.
                Then flow directions are following this convention:

                Flow direction is assumed like this meaning that the ORIGIN IS THE UPPER LEFT CORNER
                          64             y-axis
                      32      128          |
                  16      -1       1       |
                       8       2          \|/
                           4               V
                     x-axis ------>

                Sinks are marked by -1.

            Examples
            --------
            >>> # Create some data
            >>> fd = np.ma.array([[  2,   1,   1,   2,   4,   4,   8,  8,  8], 
            ...                   [  1,   2,   1,   1,   2,   4,   4,  4,  8], 
            ...                   [128,   1, 128,   1,   1,   2,   4,  4,  4], 
            ...                   [  1, 128,  64, 128, 128,   2,   4,  4,  8], 
            ...                   [128,  64,  64,  64,   1,   2,   4,  4,  4], 
            ...                   [ 64, 128,  64,  32,   1,   1,   2,  2,  4], 
            ...                   [128,  64,  64,  64,   1,   1,   1,  1,  1], 
            ...                   [128,  64, 128,  64,  32,   1, 128, 64, 64], 
            ...                   [ 64, 128,  64,  64,  64, 128,  64, 64, 32]])
            >>> sinks = np.array([[6], [8]])
            >>> river_network(fdir=fd, do_fa=True, do_co=False, sinks=sinks).fa
            masked_array(data =
             [[1.0 1.0 2.0 3.0 1.0 1.0 1.0 1.0 1.0]
             [1.0 4.0 1.0 32.0 37.0 3.0 2.0 2.0 1.0]
             [1.0 1.0 30.0 1.0 4.0 46.0 3.0 4.0 1.0]
             [1.0 5.0 19.0 2.0 1.0 1.0 50.0 5.0 2.0]
             [2.0 1.0 18.0 1.0 1.0 2.0 52.0 8.0 1.0]
             [1.0 6.0 2.0 9.0 1.0 2.0 57.0 9.0 2.0]
             [1.0 4.0 1.0 8.0 1.0 2.0 3.0 68.0 81.0]
             [2.0 1.0 3.0 2.0 2.0 1.0 4.0 3.0 1.0]
             [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0]],
                         mask =
             [[False False False False False False False False False]
             [False False False False False False False False False]
             [False False False False False False False False False]
             [False False False False False False False False False]
             [False False False False False False False False False]
             [False False False False False False False False False]
             [False False False False False False False False False]
             [False False False False False False False False False]
             [False False False False False False False False False]],
                   fill_value = -9.0)
            <BLANKLINE>
            >>> river_network(fdir=fd, do_fa=False, do_co=True, sinks=sinks).co
            masked_array(data =
             [[1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0]
             [1.0 2.0 1.0 3.0 3.0 2.0 1.0 1.0 1.0]
             [1.0 1.0 3.0 1.0 2.0 3.0 1.0 2.0 1.0]
             [1.0 2.0 3.0 1.0 1.0 1.0 3.0 2.0 1.0]
             [1.0 1.0 3.0 1.0 1.0 1.0 3.0 2.0 1.0]
             [1.0 2.0 1.0 2.0 1.0 1.0 3.0 2.0 1.0]
             [1.0 2.0 1.0 2.0 1.0 1.0 1.0 3.0 3.0]
             [1.0 1.0 2.0 1.0 1.0 1.0 2.0 2.0 1.0]
             [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0]],
                         mask =
             [[False False False False False False False False False]
             [False False False False False False False False False]
             [False False False False False False False False False]
             [False False False False False False False False False]
             [False False False False False False False False False]
             [False False False False False False False False False]
             [False False False False False False False False False]
             [False False False False False False False False False]
             [False False False False False False False False False]],
                   fill_value = -9.0)
            <BLANKLINE>
                   
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
            Written,  ST, Dec 2015
            Modified                
        """
        # initialize all arrays
        self.dem = None # digital elevation model
        self.fdir = None # flow direction
        self.sinks = None # sinks
        self.co = None # channel order
        self.fa = None # flow accumulation
        # consistency check
        if fdir is None  and dem is None:
            raise ValueError('***ERROR: specify either dem or fdir to create a river_network object')
        if fdir is None:
            if print_info:
                print('calculating flow direction')
            # create flow direction if necessary
            self.dem = dem
            self.fdir = self.flow_direction()
            if print_info:
                print('calculating flow direction... ok')
        else:
            self.fdir = fdir
        # assure that fdir is masked array
        if not np.ma.isMaskedArray(self.fdir):
            self.fdir = np.ma.array(self.fdir)
            self.fdir.mask = False
        # assign flow accumulation
        if not fa is None:
            self.fa = fa
        # assign channel order
        if not co is None:
            self.co = co
        # assign sinks
        if sinks is None and fa is None:
            raise ValueError('***ERROR: for initializing river network either the location of the sinks or flow accumulation has to be given')
        elif sinks is None and not fa is None:
            self.sinks = self._get_sinks()
        elif not sinks is None and fa is None:
            self.sinks = sinks
        # get channel order and flow accumulation
        if do_co and do_fa:
            self.co = np.ma.masked_array(np.zeros(self.fdir.shape) + missing_value,
                                         mask=np.ones(self.fdir.shape), fill_value=missing_value)
            self.fa = np.ma.masked_array(np.zeros(self.fdir.shape) + missing_value,
                                         mask=np.ones(self.fdir.shape), fill_value=missing_value)
            for ii in np.arange(self.sinks[0].shape[0]):
                self.co, self.fa = self.network_properties(self.fdir, self.sinks[0][ii], self.sinks[1][ii],
                                                           do_co=do_co, co=self.co,
                                                           do_fa=do_fa, fa=self.fa,
                                                           missing_value=missing_value,
                                                           print_info=print_info)
        elif do_co and not do_fa:
            self.co = np.ma.masked_array(np.zeros(self.fdir.shape) + missing_value,
                                         mask=np.ones(self.fdir.shape), fill_value=missing_value)
            for ii in np.arange(self.sinks[0].shape[0]):
                self.co = self.network_properties(self.fdir, self.sinks[0][ii], self.sinks[1][ii],
                                                  do_co=do_co, co=self.co,
                                                  do_fa=do_fa,
                                                  missing_value=missing_value,
                                                  print_info=print_info)
        elif not do_co and do_fa:
            self.fa = np.ma.masked_array(np.zeros(self.fdir.shape) + missing_value,
                                         mask=np.ones(self.fdir.shape), fill_value=missing_value)
            for ii in np.arange(self.sinks[0].shape[0]):
                self.fa = self.network_properties(self.fdir, self.sinks[0][ii], self.sinks[1][ii],
                                                  do_co=do_co,
                                                  do_fa=do_fa, fa=self.fa,
                                                  missing_value=missing_value,
                                                  print_info=print_info)
        

    def flow_direction(self, print_info=False):
        """
            Calculates flow direction from a DEM.

            Definition
            ----------
            def flow_direction(self, print_info=False)


            Input
            -----
            self         self - river_network object containing a dem array

            Optional Input Parameters
            -------------------------
            print_info - flag for printing additional information

            Options
            -------

            Output
            ------
            fd           array containing flow direction with the following convention

            Restrictions
            ------------
                The origin of the flow direction field is assumed to be located in the upper left corner.
                Then flow directions are following this convention:

                flow direction is assumed like this meaning that the ORIGIN IS THE UPPER LEFT CORNER
                          64             y-axis
                      32      128          |
                  16      -1       1       |
                       8       2          \|/
                           4               V
                     x-axis ------>

                Sinks are marked by -1.

            Examples
            --------
            >>> # Create some data
            >>> dem = np.ma.array([[ 30,   2,   1],
            ...                    [  5,  10,  25],
            ...                    [ 15,  23,  24]])
            >>> dem.mask = np.zeros(dem.shape, dtype='bool')
            >>> river_network(dem=dem, do_fa=False, do_co=False, sinks=np.array([[0], [2]])).fdir
            masked_array(data =
             [[1.0 1.0 -1.0]
             [128.0 128.0 64.0]
             [64.0 32.0 32.0]],
                         mask =
             [[False False False]
             [False False False]
             [False False False]],
                   fill_value = 1e+20)
                              
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

            Copyright 2009-2015 Stephan Thober, David Schaefer


            History
            -------
            Written,  ST & DS, Dec 2015
            Modified                
        """
        # global variable used: correct_direction
        fd = np.zeros(self.dem.shape)
        #
        for ii in np.arange(self.dem.shape[0]):
            for jj in np.arange(self.dem.shape[1]):
                if self.dem.mask[ii, jj]:
                    continue
                if print_info:
                    print('processing cell: ', ii, jj)
                # get mask of neighbors and y and x locations
                neighbors, yy, xx = self._get_neighbors(self.dem, ii, jj)
                # get position of cell with steepest gradient
                pos_min = np.ma.argmin(self.dem[yy, xx] - self.dem[ii, jj])
                fd[ii, jj] = local_flow_direction[neighbors][pos_min]
        return fd


    def network_properties(self, fd, yy, xx, print_info=False, do_co=True, co=None, do_fa=True, fa=None,
                           missing_value=-9999.):
        """
            Calculates channel order number and flow accumulation starting from one sink in a flow direction map

            channel order is calculated following Strahler 1952. It is ONE for headwater. If channels join, the
            channel order of the resulting stream is the highest one of the inflowing streams, if two or more than
            two inflowing streams have the highest channel order, the channel order of the resulting stream is one
            higher than the highest channel order of the inflowing streams.

            Definition
            ----------
            def network_properties(self, fd, yy, xx, print_info=False, do_co=True, co=None, do_fa=True, fa=None, missing_value=-9999.):


            Input
            -----
            self          self - river_network object
            fd            flow direction field, basically river_network.fd
            yy            row coordinate of sink
            xx            column coordinate of sink


            Optional Input Parameters
            -------------------------
            print_info    write additional info on std_out
            do_co         calculate channel order
            co            given channel order field
            do_fa         calculate flow accumulation
            fa            given flow accumulation
            missing_value floating point value for masking
            
            Options
            -------
            print_info   True: write additional information
                         False: do not write additional information (default)
            do_co        True: calculate channel order (default)
                         False: do not channel order
            co           None: no channel order field specified, will be created (default)
            do_fa        True: calculate flow accumulation (default)
                         False: do not flow accumulation
            fa           None: no flow accumulation field specified, will be created (default)

            Output
            ------
            Depending on options:
                co, fa if do_co=True and do_fa=True
                co if do_co=True and not do_fa=True
                fa if not do_co=True and do_fa=True


            Restrictions
            ------------
                The origin of the flow direction field is assumed to be located in the upper left corner.
                Then flow directions are following this convention:

                flow direction is assumed like this meaning that the ORIGIN IS THE UPPER LEFT CORNER
                          64             y-axis
                      32      128          |
                  16      -1       1       |
                       8       2          \|/
                           4               V
                     x-axis ------>

                Sinks are marked by -1.

            Examples
            --------

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

            Copyright 2009-2015 Stephan Thober, David Schaefer


            History
            -------
            Written,  ST & DS, Dec 2015
            Modified                
        """
        if co is None and do_co:
            co = np.ma.masked_array(np.zeros(fd.shape) + missing_value,
                                    mask=np.ones(fd.shape), fill_value=missing_value)
        if fa is None and do_fa:
            fa = np.ma.masked_array(np.zeros(fd.shape) + missing_value,
                                    mask=np.ones(fd.shape), fill_value=missing_value)
        # flow direction stack to emulate recursion
        if not do_co and not do_fa:
            raise ValueERROR('***ERROR: neither fa nor co calculated')
        fd_stack = [[yy, xx]] # start at initial sink
        while fd_stack:        
            if print_info:
                print('current flow accumulation stack: ', fd_stack)
            upstream = self._get_upstream(fd, fd_stack[-1])
            if print_info:
                print('upstream locations: ', upstream)
            if do_co:
                # use co for identifying upstream cells
                ext = [l for l in upstream if co.data[l[0],l[1]] == missing_value]
            else:
                # use fa for identifying upstream cells
                ext = [l for l in upstream if fa.data[l[0],l[1]] == missing_value]
            if ext:
                fd_stack.extend(ext)
                continue
            # all upstream cells are available
            # note that headwaters dont have an upstream cell
            cell = fd_stack.pop() # save cell
            if do_co:
                co_upstream = [co[loc[0], loc[1]] for loc in upstream]
                co_max = np.amax(co_upstream + [1])
                # if two streams of equal co merge, increment channel order   
                if len(np.where(co_upstream == co_max)[0]) > 1:
                    co_max += 1
                co.data[cell[0], cell[1]] = co_max
                co.mask[cell[0], cell[1]] = False
                if print_info:
                    print('co (channel order) of upstream: ', co_upstream)
            if do_fa:
                fa_upstream = [fa[loc[0], loc[1]] for loc in upstream]        
                fa.data[cell[0], cell[1]] = np.sum(fa_upstream) + 1
                fa.mask[cell[0], cell[1]] = False
                if print_info:
                    print('sum of upstream: ', np.sum(fa_upstream))
        if do_co and do_fa:
            return co, fa
        elif do_co and not do_fa:
            return co
        elif not do_co and do_fa:
            return fa


    def _get_sinks(self):
        # set sinks to maximum flow accumulation
        return np.ma.where(self.fa == np.amax(self.fa))


    def _get_neighbors(self, arr, yy_loc, xx_loc):
        # global variables used: yy_offset, xx_offset
        #
        yy_ind = yy_offset + yy_loc
        xx_ind = xx_offset + xx_loc
        # create mask for valid neighbors
        neighbors = ((yy_ind >= 0) &
                     (yy_ind < arr.shape[0]) &
                     (xx_ind >= 0) &
                     (xx_ind < arr.shape[1]))
        return neighbors, yy_ind[neighbors], xx_ind[neighbors]


    def _get_upstream(self, fd, loc):
        # global variable used: inflow_direction
        #
        # get mask of neighbors and y and x locations
        neighbors, yy, xx = self._get_neighbors(fd, loc[0], loc[1])
        # mask inflowing cells
        upstream_mask = (fd.data[yy, xx] == inflow_direction[neighbors])
        yy_upstream = yy[upstream_mask]
        xx_upstream = xx[upstream_mask]
        return [[yy_upstream[ii], xx_upstream[ii]] for ii in np.arange(np.sum(upstream_mask))]


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
