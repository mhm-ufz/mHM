#!/usr/bin/env python

#import os
#import sys

#sys.path.append( os.environ['HOME'] + '/Source/lib/chs-svn/PYTHON_chs_lib/' )
#sys.path.append( './sas/' )

#import ufz
import numpy as np
from sas import SAS
#from sas.plot_fun import plot_basin
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature         

#-- paths ---------------------------------------------------------------------

project_path = '../'

#-- initialize mhm class ------------------------------------------------------ 

mhm = SAS(model_path = project_path, file_name = 'mhm.nml', 
                 rel_path_name = True)

#-- import mhm data (input and output) ---------------------------------------- 

mhm.import_data('lat_lon')
mhm.import_data('slope')
mhm.import_data('soil_class')
mhm.import_data('lai_class')
mhm.import_data('dem')
mhm.import_data('pre')

mhm.import_data('states_and_fluxes', 'all')
mhm.combine_variables(['SWC'])

#-- postprocessing mhm --------------------------------------------------------

mhm.get_p(p_type = 'forward', time_series = 'all')

#-- plotting ------------------------------------------------------------------

#plot_basin(mhm.lon, mhm.lat, np.mean(mhm.recharge, axis = 0))
#plot_basin(mhm.lon, mhm.lat, mhm.mean_tt_opt)
