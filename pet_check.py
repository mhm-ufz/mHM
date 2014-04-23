#
# script for checking pet results
#
# written: Matthias Zink, Feb 2014
#
import numpy as np
import math  as math

deg2rad  = math.pi /180

julian   = 250
tavg     = 18.106281280517578
tmax     = 27.007287979125977
tmin     = 10.096543312072754
latitude = 48.108765567538860
netrad   = 116.85641364103674
eabs     = 1359.7239990234374 / 1000.

aer_res       = 100.0
bulksurf_res  = 100.0

# mhm werte
pet_mhm_hs  = 3.9929043540247986
pet_mhm_pt  = 3.4734362578081401
pet_mhm_pm  = 0
slope_sat_mhm = 0.13053207280816087
# konstanten
rho     = 1.225 
cp      = 1005.0
psychro = 0.0646

##############################
lat        =  latitude * deg2rad
dr         = 1. + 0.033 * np.cos(2. * math.pi * julian / 365.)
delta      = 0.4093     * np.sin(2. * math.pi * julian / 365. - 1.39)
omega      = np.arccos( -np.tan(lat) * np.tan(delta) )
#            DaySecs / PI_D    / SpecHeatET_dp * SolarConst_dp
conversion = 86400.  / math.pi /       2.45e06 * 1367.0
Ra         = conversion * dr *(omega * np.sin(lat) * np.sin(delta) + 
                                       np.cos(lat) * np.cos(delta) * np.sin(omega)) 

PET_HS     = 2.3*10**(-3) * Ra * (tavg + 17.8) * np.sqrt( tmax - tmin )

       # delta1       delta2
delta      = 0.04145 * np.exp(0.06088 * tavg)
e0_abs     = 0.6108 * np.exp(17.27 * tavg / (tavg + 237.3))
slop_sat   = 4098. * e0_abs / (tavg + 237.3)**2

PET_PT     = 1.26 * slop_sat  / (0.0646 + slop_sat ) * netrad *  86400. / 2.45e06


PET_PM     = (
  (slop_sat * netrad + rho * cp * (e0_abs - eabs)/aer_res) /
  (slop_sat + psychro * (1 + bulksurf_res/aer_res))
            ) * 86400. / 2.45e06




print 'HarSam:  ', PET_HS, 'mMM: ', pet_mhm_hs
print 'PrieTay: ', PET_PT, 'mMM: ', pet_mhm_pt
print 'PenMon:  ', PET_PM, 'mMM: ', pet_mhm_pm
