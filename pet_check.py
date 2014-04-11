#
# script for checking pet results
#
# written: Matthias Zink, Feb 2014
#
import numpy as np
import math  as math

deg2rad  = math.pi /180

tavg     = 19.3178043365479
tmax     = 28.1099205017090
tmin     = 11.5728397369385
julian   = 250
latitude = 51.7410582896352
netrad   = 95.3826454569384

pet_mhm_hs  = 3.86317306435648


lat        =  latitude * deg2rad

dr         = 1. + 0.033 * np.cos(2. * math.pi * julian / 365.)
delta      = 0.4093     * np.sin(2. * math.pi * julian / 365. - 1.39)
omega      = np.arccos( -np.tan(lat) * np.tan(delta) )

conversion = 0.408 * 24 * 60 * 0.082 / math.pi 
#            DaySecs / PI_D    / SpecHeatET_dp * SolarConst_dp
conversion = 86400.  / math.pi /       2.45e06 * 1367.0


Ra         = conversion * dr *(omega * np.sin(lat) * np.sin(delta) + 
                                       np.cos(lat) * np.cos(delta) * np.sin(omega)) 
PET_HS     = 2.3*10**(-3) * Ra * (tavg + 17.8) * np.sqrt( tmax - tmin )

       # delta1       delta2
delta      = 0.04145 * np.exp(0.06088 * tavg)
PET_PT     = 1.26 * delta / (0.0646 + delta) * netrad *  86400. / 2.45e06
