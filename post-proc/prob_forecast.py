#!/usr/bin/env python
"""

This script creates probabilistic forecast from mHM simulations given the method described in Section 2.3.2 of 

Woldemeskel et al. Hydrol Earth Syst Sci, 2018 vol. 22 (12) pp. 6257-6278.
"Evaluating post-processing approaches for monthly and seasonal streamflow forecasts."
https://www.hydrol-earth-syst-sci.net/22/6257/2018/

It fits an AR(1) model to the additive residuals of streamflow and uses it to forecast.

License
-------
This Script can be redistributed and/or modified
under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This script is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You can find the license under <http://www.gnu.org/licenses/>.

Copyright 2019-2019 Stephan Thober

History
-------
Written,  Stephan Thober - Aug 2019
Modified,
         
"""

# IMPORTS
import argparse
import textwrap
from ufz.netcdf4 import NcDataset


# FUNCTIONS
def parser():

    discharge_file = "test_domain/output_b1/discharge.nc"

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
          description:
            This script creates probabilistic forecast from mHM simulations given the method described in Section 2.3.2 of 

            Woldemeskel et al. Hydrol Earth Syst Sci, 2018 vol. 22 (12) pp. 6257-6278.
            "Evaluating post-processing approaches for monthly and seasonal streamflow forecasts."
            https://www.hydrol-earth-syst-sci.net/22/6257/2018/

            It fits an AR(1) model to the additive residuals of streamflow and uses it to forecast.

            License
            -------
            This Script can be redistributed and/or modified
            under the terms of the GNU Lesser General Public License as published by
            the Free Software Foundation, either version 3 of the License, or
            (at your option) any later version.

            This script is distributed in the hope that it will be useful,
            but WITHOUT ANY WARRANTY; without even the implied warranty of
            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
            GNU Lesser General Public License for more details.

            You can find the license under <http://www.gnu.org/licenses/>.

            Copyright 2019-2019 Stephan Thober
            

          Example:
            python prob_forecast -f test_basin/output/discharge

          '''))

    parser.add_argument('-f', '--discharge_file', action='store',
                    default=discharge_file, dest='discharge_file', metavar='discharge_file',
                    help='Name of the input discharge file (default: test_domain/output_b1/discharge.nc)')
    
    # evaluate args
    args = parser.parse_args()
    discharge_file = args.discharge_file


    return discharge_file


def read_data(discharge_file):
    
    ncin = NcDataset(discharge_file, 'r')
    vnames = list(ncin.variables.keys())
    vnames.remove('time')

    for vv in vnames:
        if vv[:4] == 'Qobs':
            runoff_obs = ncin.variables[vv][:]
        elif vv[:4] == 'Qsim':
            runoff_sim = ncin.variables[vv][:]
        else:
            raise ValueError('Netcdf file {} does not seem to be a discharge file of mHM'.format(discharge_file))
    ncin.close()

    return runoff_obs, runoff_sim


if __name__ == '__main__':

    discharge_file = parser()

    runoff_obs, runoff_sim = read_data(discharge_file)

    print('Done!')
