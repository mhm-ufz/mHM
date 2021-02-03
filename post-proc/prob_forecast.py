#!/usr/bin/env python
"""

This script creates probabilistic forecast from mHM simulations given the method described in Section 2.3.2 of 

Woldemeskel et al. Hydrol Earth Syst Sci, 2018 vol. 22 (12) pp. 6257-6278.
"Evaluating post-processing approaches for monthly and seasonal streamflow forecasts."
https://www.hydrol-earth-syst-sci.net/22/6257/2018/

It fits an AR(1) model to the additive residuals of streamflow and uses it to forecast.

It removes missing values from the simulated and observed data. This is not clean because jumps might occur if wholes are filled. If period is long enough, this might not be a problem.

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
import numpy as np
from scipy.stats import boxcox
from shutil import copyfile
from time import asctime


# FUNCTIONS
def parser():

    discharge_file = "test_domain/output_b1/discharge.nc"
    out_file = "test_domain/output_b1/prob_discharge.nc"
    lmbda = 0.2
    n_prob = 100

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

          Note:
            A default lambda value of 0.2 is used, which is hard-coded in mHM. It should always be the same that is used in mHM for calbration (opti_function 32). The script will raise a ValueError if another value is used.

          '''))

    parser.add_argument('-f', '--discharge_file', action='store',
                    default=discharge_file, dest='discharge_file', metavar='discharge_file',
                    help='Name of the input discharge file (default: test_domain/output_b1/discharge.nc)')
    parser.add_argument('-l', '--lambda', action='store',
                    default=lmbda, dest='lmbda', metavar='lmbda',
                    help='Lambda parameter for boxcox transformation (default: 0.2)')
    parser.add_argument('-n', '--n_prob', action='store',
                    default=n_prob, dest='n_prob', metavar='n_prob',
                    help='Number of probabilistic forecasts to generate (default: 100)')
    parser.add_argument('-o', '--out_file', action='store',
                    default=out_file, dest='out_file', metavar='out_file',
                    help='Name of the output file (default: test_domain/output_b1/prob_discharge.nc)')
    
    # evaluate args
    args = parser.parse_args()
    discharge_file = args.discharge_file
    out_file = args.out_file
    lmbda = args.lmbda
    n_prob = args.n_prob

    return discharge_file, out_file, lmbda, n_prob


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


def boxcox_inv(x, lmbda=0.2):

    if (lmbda != 0.2):
        raise ValueError('Inversion of boxcox only implemented for lambda value of 0.2')

    return (x * lmbda + 1.)**(1./lmbda)


def calculate_param(runoff_obs, runoff_sim):
    n_sample = runoff_sim.shape[0]
    z_f = boxcox(runoff_sim, lmbda=lmbda)
    z_q = boxcox(runoff_obs, lmbda=lmbda)
    eta = z_q - z_f

    eta_mean = eta.mean()
    eta_std = eta.std()

    nu = (eta - eta_mean) / eta_std
    rho = np.corrcoef(nu[1:], nu[:-1])[0, 1]
    sigma_y = np.std(nu[1:] - nu[:-1])

    return z_f, eta_mean, eta_std, rho, sigma_y, n_sample


def sample_forecasts(n_prob, z_f, eta_mean, eta_std, rho, sigma_y, n_sample):
    runoff_prob = np.zeros((n_sample, n_prob))
    for nn in np.arange(n_prob):
        y_hat = np.random.normal(0., sigma_y, size=n_sample)
        nu_hat = y_hat
        for ii in np.arange(n_sample - 1) + 1:
            nu_hat[ii] = nu_hat[ii - 1] * rho + y_hat[ii]
            
        eta_hat = nu_hat * eta_std + eta_mean

        runoff_prob[:, nn] = boxcox_inv(z_f + eta_hat, lmbda=lmbda)
    return runoff_prob


def write_file(runoff_prob, eta_mean, eta_std, rho, sigma_y, discharge_file, out_file):
    copyfile(discharge_file, out_file)

    ncout = NcDataset(out_file, 'a')

    ncout.createDimension("n_prob", n_prob)

    var = ncout.createVariable("runoff_prob", "f8", ("n_prob", "time"))
    var[:] = runoff_prob.transpose()
    var.createAttribute("long_name", "probabilistic forecast based on Woldemeskel et al. (2018)")
    var.createAttribute("units", "m3 s-1")

    ncout.createAttribute("Description", "Created probabilistic forecast given the output of mHM in Qsim_XXX and Qobs_xxx variable and the statistical model described in Woldemeskel et al. (2018). Lambda parameter of the boxcox transformation is set to the {:.2f}".format(lmbda))
    ncout.createAttribute("lambda", lmbda)
    ncout.createAttribute("script", "prob_forecast.py")
    ncout.createAttribute("paper", "Woldemeskel et al. Hydrol Earth Syst Sci, 2018 vol. 22 (12) pp. 6257-6278. 'Evaluating post-processing approaches for monthly and seasonal streamflow forecasts.', https://www.hydrol-earth-syst-sci.net/22/6257/2018/")
    ncout.createAttribute("created", asctime())
    ncout.createAttribute("eta_mean", "mean of boxcox transformed error: {:f}".format(eta_mean))
    ncout.createAttribute("eta_std", "standard deviation of boxcox transformed error: {:f}".format(eta_std))
    ncout.createAttribute("rho", "lag-1 auto-correlation: {:f}".format(rho))
    ncout.createAttribute("sigma_y", "lag_1 standard deviation: {:f}".format(sigma_y))
    ncout.close()


if __name__ == '__main__':

    discharge_file, out_file, lmbda, n_prob = parser()

    runoff_obs, runoff_sim = read_data(discharge_file)

    # remove missing values
    mask = np.logical_or(runoff_obs < 0., runoff_sim < 0.)
    if np.sum(mask) > 0:
        print('***WARNING: missing values are removed from observed and simulated data')
    runoff_obs = runoff_obs[~mask]
    runoff_sim = runoff_sim[~mask]

    # using nomenclature of Woldemeskel et al. 2018
    z_f, eta_mean, eta_std, rho, sigma_y, n_sample = calculate_param(runoff_obs, runoff_sim)

    # sample forecast
    runoff_prob = sample_forecasts(n_prob, z_f, eta_mean, eta_std, rho, sigma_y, n_sample)

    # add missing values
    if np.sum(mask) > 0.:
        tmp = np.zeros(mask.shape + (n_prob,)) - 9999.
        tmp[~mask] = runoff_prob
        runoff_prob = tmp

    # write to file
    write_file(runoff_prob, eta_mean, eta_std, rho, sigma_y, discharge_file, out_file)
    
    print('Done!')
