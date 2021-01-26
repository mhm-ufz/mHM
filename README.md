# The mesoscale Hydrological Model -- mHM

- The current release is **[mHM v5.10][1]**.
- The latest mHM release notes can be found in the file [RELEASES][3] or [online][4].
- General information can be found on the [mHM website](http://www.ufz.de/mhm/).
- The mHM comes with a [LICENSE][6] agreement, this includes also the GNU Lesser General Public License.
- There is a list of [publications using mHM][7].

**Please note:** The [GitLab repository](https://git.ufz.de/mhm/mhm) grants read access to the code.
If you like to contribute to the code, please contact [mhm-admin@ufz.de](mailto:mhm-admin@ufz.de).

## Documentation

The online documentation for mHM can be found here (pdf versions are provided there as well):
- stable: https://mhm.pages.ufz.de/mhm
- latest: https://mhm.pages.ufz.de/mhm/latest

## Cite as

Please refer to the main model by citing Samaniego et al. (2010) and Kumar et al. (2013):

> Samaniego L., R. Kumar, S. Attinger (2010): Multiscale parameter regionalization of a grid-based hydrologic model at the mesoscale. Water Resour. Res., 46,W05523, doi:10.1029/2008WR007327, http://onlinelibrary.wiley.com/doi/10.1029/2008WR007327/abstract

> Kumar, R., L. Samaniego, and S. Attinger (2013): Implications of distributed hydrologic model parameterization on water fluxes at multiple scales and locations, Water Resour. Res., 49, doi:10.1029/2012WR012195, http://onlinelibrary.wiley.com/doi/10.1029/2012WR012195/abstract

The model code can be generally cited as:

> **mHM:** Luis Samaniego et al., mesoscale Hydrologic Model. Zenodo. doi:10.5281/zenodo.1069202, https://doi.org/10.5281/zenodo.1069202

To cite a certain Version, have a look at the [Zenodo site][10].

## Install

mHM can be compiled with cmake. See more details under [cmake manual][9].
Please see the file [DEPENDENCIES][8] for external software required to run mHM.
See also the [documentation][5] for detailed instructions to setup mHM.


## Quick start

1. Compile mHM
2. Run mHM on the test domains with the command `./mhm`, which uses settings from [mhm.nml](mhm.nml).
3. Explore the results in the [output directory](test_domain/), e.g. by using the NetCDF viewer `ncview`.


[1]: https://git.ufz.de/mhm/mhm/tree/5.10
[3]: doc/RELEASES.md
[4]: https://git.ufz.de/mhm/mhm/tags/
[5]: https://mhm.pages.ufz.de/mhm
[6]: LICENSE
[7]: doc/mhm_papers.md
[8]: doc/DEPENDENCIES.md
[9]: doc/INSTALL.md
[10]: https://zenodo.org/record/3239055
