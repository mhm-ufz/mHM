# The mesoscale Hydrological Model -- mHM

- The current release is **[mHM v5.8][1]** ([Download as zip][2]).
- The latest mHM release notes can be found in the file [RELEASES][3] or [online][4].
- All the details can be found in the [user manual][5].
- The mHM comes with a [LICENSE][6] agreement, this includes also the GNU Lesser General Public License.
- There is a list of [publications using mHM][7].

**Please note:** This Git repository grants read access to the code of the release versions only.
For the time being, please develop your code using the [SVN repository](https://svn.ufz.de/mhm/),
access can be granted through the [registration form](http://www.ufz.de/index.php?en=40116) on the [mHM website](http://www.ufz.de/mhm/).
In the course of 2018, we are planning to migrate our development repository from SVN to [GitLab](https://git.ufz.de/mhm).
The [repository on GitHub.com](https://github.com/mhm-ufz/mhm) will be always a mirror of the releases at [GitLab](https://git.ufz.de/mhm).

## Cite as

Please refer to the main model by citing Samaniego et al. (2010) and Kumar et al. (2013):

- Samaniego L., R. Kumar, S. Attinger (2010): Multiscale parameter regionalization of a grid-based hydrologic model at the mesoscale. Water Resour. Res., 46,W05523, doi:10.1029/2008WR007327, http://onlinelibrary.wiley.com/doi/10.1029/2008WR007327/abstract
- Kumar, R., L. Samaniego, and S. Attinger (2013): Implications of distributed hydrologic model parameterization on water fluxes at multiple scales and locations, Water Resour. Res., 49, doi:10.1029/2012WR012195, http://onlinelibrary.wiley.com/doi/10.1029/2012WR012195/abstract

The model code can be cited as:

- **mHM:** Luis Samaniego et al. (2017), mesoscale Hydrologic Model, doi:10.5281/zenodo.1069202, https://doi.org/10.5281/zenodo.1069202
- **mHM v5.8:** Luis Samaniego et al. (2017), mesoscale Hydrologic Model (Version 5.8), doi:10.5281/zenodo.1069203, https://doi.org/10.5281/zenodo.1069203


## Install

Please see the file [DEPENDENCIES][8] for external software required to run mHM.
See also the [users manual][5] for detailed instructions to setup mHM.


## Quick start

1. Compile mHM with the `make` command, which uses settings from [Makefile](Makefile).
2. Run mHM on a test basin with the command `./mhm`, which uses settings from [mhm.nml](mhm.nml).
3. Explore the results in the [output directory](test_basin/), e.g. by using the NetCDF viewer `ncview`.


[1]: https://git.ufz.de/mhm/mhm/tags/v5.8
[2]: https://git.ufz.de/mhm/mhm/repository/v5.8/archive.zip
[3]: doc/RELEASES.md
[4]: https://git.ufz.de/mhm/mhm/tags/
[5]: doc/mhm_manual_v5.8.pdf
[6]: LICENSE
[7]: doc/mhm_papers.md
[8]: doc/DEPENDENCIES.md
