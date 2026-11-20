# The mesoscale Hydrological Model — mHM <!--! {#mhm_readme} -->

- The latest mHM release can be found [here][0].
- The changelog can be found in the file [RELEASES][3].
- General information can be found on the [mHM website](https://mhm-ufz.org/).
- The mHM comes with a [LICENSE][6] agreement, this includes also the GNU Lesser General Public License.
- There is a list of [publications using mHM][7].

**Please note:** The [GitLab repository](https://git.ufz.de/mhm/mhm) grants read access to the code.
If you like to contribute to the code, please contact
<a href="mailto:mhm-admin@ufz.de">the mHM administrators</a>.

<!--! [TOC] -->

## Quick start

Use CMake 3.18 or newer, C and Fortran compilers, and NetCDF support.
CMake uses a local `forces/` checkout when present, otherwise it downloads
FORCES. From the repository root, build and run an example domain:

```bash
cmake -S . -B release -DCMAKE_BUILD_TYPE=Release
cmake --build release --parallel 4
./release/app/mhm examples/domain_01
```

The run uses `mhm.nml`, `parameters.nml`, and `outputs.nml` within the supplied
directory. Use `./release/app/mhm --help` for CLI options. Installation,
additional run configurations, and validation are described in the
[build and run guide](doc/v6-getting-started.md); dependency installation is
covered in [Compilation](doc/INSTALL.md).

## Documentation

The online documentation includes HTML and PDF versions:

- [Latest v6 development documentation](https://mhm.readthedocs.io/latest/) (`main`)
- [Stable release documentation](https://mhm.readthedocs.io/stable/)
- [Legacy v5 development documentation](https://mhm.readthedocs.io/v5-dev/) (`v5.x`)

| Guide | Contents |
| --- | --- |
| [Build and run mHM v6](doc/v6-getting-started.md) | Requirements, configuration, example runs, and layered output. |
| [Developing mHM v6](doc/v6-development-guide.md) | Component architecture, extension workflows, and coverage reports. |
| [Documentation workflow](doc/documentation-guide.md) | Building and publishing the HTML/PDF manual. |
| [Python package](pybind/README.md) | Forward execution, CLI, and example downloads. |
| [Namelist reference](doc/nml/namelists.md) | Configuration, process parameters, and output settings. |
| [Parameter migration](doc/parameter-migration-v6.md) | Converting legacy parameter namelists to v6. |
| [Model overview](doc/model-overview.md) | Scientific formulation, regionalization, and calibration. |
| [mRM overview](doc/mrm-overview.md) | Routing equations and parameterizations. |

Repository workflows and coding conventions are recorded in [AGENTS.md](AGENTS.md).

## Cite as

Please refer to the main model by citing Samaniego et al. (2010) and Kumar et al. (2013):

> Samaniego L., R. Kumar, S. Attinger (2010): Multiscale parameter regionalization of a grid-based hydrologic model at the mesoscale. Water Resour. Res., 46,W05523, doi:10.1029/2008WR007327, http://onlinelibrary.wiley.com/doi/10.1029/2008WR007327/abstract

> Kumar, R., L. Samaniego, and S. Attinger (2013): Implications of distributed hydrologic model parameterization on water fluxes at multiple scales and locations, Water Resour. Res., 49, doi:10.1029/2012WR012195, http://onlinelibrary.wiley.com/doi/10.1029/2012WR012195/abstract

The model code can be generally cited as:

> **mHM:** Luis Samaniego et al., mesoscale Hydrologic Model. Zenodo. doi:10.5281/zenodo.1069202, https://doi.org/10.5281/zenodo.1069202

To cite a certain version, have a look at the [Zenodo site][10].

## License overview

LGPLv3 (c) 2005-2026 mHM-Developers. See [LICENSE.md](LICENSE.md).

[0]: https://git.ufz.de/mhm/mhm/-/releases
[3]: doc/RELEASES.md
[6]: LICENSE.md
[7]: https://mhm-ufz.org/about/publications/
[10]: https://zenodo.org/record/3239055
