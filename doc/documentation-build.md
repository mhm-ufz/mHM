# Building documentation and coverage {#documentation_build}

These commands apply to the legacy `v5.x` branch. Run them from the repository root.

## HTML and PDF manual

Use Doxygen 1.18.0, Graphviz, and a Fortran compiler. The PDF also needs a recent
TeX Live installation with `latexmk`, `makeindex`, `newunicodechar`, `tabularray`,
and `xltabular`. The custom LaTeX header and table stylesheet preserve the manual's
cover, scientific chapters, and API reference.

```bash
export FC=gfortran
export PROJECT_NUMBER="$(cat version.txt)"
doxygen doc/doxygen.config
latexmk -norc -pdf -cd -interaction=nonstopmode -halt-on-error -file-line-error latex/refman.tex
```

Open `html/index.html` or `latex/refman.pdf`. Each build documents the checked-out
commit once; it does not check out a release tag or create stable/latest copies.

The [Read the Docs v5-dev manual](https://mhm.readthedocs.io/v5-dev/) tracks
`v5.x`; `latest` tracks v6 on `main`, and `stable` selects a release.
Read the Docs uses `.readthedocs.yaml`, `doc/rtd-environment.yml`, and
`doc/readthedocs_build.sh`. The script exposes `environment`, `doxygen`, `pdf`,
and `publish` phases, or runs all phases without an argument. HTML and the separate
PDF download are staged under `${READTHEDOCS_OUTPUT:-_readthedocs}`.

RTD builds use `doc/doxygen-rtd.config` to enable doxygen-awesome's hosted search
integration. The ordinary Doxygen configuration keeps the local search index for
GitLab artifacts and offline use. On RTD, one dropdown offers documentation
versions and a PDF download for the current version. It appears in the navigation
sidebar, or inside the navigation menu on phones, and uses RTD's Addons data.
Disable the floating flyout in RTD under Settings > Addons > Flyout menu.
The dropdown receives version and download data through the Addons event even
when the flyout is disabled. Integrated search and the dropdown follow the
documentation's light/dark theme.

GitLab loads `GCC/13.3.0`, `texlive/20250308`, and `chs-fortran-conda/0.6.0`.
It builds the PDF with Tectonic, runs `makeindex`, then runs Tectonic again to
include the index and resolve cross-references. Read the Docs continues to use
`latexmk`. The GitLab Doxygen step sets `CAIRO_DEBUG_PDF=1` to generate Graphviz
PDFs without Cairo compression, avoiding Tectonic 0.15's PDF import failures.
The documentation job retains `html/`, `mhm_doc.pdf`, and Doxygen logs as job
artifacts. This legacy branch has no GitLab Pages deployment, so it cannot replace
the main branch site.

## Fortran coverage

Coverage requires GNU Fortran, NetCDF, serial pFUnit, LCOV, `lcov_cobertura`, and
the check-case Python dependencies (`numpy`, `pandas`, `xarray`, `netCDF4`, `pexpect`).
GitLab uses `chs-fortran-conda/0.6.0`, including LCOV 2.5. GitHub builds serial
pFUnit 4.12.0 and publishes the resulting report to Coveralls.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DBUILD_MHM_DRIVER=ON \
  -DBUILD_TESTING=ON -DCMAKE_WITH_COVERAGE=ON -DCMAKE_WITH_OpenMP=OFF
cmake --build build --parallel 4
cmake --build build --target mhm_coverage_CI
```

Set `CMAKE_PREFIX_PATH` to your pFUnit installation if CMake cannot find it.
The coverage target runs the pFUnit tests and the legacy numerical check cases.
It keeps the check script's existing `case_11` exclusion and stages writable cases
under `build/coverage-runs/`, leaving source inputs and reference outputs intact.
Missing pFUnit is an error for coverage builds.

Open `build/mhm_coverage_CI/index.html`. CI normalizes checkout prefixes with LCOV
before converting `coverage/coverage.info` to `coverage/coverage.xml` with
`lcov_cobertura`; no custom Python collection script is needed. GitLab exposes the
HTML report and Cobertura data as job artifacts using the instance's default
retention. GitHub retains its coverage artifacts for 14 days.
