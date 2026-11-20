# Building and publishing the documentation <!--! {#mhm_documentation} -->

Generate the HTML/PDF manual and maintain the GitLab and Read the Docs
publishing workflows. Run all commands from the repository root. For component
documentation and coverage reports, see the [developer guide](v6-development-guide.md).

<!--! [TOC] -->

## Generate HTML and PDF

```bash
PROJECT_NUMBER="$(cat version.txt)" doxygen doc/doxygen.config
make -C latex
```

Open `html/index.html` for the manual or `latex/refman.pdf` for its PDF.
Documentation builds need Doxygen 1.18.0, pinned in `doc/rtd-environment.yml`,
Graphviz, `awk`, and a Fortran compiler for `doc/pre-proc.sh`. Export `FC` to select
a compiler. The preprocessor otherwise uses `gfortran` when available.
System TeX Live supplies the PDF tools, including `makeindex`, and the
`newunicodechar`, `xltabular`, and `tabularray` packages. The Unicode mapping
supports pdfLaTeX and XeTeX, including Tectonic.

The extra stylesheet `doc/mhm_api_tables.sty` uses `longtable` for API parameter
tables and `xltabular` for simple documentation tables to reduce PDF build time.
Merged cells, captions, and labels use Doxygen's standard renderer. The manual
retains its source listings, bibliography, and indices.

The vendored theme is [doxygen-awesome v2.5.0](https://github.com/jothepro/doxygen-awesome-css/releases/tag/v2.5.0),
copied unchanged under `doc/doxygen-awesome-css/`. The HTML header displays the
project version from `version.txt`. Local and GitLab builds use Doxygen's local search.
Email obfuscation is disabled because Doxygen 1.18 otherwise emits malformed
`mailto:` links in the generated HTML.

## Maintain pages and navigation

Write standalone prose in Markdown. The [model overview](model-overview.md),
[mRM overview](mrm-overview.md), and [namelist index](nml/namelists.md) retain
Doxygen commands where needed for citations, PDF figures,
and subpage relationships. Use ordinary Markdown links so repository readers
can follow the guides. Add `<!--! [TOC] -->` to every handwritten Markdown page
for an automatic Doxygen table of contents, hidden in repository views.
Generated namelist references already include their own TOC commands.
`TOC_INCLUDE_HEADINGS=1` keeps the automatically included heading overview short.
Doxygen's `<!--! ... -->` comments hide
manual metadata in repository views. For format-specific figures, put the
repository's `<img>` element between hidden `\cond` and `\endcond` commands,
then keep the manual's HTML/PDF image commands in the same hidden comment.
For SVGs with theme variants, use an `<img>` and caption in a Doxygen `htmlonly` block
and copy both variants through `HTML_EXTRA_FILES`. This allows the theme toggle
to replace the displayed image; Doxygen's default SVG objects retain their content.
Set `width="500"` on the repository image to match the manual's size limit.
This preserves vector PDF figures and existing HTML styling. Namelist entries link
only the page title, leaving descriptions and literal namelist names outside the link.
A final relationship block records `\subpage` declarations, preserving the PDF
hierarchy. Its hidden HTML wrapper and LaTeX conditional suppress duplicate links
in both formats. The navigation layout lists the individual namelist
pages under their configuration, parameter, or output category; keep it aligned
with the overview when adding a namelist. Their explicit page URLs let the tree
select pages opened directly or at a particular field anchor.

In Markdown pages, write inline math as `$...$` on a single line.
Put display math between `$$` lines, with blank lines around the block. For numbered
equations or aligned rows, put `\begin{equation}` or `\begin{align}` (or its
starred form) directly after the opening `$$` and the matching `\end` directly
before the closing `$$`. GitLab and GitHub render the math from the Markdown.
Doxygen runs `doc/markdown-filter.awk` on Markdown inputs to translate the
delimiters and preserve `%` in single-line inline code such as
`exchange%config%processes`. Put code examples in fenced blocks and escape
literal dollar signs in prose as `\$`. Use Doxygen's math commands in Fortran API
comments and escape `%` there as `\%` in prose or `%%` inside backticks.

The cell illustration and MPR diagram use SVGs converted from their original
PDFs, preserving vector geometry and labels. In the dark MPR diagram, the
coloured layers are vector polygons fitted to the original raster geometry.
The hierarchy PDF contains a raster;
its SVG keeps that terrain image and replaces the blurry annotations with larger
vector labels using the original notation. Mathematical labels are outlined TeX
glyphs, so they require no browser fonts. Light variants include an opaque white
canvas so labels remain readable on dark repository pages when the renderer
ignores the `<picture>` theme selection. Each figure also has a `*-dark.svg`
variant with adapted neutral surfaces, a transparent canvas, and light annotations.
MPR variables on the coloured levels and the cell's tree label use black lettering.
Keep labels as vectors
when editing these assets so they stay sharp at the 500-pixel display width and
when zoomed. The original PNGs and PDF figures remain unchanged. SVG descriptions
record their source, adaptations, and CC BY 4.0 licensing.
Doxygen selects the variants through `html.dark-mode` CSS
rules and copies them with `HTML_EXTRA_FILES`. Repository views use `<picture>`
with `prefers-color-scheme`; their selection follows the browser's colour preference.
Keep the 500-pixel image width and test both the system preference and the
manual's explicit theme toggle when changing these figures.
Keep structural directory/group declarations in
`doc/directories.dox`; module headers declare their group with `\ingroup`.
The generated files under `doc/nml/` are maintained through nml-tools and their
schemas, as described in the [build and run guide](v6-getting-started.md).

The navigation starts with **mHM**, **Guides**, **Configuration**, **About**,
**Model components**, and **API reference**. The scientific model overview and
mHM logo form the landing page. Guides starts with the README's project information
and quick start. Configuration contains the namelist reference, its generated
subpages, and parameter migration guidance. Component and API indices remain
automatic, including Python modules. Maintain the menu in `doc/doxygen-layout.xml`
using stable page IDs. Groups use `url="[none]"` to expand or collapse without
opening a placeholder page. Local CSS in `doc/html_files/html_extra.css` hides
the redundant outer project row and removes its indentation while preserving
Doxygen's tree structure. It also gives named member groups theme-aware colors.
Check desktop navigation, the mobile menu, and both color themes after changes.

Give each standalone Markdown page an explicit first-heading ID inside a
Doxygen comment, for example `# Project information <!--! {#mhm_readme} -->`.
Doxygen reads the ID; GitHub and GitLab hide the comment. This controls its
generated filename,
`mhm_readme.html`. The main-page selection overrides the filename: the model
overview's `model_overview` ID resolves to `index.html` because
`USE_MDFILE_AS_MAINPAGE = doc/model-overview.md`. Use `@ref` with these IDs in the layout, and relative
Markdown links in prose so links also work in the repository.
`IMPLICIT_DIR_DOCS = NO` lets the Python README become a standalone page;
directory descriptions remain explicit in `doc/directories.dox`.
`MARKDOWN_ID_STYLE = GITHUB` controls automatic heading anchors such as
`mhm_readme.html#quick-start`; it does not choose page filenames. Preserve explicit
section anchors when converting existing documentation; `<a name="description"></a>`
preserves an old anchor in both renderers without adding visible heading markup.
Doxygen continues to name API pages and member anchors.

The scientific introduction remains at `index.html`, including existing section
links such as `index.html#description`. The README's project information and
quick start live at `mhm_readme.html`. The mRM and namelist pages retain
`m_r_m.html` and `nml.html`.
License, repository instructions, and Python guidance use `license.html`,
`repository_instructions.html`, and `mhm_python.html`. Update links to their
former automatically generated names; no compatibility redirects are generated.

## GitLab publishing

Each pipeline documents its checked-out commit once for HTML and LaTeX. The job
loads `GCC/13.3.0`, `texlive/20250308`, and `chs-fortran-conda/0.6.0`, taking its
Fortran compiler from the Conda module's `FC` setting. Graphviz rendering is
limited to four workers. CI sets `CAIRO_DEBUG_PDF=1` during generation to avoid
Cairo PDF object streams that can fail in Tectonic's reader.

PDF compilation uses Tectonic with `-Z continue-on-errors`, TeX Live's
`makeindex`, and a second Tectonic compilation to include the index.
Recoverable errors remain in the logs, and unreadable graph PDFs may be omitted.
Tectonic chooses its own reruns; `-r 10` would force ten additional passes rather
than set a maximum. The job retains `html/`, `mhm_doc.pdf`, `doxygen_log.txt`,
`doxygen_warn.txt`, `doxygen_latex.txt`, and `latex/refman.log` as artifacts,
including on failure. GitLab Pages publishes the `main` branch HTML at its root;
Read the Docs provides versioned documentation. Coverage reports remain
separate CI artifacts.

## Read the Docs

[Read the Docs](https://mhm.readthedocs.io/) provides versioned HTML/PDF manuals
and hosted search. Its navigation dropdown offers version selection and the
current version's PDF download on desktop and mobile, preserving the current
page and anchor when switching versions. Disable the floating flyout under
**Settings > Addons > Flyout menu**; the dropdown still receives its Addons data.

Environment creation, Doxygen, PDF compilation, and publishing run as separate
commands so the build log identifies the current phase:

```bash
bash doc/readthedocs_build.sh environment
bash doc/readthedocs_build.sh doxygen
bash doc/readthedocs_build.sh pdf
bash doc/readthedocs_build.sh publish
```

With no argument, `doc/readthedocs_build.sh` runs all phases, including creating
the Conda environment. Set `MHM_DOCS_ENV` to override its default name `mhm-docs`.
PDF generation uses system `latexmk -pdf` to schedule compiler, bibliography,
and index passes. The complete log is `doxygen_latex_make.log`; the last 80 lines
are printed on failure. The script stages HTML and `pdf/mhm.pdf` under
`$READTHEDOCS_OUTPUT` (or `_readthedocs` locally).
