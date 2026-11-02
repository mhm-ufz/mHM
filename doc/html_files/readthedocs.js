// RTD's Doxygen configuration generates serverUrl in search/search.js.
// Initialize before Doxygen constructs its SearchBox. Local/GitLab builds
// have no serverUrl and retain their built-in search unchanged.
if (typeof serverUrl !== "undefined") {
  DoxygenAwesomeReadtheDocsSearch.init("leftAlign");
}

(() => {
  function addNavigation(config) {
    const sidebar = document.getElementById("nav-tree");
    const menu = document.getElementById("main-menu");
    const current = config.versions?.current;
    if (!sidebar || !menu || !current || document.getElementById("mhm-rtd-sidebar")) return;

    const versions = [current, ...config.versions.active.filter(v => v.slug !== current.slug)];

    function createControls(id) {
      const controls = document.createElement("nav");
      controls.className = "mhm-rtd-controls";
      controls.setAttribute("aria-label", "Documentation versions and downloads");

      const select = document.createElement("select");
      select.id = id;
      select.setAttribute("aria-label", "Documentation version or PDF download");
      const versionGroup = document.createElement("optgroup");
      versionGroup.label = "Versions";
      for (const version of versions) {
        const option = document.createElement("option");
        option.textContent = version.verbose_name || version.slug;
        option.value = version.urls.documentation;
        option.selected = version.slug === current.slug;
        versionGroup.appendChild(option);
      }
      select.appendChild(versionGroup);

      if (current.downloads?.pdf) {
        const downloads = document.createElement("optgroup");
        downloads.label = "Downloads";
        const pdf = document.createElement("option");
        pdf.value = current.downloads.pdf;
        pdf.dataset.download = "pdf";
        pdf.textContent = "Download PDF (current version)";
        downloads.appendChild(pdf);
        select.appendChild(downloads);
      }

      select.addEventListener("change", () => {
        if (select.selectedOptions[0].dataset.download === "pdf") {
          select.value = current.urls.documentation;
          window.location.assign(current.downloads.pdf);
          return;
        }
        // Keep the current page when switching versions.
        const filename = config.readthedocs?.resolver?.filename ||
          document.querySelector('meta[name="readthedocs-resolver-filename"]')?.content || "";
        const target = new URL(filename.replace(/^\/+/, ""), select.value);
        target.hash = window.location.hash;
        window.location.assign(target.href);
      });
      controls.appendChild(select);
      return controls;
    }

    const desktop = createControls("mhm-rtd-version-desktop");
    desktop.id = "mhm-rtd-sidebar";
    const mobile = document.createElement("li");
    mobile.id = "mhm-rtd-mobile";
    mobile.appendChild(createControls("mhm-rtd-version-mobile"));
    sidebar.prepend(desktop);
    // Keep Doxygen's search item last: the theme styles it specially.
    menu.insertBefore(mobile, document.getElementById("searchBoxPos2"));
  }

  document.addEventListener("readthedocs-addons-data-ready", event => {
    const config = event.detail.data();
    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", () => addNavigation(config), { once: true });
    } else {
      addNavigation(config);
    }
  });
})();
