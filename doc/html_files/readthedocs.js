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

      const label = document.createElement("label");
      label.htmlFor = id;
      label.textContent = "Version";
      const select = document.createElement("select");
      select.id = id;
      for (const version of versions) {
        const option = document.createElement("option");
        option.textContent = version.verbose_name || version.slug;
        option.value = version.urls.documentation;
        option.selected = version.slug === current.slug;
        select.appendChild(option);
      }
      select.addEventListener("change", () => {
        // Like RTD's flyout, keep the current page when switching versions.
        const filename = config.readthedocs?.resolver?.filename ||
          document.querySelector('meta[name="readthedocs-resolver-filename"]')?.content || "";
        const target = new URL(filename.replace(/^\/+/, ""), select.value);
        target.hash = window.location.hash;
        window.location.assign(target.href);
      });
      controls.append(label, select);

      if (current.downloads?.pdf) {
        const pdf = document.createElement("a");
        pdf.href = current.downloads.pdf;
        pdf.textContent = "Download PDF";
        controls.appendChild(pdf);
      }
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

    // Hide the floating widget only after both replacements are available.
    document.documentElement.classList.add("mhm-rtd-navigation-ready");
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
