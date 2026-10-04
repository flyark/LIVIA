// LIVIA's build: one string every page shows at its foot and every generated script, report and settings file carries, so an
// export can be traced to the page that made it. Change BUILD whenever the pages change (the date of the change).
(function () {
  window.LIVIA_BUILD = '2026-10-03';
  function stamp() {
    if (document.getElementById('livia-build')) return;
    const d = document.createElement('div'); d.id = 'livia-build';
    d.style.cssText = 'max-width:1100px;margin:24px auto 18px;padding:0 20px;font:12px/1.5 system-ui,sans-serif;color:#5F6771;text-align:center';
    d.innerHTML = 'LIVIA build ' + window.LIVIA_BUILD + ' · <a href="https://github.com/flyark/LIVIA/commits/main" target="_blank" rel="noopener" style="color:#2471A3">changes</a> · Kim &amp; Perrimon (2026), <a href="https://doi.org/10.64898/2026.05.01.721633" target="_blank" rel="noopener" style="color:#2471A3">bioRxiv</a>';
    document.body.appendChild(d);
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', stamp); else stamp();
})();
