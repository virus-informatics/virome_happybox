(function () {
  /* __TEMPLATE_VARS_START__ */
  const divId = "";
  const cdsMap = {};
  const autoVis = {};
  const allIndices = [];
  const cdsIds = [];
  const cbClass = "";
  /* __TEMPLATE_VARS_END__ */

  window["restyleLabels_" + divId] = function (mode) {
    const vis = allIndices.map(function (i) {
      if (mode === "all") return true;
      if (mode === "off") return "legendonly";
      return autoVis[String(i)] ? true : "legendonly";
    });
    Plotly.restyle(divId, { visible: vis }, allIndices);
    document.querySelectorAll("." + cbClass).forEach(function (cb) {
      const idx = parseInt(cb.dataset.idx);
      const cdsId = cdsIds[idx];
      const indices = cdsMap[cdsId] || [];
      if (mode === "all") cb.checked = true;
      else if (mode === "off") cb.checked = false;
      else cb.checked = indices.some(function (i) { return autoVis[String(i)]; });
    });
  };

  window["toggleCds_" + divId] = function (cb) {
    const idx = parseInt(cb.dataset.idx);
    const cdsId = cdsIds[idx];
    const indices = cdsMap[cdsId] || [];
    const vis = indices.map(function () { return cb.checked ? true : "legendonly"; });
    Plotly.restyle(divId, { visible: vis }, indices);
  };
})();
