"""python diff_html.py <box.npz> <minuend> <subtrahend> <out.html> "<title>": interactive Mach difference map
(<minuend> - <subtrahend>) over the refinement box, in the same format as figures/2026-10-29_srp_SA1_diff-contours-Mach_
fine-vs-veryFine.html: 0.02 mm sampling, colour clipped at +/-0.05, hover shows the true difference and the minuend Mach,
dashed standard lines, wheel zoom and right/middle-drag pan. <box.npz> comes from resample_box.py (0.01 mm)."""
import json
import sys
import numpy as np
z = np.load(sys.argv[1]); a, b, out, title = sys.argv[2:6]
x, r = z["x_mm"][::2], z["r_mm"][::2]
ma = z[f"{a}_Mach"][::2, ::2].astype(float)
d = ma - z[f"{b}_Mach"][::2, ::2].astype(float)
grid = lambda m, n: [[None if not np.isfinite(v) else round(float(v), n) for v in row] for row in m]
CLIP = 0.05
LINES = [0.01, 1.3, 2.3, 3.0, 6.0, 9.0]
trace = {"type": "heatmap", "x": [round(float(v), 4) for v in x], "y": [round(float(v), 4) for v in r], "z": grid(d, 4),
         "customdata": grid(ma, 3), "zmin": -CLIP, "zmax": CLIP, "colorscale": "RdBu", "reversescale": False,
         "colorbar": {"title": {"text": f"dMach<br>(clipped ±{CLIP:g})"}},
         "hovertemplate": f"x=%{{x:.2f}} mm<br>r=%{{y:.2f}} mm<br>dMach=%{{z:.4f}}<br>{a} Mach=%{{customdata:.3f}}<extra></extra>"}
m = np.isfinite(d)
stats = (f"mean |dM| {np.abs(d[m]).mean():.4f}, p99 {np.percentile(np.abs(d[m]), 99):.4f}, "
         f"max {np.abs(d[m]).max():.3f}")
layout = {"title": {"text": f"{title}<br><sup>{stats}. Wheel zooms, left drag box-zooms, right or middle drag pans, "
                            "double-click resets; colour clipped, hover shows the true value; dashed: standard lines "
                            "r = 0.01, 1.3, 2.3, 3, 6, 9 mm</sup>"},
          "xaxis": {"title": "x [mm] (freestream from the right)", "range": [float(x[0]), float(x[-1])], "constrain": "domain"},
          "yaxis": {"title": "r [mm]", "range": [float(r[0]), float(r[-1])], "scaleanchor": "x", "constrain": "domain"},
          "shapes": [{"type": "line", "xref": "x", "yref": "y", "x0": 1, "x1": 20, "y0": rl, "y1": rl,
                      "line": {"color": "black", "width": 0.6, "dash": "dash"}, "opacity": 0.35} for rl in LINES],
          "height": 700, "template": "plotly_white", "margin": {"l": 60, "r": 40, "t": 90, "b": 60}}
PAN = """<script>/* Mouse add-ons: wheel zooms (scrollZoom); right- or middle-button drag pans the subplot under the cursor
   (left drag still box-zooms, toolbar unchanged). */
(function () {
  var gd = document.getElementById("p"), st = null, raf = 0, last = null;
  gd.addEventListener("contextmenu", function (e) { e.preventDefault(); });
  gd.addEventListener("mousedown", function (e) {
    if (e.button !== 1 && e.button !== 2) return;
    var fl = gd._fullLayout, bb = gd.getBoundingClientRect(), px = e.clientX - bb.left, py = e.clientY - bb.top, hit = null;
    Object.keys(fl._plots).forEach(function (k) {
      var p = fl._plots[k], xa = p.xaxis, ya = p.yaxis;
      if (!hit && px >= xa._offset && px <= xa._offset + xa._length && py >= ya._offset && py <= ya._offset + ya._length) hit = xa;
    });
    if (!hit) return;
    var axes = [hit];
    Object.keys(fl._plots).forEach(function (k) {
      var ya = fl._plots[k].yaxis;
      if (fl._plots[k].xaxis === hit && axes.indexOf(ya) < 0) axes.push(ya);
    });
    st = { x0: e.clientX, y0: e.clientY, axes: axes.map(function (a) {
      return { name: a._name, isX: a._id.charAt(0) === "x", r: a.range.slice(), len: a._length };
    }) };
    e.preventDefault(); e.stopPropagation();
  }, true);
  window.addEventListener("mousemove", function (e) {
    if (!st) return;
    last = e;
    if (raf) return;
    raf = requestAnimationFrame(function () {
      raf = 0;
      var dx = last.clientX - st.x0, dy = last.clientY - st.y0, upd = {};
      st.axes.forEach(function (a) {
        var span = a.r[1] - a.r[0], d = (a.isX ? -dx : dy) * span / a.len;
        upd[a.name + ".range"] = [a.r[0] + d, a.r[1] + d];
      });
      Plotly.relayout(gd, upd);
    });
  });
  window.addEventListener("mouseup", function () { st = null; });
})();
</script>"""
html = (f'<!doctype html><html><head><meta charset="utf-8"><title>dMach {a} vs {b}</title>\n'
        '<script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>\n'
        '<style>body{margin:0;font-family:sans-serif;background:#fff}</style></head><body><div id="p" style="width:100%"></div>\n'
        f'<script>Plotly.newPlot("p", [{json.dumps(trace, separators=(",", ":"))}], {json.dumps(layout)}, '
        '{responsive:true, displaylogo:false, scrollZoom:true,\n toImageButtonOptions:{format:"png", scale:2}});</script>'
        + PAN + "\n</body></html>\n")
open(out, "w", encoding="utf-8").write(html)
print("saved", out, f"{len(html) / 1e6:.1f} MB;", stats)
