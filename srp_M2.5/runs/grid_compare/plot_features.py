"""python plot_features.py <feature_positions_all.csv> <machdisk_fine_sst.csv> <out_prefix>
Mach disk x and r, bow shock x and stagnation x (all on the axis except r) against mesh for the SST cases
(user, 2026-11-04). Writes <out_prefix>.png (matplotlib, 2x2) and <out_prefix>.xlsx (data, native Excel charts, notes).
Mach disk x and r come from the 0.005 mm box (machdisk_fine.py); r is the subsonic-behind edge. Bow shock and
stagnation come from the 0.025 mm box (feature_positions.py). Mesh axis: ordered categories.
"""
import math
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpyxl import Workbook
from openpyxl.chart import LineChart, Reference
from openpyxl.chart.layout import Layout, ManualLayout
from openpyxl.drawing.image import Image
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter

MESHES = ["Coarse", "Medium", "Fine", "FineJet", "FineJet2", "VeryFine"]
SERIES = [("SST", "p0j7.4945", "SST, p0j 7.4945 bar (total conditions)", "#1f77b4", "o", "-"),
          ("SST", "massflow_rhoj8.681", "SST, MASS_FLOW (rho 8.681, V 99.11)", "#2ca02c", "s", "--"),
          ("SST", "p0j8.015", "SST, p0j 8.015 bar (old inlet)", "#ff7f0e", "^", "-"),
          ("SSTsarkar", "p0j8.015", "SST-Sarkar, p0j 8.015 bar (old inlet)", "#d62728", "v", "-")]
QUANT = [("mach_disk_x_fine", "Mach disk x on the axis (mm)", "Mach disk x"),
         ("mach_disk_subsonic_edge_r_fine", "Mach disk radius r (mm)", "Mach disk r"),
         ("bow_x_axis", "Bow shock x on the axis (mm)", "Bow shock x"),
         ("stagnation_x", "Stagnation point x on the axis (mm)", "Stagnation x")]

f = pd.read_csv(sys.argv[1], comment="#")
m = pd.read_csv(sys.argv[2], comment="#")
d = f.merge(m, on=["model", "inlet", "mesh", "stage"], how="inner")
d = d[d.model.isin(["SST", "SSTsarkar"])].copy()
d["mesh"] = pd.Categorical(d.mesh, MESHES, ordered=True)
d = d.sort_values(["model", "inlet", "mesh"])
out = sys.argv[3]


def table(col):
    t = pd.DataFrame(index=MESHES)
    for mo, inl, lab, *_ in SERIES:
        s = d[(d.model == mo) & (d.inlet == inl)].set_index("mesh")[col]
        t[lab] = [s.get(k, float("nan")) for k in MESHES]
    return t


# ---- PNG ----
fig, axs = plt.subplots(2, 2, figsize=(12, 8.5), constrained_layout=True)
for ax, (col, ylab, title) in zip(axs.flat, QUANT):
    t = table(col)
    for (mo, inl, lab, c, mk, ls) in SERIES:
        y = t[lab]; ok = y.notna()
        ax.plot([i for i, v in enumerate(ok) if v], y[ok], color=c, marker=mk, ls=ls, label=lab)
    ax.set_xticks(range(len(MESHES)), MESHES)
    ax.set_ylabel(ylab); ax.set_title(title); ax.grid(alpha=0.3)
axs[0, 0].legend(fontsize=8)
fig.suptitle("srp_M2.5 SST final solutions: feature positions vs mesh (x from CFD origin; nozzle exit x = 1.875 mm)")
fig.savefig(out + ".png", dpi=150)
print("saved", out + ".png")

# ---- Excel ----
AR = Font(name="Arial", size=10); ARB = Font(name="Arial", size=10, bold=True); HF = PatternFill("solid", fgColor="DDEBF7")
wb = Workbook()
ws = wb.active; ws.title = "Charts"
pd_ws = wb.create_sheet("Plot data")
row = 1
anchors = ["A1", "K1", "A24", "K24"]
for (col, ylab, title), anc in zip(QUANT, anchors):
    t = table(col)
    pd_ws.cell(row, 1, title).font = ARB
    hdr = row + 1
    pd_ws.cell(hdr, 1, "Mesh").font = ARB
    for j, lab in enumerate(t.columns, 2):
        c = pd_ws.cell(hdr, j, lab); c.font = ARB; c.fill = HF; c.alignment = Alignment(wrap_text=True)
    for i, mesh in enumerate(MESHES):
        pd_ws.cell(hdr + 1 + i, 1, mesh).font = AR
        for j, lab in enumerate(t.columns, 2):
            v = t.loc[mesh, lab]
            c = pd_ws.cell(hdr + 1 + i, j, None if pd.isna(v) else round(float(v), 4))
            c.font = AR; c.number_format = "0.000"
    ch = LineChart(); ch.title = title; ch.y_axis.title = ylab; ch.x_axis.title = "Mesh"
    ch.add_data(Reference(pd_ws, min_col=2, max_col=1 + len(SERIES), min_row=hdr, max_row=hdr + len(MESHES)),
                titles_from_data=True)
    ch.set_categories(Reference(pd_ws, min_col=1, min_row=hdr + 1, max_row=hdr + len(MESHES)))
    ch.display_blanks = "span"
    for s, (*_, colr, mk, ls) in zip(ch.series, SERIES):
        s.graphicalProperties.line.solidFill = colr[1:]
        s.graphicalProperties.line.width = 22000
        if ls == "--":
            s.graphicalProperties.line.dashStyle = "dash"
        s.marker.symbol = {"o": "circle", "s": "square", "^": "triangle", "v": "diamond"}[mk]
        s.marker.size = 7
        s.marker.graphicalProperties.solidFill = colr[1:]
        s.marker.graphicalProperties.line.solidFill = colr[1:]
        s.smooth = False
    lo, hi = float(t.min().min()), float(t.max().max())
    step = next(s for s in (0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.5, 1.0) if (hi - lo) / s <= 7)
    ch.y_axis.scaling.min = math.floor(lo / step + 1e-9) * step
    ch.y_axis.scaling.max = math.ceil(hi / step + 1e-9) * step
    ch.y_axis.majorUnit = step
    ch.y_axis.number_format = "0.000" if col.startswith("mach_disk") else "0.00"
    ch.x_axis.delete = False; ch.y_axis.delete = False
    ch.x_axis.tickLblPos = "low"
    ch.legend.position = "b"; ch.legend.overlay = False
    ch.layout = Layout(manualLayout=ManualLayout(x=0.1, y=0.08, w=0.87, h=0.6, xMode="edge", yMode="edge"))
    ch.width, ch.height = 16, 11
    ws.add_chart(ch, anc)
    row = hdr + len(MESHES) + 3
pd_ws.column_dimensions["A"].width = 12
for j in range(2, 2 + len(SERIES)):
    pd_ws.column_dimensions[get_column_letter(j)].width = 22

img_ws = wb.create_sheet("Figure (PNG)")
img = Image(out + ".png"); img.width, img.height = img.width * 0.55, img.height * 0.55
img_ws.add_image(img, "A1")

dw = wb.create_sheet("All SST features")
cols = list(d.columns)
for j, c in enumerate(cols, 1):
    x = dw.cell(1, j, c); x.font = ARB; x.fill = HF
for i, rec in enumerate(d.itertuples(index=False), 2):
    for j, v in enumerate(rec, 1):
        c = dw.cell(i, j, None if (not isinstance(v, str) and pd.isna(v)) else (v if isinstance(v, str) else float(v)))
        c.font = AR
        if isinstance(v, float):
            c.number_format = "0.0000"
dw.freeze_panes = "E2"

nw = wb.create_sheet("Notes")
notes = [
    "srp_M2.5 feature positions vs mesh, SST cases (user request 2026-11-04). All values in mm; x from the CFD origin "
    "(nozzle exit at x = 1.875 mm), r radial. Jet flows +x, freestream -x.",
    "Source: final solutions (coarse t4_cap3 where it exists, else t3_cap3) on M3, "
    "runs/grid_compare/features/ (feature_positions_all.csv, machdisk_fine_sst.csv). Scripts: cfd-tools "
    "srp_M2.5/runs/grid_compare/ (feature_positions.py, machdisk_fine.py, plot_features.py).",
    "Mach disk x: largest positive d(rho)/dx peak on the axis, parabolic sub-cell fit, 0.005 mm box (x 4-5, r 0-0.5 mm).",
    "Mach disk r: r where the minimum Mach within 0.3 mm behind the tracked shock reaches 1 (end of the subsonic-behind "
    "part), same 0.005 mm box. 0 = no Mach disk (flow behind the shock supersonic on the axis; all coarse meshes). "
    "Sampling finer than the mesh only interpolates, so this is mesh-limited.",
    "Bow shock x: d(rho)/dx minimum on the axis near the point where Mach first drops below 97 % of the upstream value, "
    "0.025 mm box.",
    "Stagnation x: first u_x zero crossing (+ to -) on the axis after the nozzle exit, 0.025 mm box. Coarse meshes have "
    "a recirculation pocket with further crossings (column axis_u0_x in 'All SST features').",
    "Mesh axis: ordered categories; FineJet and FineJet2 are local refinements of Fine, VeryFine is a global refinement. "
    "MASS_FLOW exists on Coarse and FineJet only; p0j 7.4945 bar has no VeryFine.",
    "Convergence of each solution is judged separately (line profiles); see context/m3-environment.md.",
]
nw.column_dimensions["A"].width = 140
for i, s in enumerate(notes, 1):
    c = nw.cell(i, 1, s); c.font = AR; c.alignment = Alignment(wrap_text=True, vertical="top")
for sh in (ws, pd_ws, img_ws, dw, nw):
    sh.sheet_view.showGridLines = sh is not ws
wb.save(out + ".xlsx")
print("saved", out + ".xlsx")
