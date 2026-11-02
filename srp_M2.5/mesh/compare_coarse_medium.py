"""Coarse (v4) vs meshMedium: zoomed mesh plots, cell-size maps and cell counts by region.
Usage (compute/desktop node): python compare_coarse_medium.py <out.png>"""
import sys
import numpy as np, gmsh, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Rectangle

def load(path):
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.open(path)
    nt, xyz, _ = gmsh.model.mesh.getNodes(); xyz = xyz.reshape(-1, 3)[:, :2] * 1e3
    idx = np.zeros(int(nt.max()) + 1, int); idx[nt.astype(int)] = np.arange(len(nt))
    polys = []
    for et, _, en in zip(*gmsh.model.mesh.getElements(2)):
        n = gmsh.model.mesh.getElementProperties(et)[3]
        polys += list(xyz[idx[en.astype(int).reshape(-1, n)]])
    gmsh.finalize()
    c = np.array([p.mean(0) for p in polys])
    a = np.array([0.5 * abs(np.dot(p[:, 0], np.roll(p[:, 1], -1)) - np.dot(p[:, 1], np.roll(p[:, 0], -1))) for p in polys])
    return polys, c, np.sqrt(a)

meshes = {"coarse (meshCoarse, v4)": load("srp_M2.5_coarse_v4.msh"), "medium (meshMedium)": load("srp_M2.5_meshMedium.msh")}
regions = {"refined box x 0.5-20, r 0-9.5 mm": (0.5, 20, 0, 9.5), "rest of plume box x 0-30, r 0-9 mm": None,
           "nozzle + plenum x<0.5, r<3": (-8.2, 0.5, 0, 3), "everything else": "rest"}
print(f"{'region':40s}" + "".join(f"{k.split()[0]:>12s}" for k in meshes) + "    ratio")
rows = {}
for name, (polys, c, s) in meshes.items():
    inbox = (c[:, 0] >= 0.5) & (c[:, 0] < 20) & (c[:, 1] < 9.5)
    plume = (c[:, 0] >= 0) & (c[:, 0] < 30) & (c[:, 1] < 9) & ~inbox
    noz = (c[:, 0] < 0.5) & (c[:, 0] >= -8.2) & (c[:, 1] < 3)
    rest = ~(inbox | plume | noz)
    rows[name] = [inbox.sum(), plume.sum(), noz.sum(), rest.sum(), len(s), np.median(s[inbox])]
labels = list(regions) + ["total", "median size in box (mm)"]
for i, lab in enumerate(labels):
    v = [rows[k][i] for k in meshes]
    if i < 5: print(f"{lab:40s}" + "".join(f"{x:12d}" for x in v) + f"{v[1]/v[0]:9.2f}")
    else: print(f"{lab:40s}" + "".join(f"{x:12.4f}" for x in v))

fig, ax = plt.subplots(3, 2, figsize=(15, 13), constrained_layout=True)
views = [((-1.5, 24), (0, 12), "jet, Mach disk and bow shock"), ((0.0, 3.5), (1.5, 4.0), "nozzle lip (detail)")]
for j, (name, (polys, c, s)) in enumerate(meshes.items()):
    for i, (xl, yl, title) in enumerate(views):
        a = ax[i, j]
        m = (c[:, 0] > xl[0] - 1) & (c[:, 0] < xl[1] + 1) & (c[:, 1] > yl[0] - 1) & (c[:, 1] < yl[1] + 1)
        a.add_collection(PolyCollection([polys[k] for k in np.flatnonzero(m)], facecolors="none", edgecolors="k", linewidths=0.08))
        a.add_patch(Rectangle((0.5, 0), 19.5, 9.5, fill=False, ec="r", lw=1.2, ls="--"))
        a.set_xlim(xl); a.set_ylim(yl); a.set_aspect("equal"); a.set_xlabel("x [mm]"); a.set_ylabel("r [mm]")
        a.set_title(f"{name}: {title}  ({len(s):,} cells total)", fontsize=10)
    a = ax[2, j]
    m = (c[:, 0] > -3) & (c[:, 0] < 35) & (c[:, 1] < 15)
    sc = a.scatter(c[m, 0], c[m, 1], c=s[m], s=0.3, cmap="viridis", vmin=0.0, vmax=0.15, rasterized=True)
    a.add_patch(Rectangle((0.5, 0), 19.5, 9.5, fill=False, ec="r", lw=1.2, ls="--"))
    a.set_xlim(-3, 35); a.set_ylim(0, 15); a.set_aspect("equal"); a.set_xlabel("x [mm]"); a.set_ylabel("r [mm]")
    a.set_title(f"{name}: cell size (sqrt area) [mm]", fontsize=10)
    fig.colorbar(sc, ax=a, shrink=0.8)
fig.suptitle("srp_M2.5: meshCoarse vs meshMedium (red box: 0.025 mm refinement, x 0.5-20 mm, r 0-9.5 mm)")
fig.savefig(sys.argv[1], dpi=170)
print("saved", sys.argv[1])
