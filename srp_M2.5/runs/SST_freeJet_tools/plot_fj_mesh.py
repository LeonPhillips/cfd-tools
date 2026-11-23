"""python plot_fj_mesh.py <mesh.msh> <boundaries.vtk> <out.png>
Free-jet mesh views and the boundary-marker map (user, 2026-11-23): whole domain with markers coloured, the plume,
the near field, the nozzle and lip, and the lip corner. Run on a compute node (fine mesh ~1M cells).
"""
import sys
import numpy as np
import gmsh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.open(sys.argv[1])
tags, coords, _ = gmsh.model.mesh.getNodes()
pos = {t: i for i, t in enumerate(tags)}
xy = coords.reshape(-1, 3)[:, :2] * 1000
segs = []
for etype, _, enodes in zip(*gmsh.model.mesh.getElements(2)):
    n = gmsh.model.mesh.getElementProperties(etype)[3]
    idx = np.vectorize(pos.get)(enodes.reshape(-1, n))
    for i in range(n):
        segs.append(np.stack([xy[idx[:, i]], xy[idx[:, (i + 1) % n]]], 1))
gmsh.finalize()
segs = np.concatenate(segs)

# boundary markers from the legacy VTK written by mesh_freejet.py
txt = open(sys.argv[2]).read().split("\n")
names = dict(kv.split("=") for kv in txt[1].split(": ")[1].split())
i = 4; npts = int(txt[i].split()[1]); pts = np.array([list(map(float, l.split()[:2])) for l in txt[i + 1:i + 1 + npts]]) * 1000
i += 1 + npts; ncell = int(txt[i].split()[1])
cells = np.array([list(map(int, l.split()[1:3])) for l in txt[i + 1:i + 1 + ncell]])
ids = np.array(list(map(int, txt[i + 1 + ncell + 1 + ncell + 3:i + 1 + ncell + 1 + ncell + 3 + ncell])))

views = [("Whole domain, boundary markers", (-320, 250, -2, 125)), ("Plume, x -5..110 mm", (-5, 110, 0, 20)),
         ("Near field", (0, 20, 0, 6)), ("Nozzle and lip", (-8.5, 4, 0, 3.4)), ("Lip corner", (1.2, 2.6, 2.0, 3.0)),
         ("Body and base", (-65, 5, 0, 15))]
fig, axes = plt.subplots(3, 2, figsize=(15, 14))
for ax, (title, (x0, x1, y0, y1)) in zip(axes.flat, views):
    m = ((segs[:, :, 0].max(1) > x0) & (segs[:, :, 0].min(1) < x1) & (segs[:, :, 1].max(1) > y0) & (segs[:, :, 1].min(1) < y1))
    ax.add_collection(LineCollection(segs[m], lw=0.15 if "Whole" in title else 0.25, colors="0.4"))
    if "markers" in title:
        for k, name in names.items():
            sel = cells[ids == int(k)]
            ax.add_collection(LineCollection(pts[sel], lw=3, colors=f"C{int(k)}", label=name))
        ax.legend(loc="upper left", fontsize=8)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    if "Whole" not in title:
        ax.set_aspect("equal")
    ax.set_title(title); ax.set_xlabel("x [mm]"); ax.set_ylabel("r [mm]")
fig.suptitle(sys.argv[1])
fig.tight_layout()
fig.savefig(sys.argv[3], dpi=110)
print("wrote", sys.argv[3])
