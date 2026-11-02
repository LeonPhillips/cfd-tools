"""Plot zoomed views of a 2D Gmsh mesh. Usage: python plot_mesh.py mesh.msh out.png"""
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
xy = dict(zip(tags, coords.reshape(-1, 3)[:, :2] * 1000))  # mm
segs = []
for etype, _, enodes in zip(*gmsh.model.mesh.getElements(2)):
    n = gmsh.model.mesh.getElementProperties(etype)[3]
    for el in enodes.reshape(-1, n):
        p = [xy[t] for t in el]
        segs += [(p[i], p[(i + 1) % n]) for i in range(n)]
gmsh.finalize()
segs = np.array(segs)

views = [("Whole domain", (-1200, 80, 0, 130)), ("Model region", (-65, 78, 0, 47)),
         ("Nozzle and lip", (-8.5, 2.5, 0, 3.2)), ("Lip corner", (0.3, 2.0, 2.2, 3.1)),
         ("Shoulder fillet (edge 5)", (0.2, 0.6, 8.75, 9.1)), ("Wake / base", (-80, -40, 0, 20)),
         ("Base/axis corner (stretched)", (-59.8, -57.9, 0, 0.12)),
         ("Plenum wall/jet-inlet corner", (-8.2, -7.9, 2.45, 2.9))]
fig, axes = plt.subplots(4, 2, figsize=(14, 17))
for ax, (title, (x0, x1, y0, y1)) in zip(axes.flat, views):
    m = ((segs[:, :, 0].max(1) > x0) & (segs[:, :, 0].min(1) < x1) &
         (segs[:, :, 1].max(1) > y0) & (segs[:, :, 1].min(1) < y1))
    ax.add_collection(LineCollection(segs[m], lw=0.25, colors="k"))
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    if title not in ("Whole domain", "Base/axis corner (stretched)"):
        ax.set_aspect("equal")
    ax.set_title(title); ax.set_xlabel("x [mm]"); ax.set_ylabel("r [mm]")
fig.suptitle(sys.argv[1])
fig.tight_layout()
fig.savefig(sys.argv[2], dpi=110)
print("saved", sys.argv[2])
