"""Count cells by region and report typical size. Usage: python cell_census.py mesh.msh"""
import sys
import numpy as np
import gmsh

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.open(sys.argv[1])
ntags, coords, _ = gmsh.model.mesh.getNodes()
idx = {t: i for i, t in enumerate(ntags)}
xy = coords.reshape(-1, 3)[:, :2] * 1000
cen, size = [], []
for etype, _, enodes in zip(*gmsh.model.mesh.getElements(2)):
    n = gmsh.model.mesh.getElementProperties(etype)[3]
    pts = xy[[idx[t] for t in enodes]].reshape(-1, n, 2)
    cen.append(pts.mean(1))
    x, y = pts[..., 0], pts[..., 1]
    size.append(np.sqrt(0.5 * np.abs((x * np.roll(y, -1, 1) - y * np.roll(x, -1, 1)).sum(1))))
gmsh.finalize()
cen, size = np.vstack(cen), np.concatenate(size)
regions = [("nozzle/plenum x[-8.2,2] r<3", -8.2, 2, 0, 3), ("plume box x[2,30] r<9", 2, 30, 0, 9), ("jet box x[-10,75] r<46", -10, 75, 0, 46),
           ("body box x[-60,-10] r<20", -60, -10, 0, 20), ("wake x[-300,-60]", -300, -60, 0, 200),
           ("downstream x<-300", -2000, -300, 0, 200)]
print(f"total {len(size)}")
for name, x0, x1, y0, y1 in regions:
    m = (cen[:, 0] >= x0) & (cen[:, 0] < x1) & (cen[:, 1] >= y0) & (cen[:, 1] < y1)
    print(f"{name:32s} {m.sum():8d} cells  median size {np.median(size[m]) if m.any() else 0:.4f} mm")
