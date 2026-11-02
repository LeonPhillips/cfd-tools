"""List where the worst-shaped and most abruptly-sized cells are. Usage: python find_bad_cells.py mesh.msh"""
import sys
from collections import defaultdict
import numpy as np
import gmsh

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.open(sys.argv[1])
ntags, coords, _ = gmsh.model.mesh.getNodes()
P = dict(zip(ntags, coords.reshape(-1, 3)[:, :2]))
cells = []
for etype, etags, enodes in zip(*gmsh.model.mesh.getElements(2)):
    n = gmsh.model.mesh.getElementProperties(etype)[3]
    cells += [(t, list(e)) for t, e in zip(etags, enodes.reshape(-1, n))]
gmsh.finalize()


def area(nodes):
    x, y = np.array([P[t] for t in nodes]).T
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def min_angle(nodes):
    p = np.array([P[t] for t in nodes]); n = len(p); a = 180.0
    for i in range(n):
        u, v = p[i - 1] - p[i], p[(i + 1) % n] - p[i]
        c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        a = min(a, np.degrees(np.arccos(np.clip(c, -1, 1))))
    return a


A = {t: area(nd) for t, nd in cells}
edge_cells = defaultdict(list)
for t, nd in cells:
    for i in range(len(nd)):
        edge_cells[tuple(sorted((nd[i], nd[(i + 1) % len(nd)])))].append(t)
ratio = defaultdict(float)
for cs in edge_cells.values():
    if len(cs) == 2:
        r = max(A[cs[0]], A[cs[1]]) / min(A[cs[0]], A[cs[1]])
        for c in cs:
            ratio[c] = max(ratio[c], r)

cent = {t: np.mean([P[k] for k in nd], axis=0) * 1000 for t, nd in cells}
ang = {t: min_angle(nd) for t, nd in cells}
print(f"cells: {len(cells)}  min angle: {min(ang.values()):.3f} deg  "
      f"<10 deg: {sum(a < 10 for a in ang.values())}  <20 deg: {sum(a < 20 for a in ang.values())}")
print(f"max neighbour area ratio: {max(ratio.values()):.1f}  >10: {sum(r > 10 for r in ratio.values())}  "
      f">100: {sum(r > 100 for r in ratio.values())}")
print("worst angles (deg) at x,y [mm]:")
for t in sorted(ang, key=ang.get)[:12]:
    print(f"  {ang[t]:7.3f}  ({cent[t][0]:9.4f}, {cent[t][1]:8.4f})  nverts={len(dict(cells)[t]) if False else ''}")
print("worst neighbour area ratios at x,y [mm]:")
for t in sorted(ratio, key=ratio.get, reverse=True)[:12]:
    print(f"  {ratio[t]:10.1f}  ({cent[t][0]:9.4f}, {cent[t][1]:8.4f})")
