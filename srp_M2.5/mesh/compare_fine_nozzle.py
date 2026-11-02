"""meshFine vs meshFineNozzle: nozzle views, cell-size maps and cell counts by region.
Usage (compute/desktop node): python compare_fine_nozzle.py <out.png>"""
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


meshes = {"fine (meshFine)": load("srp_M2.5_meshFine.msh"), "fine + nozzle (meshFineNozzle)": load("srp_M2.5_meshFineNozzle.msh")}
regions = {"nozzle + plenum x -8.2-0.5, r<3": lambda c: (c[:, 0] >= -8.2) & (c[:, 0] < 0.5) & (c[:, 1] < 3),
           "jet box x 0.5-20, r<9.5": lambda c: (c[:, 0] >= 0.5) & (c[:, 0] < 20) & (c[:, 1] < 9.5)}
print(f"{'region':36s}{'fine':>10s}{'fineNozzle':>12s}{'ratio':>8s}{'median mm fine':>16s}{'fineNozzle':>12s}")
for lab, f in regions.items():
    v = [f(c).sum() for _, c, _ in meshes.values()]; md = [np.median(s[f(c)]) for _, c, s in meshes.values()]
    print(f"{lab:36s}{v[0]:10d}{v[1]:12d}{v[1]/v[0]:8.2f}{md[0]:16.4f}{md[1]:12.4f}")
tot = [len(s) for _, _, s in meshes.values()]
print(f"{'total':36s}{tot[0]:10d}{tot[1]:12d}{tot[1]/tot[0]:8.2f}")
fig, ax = plt.subplots(3, 2, figsize=(16, 13), constrained_layout=True)
views = [((-8.6, 2.6), (0, 3.4), "nozzle + plenum"), ((-1.0, 2.6), (1.4, 3.2), "throat-to-exit and lip (detail)")]
for j, (name, (polys, c, s)) in enumerate(meshes.items()):
    for i, (xl, yl, title) in enumerate(views):
        a = ax[i, j]
        m = (c[:, 0] > xl[0] - 0.5) & (c[:, 0] < xl[1] + 0.5) & (c[:, 1] > yl[0] - 0.5) & (c[:, 1] < yl[1] + 0.5)
        a.add_collection(PolyCollection([polys[k] for k in np.flatnonzero(m)], facecolors="none", edgecolors="k",
                                        linewidths=0.06))
        a.add_patch(Rectangle((-8.2, 0), 8.7, 3, fill=False, ec="r", lw=1.2, ls="--"))
        a.set_xlim(xl); a.set_ylim(yl); a.set_aspect("equal"); a.set_xlabel("x [mm]"); a.set_ylabel("r [mm]")
        a.set_title(f"{name}: {title}  ({len(s):,} cells)", fontsize=10)
    a = ax[2, j]
    m = (c[:, 0] > -9) & (c[:, 0] < 4) & (c[:, 1] < 4)
    sc = a.scatter(c[m, 0], c[m, 1], c=s[m], s=0.3, cmap="viridis", vmin=0.0, vmax=0.06, rasterized=True)
    a.add_patch(Rectangle((-8.2, 0), 8.7, 3, fill=False, ec="r", lw=1.2, ls="--"))
    a.set_xlim(-9, 4); a.set_ylim(0, 4); a.set_aspect("equal"); a.set_xlabel("x [mm]"); a.set_ylabel("r [mm]")
    a.set_title(f"{name}: cell size (sqrt area) [mm]", fontsize=10)
    fig.colorbar(sc, ax=a, shrink=0.8)
fig.suptitle("srp_M2.5: meshFine vs meshFineNozzle (red box: nozzle + plenum x -8.2..0.5 mm, r 0..3 mm, 0.05 -> 0.0177 mm)")
fig.savefig(sys.argv[1], dpi=150)
print("saved", sys.argv[1])
