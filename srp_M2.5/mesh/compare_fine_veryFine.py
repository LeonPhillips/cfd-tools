"""meshFine vs meshVeryFine: jet-box views, cell-size maps and cell counts by region.
Usage (compute/desktop node): python compare_fine_veryFine.py <out.png>"""
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


meshes = {"fine (meshFine)": load("srp_M2.5_meshFine.msh"), "very fine (meshVeryFine)": load("srp_M2.5_meshVeryFine.msh")}
regions = {"nozzle + plenum x -8.2-0.5, r<3": lambda c: (c[:, 0] >= -8.2) & (c[:, 0] < 0.5) & (c[:, 1] < 3),
           "jet box x 0.5-20, r<9.5": lambda c: (c[:, 0] >= 0.5) & (c[:, 0] < 20) & (c[:, 1] < 9.5)}
print(f"{'region':36s}{'fine':>10s}{'veryFine':>12s}{'ratio':>8s}{'median mm fine':>16s}{'veryFine':>12s}")
for lab, f in regions.items():
    v = [f(c).sum() for _, c, _ in meshes.values()]; md = [np.median(s[f(c)]) for _, c, s in meshes.values()]
    print(f"{lab:36s}{v[0]:10d}{v[1]:12d}{v[1]/v[0]:8.2f}{md[0]:16.4f}{md[1]:12.4f}")
tot = [len(s) for _, _, s in meshes.values()]
print(f"{'total':36s}{tot[0]:10d}{tot[1]:12d}{tot[1]/tot[0]:8.2f}")
fig, ax = plt.subplots(3, 2, figsize=(16, 13), constrained_layout=True)
views = [((-1.5, 24), (0, 12), "jet, Mach disk and bow shock"), ((0.0, 3.5), (1.5, 4.0), "nozzle lip (detail)")]
for j, (name, (polys, c, s)) in enumerate(meshes.items()):
    for i, (xl, yl, title) in enumerate(views):
        a = ax[i, j]
        m = (c[:, 0] > xl[0] - 0.5) & (c[:, 0] < xl[1] + 0.5) & (c[:, 1] > yl[0] - 0.5) & (c[:, 1] < yl[1] + 0.5)
        a.add_collection(PolyCollection([polys[k] for k in np.flatnonzero(m)], facecolors="none", edgecolors="k",
                                        linewidths=0.05 if i == 0 else 0.08))
        a.add_patch(Rectangle((0.5, 0), 19.5, 9.5, fill=False, ec="r", lw=1.2, ls="--"))
        a.set_xlim(xl); a.set_ylim(yl); a.set_aspect("equal"); a.set_xlabel("x [mm]"); a.set_ylabel("r [mm]")
        a.set_title(f"{name}: {title}  ({len(s):,} cells)", fontsize=10)
    a = ax[2, j]
    m = (c[:, 0] > -3) & (c[:, 0] < 35) & (c[:, 1] < 15)
    sc = a.scatter(c[m, 0], c[m, 1], c=s[m], s=0.3, cmap="viridis", vmin=0.0, vmax=0.15, rasterized=True)
    a.add_patch(Rectangle((0.5, 0), 19.5, 9.5, fill=False, ec="r", lw=1.2, ls="--"))
    a.set_xlim(-3, 35); a.set_ylim(0, 15); a.set_aspect("equal"); a.set_xlabel("x [mm]"); a.set_ylabel("r [mm]")
    a.set_title(f"{name}: cell size (sqrt area) [mm]", fontsize=10)
    fig.colorbar(sc, ax=a, shrink=0.8)
fig.suptitle("srp_M2.5: meshFine vs meshVeryFine (red box: x 0.5-20 mm, r 0-9.5 mm at 0.0177 -> 0.0125 mm)")
fig.savefig(sys.argv[1], dpi=150)
print("saved", sys.argv[1])
