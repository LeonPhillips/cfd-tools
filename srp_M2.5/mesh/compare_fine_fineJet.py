"""meshFine vs meshFineJet vs meshVeryFine: jet-region views, cell-size maps and cell counts by region.
Usage (compute/desktop node): python compare_fine_fineJet.py <out.png>"""
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


meshes = {"fine (meshFine)": load("srp_M2.5_meshFine.msh"), "fine + jet box (meshFineJet)": load("srp_M2.5_meshFineJet.msh"),
          "very fine (meshVeryFine)": load("srp_M2.5_meshVeryFine.msh")}
JB = (1.5, 8.5, 0.0, 3.6)
regions = {"jet box x 1.5-8.5, r<3.6": lambda c: (c[:, 0] >= JB[0]) & (c[:, 0] < JB[1]) & (c[:, 1] < JB[3]),
           "rest of box x 0.5-20, r<9.5": lambda c: (c[:, 0] >= 0.5) & (c[:, 0] < 20) & (c[:, 1] < 9.5)
           & ~((c[:, 0] >= JB[0]) & (c[:, 0] < JB[1]) & (c[:, 1] < JB[3])),
           "nozzle + plenum x -8.2-0.5, r<3": lambda c: (c[:, 0] >= -8.2) & (c[:, 0] < 0.5) & (c[:, 1] < 3)}
names = list(meshes)
print(f"{'region':34s}" + "".join(f"{n.split(' (')[0]:>16s}" for n in names) + "   median size mm (fine / fineJet / veryFine)")
for lab, f in regions.items():
    v = [f(c).sum() for _, c, _ in meshes.values()]; md = [np.median(s[f(c)]) for _, c, s in meshes.values()]
    print(f"{lab:34s}" + "".join(f"{x:16d}" for x in v) + "   " + " / ".join(f"{m:.4f}" for m in md))
tot = [len(s) for _, _, s in meshes.values()]
print(f"{'total':34s}" + "".join(f"{x:16d}" for x in tot))
fig, ax = plt.subplots(3, 3, figsize=(24, 13), constrained_layout=True)
views = [((0.5, 10), (0, 4.5), "jet region"), ((1.2, 3.2), (1.8, 3.2), "nozzle lip and start of the shear layer (detail)")]
for j, (name, (polys, c, s)) in enumerate(meshes.items()):
    for i, (xl, yl, title) in enumerate(views):
        a = ax[i, j]
        m = (c[:, 0] > xl[0] - 0.3) & (c[:, 0] < xl[1] + 0.3) & (c[:, 1] > yl[0] - 0.3) & (c[:, 1] < yl[1] + 0.3)
        a.add_collection(PolyCollection([polys[k] for k in np.flatnonzero(m)], facecolors="none", edgecolors="k",
                                        linewidths=0.05 if i == 0 else 0.1))
        a.add_patch(Rectangle((JB[0], 0), JB[1] - JB[0], JB[3], fill=False, ec="r", lw=1.2, ls="--"))
        a.set_xlim(xl); a.set_ylim(yl); a.set_aspect("equal"); a.set_xlabel("x [mm]"); a.set_ylabel("r [mm]")
        a.set_title(f"{name}: {title}  ({len(s):,} cells)", fontsize=10)
    a = ax[2, j]
    m = (c[:, 0] > -1) & (c[:, 0] < 21) & (c[:, 1] < 10)
    sc = a.scatter(c[m, 0], c[m, 1], c=s[m], s=0.15, cmap="viridis", vmin=0.01, vmax=0.03, rasterized=True)
    a.add_patch(Rectangle((JB[0], 0), JB[1] - JB[0], JB[3], fill=False, ec="r", lw=1.2, ls="--"))
    a.set_xlim(-1, 21); a.set_ylim(0, 10); a.set_aspect("equal"); a.set_xlabel("x [mm]"); a.set_ylabel("r [mm]")
    a.set_title(f"{name}: cell size (sqrt area) [mm]", fontsize=10)
    fig.colorbar(sc, ax=a, shrink=0.8)
fig.suptitle("srp_M2.5: meshFine / meshFineJet / meshVeryFine (red box: x 1.5-8.5 mm, r 0-3.6 mm at 0.0125 mm in meshFineJet)")
fig.savefig(sys.argv[1], dpi=130)
print("saved", sys.argv[1])
