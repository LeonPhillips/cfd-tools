import sys
import numpy as np
import gmsh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

step, out = sys.argv[1], sys.argv[2]
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setString("Geometry.OCCTargetUnit", "M")
gmsh.model.occ.importShapes(step)
gmsh.model.occ.synchronize()

# Proposed marker for each curve (from the boundary-condition screenshots)
marker = {4: "axis", 9: "axis", 5: "inlet_freestream", 3: "inlet_jet", 6: "wall_slip", 7: "wall_slip", 8: "outlet"}
for c in [1, 2] + list(range(10, 20)):  # v4 geometry (su2_shortDiff_wFillet_v4.STEP)
    marker[c] = "wall_noslip"
colour = {"axis": "tab:gray", "inlet_freestream": "tab:blue", "inlet_jet": "tab:cyan",
          "wall_slip": "tab:green", "wall_noslip": "tab:red", "outlet": "tab:purple"}

curves = {}
for _, tag in gmsh.model.getEntities(1):
    lo, hi = gmsh.model.getParametrizationBounds(1, tag)
    t = np.linspace(lo[0], hi[0], 200)
    xyz = np.array(gmsh.model.getValue(1, tag, t)).reshape(-1, 3)
    curves[tag] = xyz * 1000  # plot in mm
gmsh.finalize()

views = [("Whole domain", None), ("Model region", (-65, 80, -2, 50)),
         ("Nozzle and lip", (-9, 3, -0.3, 3.4))]
fig, axes = plt.subplots(3, 1, figsize=(11, 13))
for ax, (title, lim) in zip(axes, views):
    for tag, xyz in curves.items():
        m = marker[tag]
        ax.plot(xyz[:, 0], xyz[:, 1], color=colour[m], lw=2)
        mid = xyz[len(xyz) // 2]
        if lim is None or (lim[0] < mid[0] < lim[1] and lim[2] < mid[1] < lim[3]):
            ax.annotate(str(tag), mid[:2], fontsize=8, ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec=colour[m]))
    if lim:
        ax.set_xlim(lim[0], lim[1]); ax.set_ylim(lim[2], lim[3])
    ax.set_aspect("equal" if lim else "auto")
    ax.set_title(title); ax.set_xlabel("x [mm]"); ax.set_ylabel("r [mm]")
    ax.grid(alpha=0.3)
handles = [plt.Line2D([], [], color=c, lw=3, label=m) for m, c in colour.items()]
axes[0].legend(handles=handles, loc="upper right", fontsize=8)
fig.suptitle("srp_M2.5 geometry: curve numbers and proposed SU2 markers")
fig.tight_layout()
fig.savefig(out, dpi=110)
print("saved", out)
