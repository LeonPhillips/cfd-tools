"""python plot_conv_lines.py <stage folder> [previous stage folder]: <stage>/conv_lines.png (user, 2026-11-04).
Six panels, one per convergence line, on a single Mach axis: Mach along the line for every snapshot of the stage
(line_history.npz, conv_update.py), coloured by iteration on one colour scale shared by all panels (light = early,
dark = late); the latest snapshot is black and thicker; the previous stage's last snapshot (the restart) is grey dashed.
The mean / p99 / max |dMach| per 1000 iterations against the snapshot before are printed in each title (line_change.csv).
Lines: r = 0.01, 1.3, 2.3, 3, 5 mm (x 1-20 mm) and the jet shear layer (2.0, 2.34) -> (7.84, 4.70) mm (snap.py).
No pass threshold is set: report only.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

folder = sys.argv[1].rstrip("/")
prev = sys.argv[2].rstrip("/") if len(sys.argv) > 2 and sys.argv[2] != "-" else None
hf = f"{folder}/line_history.npz"
if not os.path.exists(hf):
    sys.exit("no line history yet")
h = np.load(hf)
its, LM = h["iter"], h["line_Mach"]
names, L = [str(s) for s in h["line_names"]], h["lines_mm"]
restart = None
if prev and os.path.exists(f"{prev}/line_history.npz"):
    with np.load(f"{prev}/line_history.npz") as hp:
        if hp["line_Mach"].shape[1:] == LM.shape[1:]:
            restart = (int(hp["iter"][-1]), hp["line_Mach"][-1])
lc = f"{folder}/line_change.csv"
last = None
if os.path.exists(lc) and len(open(lc).readlines()) > 1:
    tab = np.atleast_2d(np.genfromtxt(lc, delimiter=",", skip_header=1))
    if tab.size:
        last = tab[-1]

cmap = plt.cm.viridis_r
norm = Normalize(vmin=its[0], vmax=max(its[-1], its[0] + 1))
fig, axs = plt.subplots(2, 3, figsize=(19, 9.5), constrained_layout=True)
for i, (a, n) in enumerate(zip(axs.flat, names)):
    x0, r0, x1, r1 = L[i]
    xs = np.linspace(x0, x1, LM.shape[2])
    if restart is not None:
        a.plot(xs, restart[1][i], color="0.55", lw=1.0, ls="--")
    for k in range(len(its) - 1):
        a.plot(xs, LM[k, i], color=cmap(norm(its[k])), lw=0.8)
    a.plot(xs, LM[-1, i], color="k", lw=1.6)
    t = f"{n}: (x, r) ({x0:g}, {r0:g}) -> ({x1:g}, {r1:g}) mm"
    if last is not None:
        t += f"\nlatest |dMach| per 1000 it: mean {last[2 + 3 * i]:.2e}, p99 {last[3 + 3 * i]:.2e}, max {last[4 + 3 * i]:.2e}"
    a.set_title(t, fontsize=9)
    a.set_xlabel("x [mm]"); a.set_ylabel("Mach"); a.grid(alpha=0.3)

sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
fig.colorbar(sm, ax=axs, shrink=0.85, pad=0.01, label="snapshot iteration")
hand = [Line2D([], [], color="k", lw=1.6, label=f"latest (iteration {int(its[-1])})"),
        Line2D([], [], color=cmap(0.5), lw=0.8, label=f"earlier snapshots, every 1000 it ({len(its) - 1})")]
if restart is not None:
    hand.append(Line2D([], [], color="0.55", lw=1.0, ls="--", label=f"previous stage, last snapshot ({restart[0]})"))
axs.flat[0].legend(handles=hand, fontsize=8, loc="best")
run = ([q for q in folder.split("/") if "_mesh" in q] or ["run"])[-1]
fig.suptitle(f"{run} / {folder.split('/')[-1]}: Mach along the convergence lines, all {len(its)} snapshots"
             "   (no pass threshold set; report only)", fontsize=11)
fig.savefig(f"{folder}/conv_lines.png", dpi=100)
print("saved", f"{folder}/conv_lines.png")
