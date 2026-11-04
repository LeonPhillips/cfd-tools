"""python plot_conv_lines.py <stage folder> [previous stage folder]: <stage>/conv_lines.png (user, 2026-11-04).
Six panels, one per convergence line: Mach along the line for the latest snapshot and the one before (the two kept
snapshots), on a single Mach axis. The mean / p99 / max |dMach| per 1000 iterations of each line are printed in its
title and logged in line_change.csv (conv_update.py).
Lines: r = 0.01, 1.3, 2.3, 3, 5 mm (x 1-20 mm) and the jet shear layer (2.0, 2.34) -> (7.84, 4.70) mm (snap.py).
No pass threshold is set: report only.
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

folder = sys.argv[1].rstrip("/")
prev = sys.argv[2].rstrip("/") if len(sys.argv) > 2 and sys.argv[2] != "-" else None
key = lambda f: int(f.rsplit("_", 1)[-1][:-4])
snaps = sorted(glob.glob(f"{folder}/snaps/flow_*.npz"), key=key)
if len(snaps) < 2 and prev:
    old = sorted(glob.glob(f"{prev}/snaps/flow_*.npz"), key=key)[-1:]
    if old and "line_Mach" in np.load(old[0]).files:
        snaps = old + snaps
if not snaps:
    sys.exit("no snapshots yet")
zb = np.load(snaps[-1])
names, L = [str(s) for s in zb["line_names"]], zb["lines_mm"]
za = np.load(snaps[-2]) if len(snaps) > 1 else None
lc = f"{folder}/line_change.csv"
last = None
if os.path.exists(lc) and os.path.getsize(lc) > 80:
    tab = np.atleast_2d(np.genfromtxt(lc, delimiter=",", skip_header=1))
    if tab.size:
        last = tab[-1]
cols = plt.cm.viridis(np.linspace(0, 0.9, len(names)))

fig, axs = plt.subplots(2, 3, figsize=(18, 9))
for i, (a, n, c) in enumerate(zip(axs.flat, names, cols)):
    x0, r0, x1, r1 = L[i]
    xs = np.linspace(x0, x1, zb["line_Mach"].shape[1])
    a.plot(xs, zb["line_Mach"][i], color=c, lw=1.3, label=f"iteration {key(snaps[-1])}")
    if za is not None:
        lab = f"{'prev. stage ' if prev and snaps[-2].startswith(prev) else ''}iteration {key(snaps[-2])}"
        a.plot(xs, za["line_Mach"][i], color="0.4", lw=1.0, ls="--", label=lab)
    t = f"{n}: (x, r) ({x0:g}, {r0:g}) -> ({x1:g}, {r1:g}) mm"
    if last is not None:
        t += f"\n|dMach| per 1000 it: mean {last[2 + 3 * i]:.2e}, p99 {last[3 + 3 * i]:.2e}, max {last[4 + 3 * i]:.2e}"
    a.set_title(t, fontsize=9)
    a.set_xlabel("x [mm]"); a.set_ylabel("Mach"); a.grid(alpha=0.3); a.legend(fontsize=8, loc="best")

run = ([q for q in folder.split("/") if "_mesh" in q] or ["run"])[-1]
fig.suptitle(f"{run} / {folder.split('/')[-1]}: Mach along the convergence lines, latest snapshot vs the one before"
             "   (no pass threshold set; report only)", fontsize=11)
fig.tight_layout()
fig.savefig(f"{folder}/conv_lines.png", dpi=100)
print("saved", f"{folder}/conv_lines.png")
