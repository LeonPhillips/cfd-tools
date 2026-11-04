"""python plot_grid_features.py <out.png> "<title>" <label>=<stage folder> ...   (coarsest first, finest last)
Grid-study version of plot_conv_features.py (user, 2026-11-04): the same 2x2 (Mach disk x, Mach disk r, bow shock x,
stagnation x, all on the axis except the Mach disk r), with the mesh on the x-axis instead of the iteration. Each value
is the last row of that stage's features.csv (the live monitor, conv_update.py), i.e. the final snapshot.
Definitions as in grid_compare/feature_positions.py (0.025 mm box, sub-cell fit): Mach disk = the normal jet shock at
~8.6 mm, r = its subsonic-behind edge; stagnation = first axis u = 0 crossing. x from the CFD origin (nozzle exit
1.875 mm). All features.csv columns of the chosen rows are written to <out>.csv. Report only.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

out, title = sys.argv[1], sys.argv[2]
runs = [a.split("=", 1) for a in sys.argv[3:]]
lab, rows, head = [], [], None
for k, f in runs:
    lines = [l.rstrip("\n") for l in open(os.path.join(f, "features.csv")) if l.strip()]
    h = lines[0].split(",")
    if head is None:
        head = h
    elif h != head:
        sys.exit(f"{k}: features.csv columns differ")
    if len(lines) < 2:
        sys.exit(f"{k}: no features yet")
    lab.append(k); rows.append(lines[-1].split(",", len(head) - 1))
val = {c: np.array([float(r[i]) if r[i] not in ("", "nan") else np.nan for r in rows]) for i, c in enumerate(head[:-1])}
with open(os.path.splitext(out)[0] + ".csv", "w") as fh:
    fh.write("mesh," + ",".join(head) + "\n")
    for k, r in zip(lab, rows):
        fh.write(k + "," + ",".join(r) + "\n")

Q = [("mach_disk_x", "Mach disk x on the axis (mm)", "Mach disk x"),
     ("mach_disk_subsonic_edge_r", "Mach disk radius r (mm)", "Mach disk r"),
     ("bow_x_axis", "Bow shock x on the axis (mm)", "Bow shock x"),
     ("stagnation_x", "Stagnation point x on the axis (mm)", "Stagnation x")]
xs = np.arange(len(lab))
fig, axs = plt.subplots(2, 2, figsize=(12, 8.5), constrained_layout=True)
for a, (k, ylab, t) in zip(axs.flat, Q):
    v = val[k]
    a.plot(xs, v, "o-", color="#1f77b4", ms=6)
    for x, y in zip(xs, v):
        a.annotate(f"{y:.3f}", (x, y), (0, 7), textcoords="offset points", ha="center", fontsize=8)
    a.set_title(f"{t}: {lab[-1]} {v[-1]:.3f} mm, {lab[-2]} {v[-2] - v[-1]:+.3f} mm from it" if len(v) > 1 else t)
    a.set_xticks(xs, lab); a.set_xlabel("mesh"); a.set_ylabel(ylab); a.grid(alpha=0.3)
    a.margins(y=0.15)
its = ", ".join(f"{k} {int(val['iter'][j])}" for j, k in enumerate(lab))
fig.suptitle(f"{title}: feature positions vs mesh (x from CFD origin; nozzle exit x = 1.875 mm)\n"
             f"final snapshot iterations: {its}", fontsize=10)
fig.savefig(out, dpi=110)
print("saved", out)
for k, j in zip(lab, range(len(lab))):
    print(f"  {k:>10}: " + "  ".join(f"{c} {val[c][j]:.3f}" for c, _, _ in Q))
