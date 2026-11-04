"""python plot_conv_features.py <stage folder>: <stage>/conv_features.png from features.csv (conv_update.py; user,
2026-11-04). A 2x2 figure like figures/srp_M2.5/comparisons/feature-positions/..._vs-mesh.png, against iteration:
Mach disk x (axis), Mach disk r, bow shock x (axis), stagnation x (axis). One point per snapshot (every 1000 iterations).
Definitions as in grid_compare/feature_positions.py: Mach disk = the normal jet shock at ~8.6 mm, r = its subsonic-behind
edge; stagnation = first axis u = 0 crossing. x in mm from the CFD origin (nozzle exit 1.875 mm). 0.025 mm box with a
sub-cell fit. Report only.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

folder = sys.argv[1].rstrip("/")
f = f"{folder}/features.csv"
if not os.path.exists(f) or os.path.getsize(f) < 60:
    sys.exit("no features yet")
head = open(f).readline().strip().split(",")
rows = [l.rstrip("\n").split(",", len(head) - 1) for l in open(f).readlines()[1:] if l.strip()]
it = np.array([float(r[0]) for r in rows])
col = {h: np.array([float(r[i]) if r[i] not in ("", "nan") else np.nan for r in rows]) for i, h in enumerate(head[1:-1], 1)}
o = np.argsort(it); it = it[o]; col = {k: v[o] for k, v in col.items()}

Q = [("mach_disk_x", "Mach disk x on the axis (mm)", "Mach disk x"),
     ("mach_disk_subsonic_edge_r", "Mach disk radius r (mm)", "Mach disk r"),
     ("bow_x_axis", "Bow shock x on the axis (mm)", "Bow shock x"),
     ("stagnation_x", "Stagnation point x on the axis (mm)", "Stagnation x")]
fig, axs = plt.subplots(2, 2, figsize=(12, 8.5), constrained_layout=True)
for a, (k, ylab, title) in zip(axs.flat, Q):
    v = col[k]
    a.plot(it, v, "o-", color="#1f77b4", ms=5)
    d = f", last step {v[-1] - v[-2]:+.3f} mm" if len(v) > 1 else ""
    a.set_title(f"{title}: now {v[-1]:.3f} mm{d}")
    a.set_xlabel("iteration"); a.set_ylabel(ylab); a.grid(alpha=0.3)
run = ([q for q in folder.split("/") if "_mesh" in q] or ["run"])[-1]
fig.suptitle(f"{run} / {folder.split('/')[-1]}: feature positions vs iteration "
             "(x from CFD origin; nozzle exit x = 1.875 mm)")
fig.savefig(f"{folder}/conv_features.png", dpi=110)
print("saved", f"{folder}/conv_features.png")
