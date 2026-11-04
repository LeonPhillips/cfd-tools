"""python plot_conv_features.py <stage folder>: <stage>/conv_features.png from features.csv (conv_update.py; user,
2026-11-04). One row per snapshot (every 1000 iterations). Definitions as in grid_compare/feature_positions.py:
Mach disk = the normal jet shock at ~8.6 mm, r = its subsonic-behind edge; shock crossing = barrel shocks on the axis
(~4.5 mm); bow shock on the axis, at 1 D and 2 D (D = 4.6 mm); stagnation = first axis u = 0 crossing. x in mm from the
CFD origin (nozzle exit 1.875 mm). Box resolution 0.025 mm; positions use a sub-cell fit. Report only.
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
flags = rows[o[-1]][-1]

AX = [("mach_disk_x", "Mach disk x (axis)", "tab:blue"), ("shock_crossing_x", "shock crossing x (axis)", "tab:purple"),
      ("stagnation_x", "stagnation x (axis)", "tab:green"), ("bow_x_axis", "bow shock x (axis)", "tab:red"),
      ("bow_x_1D", "bow shock x (r = 1 D)", "tab:orange"), ("bow_x_2D", "bow shock x (r = 2 D)", "tab:brown")]
fig, ax = plt.subplots(2, 2, figsize=(17, 10))
for k, lab, c in AX:
    ax[0, 0].plot(it, col[k] - col[k][-1], "o-", ms=3, color=c, label=f"{lab}: now {col[k][-1]:.3f} mm")
ax[0, 0].axhline(0, color="k", lw=0.5)
ax[0, 0].set_ylabel("position - latest value [mm]")
ax[0, 0].set_title("Feature positions, each relative to its latest value (box cell = 0.025 mm)", fontsize=10)

ax[0, 1].plot(it, col["mach_disk_subsonic_edge_r"], "o-", ms=3, color="tab:blue",
              label=f"subsonic-behind edge (reported): now {col['mach_disk_subsonic_edge_r'][-1]:.3f} mm")
ax[0, 1].plot(it, col["mach_disk_edge_r"], "s--", ms=3, color="0.5",
              label=f"50 % of axis strength (fragile): now {col['mach_disk_edge_r'][-1]:.3f} mm")
ax[0, 1].set_ylabel("Mach disk radius r [mm]"); ax[0, 1].set_title("Mach disk radius (triple-point estimate)", fontsize=10)

if len(it) > 1:
    dit = np.diff(it)
    for k, lab, c in AX + [("mach_disk_subsonic_edge_r", "Mach disk r", "k")]:
        ax[1, 0].semilogy(it[1:], np.maximum(np.abs(np.diff(col[k])) * 1000 / dit, 1e-5), "o-", ms=3, color=c, label=lab)
    ax[1, 0].axhline(0.025, color="k", ls=":", lw=1, label="one box cell (0.025 mm)")
ax[1, 0].set_ylabel("|change| per 1000 iterations [mm] (floor 1e-5)")
ax[1, 0].set_title("Feature movement between consecutive snapshots", fontsize=10)

ax[1, 1].axis("off")
txt = ["latest snapshot: iteration %d" % it[-1], ""]
for k, lab, _ in AX + [("mach_disk_subsonic_edge_r", "Mach disk r", None)]:
    v = col[k]
    d1 = v[-1] - v[-2] if len(v) > 1 else np.nan
    d5 = v[-1] - v[-6] if len(v) > 5 else np.nan
    txt.append(f"{lab:<24s} {v[-1]:8.3f} mm   last step {d1:+.3f}   last 5 steps {d5:+.3f}")
txt += ["", f"flags: {flags or 'none'}", "", "No pass threshold set: report only (convergence.md)."]
ax[1, 1].text(0.0, 1.0, "\n".join(txt), va="top", family="monospace", fontsize=9, transform=ax[1, 1].transAxes)
for a in (ax[0, 0], ax[0, 1], ax[1, 0]):
    a.grid(alpha=0.3, which="both"); a.set_xlabel("iteration"); a.legend(fontsize=8)
run = ([q for q in folder.split("/") if "_mesh" in q] or ["run"])[-1]
fig.suptitle(f"{run} / {folder.split('/')[-1]}: key-feature positions every 1000 iterations", fontsize=11)
fig.tight_layout()
fig.savefig(f"{folder}/conv_features.png", dpi=100)
print("saved", f"{folder}/conv_features.png")
