"""python plot_grid_lines.py <out.png> "<title>" <label>=<stage folder> ...   (coarsest first, finest last)
Grid-study version of plot_conv_lines.py (user, 2026-11-04): the same six panels (r = 0.01, 1.3, 2.3, 3, 5 mm and the
jet shear layer; Mach on one left axis), but one curve per mesh instead of per snapshot. Each curve is the last
snapshot in that stage's line_history.npz (written by the live monitor, conv_update.py), i.e. the final solution.
Light to dark = coarse to fine; the finest is black and thicker. Titles: mean / p99 / max |dMach| of the second-finest
mesh against the finest. The same table for every mesh is printed and written to <out>.csv. No pass threshold: report only.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

out, title = sys.argv[1], sys.argv[2]
runs = [a.split("=", 1) for a in sys.argv[3:]]
if len(runs) < 2:
    sys.exit("need at least two <label>=<stage folder>")
lab, LM, its, L, names = [], [], [], None, None
for k, f in runs:
    h = np.load(os.path.join(f, "line_history.npz"))
    if L is None:
        L, names = h["lines_mm"], [str(s) for s in h["line_names"]]
    elif not np.allclose(h["lines_mm"], L) or h["line_Mach"].shape[1:] != LM[0].shape:
        sys.exit(f"{k}: different lines from {lab[0]}")
    lab.append(k); LM.append(h["line_Mach"][-1]); its.append(int(h["iter"][-1]))
LM = np.array(LM)
n = len(lab)


def stats(a, b):
    d = np.abs(a - b); d = d[np.isfinite(d)]
    return (d.mean(), np.percentile(d, 99), d.max()) if d.size else (np.nan,) * 3


rows = [[lab[j]] + [v for i in range(len(names)) for v in stats(LM[j, i], LM[-1, i])] for j in range(n - 1)]
with open(os.path.splitext(out)[0] + ".csv", "w") as fh:
    fh.write("mesh_vs_" + lab[-1] + "," + ",".join(f"{s}_{q}" for q in names for s in ("mean", "p99", "max")) + "\n")
    for r in rows:
        fh.write(r[0] + "," + ",".join(f"{v:.4e}" for v in r[1:]) + "\n")
print(f"|dMach| against {lab[-1]} (mean / p99 / max):")
for r in rows:
    print(f"  {r[0]:>10}: " + "  ".join(f"{q} {r[1 + 3 * i]:.3f}/{r[2 + 3 * i]:.3f}/{r[3 + 3 * i]:.3f}" for i, q in enumerate(names)))

cols = [plt.cm.viridis_r(0.15 + 0.7 * j / max(n - 2, 1)) for j in range(n - 1)] + ["k"]
fig, axs = plt.subplots(2, 3, figsize=(19, 9.5), constrained_layout=True)
for i, (a, q) in enumerate(zip(axs.flat, names)):
    x0, r0, x1, r1 = L[i]
    xs = np.linspace(x0, x1, LM.shape[2])
    for j in range(n):
        a.plot(xs, LM[j, i], color=cols[j], lw=1.6 if j == n - 1 else 1.0)
    m, p, mx = stats(LM[-2, i], LM[-1, i])
    a.set_title(f"{q}: (x, r) ({x0:g}, {r0:g}) -> ({x1:g}, {r1:g}) mm\n"
                f"|dMach| {lab[-2]} vs {lab[-1]}: mean {m:.2e}, p99 {p:.2e}, max {mx:.2e}", fontsize=9)
    a.set_xlabel("x [mm]"); a.set_ylabel("Mach"); a.grid(alpha=0.3)
hand = [Line2D([], [], color=cols[j], lw=1.6 if j == n - 1 else 1.0, label=f"{lab[j]} (iteration {its[j]})") for j in range(n)]
axs.flat[0].legend(handles=hand, fontsize=8, loc="best")
fig.suptitle(f"{title}: Mach along the convergence lines, one curve per mesh (no pass threshold set; report only)", fontsize=11)
fig.savefig(out, dpi=100)
print("saved", out)
