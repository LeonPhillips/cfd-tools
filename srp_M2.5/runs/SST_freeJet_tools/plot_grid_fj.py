"""python plot_grid_fj.py <out prefix> "<title>" <label>=<stage folder> ...   (coarsest first, finest last)
Free-jet grid study (user, 2026-11-24): the final solution of each mesh, from the live monitor's small files
(line_history.npz and features.csv of each mesh's t4_cap3; conv_fj.py). Writes
    <prefix>_lines.png      Mach on the five monitor lines, one curve per mesh (coarse light, finest black). Titles:
                            mean / p99 / max |dMach| of the second-finest mesh against the finest.
    <prefix>_lines.csv      that table for every mesh against the finest
    <prefix>_features.png   the x of the first four axis-Mach peaks and the supersonic core length, against mesh
    <prefix>_features.csv   the last features.csv row of each mesh
No pass threshold is set: report only.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

pre, title = sys.argv[1], sys.argv[2]
runs = [a.split("=", 1) for a in sys.argv[3:]]
if len(runs) < 2:
    sys.exit("need at least two <label>=<stage folder>")

# ---- lines -----------------------------------------------------------------------------------------------------------
lab, LM, its, L, names = [], [], [], None, None
for k, f in runs:
    h = np.load(os.path.join(f, "line_history.npz"), allow_pickle=True)
    if L is None:
        L, names = h["lines_mm"], [str(s) for s in h["line_names"]]
    elif not np.allclose(h["lines_mm"], L) or h["mach"].shape[1:] != LM[0].shape:
        sys.exit(f"{k}: different lines from {lab[0]}")
    lab.append(k); LM.append(h["mach"][-1]); its.append(int(h["iterations"][-1]))
LM = np.array(LM)
n = len(lab)


def stats(a, b):
    d = np.abs(a - b); d = d[np.isfinite(d)]
    return (d.mean(), np.percentile(d, 99), d.max()) if d.size else (np.nan,) * 3


rows = [[lab[j]] + [v for i in range(len(names)) for v in stats(LM[j, i], LM[-1, i])] for j in range(n - 1)]
with open(pre + "_lines.csv", "w") as fh:
    fh.write("mesh_vs_" + lab[-1] + "," + ",".join(f"{s}_{q}" for q in names for s in ("mean", "p99", "max")) + "\n")
    for r in rows:
        fh.write(r[0] + "," + ",".join(f"{v:.4e}" for v in r[1:]) + "\n")
print(f"|dMach| against {lab[-1]} (mean / p99 / max):")
for r in rows:
    print(f"  {r[0]:>8}: " + "  ".join(f"{q} {r[1 + 3 * i]:.3f}/{r[2 + 3 * i]:.3f}/{r[3 + 3 * i]:.3f}"
                                       for i, q in enumerate(names)))

cols = [plt.cm.viridis_r(0.15 + 0.7 * j / max(n - 2, 1)) for j in range(n - 1)] + ["k"]
fig, axs = plt.subplots(2, 3, figsize=(19, 9.5), constrained_layout=True)
for i, a in enumerate(axs.flat):
    if i >= len(names):
        a.axis("off"); continue
    x0, r0, x1, r1 = L[i]
    xs = np.linspace(x0, x1, LM.shape[2])
    for j in range(n):
        a.plot(xs, LM[j, i], color=cols[j], lw=1.6 if j == n - 1 else 1.0)
    m, p, mx = stats(LM[-2, i], LM[-1, i])
    a.set_title(f"{names[i]}: (x, r) ({x0:g}, {r0:g}) -> ({x1:g}, {r1:g}) mm\n"
                f"|dMach| {lab[-2]} vs {lab[-1]}: mean {m:.2e}, p99 {p:.2e}, max {mx:.2e}", fontsize=9)
    a.set_xlabel("x [mm]"); a.set_ylabel("Mach"); a.grid(alpha=0.3)
hand = [Line2D([], [], color=cols[j], lw=1.6 if j == n - 1 else 1.0, label=f"{lab[j]} (t4 iteration {its[j]})")
        for j in range(n)]
axs.flat[0].legend(handles=hand, fontsize=8, loc="best")
fig.suptitle(f"{title}: Mach along the monitor lines, final solution of each mesh (no pass threshold set)", fontsize=11)
fig.savefig(pre + "_lines.png", dpi=100)
print("saved", pre + "_lines.png")

# ---- features --------------------------------------------------------------------------------------------------------
head, frows = None, []
for k, f in runs:
    ls = [l.strip() for l in open(os.path.join(f, "features.csv")) if l.strip()]
    if head is None:
        head = ls[0].split(",")
    elif ls[0].split(",") != head:
        sys.exit(f"{k}: features.csv columns differ")
    frows.append(ls[-1].split(","))
val = {c: np.array([float(r[i]) for r in frows]) for i, c in enumerate(head)}
with open(pre + "_features.csv", "w") as fh:
    fh.write("mesh," + ",".join(head) + "\n")
    for k, r in zip(lab, frows):
        fh.write(k + "," + ",".join(r) + "\n")
Q = [c for c in head if c != "iteration"]
fig, axs = plt.subplots(1, len(Q), figsize=(4 * len(Q), 4.2), constrained_layout=True)
xs = np.arange(n)
for a, c in zip(axs, Q):
    v = val[c]
    a.plot(xs, v, "o-", color="#1f77b4", ms=6)
    for x, y in zip(xs, v):
        a.annotate(f"{y:.3f}", (x, y), (0, 7), textcoords="offset points", ha="center", fontsize=8)
    t = c.replace("_x_mm", " x").replace("_mm", "").replace("_", " ")
    a.set_title(f"{t}: {lab[-1]} {v[-1]:.3f} mm\n{lab[-2]} {v[-2] - v[-1]:+.3f} mm from it", fontsize=9)
    a.set_xticks(xs, lab); a.set_ylabel("x [mm]"); a.grid(alpha=0.3); a.margins(y=0.2)
fig.suptitle(f"{title}: axis features vs mesh (Mach peaks = end of each shock cell; core = last x with Mach >= 1; "
             "x from the CFD origin, nozzle exit 1.876 mm)", fontsize=10)
fig.savefig(pre + "_features.png", dpi=110)
print("saved", pre + "_features.png")
for j, k in enumerate(lab):
    print(f"  {k:>8}: " + "  ".join(f"{c} {val[c][j]:.3f}" for c in Q))
