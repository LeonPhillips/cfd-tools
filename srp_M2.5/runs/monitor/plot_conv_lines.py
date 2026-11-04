"""python plot_conv_lines.py <stage folder> [previous stage folder]: <stage>/conv_lines.png (user, 2026-11-04).
Top: mean and p99 |dMach| per 1000 iterations on each convergence line, from line_change.csv (conv_update.py).
Below: Mach along each line, latest snapshot and the one before (from the two kept snapshots), and their difference.
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
cols = plt.cm.viridis(np.linspace(0, 0.9, len(names)))

fig = plt.figure(figsize=(19, 11))
gs = fig.add_gridspec(3, 6, height_ratios=[1.15, 1, 1])
a_mean, a_p99 = fig.add_subplot(gs[0, :3]), fig.add_subplot(gs[0, 3:])
lc = f"{folder}/line_change.csv"
tab = np.atleast_2d(np.genfromtxt(lc, delimiter=",", skip_header=1)) if os.path.exists(lc) and os.path.getsize(lc) > 80 else None
if tab is not None and tab.size:
    for i, (n, c) in enumerate(zip(names, cols)):
        a_mean.semilogy(tab[:, 0], tab[:, 2 + 3 * i], "o-", ms=3, color=c, label=f"{n}  (now {tab[-1, 2 + 3 * i]:.2e})")
        a_p99.semilogy(tab[:, 0], tab[:, 3 + 3 * i], "s-", ms=3, color=c, label=f"{n}  (now {tab[-1, 3 + 3 * i]:.2e})")
    for a, t in ((a_mean, "mean"), (a_p99, "p99")):
        a.set_title(f"{t} |dMach| per 1000 iterations along each line", fontsize=10)
        a.set_xlabel("iteration"); a.set_ylabel(f"{t} |dMach| / 1000 it"); a.grid(alpha=0.3, which="both")
        a.legend(fontsize=8)
else:
    a_mean.text(0.5, 0.5, "no snapshot pair yet", ha="center", transform=a_mean.transAxes)

for i, (n, c) in enumerate(zip(names, cols)):
    a = fig.add_subplot(gs[1 + i // 3, 2 * (i % 3):2 * (i % 3) + 2])
    x0, r0, x1, r1 = L[i]
    xs = np.linspace(x0, x1, zb["line_Mach"].shape[1])
    a.plot(xs, zb["line_Mach"][i], color=c, lw=1.2, label=f"iteration {key(snaps[-1])}")
    if za is not None:
        a.plot(xs, za["line_Mach"][i], color="0.5", lw=0.9, ls="--", label=f"{'prev. stage ' if prev and snaps[-2].startswith(prev) else ''}{key(snaps[-2])}")
        b = a.twinx()
        d = zb["line_Mach"][i].astype(float) - za["line_Mach"][i].astype(float)
        b.plot(xs, d, color="tab:red", lw=0.7, alpha=0.8)
        lim = max(np.nanmax(np.abs(d)) if np.isfinite(d).any() else 0, 1e-6)
        b.set_ylim(-1.1 * lim, 1.1 * lim); b.tick_params(axis="y", colors="tab:red", labelsize=7)
        b.set_ylabel("dMach (latest - previous)", color="tab:red", fontsize=8)
    a.set_title(f"{n}: (x, r) ({x0:g}, {r0:g}) -> ({x1:g}, {r1:g}) mm", fontsize=9)
    a.set_xlabel("x [mm]"); a.set_ylabel("Mach"); a.grid(alpha=0.3); a.legend(fontsize=7, loc="upper left")

run = ([q for q in folder.split("/") if "_mesh" in q] or ["run"])[-1]
fig.suptitle(f"{run} / {folder.split('/')[-1]}: Mach along the convergence lines, latest snapshot iteration {key(snaps[-1])}"
             "   (no pass threshold set; report only)", fontsize=11)
fig.tight_layout()
fig.savefig(f"{folder}/conv_lines.png", dpi=100)
print("saved", f"{folder}/conv_lines.png")
