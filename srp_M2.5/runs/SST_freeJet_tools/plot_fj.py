"""python plot_fj.py <run folder> "<title>"
Free-jet convergence figures for a whole chain (user, 2026-11-23). The stages in <run>/2_jet_staged share one
cumulative iteration axis (each is offset by the ITER of the stages before it). Writes in <run>:
    convergence_chain.png    residuals, CFL, boundary mass flows, mass imbalance (information only)
    conv_lines_chain.png     Mach on each monitor line for every snapshot, coloured by iteration, latest in black
    conv_features_chain.png  axis shock-cell peaks, supersonic core length, and mean |dMach| per 1000 iterations
Prints one status line. Reads only history.csv, line_history.npz, line_change.csv and features.csv. Stops nothing.
"""
import csv
import glob
import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

run = sys.argv[1].rstrip("/")
title = sys.argv[2] if len(sys.argv) > 2 else os.path.basename(run)
stages = sorted(d for d in glob.glob(f"{run}/2_jet_staged/t*") if os.path.isdir(d))


def iters(s):
    t = open(f"{s}/{os.path.basename(s)}.cfg").read()
    return int(re.search(r"^ITER= *(\d+)", t, re.M).group(1))


def table(path):
    rows = list(csv.reader(open(path)))
    head = [c.strip().strip('"') for c in rows[0]]
    data = np.array([[float(v) if v.strip() else np.nan for v in r[:len(head)]] for r in rows[1:] if len(r) >= len(head)])
    return head, data


H, LH, LC, FE, bounds, off = [], [], [], [], [], 0
for s in stages:
    name = os.path.basename(s)
    if os.path.exists(f"{s}/history.csv"):
        h, d = table(f"{s}/history.csv")
        if len(d):
            H.append((name, off, h, d))
    if os.path.exists(f"{s}/line_history.npz"):
        z = np.load(f"{s}/line_history.npz", allow_pickle=True)
        LH.append((z["iterations"] + off, z["mach"], z["lines_mm"], z["line_names"]))
    for lst, fn in ((LC, "line_change.csv"), (FE, "features.csv")):
        if os.path.exists(f"{s}/{fn}") and os.path.getsize(f"{s}/{fn}") > 0:
            h, d = table(f"{s}/{fn}")
            if len(d):
                d = d.copy(); d[:, 0] += off
                lst.append((h, d))
    bounds.append((off, name))
    off += iters(s)
if not H:
    sys.exit("no history yet")


def col(h, d, name):
    return d[:, h.index(name)] if name in h else np.full(len(d), np.nan)


def stage_lines(ax):
    for o, n in bounds[1:]:
        ax.axvline(o, color="0.6", lw=0.6, ls=":")


# ---- convergence_chain.png ---------------------------------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(14, 9))
for name, o, h, d in H:
    it = col(h, d, "Inner_Iter") + o
    for k, c in (("rms[Rho]", "C0"), ("rms[RhoU]", "C1"), ("rms[RhoE]", "C2"), ("rms[k]", "C3"), ("rms[w]", "C4")):
        ax[0, 0].plot(it, col(h, d, k), c, lw=0.5, label=k if name == H[0][0] else None)
    ax[0, 1].plot(it, col(h, d, "AvgCFL"), "k", lw=0.5)
    flows = {m: np.abs(col(h, d, f"Avg_Massflow({m})")) for m in ("inlet_jet", "ambient_side", "ambient_rear", "outlet")}
    for m, c in zip(flows, ("C0", "C1", "C2", "k")):
        ax[1, 0].plot(it, flows[m], c, lw=0.7, label=m if name == H[0][0] else None)
    inflow = flows["inlet_jet"] + flows["ambient_side"] + flows["ambient_rear"]
    ax[1, 1].plot(it, 100 * (inflow - flows["outlet"]) / flows["outlet"], "k", lw=0.5)
ax[0, 0].set_ylabel("log10 rms residual"); ax[0, 0].legend(fontsize=8)
ax[0, 1].set_ylabel("average CFL")
ax[1, 0].set_ylabel("|mass flow| [kg/s per rad]"); ax[1, 0].set_yscale("log"); ax[1, 0].legend(fontsize=8)
ax[1, 1].set_ylabel("(in - out) / out [%]  (sign of each flow ignored)")
for a in ax.flat:
    a.set_xlabel("cumulative iteration"); stage_lines(a); a.grid(alpha=0.3)
fig.suptitle(f"{title}: residuals, CFL and mass flows (information only)")
fig.tight_layout(); fig.savefig(f"{run}/convergence_chain.png", dpi=110); plt.close(fig)

# ---- conv_lines_chain.png -----------------------------------------------------------------------------------------
last_change = {}
if LC:
    h, d = LC[-1]
    for j in range(len(h) // 3):
        last_change[j] = d[-1, 1 + 3 * j]
if LH:
    its = np.concatenate([a[0] for a in LH]); mach = np.concatenate([a[1] for a in LH])
    L, names = LH[-1][2], LH[-1][3]
    fig, ax = plt.subplots(len(L), 1, figsize=(14, 3.2 * len(L)))
    cmap = plt.get_cmap("viridis")
    for j, (x0, r0, x1, r1) in enumerate(L):
        s = np.linspace(0, 1, mach.shape[2])
        xs = x0 + s * (x1 - x0)
        for k in range(len(its)):
            last = k == len(its) - 1
            ax[j].plot(xs, mach[k, j], color="k" if last else cmap(k / max(1, len(its) - 1)), lw=1.0 if last else 0.5)
        extra = f"; latest mean |dMach| per 1000 it {last_change[j]:.2e}" if j in last_change else ""
        ax[j].set_title(f"{names[j]}: ({x0:g}, {r0:g}) to ({x1:g}, {r1:g}) mm{extra}", fontsize=10)
        ax[j].set_xlabel("x [mm]"); ax[j].set_ylabel("Mach"); ax[j].grid(alpha=0.3)
    fig.suptitle(f"{title}: Mach on the monitor lines, iterations {its[0]}-{its[-1]} (latest black)")
    fig.tight_layout(); fig.savefig(f"{run}/conv_lines_chain.png", dpi=100); plt.close(fig)

# ---- conv_features_chain.png --------------------------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
for h, d in FE:
    for i in range(4):
        ax[0].plot(d[:, 0], d[:, 1 + i], "o-", ms=2, lw=0.8, color=f"C{i}")
    ax[1].plot(d[:, 0], d[:, 5], "ko-", ms=2, lw=0.8)
for h, d in LC:
    for j in range(len(h) // 3):
        ax[2].semilogy(d[:, 0], d[:, 1 + 3 * j], "o-", ms=2, lw=0.8, color=f"C{j}")
if LH:
    ax[2].legend([str(n) for n in LH[-1][3]], fontsize=7)
ax[0].set_ylabel("x of axis-Mach peak [mm] (peaks 1-4)"); ax[1].set_ylabel("supersonic core length [mm]")
ax[2].set_ylabel("mean |dMach| per 1000 iterations")
for a in ax:
    a.set_xlabel("cumulative iteration"); stage_lines(a); a.grid(alpha=0.3)
fig.suptitle(f"{title}: flow features and line changes")
fig.tight_layout(); fig.savefig(f"{run}/conv_features_chain.png", dpi=110); plt.close(fig)

name, o, h, d = H[-1]
msg = f"{name} it {int(col(h, d, 'Inner_Iter')[-1])} (cum {int(col(h, d, 'Inner_Iter')[-1] + o)}): rms[Rho] {col(h, d, 'rms[Rho]')[-1]:.3f}, CFL {col(h, d, 'AvgCFL')[-1]:.3g}, jet {abs(col(h, d, 'Avg_Massflow(inlet_jet)')[-1]):.5g} kg/s"
if FE:
    msg += f", peaks {FE[-1][1][-1, 1:5].round(3).tolist()} mm, core {FE[-1][1][-1, 5]:.2f} mm"
print(msg)
