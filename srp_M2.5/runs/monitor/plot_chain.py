"""python plot_chain.py <run folder> [title tag]: the three monitor figures for a whole staged chain (user, 2026-11-05).
Joins the stages in <run>/2_jet_staged/t* (sorted: t1, t2, t3) on one cumulative iteration axis: stage k's iterations are
offset by the ITER of the stages before it (from each stage's .cfg), so with 5000 / 25000 / 20000, t2 iteration 20000 is
25000 and t3 iteration 20000 is 50000. Stage boundaries are dashed vertical lines.
Writes, in <run>:
- convergence_chain.png: residuals | probe Mach | box change per 1000 it
                         mass flows | average CFL | Mach difference map (latest two snapshots of the last stage)
  The probes are re-placed at the start of each stage (probes_md.py), so a step in probe Mach at a boundary is partly the
  probe moving; each stage's probe x on the axis is in the legend.
- conv_lines_chain.png: Mach on the six lines for every snapshot of every stage, coloured by cumulative iteration; the
  last snapshot of each earlier stage is dashed, the latest is black.
- conv_features_chain.png: Mach disk x, Mach disk r, bow shock x, stagnation x against cumulative iteration, one colour
  per stage.
Reads only the monitor's own files (history.csv, box_change.csv, line_change.csv, features.csv, line_history.npz, snaps/).
Report only; no pass threshold set.
"""
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

MAP = (0.0, 16.0, 0.0, 6.0)
run = sys.argv[1].rstrip("/")
tag = sys.argv[2] if len(sys.argv) > 2 else os.path.basename(run)
stages = sorted(d for d in glob.glob(f"{run}/2_jet_staged/t*") if os.path.isdir(d))
stages = [d for d in stages if os.path.isfile(f"{d}/history.csv")]   # skip stages set up but not started yet
SC = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]   # one colour per stage


def cfg_of(s):
    return open(glob.glob(f"{s}/*.cfg")[0]).read()


def iters(s):
    return int(re.search(r"^ITER=\s*(\d+)", cfg_of(s), re.M).group(1))


def probes(s):
    m = re.search(r"^CUSTOM_OUTPUTS=(.*)$", cfg_of(s), re.M)
    return {n: (float(x) * 1e3, float(r) * 1e3) for n, x, r in
            re.findall(r"(\w+)\s*:\s*Probe\{[^}]*\}\[\s*([-0-9.eE]+)\s*,\s*([-0-9.eE]+)\s*\]", m.group(1) if m else "")}


def table(path):
    """(header, float array) of a CSV; empty array if missing."""
    if not os.path.exists(path):
        return [], np.zeros((0, 0))
    head = [c.strip().strip('"') for c in open(path).readline().split(",")]
    a = np.genfromtxt(path, delimiter=",", skip_header=1, usecols=range(len(head) - (head[-1] == "flags")))
    return head, np.atleast_2d(a)


off, o = [], 0
for s in stages:
    off.append(o); o += iters(s)
name = [os.path.basename(s) for s in stages]
nice = [n.replace("_cap", " cap ") for n in name]


def bounds(a, top=True):
    for k in range(1, len(stages)):
        a.axvline(off[k], color="0.4", ls="--", lw=0.8)
    if top:
        for k in range(len(stages)):
            a.text(off[k] + 0.005 * o, 0.985, name[k].split("_")[0], transform=a.get_xaxis_transform(), fontsize=8,
                   va="top", color=SC[k], fontweight="bold")


# ---------- convergence_chain.png
H = []
for k, s in enumerate(stages):
    head = [c.strip().strip('"') for c in open(f"{s}/history.csv").readline().split(",")]
    a = np.loadtxt(f"{s}/history.csv", delimiter=",", skiprows=1)
    H.append(({h: a[:, i] for i, h in enumerate(head)}, a[:, head.index("Inner_Iter")] + off[k]))
fig, ax = plt.subplots(2, 3, figsize=(19, 9.5))
RES = [r for r in ["rms[Rho]", "rms[RhoU]", "rms[RhoE]", "rms[nu]", "rms[k]", "rms[w]"] if r in H[0][0]]
for j, r in enumerate(RES):
    for k, (c, it) in enumerate(H):
        ax[0, 0].plot(it, c[r], lw=0.6, color=f"C{j}", label=r if k == 0 else None)
ax[0, 0].set_ylabel("log10 residual"); ax[0, 0].legend(fontsize=8)
ax[0, 0].set_title("residuals (information only)", fontsize=9)

PN = sorted({p for s in stages for p in probes(s)}, key=lambda p: (p[2:], p[:2]))
for j, p in enumerate(PN):
    for k, (c, it) in enumerate(H):
        if p in c:
            ax[0, 1].plot(it, c[p], lw=0.8, color=f"C{j}", label=p if k == 0 else None)
xs = "; ".join(f"{nice[k]} up0/dn0 x {probes(s).get('up0', (np.nan,))[0]:.2f}/{probes(s).get('dn0', (np.nan,))[0]:.2f}"
               for k, s in enumerate(stages))
ax[0, 1].set_ylabel("probe Mach"); ax[0, 1].legend(fontsize=7, ncol=2)
ax[0, 1].set_title(f"Mach at the probes (re-placed 1 mm either side of the Mach disk at each stage start)\n{xs} mm",
                   fontsize=8)

for k, s in enumerate(stages):
    h, p = table(f"{s}/box_change.csv")
    if p.size:
        i = p[:, 0] + off[k]
        ax[0, 2].semilogy(i, p[:, 2], "o-", ms=2.5, color="C0", label="mean |dMach|" if k == 0 else None)
        ax[0, 2].semilogy(i, p[:, 3], "s-", ms=2.5, color="C1", label="p99 |dMach|" if k == 0 else None)
        ax[0, 2].semilogy(i, p[:, 5], "D-", ms=2.5, color="C2", label="mean |dMach| of the top 1 %" if k == 0 else None)
        ax[0, 2].semilogy(i, p[:, 4], "^--", ms=2.5, color="C3", label="mean |drho|/rho [%]" if k == 0 else None)
        last = p[-1]
ax[0, 2].legend(fontsize=8); ax[0, 2].set_ylabel("change per 1000 iterations")
ax[0, 2].set_title(f"box x {MAP[0]:g}-{MAP[1]:g}, r {MAP[2]:g}-{MAP[3]:g} mm: change per 1000 iterations\n"
                   f"latest: mean |dM| {last[2]:.2e}, p99 {last[3]:.2e}, top-1 % {last[5]:.2e}, "
                   f"|drho|/rho {last[4]:.2e} %", fontsize=9)

FL = [h for h in H[0][0] if h.startswith("Avg_Massflow(")]
for j, f in enumerate(FL):
    for k, (c, it) in enumerate(H):
        ax[1, 0].plot(it, c[f], lw=0.8, color=f"C{j}", label=f.replace("Avg_Massflow", "") if k == 0 else None)
ax[1, 0].set_ylabel("mass flow [kg/s, 360 deg]"); ax[1, 0].legend(fontsize=8, loc="center left")
if "Avg_Massflow(inlet_jet)" in H[0][0]:   # the jet is ~1 % of the freestream: its own axis
    tw = ax[1, 0].twinx()
    for k, (c, it) in enumerate(H):
        tw.plot(it, c["Avg_Massflow(inlet_jet)"], lw=0.9, color="k")
    tw.set_ylabel("jet mass flow [kg/s] (black)")
    ax[1, 0].set_title(f"mass flows; jet (black, right axis) now {H[-1][0]['Avg_Massflow(inlet_jet)'][-1]:.6f} kg/s",
                       fontsize=9)
cl = "Avg CFL" if "Avg CFL" in H[0][0] else "AvgCFL"
for k, (c, it) in enumerate(H):
    ax[1, 1].semilogy(it, c[cl], lw=0.8, color=SC[k])
ax[1, 1].set_ylabel("average CFL"); ax[1, 1].set_title(f"CFL now {H[-1][0][cl][-1]:.3g}", fontsize=9)
for a in list(ax[0, :]) + list(ax[1, :2]):
    a.grid(alpha=0.3); a.set_xlabel("cumulative iteration"); bounds(a)

key = lambda f: int(f.rsplit("_", 1)[-1][:-4])
sn = sorted(glob.glob(f"{stages[-1]}/snaps/flow_*.npz"), key=key)
if len(sn) >= 2:
    za, zb = np.load(sn[-2]), np.load(sn[-1])
    ix = (zb["x_mm"] >= MAP[0]) & (zb["x_mm"] <= MAP[1]); ir = (zb["r_mm"] >= MAP[2]) & (zb["r_mm"] <= MAP[3])
    xm, rm = zb["x_mm"][ix], zb["r_mm"][ir]
    cm = lambda a: a[np.ix_(ir, ix)].astype(float)
    d = cm(zb["Mach"]) - cm(za["Mach"])
    lim = max(np.nanpercentile(np.abs(d), 99.5), 1e-6)
    im = ax[1, 2].imshow(np.ma.masked_invalid(d), origin="lower", cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="equal",
                         extent=[xm[0], xm[-1], rm[0], rm[-1]], interpolation="nearest")
    fig.colorbar(im, ax=ax[1, 2], shrink=0.8, label="dMach (clipped at p99.5)")
    ax[1, 2].contour(xm, rm, np.ma.masked_invalid(cm(zb["Mach"])), levels=[1.0, 2.0, 2.5, 3.0], colors="k",
                     linewidths=0.4, alpha=0.5)
    g = lambda f: key(f) + off[-1]
    ax[1, 2].set_title(f"Mach, snapshot {g(sn[-1])} minus {g(sn[-2])} (cumulative)   max |dM| {np.nanmax(np.abs(d)):.2e}\n"
                       "(grey: current Mach 1, 2, 2.5, 3)", fontsize=9)
    ax[1, 2].set_xlabel("x [mm] (freestream from the right)"); ax[1, 2].set_ylabel("r [mm]")
c = H[-1][0]
mi = c["Avg_Massflow(inlet_freestream)"][-1] + c.get("Avg_Massflow(inlet_jet)", np.zeros(1))[-1]
imb = 100 * abs(mi + c["Avg_Massflow(outlet)"][-1]) / mi
fig.suptitle(f"{tag}: {' + '.join(nice)}   cumulative iteration {int(H[-1][1][-1])}   rms[Rho] {c['rms[Rho]'][-1]:.2f}   "
             f"mass imbalance {imb:.3f} %   jet {c.get('Avg_Massflow(inlet_jet)', [np.nan])[-1]:.6f} kg/s   "
             f"CD {c['CD'][-1]:.5f} (info only)", fontsize=10)
fig.tight_layout()
fig.savefig(f"{run}/convergence_chain.png", dpi=110); plt.close(fig)
print("saved", f"{run}/convergence_chain.png")

# ---------- conv_lines_chain.png
IT, LM, stage_end = [], [], []
for k, s in enumerate(stages):
    with np.load(f"{s}/line_history.npz") as h:
        IT += list(h["iter"] + off[k]); LM += list(h["line_Mach"])
        names, L = [str(q) for q in h["line_names"]], h["lines_mm"]
    stage_end.append(len(IT) - 1)
IT, LM = np.array(IT), np.array(LM)
_, lc = table(f"{stages[-1]}/line_change.csv")
cmap = plt.cm.viridis_r
norm = Normalize(vmin=IT[0], vmax=IT[-1])
fig, axs = plt.subplots(2, 3, figsize=(19, 9.5), constrained_layout=True)
for i, (a, n) in enumerate(zip(axs.flat, names)):
    x0, r0, x1, r1 = L[i]
    xs = np.linspace(x0, x1, LM.shape[2])
    for k in range(len(IT) - 1):
        a.plot(xs, LM[k, i], color=cmap(norm(IT[k])), lw=0.7)
    for k, e in enumerate(stage_end[:-1]):
        a.plot(xs, LM[e, i], color=SC[k], lw=1.2, ls="--")
    a.plot(xs, LM[-1, i], color="k", lw=1.6)
    t = f"{n}: (x, r) ({x0:g}, {r0:g}) -> ({x1:g}, {r1:g}) mm"
    if lc.size:
        t += f"\nlatest |dMach| per 1000 it: mean {lc[-1, 2 + 3 * i]:.2e}, p99 {lc[-1, 3 + 3 * i]:.2e}, max {lc[-1, 4 + 3 * i]:.2e}"
    a.set_title(t, fontsize=9); a.set_xlabel("x [mm]"); a.set_ylabel("Mach"); a.grid(alpha=0.3)
fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs, shrink=0.85, pad=0.01,
             label="snapshot cumulative iteration")
hand = [Line2D([], [], color="k", lw=1.6, label=f"latest (cumulative iteration {int(IT[-1])})")]
hand += [Line2D([], [], color=SC[k], lw=1.2, ls="--", label=f"end of {nice[k]} ({int(IT[e])})")
         for k, e in enumerate(stage_end[:-1])]
hand.append(Line2D([], [], color=cmap(0.5), lw=0.7, label=f"other snapshots, every 1000 it ({len(IT) - len(stages)})"))
axs.flat[0].legend(handles=hand, fontsize=8, loc="best")
fig.suptitle(f"{tag}: Mach along the convergence lines, {' + '.join(nice)}, all {len(IT)} snapshots on one cumulative "
             "iteration scale   (no pass threshold set; report only)", fontsize=11)
fig.savefig(f"{run}/conv_lines_chain.png", dpi=100); plt.close(fig)
print("saved", f"{run}/conv_lines_chain.png")

# ---------- conv_features_chain.png
Q = [("mach_disk_x", "Mach disk x on the axis (mm)", "Mach disk x"),
     ("mach_disk_subsonic_edge_r", "Mach disk radius r (mm)", "Mach disk r"),
     ("bow_x_axis", "Bow shock x on the axis (mm)", "Bow shock x"),
     ("stagnation_x", "Stagnation point x on the axis (mm)", "Stagnation x")]
F = [table(f"{s}/features.csv") for s in stages]
fig, axs = plt.subplots(2, 2, figsize=(13, 8.5), constrained_layout=True)
for a, (q, ylab, title) in zip(axs.flat, Q):
    ends = []
    for k, (h, v) in enumerate(F):
        j = h.index(q)
        a.plot(v[:, 0] + off[k], v[:, j], "o-", ms=3.5, color=SC[k], label=nice[k])
        ends.append(v[-1, j])
    step = " / ".join(f"{e:.3f}" for e in ends)
    a.set_title(f"{title}: end of each stage {step} mm", fontsize=10, pad=14)
    a.set_xlabel("cumulative iteration"); a.set_ylabel(ylab); a.grid(alpha=0.3); bounds(a, top=False)
axs.flat[0].legend(fontsize=8)
fig.suptitle(f"{tag}: feature positions vs cumulative iteration, {' + '.join(nice)} "
             "(x from CFD origin; nozzle exit x = 1.875 mm)")
fig.savefig(f"{run}/conv_features_chain.png", dpi=110); plt.close(fig)
print("saved", f"{run}/conv_features_chain.png")
