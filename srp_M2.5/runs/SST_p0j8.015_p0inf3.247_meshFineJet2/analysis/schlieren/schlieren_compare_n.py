"""python schlieren_compare_n.py <prefix> <tag>=<case.npz>=<label> ... (user, 2026-11-05)
Any number of CFD solutions against the same mean experimental image, from the <case>.npz files written by
schlieren_compare.py (each keeps its own LOS fit E = s*C + c). All must use the same experimental image.
Same panels as schlieren_compare_two.py's <prefix>_mirror_LOS.png (experiment above the axis, CFD LOS below, one shared
colour scale), but for N cases:
- <prefix>_mirror_LOS.png: all cases in a single column, in the order given
- <prefix>_mirror_LOS_<tag>.png: each case on its own
- <prefix>_features.csv: feature positions (schlieren_compare.py definitions) for the experiment and every case
Env: CASE_LABEL for the title.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

pre = sys.argv[1]
CASE_LABEL = os.environ.get("CASE_LABEL", "NPR 501 case, p0inf 3.247 bar")
S = []
for a in sys.argv[2:]:
    tag, fn, lab = a.split("=", 2)
    z = np.load(fn)
    s, k, cc = z["fit_los"]
    S.append(dict(tag=tag, lab=lab, L=s * z["los"] + k, s=s, k=k, cc=cc, blocked=z["blocked"], z=z))
z = S[0]["z"]
xcol, yrow, E, expbody, D = z["xD"], z["yD"], z["exp"], z["expbody"], float(z["D_mm"])
for d in S[1:]:
    if not np.array_equal(E, d["z"]["exp"], equal_nan=True):
        sys.exit(f"{d['tag']} uses a different experimental image")

use = ~S[0]["blocked"] & ~expbody & np.isfinite(E)
lim = np.nanpercentile(np.abs(E[use]), 99.5)    # one colour scale for every panel and file
ext = [xcol[0], xcol[-1], yrow[0], yrow[-1]]
cm = plt.get_cmap("gray").copy(); cm.set_bad("0.75")
SUP = f"d(rho)/dx schlieren, {CASE_LABEL}, D = {D} mm; each CFD scaled by its own fit; grey: model / masked"


def panel(fig, ax, d, xlabel=True):
    comp = np.where(yrow[:, None] >= 0, E, d["L"])
    comp = np.where(d["blocked"] | expbody, np.nan, comp)
    im = ax.imshow(comp, extent=ext, origin="lower" if yrow[0] < yrow[-1] else "upper", cmap=cm,
                   vmin=-lim, vmax=lim, aspect="equal", interpolation="nearest")
    ax.axhline(0, color="y", lw=0.6, ls="--")
    kw = dict(transform=ax.transAxes, ha="left", fontsize=12, color="w", fontweight="bold",
              bbox=dict(boxstyle="round", fc="k", alpha=0.6))
    ax.text(0.01, 0.97, "experiment, mean (above axis)", va="top", **kw)
    ax.text(0.01, 0.03, f"CFD {d['lab']}, line-of-sight integrated (below axis)", va="bottom", **kw)
    ax.set_ylabel("y / D")
    ax.set_title(f"CFD {d['lab']}: E = {d['s']:.3g} C + {d['k']:.3g}, correlation {d['cc']:.2f}")
    fig.colorbar(im, ax=ax, shrink=0.8).set_label("schlieren signal [exp grey levels]")
    if xlabel:
        ax.set_xlabel("x / D from nozzle exit (freestream from the right)")


fig, axs = plt.subplots(len(S), 1, figsize=(12, 8 * len(S)))
for i, (ax, d) in enumerate(zip(np.atleast_1d(axs), S)):
    panel(fig, ax, d, xlabel=i == len(S) - 1)
fig.suptitle(SUP)
fig.tight_layout(); fig.savefig(f"{pre}_mirror_LOS.png", dpi=150); plt.close(fig)
print("saved", f"{pre}_mirror_LOS.png")
for d in S:
    fig, ax = plt.subplots(figsize=(12, 8.6))
    panel(fig, ax, d)
    fig.suptitle(SUP, fontsize=10)
    fig.tight_layout(); fig.savefig(f"{pre}_mirror_LOS_{d['tag']}.png", dpi=150); plt.close(fig)
    print("saved", f"{pre}_mirror_LOS_{d['tag']}.png")


# feature positions, same definitions as schlieren_compare.py / schlieren_compare_two.py
def rowavg(A, ymm):
    iu = np.argmin(np.abs(yrow * D - ymm)); il = np.argmin(np.abs(yrow * D + ymm))
    return 0.5 * (A[iu] + A[il]) if iu != il else A[iu]


def pick(p, lo, hi, sgn):
    m = (xcol >= lo) & (xcol <= hi) & np.isfinite(p)
    return xcol[np.where(m)[0][np.argmax(sgn * p[m])]]


def edge(p, lo, hi):
    m = (xcol >= lo) & (xcol <= hi) & np.isfinite(p)
    xs, ps = xcol[m][::-1], p[m][::-1]
    return xs[np.argmax(ps < 0.5 * ps.min())]


A = [E] + [d["L"] for d in S]
rows = []
for ymm, feat, lo, hi, sgn in ((0.01, "terminal shock (max)", 0.8, 2.2, 1), (0.01, "bow shock (min)", 1.6, 4.5, -1),
                               (4.6, "bow shock (min)", 1.0, 4.5, -1), (9.2, "bow shock (min)", 0.5, 4.5, -1)):
    rows.append((ymm, feat, [pick(rowavg(a, ymm), lo, hi, sgn) for a in A]))
for ymm, lo, hi in ((0.01, 1.6, 4.5), (4.6, 1.0, 4.5), (9.2, 0.5, 4.5)):
    rows.append((ymm, "bow shock upstream edge", [edge(rowavg(a, ymm), lo, hi) for a in A]))
T = [d["tag"] for d in S]
with open(f"{pre}_features.csv", "w") as fo:
    fo.write("# " + "; ".join(f"{d['tag']} = {d['lab']}" for d in S) + "\n")
    fo.write("y_mm,feature,exp_xD," + ",".join(f"{t}_xD" for t in T) + "," + ",".join(f"{t}_minus_exp_mm" for t in T) + "\n")
    print(f"{'y mm':>6} {'feature':<24} {'exp x/D':>8} " + " ".join(f"{t + ' -exp mm':>14}" for t in T))
    for ymm, feat, v in rows:
        e = v[0]
        fo.write(f"{ymm},{feat},{e:.3f}," + ",".join(f"{x:.3f}" for x in v[1:]) + ","
                 + ",".join(f"{(x - e) * D:.2f}" for x in v[1:]) + "\n")
        print(f"{ymm:6.2f} {feat:<24} {e:8.3f} " + " ".join(f"{(x - e) * D:14.2f}" for x in v[1:]))
print("saved", f"{pre}_features.csv")
