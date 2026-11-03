"""python schlieren_compare_two.py <a.npz> <label a> <b.npz> <label b> <prefix>
Two CFD solutions against the same mean experimental image, from the <prefix>.npz files written by schlieren_compare.py
(each keeps its own LOS fit E = s*C + c). Both must use the same experimental image.
Writes <prefix>_mirror_LOS.png (one mirror panel per solution, same colour scale), <prefix>_profiles.png and
<prefix>_features.csv. Env: CASE_LABEL for the title."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fa, la, fb, lb, pre = sys.argv[1:6]
CASE_LABEL = os.environ.get("CASE_LABEL", "NPR 501 case")
S = []
for fn, lab in ((fa, la), (fb, lb)):
    z = np.load(fn)
    s, k, cc = z["fit_los"]
    S.append(dict(lab=lab, L=s * z["los"] + k, s=s, k=k, cc=cc, blocked=z["blocked"], z=z))
z = S[0]["z"]
xcol, yrow, E, expbody, D = z["xD"], z["yD"], z["exp"], z["expbody"], float(z["D_mm"])
if not np.array_equal(E, S[1]["z"]["exp"], equal_nan=True):
    sys.exit("the two files use different experimental images")

use = ~S[0]["blocked"] & ~expbody & np.isfinite(E)
lim = np.nanpercentile(np.abs(E[use]), 99.5)
ext = [xcol[0], xcol[-1], yrow[0], yrow[-1]]
cm = plt.get_cmap("gray").copy(); cm.set_bad("0.75")
fig, axs = plt.subplots(2, 1, figsize=(12, 16))
for ax, d in zip(axs, S):
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
axs[-1].set_xlabel("x / D from nozzle exit (freestream from the right)")
fig.suptitle(f"d(rho)/dx schlieren, {CASE_LABEL}, D = {D} mm; each CFD scaled by its own fit; grey: model / masked")
fig.tight_layout(); fig.savefig(f"{pre}_mirror_LOS.png", dpi=150); plt.close(fig)
print("saved", f"{pre}_mirror_LOS.png")

# profiles along y = const (both halves of the experiment averaged), as in schlieren_compare.py
lines_mm = [0.01, 1.3, 2.3, 3.0, 6.0, 9.0]


def rowavg(A, ymm):
    iu = np.argmin(np.abs(yrow * D - ymm)); il = np.argmin(np.abs(yrow * D + ymm))
    return 0.5 * (A[iu] + A[il]) if iu != il else A[iu]


fig, axs = plt.subplots(3, 2, figsize=(15, 12), sharex=True)
for ax, ymm in zip(axs.flat, lines_mm):
    bm = rowavg((S[0]["blocked"] | S[1]["blocked"] | expbody).astype(float), ymm) > 0
    for A, lab, sty in ((E, "experiment", "k-"), (S[0]["L"], f"CFD {la}", "C0-"), (S[1]["L"], f"CFD {lb}", "C3-")):
        p = rowavg(A, ymm).copy(); p[bm] = np.nan
        ax.plot(xcol, p, sty, lw=0.9 if lab != "experiment" else 1.2, label=lab)
    ax.set_title(f"y = {ymm} mm ({ymm / D:.3f} D)"); ax.grid(alpha=0.3); ax.set_ylabel("signal [grey levels]")
axs[0, 0].legend()
for ax in axs[-1]:
    ax.set_xlabel("x / D from nozzle exit")
fig.suptitle("Schlieren d(rho)/dx along lines of constant y: experiment (mean of +y and -y rows) vs both CFD (LOS, fitted)")
fig.tight_layout(); fig.savefig(f"{pre}_profiles.png", dpi=130); plt.close(fig)
print("saved", f"{pre}_profiles.png")


# feature positions, same definitions as schlieren_compare.py
def pick(p, lo, hi, sgn):
    m = (xcol >= lo) & (xcol <= hi) & np.isfinite(p)
    return xcol[np.where(m)[0][np.argmax(sgn * p[m])]]


def edge(p, lo, hi):
    m = (xcol >= lo) & (xcol <= hi) & np.isfinite(p)
    xs, ps = xcol[m][::-1], p[m][::-1]
    return xs[np.argmax(ps < 0.5 * ps.min())]


rows = []
for ymm, feat, lo, hi, sgn in ((0.01, "terminal shock (max)", 0.8, 2.2, 1), (0.01, "bow shock (min)", 1.6, 4.5, -1),
                               (4.6, "bow shock (min)", 1.0, 4.5, -1), (9.2, "bow shock (min)", 0.5, 4.5, -1)):
    rows.append((ymm, feat, *[pick(rowavg(A, ymm), lo, hi, sgn) for A in (E, S[0]["L"], S[1]["L"])]))
for ymm, lo, hi in ((0.01, 1.6, 4.5), (4.6, 1.0, 4.5), (9.2, 0.5, 4.5)):
    rows.append((ymm, "bow shock upstream edge", *[edge(rowavg(A, ymm), lo, hi) for A in (E, S[0]["L"], S[1]["L"])]))
with open(f"{pre}_features.csv", "w") as fo:
    fo.write(f"# A = {la}; B = {lb}\ny_mm,feature,exp_xD,A_xD,B_xD,A_minus_exp_mm,B_minus_exp_mm\n")
    print(f"A = {la}; B = {lb}")
    print(f"{'y mm':>6} {'feature':<24} {'exp x/D':>8} {'A x/D':>7} {'B x/D':>7} {'A-exp mm':>9} {'B-exp mm':>9}")
    for ymm, feat, e, a, b in rows:
        fo.write(f"{ymm},{feat},{e:.3f},{a:.3f},{b:.3f},{(a - e) * D:.2f},{(b - e) * D:.2f}\n")
        print(f"{ymm:6.2f} {feat:<24} {e:8.3f} {a:7.3f} {b:7.3f} {(a - e) * D:9.2f} {(b - e) * D:9.2f}")
print("saved", f"{pre}_features.csv")
