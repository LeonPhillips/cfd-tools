"""python conv_fj.py <stage folder>
Free-jet live monitor (user, 2026-11-23): folds new snaps/flow_<iter>.npz into
    line_history.npz   Mach on every line for every snapshot of the stage (small; kept)
    line_change.csv    per snapshot and line: mean, p99 and max |dMach| per 1000 iterations against the previous snapshot
    features.csv       per snapshot, from the axis line: x of the first four axis-Mach peaks (the end of each shock cell,
                       just upstream of its shock) and the supersonic core length (last x with Mach >= 1)
Only the latest two snapshot npz files are kept. Nothing here stops a run.
"""
import glob
import os
import sys
import numpy as np

D = sys.argv[1].rstrip("/")
NPEAK = 4


def features(x, m):
    ok = np.isfinite(m)
    x, m = x[ok], m[ok]
    peaks = []
    for i in range(2, len(m) - 2):   # local maxima that stand at least 0.05 above the minima either side
        if m[i] >= m[i - 1] and m[i] > m[i + 1] and m[i] >= 1.0:
            lo = max(0, i - 200); hi = min(len(m), i + 200)
            if m[i] - min(m[lo:i].min(), m[i:hi].min()) > 0.05 and (not peaks or x[i] - peaks[-1] > 1.0):
                peaks.append(x[i])
        if len(peaks) == NPEAK:
            break
    sup = np.where(m >= 1.0)[0]
    core = x[sup[-1]] if len(sup) else np.nan
    return peaks + [np.nan] * (NPEAK - len(peaks)), core


files = sorted(glob.glob(f"{D}/snaps/flow_*.npz"))
hist = f"{D}/line_history.npz"
H = dict(np.load(hist, allow_pickle=True)) if os.path.exists(hist) else None
done = set(H["iterations"].tolist()) if H is not None else set()
new = [f for f in files if int(os.path.basename(f)[5:11]) not in done]
for f in new:
    z = np.load(f, allow_pickle=True)
    it = int(z["iteration"])
    if H is None:
        H = {"iterations": np.array([it]), "mach": z["mach"][None], "lines_mm": z["lines_mm"],
             "line_names": z["line_names"]}
    else:
        H["iterations"] = np.append(H["iterations"], it)
        H["mach"] = np.concatenate([H["mach"], z["mach"][None]])
if new:
    o = np.argsort(H["iterations"])
    H["iterations"], H["mach"] = H["iterations"][o], H["mach"][o]
    np.savez_compressed(hist, **H)
    L = H["lines_mm"]
    x_axis = np.linspace(L[0][0], L[0][2], H["mach"].shape[2])
    with open(f"{D}/line_change.csv", "w") as fh:
        fh.write("iteration," + ",".join(f"L{j}_{s}" for j in range(len(L)) for s in ("mean", "p99", "max")) + "\n")
        for k in range(1, len(H["iterations"])):
            per = (H["iterations"][k] - H["iterations"][k - 1]) / 1000.0
            d = np.abs(H["mach"][k] - H["mach"][k - 1]) / per
            row = []
            for j in range(len(L)):
                v = d[j][np.isfinite(d[j])]
                row += [v.mean(), np.percentile(v, 99), v.max()] if len(v) else [np.nan] * 3
            fh.write(f"{H['iterations'][k]}," + ",".join(f"{v:.6g}" for v in row) + "\n")
    with open(f"{D}/features.csv", "w") as fh:
        fh.write("iteration," + ",".join(f"peak{i + 1}_x_mm" for i in range(NPEAK)) + ",core_length_mm\n")
        for k, it in enumerate(H["iterations"]):
            peaks, core = features(x_axis, H["mach"][k][0])
            fh.write(f"{it}," + ",".join(f"{v:.4f}" for v in peaks + [core]) + "\n")
for f in files[:-2]:
    os.remove(f)
