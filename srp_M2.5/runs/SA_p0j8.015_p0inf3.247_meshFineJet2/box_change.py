"""python box_change.py <stage folder> [previous stage folder]: for each new pair of consecutive jet-box snapshots
(<stage>/snaps/flow_*.npz from box_snap.py) append one row to <stage>/box_change.csv, then delete all but the latest two
snapshots (user, 2026-10-29: only the latest pair is kept, for the difference map in plot_conv_roi.py).
The first snapshot of a stage is compared with the latest snapshot of the previous stage (the restart it started from).
Row: iteration, iterations since the previous snapshot, mean / p99 |dMach|, mean |drho|/rho [%], and the mean |dMach| of
the points at or above p99 (the most-changing 1 %; user, 2026-10-30), all per 1000 iterations. A log written before that
column existed is rewritten with an empty (nan) value in its old rows.
"""
import glob
import os
import sys

import numpy as np

folder = sys.argv[1].rstrip("/")
prev = sys.argv[2].rstrip("/") if len(sys.argv) > 2 and sys.argv[2] != "-" else None
key = lambda f: int(f.rsplit("_", 1)[-1][:-4])
snaps = sorted(glob.glob(f"{folder}/snaps/flow_*.npz"), key=key)
out = f"{folder}/box_change.csv"
HEAD = "iter,n_iter,mean_dMach,p99_dMach,mean_drho_pct,top1_mean_dMach\n"
done = set()
if os.path.exists(out):
    lines = [l for l in open(out).read().split("\n") if l]
    if lines and lines[0] != HEAD.strip():   # old 5-column log: add the new column, nan for past rows
        open(out, "w").write(HEAD + "".join(l + ",nan\n" for l in lines[1:]))
    done = {int(l.split(",")[0]) for l in lines[1:]}
else:
    open(out, "w").write(HEAD)
pairs = list(zip(snaps[:-1], snaps[1:]))
if prev and snaps:
    old = sorted(glob.glob(f"{prev}/snaps/flow_*.npz"), key=key)
    if old:
        pairs.insert(0, (old[-1], snaps[0]))
with open(out, "a") as fh:
    for a, b in pairs:
        ib = key(b)
        if ib in done:
            continue
        za, zb = np.load(a), np.load(b)
        n = ib - (0 if a.startswith(prev or "\0") else key(a))
        dm = np.abs(zb["Mach"] - za["Mach"]); dr = 100 * np.abs(zb["Density"] - za["Density"]) / zb["Density"]
        ok = np.isfinite(dm) & np.isfinite(dr)
        k = 1000 / max(n, 1)
        p99 = np.percentile(dm[ok], 99)
        fh.write(f"{ib},{n},{dm[ok].mean() * k:.6e},{p99 * k:.6e},{dr[ok].mean() * k:.6e},"
                 f"{dm[ok][dm[ok] >= p99].mean() * k:.6e}\n")
for f in snaps[:-2]:
    os.remove(f)
