"""python conv_update.py <stage folder> [previous stage folder]: called by watch_conv.sh after snap.py (user, 2026-11-04).
For the stage's snapshots (<stage>/snaps/flow_*.npz) it appends to three small logs, one row per new snapshot:
- box_change.csv: Mach/density change over the map region x 0-16, r 0-6 mm, per 1000 iterations (same columns as the
  old jet-box log: iter, n_iter, mean_dMach, p99_dMach, mean_drho_pct, top1_mean_dMach).
- line_change.csv: per convergence line, mean / p99 / max |dMach| per 1000 iterations against the previous snapshot.
- features.csv: key-feature positions of each snapshot (feature_positions.py definitions), mm.
Then it deletes all but the latest two snapshots. The first snapshot of a stage is compared with the previous stage's
latest snapshot (the solution it restarted from), if that one has the same format.
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "grid_compare"))
import feature_positions as fp   # noqa: E402

MAP = (0.0, 16.0, 0.0, 6.0)   # x0, x1, r0, r1 mm: difference map and box_change region (user, 2026-11-04)
FEAT = ["mach_disk_x", "mach_disk_subsonic_edge_r", "mach_disk_edge_r", "shock_crossing_x", "stagnation_x",
        "bow_x_axis", "bow_x_1D", "bow_x_2D"]
folder = sys.argv[1].rstrip("/")
prev = sys.argv[2].rstrip("/") if len(sys.argv) > 2 and sys.argv[2] != "-" else None
key = lambda f: int(f.rsplit("_", 1)[-1][:-4])
snaps = sorted(glob.glob(f"{folder}/snaps/flow_*.npz"), key=key)


def done(path, head):
    """Iterations already logged in path (created with head if missing)."""
    if not os.path.exists(path):
        open(path, "w").write(head + "\n")
        return set()
    return {int(l.split(",")[0]) for l in open(path).read().split("\n")[1:] if l}


def crop(z, a):
    x, r = z["x_mm"], z["r_mm"]
    return a[np.ix_((r >= MAP[2]) & (r <= MAP[3]), (x >= MAP[0]) & (x <= MAP[1]))]


pairs = list(zip(snaps[:-1], snaps[1:]))
if prev and snaps:
    old = sorted(glob.glob(f"{prev}/snaps/flow_*.npz"), key=key)
    if old:
        with np.load(old[-1]) as zo, np.load(snaps[0]) as zn:
            same = "line_Mach" in zo.files and zo["Mach"].shape == zn["Mach"].shape
        if same:
            pairs.insert(0, (old[-1], snaps[0]))

if snaps:
    with np.load(snaps[0]) as z0:
        names = [str(s) for s in z0["line_names"]]
    tag = [n.replace("r=", "r").replace(" ", "_") for n in names]
    bhead = "iter,n_iter,mean_dMach,p99_dMach,mean_drho_pct,top1_mean_dMach"
    lhead = "iter,n_iter," + ",".join(f"{s}_{t}" for t in tag for s in ("mean", "p99", "max"))
    bdone, ldone = done(f"{folder}/box_change.csv", bhead), done(f"{folder}/line_change.csv", lhead)
    with open(f"{folder}/box_change.csv", "a") as fb, open(f"{folder}/line_change.csv", "a") as fl:
        for a, b in pairs:
            ib = key(b)
            if ib in bdone and ib in ldone:
                continue
            za, zb = np.load(a), np.load(b)
            n = ib - (0 if (prev and a.startswith(prev)) else key(a))
            k = 1000 / max(n, 1)
            if ib not in bdone:
                dm = np.abs(crop(zb, zb["Mach"]) - crop(za, za["Mach"]))
                dr = 100 * np.abs(crop(zb, zb["Density"]) - crop(za, za["Density"])) / crop(zb, zb["Density"])
                ok = np.isfinite(dm) & np.isfinite(dr)
                p99 = np.percentile(dm[ok], 99)
                fb.write(f"{ib},{n},{dm[ok].mean() * k:.6e},{p99 * k:.6e},{dr[ok].mean() * k:.6e},"
                         f"{dm[ok][dm[ok] >= p99].mean() * k:.6e}\n")
            if ib not in ldone:
                d = np.abs(zb["line_Mach"].astype(float) - za["line_Mach"].astype(float))
                vals = []
                for row in d:
                    v = row[np.isfinite(row)]
                    vals += [v.mean() * k, np.percentile(v, 99) * k, v.max() * k] if len(v) else [np.nan] * 3
                fl.write(f"{ib},{n}," + ",".join(f"{v:.6e}" for v in vals) + "\n")

    fdone = done(f"{folder}/features.csv", "iter," + ",".join(FEAT) + ",flags")
    with open(f"{folder}/features.csv", "a") as ff:
        for s in snaps:
            if key(s) in fdone:
                continue
            with np.load(s) as z:
                fp.set_grid(z["x_mm"], z["r_mm"])
                o = fp.features(z["Mach"].astype(float), z["Density"].astype(float), z["Velocity_x"].astype(float))
            ff.write(f"{key(s)}," + ",".join(f"{o.get(c, np.nan):.4f}" for c in FEAT) + f",{o['flags']}\n")

for f in snaps[:-2]:
    os.remove(f)
