"""Auto-stop test for the CT sweep (user, 2026-11-27: "Auto-stop, cap 250k"), from the live monitor's small files.

PASS when, over the last WINDOW (10,000) iterations:
  - Mach disk x, bow shock x on the axis, at 1 D and at 2 D each move by less than TOL_MM (0.01 mm) (max - min of the
    snapshots in the window, including the snapshot at its start), and
  - the line mean |dMach| per 1000 iterations is below TOL_DMACH (0.002) on every line, in every snapshot of the window.
The stagnation point is not used (it wanders in weak-jet cases, e.g. CT 0.163). Missing or NaN values fail.
The run driver (SST_SLAU2-VW_meshFine_tools/run_ct_autostop.sh) applies it only from 50,000 cumulative iterations, and
stops at 250,000 whatever the result. CD and residuals are not used (convergence.md).

  python autostop_check.py <stage folder> [<previous stage folder>]
      prints one line "PASS ..." / "FAIL ..." / "NODATA ..."; exit 0 = PASS, 1 = FAIL, 2 = not enough data.
  python autostop_check.py --replay <run folder> [min_cum=50000]
      replays the test on every 10k-iteration window of a finished chain (for checking the criterion).
"""
import csv, glob, math, os, sys

WINDOW, TOL_MM, TOL_DMACH = 10000, 0.01, 0.002
FEAT = ["mach_disk_x", "bow_x_axis", "bow_x_1D", "bow_x_2D"]


def rows(path):
    if not os.path.isfile(path):
        return [], []
    r = list(csv.reader(open(path)))
    if not r:
        return [], []
    return r[0], [x for x in r[1:] if x]


def num(v):
    try:
        f = float(v)
        return f if math.isfinite(f) else float("nan")
    except ValueError:
        return float("nan")


def stage_series(d, off=0):
    """[(cumulative iter, {feature: value}, max line mean dMach or None)] for one stage folder."""
    h, R = rows(f"{d}/features.csv")
    hl, L = rows(f"{d}/line_change.csv")
    lm = {}
    if hl:
        cols = [i for i, c in enumerate(hl) if c.startswith("mean_")]
        for x in L:
            lm[int(num(x[0]))] = max(num(x[i]) for i in cols) if all(i < len(x) for i in cols) else float("nan")
    out = []
    for x in R:
        it = int(num(x[0]))
        out.append((it + off, {k: num(x[h.index(k)]) if k in h and h.index(k) < len(x) else float("nan") for k in FEAT},
                    lm.get(it)))
    return out


def test(series):
    """series: [(iter, feats, lmax)] in order; tests the last WINDOW iterations."""
    if not series:
        return 2, "NODATA no snapshots"
    end = series[-1][0]
    win = [s for s in series if s[0] >= end - WINDOW - 1]
    if win[-1][0] - win[0][0] < WINDOW - 1001:
        return 2, f"NODATA window spans only {win[-1][0] - win[0][0]} iterations"
    msg, ok = [], True
    for k in FEAT:
        v = [s[1][k] for s in win]
        rng = max(v) - min(v) if all(math.isfinite(a) for a in v) else float("nan")
        good = math.isfinite(rng) and rng < TOL_MM
        ok &= good
        msg.append(f"{k} {rng:.4f}{'' if good else '!'}")
    lm = [s[2] for s in win[1:]]   # the window's first snapshot is the reference; its dMach belongs to the window before
    lmax = max(lm) if lm and all(a is not None and math.isfinite(a) for a in lm) else float("nan")
    good = math.isfinite(lmax) and lmax < TOL_DMACH
    ok &= good
    msg.append(f"line mean dMach max {lmax:.5f}{'' if good else '!'}")
    return (0 if ok else 1), f"{'PASS' if ok else 'FAIL'} it {end}, last {win[-1][0] - win[0][0]} it: " + ", ".join(msg) \
        + f"  (limits {TOL_MM} mm, {TOL_DMACH})"


def cfg_iter(d):
    for c in glob.glob(f"{d}/*.cfg"):
        for line in open(c):
            if line.startswith("ITER="):
                return int(line.split("=")[1].split()[0])
    return 0


if __name__ == "__main__":
    if sys.argv[1] == "--replay":
        run = sys.argv[2]
        mn = int(sys.argv[3].split("=")[1]) if len(sys.argv) > 3 else 50000
        S, off = [], 0
        for d in sorted(glob.glob(f"{run}/2_jet_staged/t*")):
            S += stage_series(d, off)
            off += cfg_iter(d)
        first = None
        for i in range(len(S)):
            if S[i][0] + 1 >= mn and round((S[i][0] + 1) / 1000) % (WINDOW // 1000) == 0:
                rc, m = test(S[: i + 1])
                print(m)
                if rc == 0 and first is None:
                    first = S[i][0]
        print(f"first PASS at cumulative iteration {first}" if first is not None else "never PASSES")
        sys.exit(0)
    d = sys.argv[1]
    prev = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] != "-" else None
    S = (stage_series(prev, -cfg_iter(prev))[-1:] if prev else []) + stage_series(d)
    rc, m = test(S)
    print(m)
    sys.exit(rc)
