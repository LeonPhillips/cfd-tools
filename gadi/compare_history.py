"""Compare two SU2 history.csv files row by row (same config, e.g. an M3 run and its Gadi rerun).

Usage: python compare_history.py <reference history.csv> <new history.csv>
Prints, for the shared iterations, the largest differences in the residuals, CFL, CD, jet mass flow and probes, plus the
mean wall time per iteration of each file (column "Time(sec)", the second one: time per iteration).
"""
import csv
import sys


def load(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    head = [h.strip().strip('"') for h in rows[0]]
    data = [[float(v) for v in r] for r in rows[1:] if len(r) == len(head)]
    return head, data


def main(ref, new):
    h1, a = load(ref)
    h2, b = load(new)
    it1, it2 = h1.index("Inner_Iter"), h2.index("Inner_Iter")
    rows_a = {int(r[it1]): r for r in a}
    shared = [(rows_a[int(r[it2])], r) for r in b if int(r[it2]) in rows_a]
    print(f"shared iterations: {len(shared)} ({int(shared[0][1][it2])}..{int(shared[-1][1][it2])})")
    cols = [c for c in h2 if c in h1 and (c.startswith("rms[") or c in ("Avg CFL", "CD", "Avg_Massflow")
                                          or c.startswith("p") and len(c) == 2)]
    print(f"{'column':>14} {'max |diff|':>12} {'ref at last':>14} {'new at last':>14}")
    for c in cols:
        i, j = h1.index(c), h2.index(c)
        d = max(abs(x[i] - y[j]) for x, y in shared)
        print(f"{c:>14} {d:12.3e} {shared[-1][0][i]:14.6g} {shared[-1][1][j]:14.6g}")
    for name, h, rows in (("ref", h1, [x for x, _ in shared]), ("new", h2, [y for _, y in shared])):
        k = [n for n, c in enumerate(h) if c == "Time(sec)"]
        if len(k) > 1:
            t = [r[k[1]] for r in rows[10:]]
            print(f"{name}: mean time per iteration {sum(t) / len(t):.4f} s (iterations 10 on)")


if __name__ == "__main__":
    main(*sys.argv[1:3])
