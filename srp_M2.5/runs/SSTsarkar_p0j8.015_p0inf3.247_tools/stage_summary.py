import sys
import numpy as np
for p in sys.argv[1:]:
    h = [c.strip().strip('"') for c in open(p).readline().split(",")]
    d = np.genfromtxt(p, delimiter=",", skip_header=1, invalid_raise=False)
    c = {k: d[:, i] for i, k in enumerate(h)}
    it, cfl, n = c["Inner_Iter"], c["Avg CFL"], len(d)
    jet = c.get("Avg_Massflow(inlet_jet)")
    print(f"== {p.split('/')[-2]}  iterations {n}  s/it {np.median(d[-500:, 3]):.4f}  CFL drops (whole stage) "
          f"{int(np.sum(np.diff(cfl) < -0.25 * cfl[:-1]))}")
    for i in [int(v) for v in np.linspace(0, n - 1, 6)]:
        s = f"  it {int(it[i]):6d}  rms[Rho] {c['rms[Rho]'][i]:7.3f}  rms[k] {c['rms[k]'][i]:7.3f}  rms[w] {c['rms[w]'][i]:6.3f}  CFL {cfl[i]:6.3f}"
        if jet is not None:
            s += f"  jet {jet[i]:.5f}"
        print(s + f"  CD {c['CD'][i]:8.4f}")
    if "p1" in c:
        print("  probes at end:", " ".join(f"{q}={c[q][-1]:.3f}" for q in ("p1", "p2", "p3", "p4", "p5", "p6", "p7")),
              "| change over last 10000:", " ".join(f"{q}={c[q][-1] - c[q][max(0, n - 10001)]:+.3f}" for q in ("p2", "p4", "p7")))
