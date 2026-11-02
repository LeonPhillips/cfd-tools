"""One-line health check of a running phase. Usage: python check_run.py <phase folder>
Prints the iteration, CD and its Cauchy value, rms[Rho] (last value and max over the last 500
iterations, to show residual spikes), CFL, and flags:
  CORRUPT   the freestream inlet mass flow has moved off its fixed value (1.86533 kg/s) by > 0.1 %,
            which is how the corrupted Phase 2 run first showed up in the history
  OSCILLATE CD over the last 2000 iterations, after removing its linear trend, has >= 6 swings and
            a band > 1 % of max(|CD|, 0.1), i.e. it is cycling rather than drifting (a steady drift
            alone does not trigger it)
  CONVERGED CD varies by < 1 % of |mean CD| over the last 10000 iterations (no 0.1 floor; needs
            >= 10000 iterations in this run) and |mass in - mass out| / mass in < 0.5 %
            (user criterion, 2026-10-27)
  DIVERGE   rms[Rho] is more than 1.5 orders above its mean over the first 50 iterations of this run
            (added for the CFL 100 continuation; earlier CFL 5 / 50 runs rose ~1.5 orders or corrupted)
  CFLCYCLE  the average CFL dropped suddenly (to < 60 % of the previous iteration: a halving or a reset)
            >= 3 times in the last 2000 iterations, checked after iteration 2000 (user, 2026-10-27)
"""
import csv
import sys

import numpy as np

rows = list(csv.reader(open(f"{sys.argv[1].rstrip('/')}/history.csv")))
head = [c.strip().strip('"') for c in rows[0]]
d = np.array([[float(v) for v in r] for r in rows[1:] if len(r) == len(head)])
c = {h: d[:, i] for i, h in enumerate(head)}
it, cd, rho = c["Inner_Iter"], c["CD"], c["rms[Rho]"]
mf = c["Avg_Massflow(inlet_freestream)"]
w5, w20 = it >= it[-1] - 500, it >= it[-1] - 2000
flags = []
if np.max(np.abs(mf[w5] / 1.86533 - 1)) > 1e-3:
    flags.append("CORRUPT")
if len(rho) > 50 and rho[-1] > rho[:50].mean() + 1.5:
    flags.append("DIVERGE")
cfl = c["Avg CFL"] if "Avg CFL" in c else c["AvgCFL"]
cw = cfl[w20]
drops = int(np.sum(cw[1:] < 0.6 * cw[:-1]))
if it[-1] >= 2000 and drops >= 3:
    flags.append("CFLCYCLE")
y = cd[w20]
resid = y - np.polyval(np.polyfit(it[w20], y, 1), it[w20]) if w20.sum() > 10 else y * 0
swings = int(np.sum(np.diff(np.sign(resid)) != 0))
band = 100 * (y.max() - y.min()) / max(abs(y.mean()), 0.1)
dband = 100 * (resid.max() - resid.min()) / max(abs(y.mean()), 0.1)
if swings >= 6 and dband > 1:
    flags.append("OSCILLATE")
mass_in = mf[-1] + c.get("Avg_Massflow(inlet_jet)", np.zeros(1))[-1]
imbalance = 100 * abs(mass_in + c["Avg_Massflow(outlet)"][-1]) / mass_in
w100 = it >= it[-1] - 10000
y10 = cd[w100]
cd_rel = 100 * (y10.max() - y10.min()) / abs(y10.mean())
if w100.sum() >= 10000 and cd_rel < 1 and imbalance < 0.5:
    flags.append("CONVERGED")
cauchy = c.get("Cauchy[CD]", c.get("Cauchy[DRAG]", np.array([np.nan])))[-1]
print(f"it {int(it[-1])}  CD {cd[-1]:.5f}  CD band(2000) {band:.2f}% (detrended {dband:.2f}%)  swings {swings}  cauchy {cauchy:.2e}  "
      f"rms[Rho] {rho[-1]:.3f} (max500 {rho[w5].max():.3f})  CFL {c['Avg CFL'][-1] if 'Avg CFL' in c else c['AvgCFL'][-1]:.2f}  "
      f"CFL drops(2000) {drops}  out {c['Avg_Massflow(outlet)'][-1]:.4f}  CD rel(10000) {cd_rel:.2f}%  imbalance {imbalance:.2f}%  flags {','.join(flags) or 'none'}")
