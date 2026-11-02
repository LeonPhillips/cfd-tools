"""Convergence plot and oscillation check for one phase of SA_p0j8.015_p0inf3.247.

Usage (on a compute/desktop node): python plot_history.py <phase folder> [window]
Writes <phase folder>/convergence.png and prints, over the last `window` iterations
(default 2000), the mean, peak-to-peak band (% of mean) and the number of zero crossings of the
detrended signal for CD and each surface mass flow. A steady, converging solution has a small
band that shrinks with time. A persistent band with many regular crossings suggests the flow is
oscillating (possibly physical unsteadiness).
"""
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

folder = sys.argv[1].rstrip("/")
window = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
rows = list(csv.reader(open(f"{folder}/history.csv")))
head = [c.strip().strip('"') for c in rows[0]]
data = np.array([[float(v) for v in r] for r in rows[1:]])
col = {h: data[:, i] for i, h in enumerate(head)}
it = col["Inner_Iter"]
flows = [h for h in head if h.startswith("Avg_Massflow(")]

print(f"{folder}: iterations {int(it[0])}-{int(it[-1])}")
print(f"rms[Rho] now {col['rms[Rho]'][-1]:.3f}, relrms[Rho] now {col['relrms[Rho]'][-1]:.3f} "
      f"(min {col['relrms[Rho]'].min():.3f})")
w = it >= it[-1] - window
for name in ["CD"] + flows:
    y = col[name][w]
    mean = y.mean()
    band = 100 * (y.max() - y.min()) / abs(mean) if mean else float("nan")
    resid = y - np.polyval(np.polyfit(it[w], y, 1), it[w])
    crossings = int(np.sum(np.diff(np.sign(resid)) != 0))
    print(f"  {name:32s} mean {mean: .5e}  band {band:7.3f} %  crossings {crossings}")

fig, ax = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
for k in ["rms[Rho]", "rms[RhoU]", "rms[RhoE]", "rms[nu]"]:
    ax[0, 0].plot(it, col[k], lw=0.8, label=k)
ax[0, 0].plot(it, col["relrms[Rho]"], "k--", lw=1, label="relrms[Rho] (drop from peak)")
ax[0, 0].set_ylabel("log10 residual"); ax[0, 0].legend(fontsize=8); ax[0, 0].grid(alpha=0.3)
ax[0, 1].plot(it, col["CD"], lw=0.8); ax[0, 1].set_ylabel("CD (ref area pi*(9 mm)^2)")
ax[0, 1].grid(alpha=0.3)
for k in flows:
    ax[1, 0].plot(it, col[k], lw=0.8, label=k.replace("Avg_Massflow", ""))
ax[1, 0].set_ylabel("mass flow [kg/s, 360 deg]"); ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=0.3)
ax[1, 1].semilogy(it, col["Avg CFL"] if "Avg CFL" in col else col["AvgCFL"], lw=0.8)
ax[1, 1].set_ylabel("average CFL"); ax[1, 1].grid(alpha=0.3)
for a in ax[1]:
    a.set_xlabel("iteration")
mass_in = col["Avg_Massflow(inlet_freestream)"][-1] + col.get("Avg_Massflow(inlet_jet)", np.zeros(1))[-1]
imb = 100 * abs(mass_in + col["Avg_Massflow(outlet)"][-1]) / mass_in if "Avg_Massflow(outlet)" in col else float("nan")
cdw = col["CD"][it >= it[-1] - 10000]
fig.suptitle(f"SA_p0j8.015_p0inf3.247 / {folder.split('/')[-1]}   iteration {int(it[-1])}   CD {col['CD'][-1]:.5f}   "
             f"CD band last 10000: {100 * (cdw.max() - cdw.min()) / abs(cdw.mean()):.2f} % (target < 1)   "
             f"mass imbalance {imb:.2f} % (target < 0.5)", fontsize=10)
fig.tight_layout()
fig.savefig(f"{folder}/convergence.png", dpi=110)
print("saved", f"{folder}/convergence.png")
