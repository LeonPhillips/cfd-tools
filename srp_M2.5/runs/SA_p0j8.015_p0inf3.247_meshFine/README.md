# SA, p0j 8.015 bar, p0inf 3.247 bar: meshFine (2026-10-28)

Same case and assumptions as `../SA_p0j8.015_p0inf3.247_meshCoarse/README_case_assumptions.md`.
Mesh: `../../mesh/srp_M2.5_meshFine.su2`, which is meshCoarse plus a 0.0177 mm box at x 0.5-20 mm, r 0-9.5 mm. It has 772,018 cells.

## 0_interp
The meshMedium m3_cap3 solution was moved onto meshFine with SU2's restart interpolation (1 iteration at CFL 1e-8).
50 of 765,737 points came out off by more than 0.5 Mach. These sat at the shocks and washed out in f1-f3.

## 2_jet_staged (run_fine.sh, tmux srp_fine, 15:40-17:22, all exit 0)
| stage | iterations | CFL cap | CD end | rms[Rho] | CFL drops |
|---|---|---|---|---|---|
| f1_cap0.5 | 5000 | 0.5 | -0.12408 | -5.04 | 0 |
| f2_cap1 | 25000 | 1 | -0.12392 | -5.94 | 0 |
| f3_cap3 | 20000 | 3 | -0.12318 | -6.51 | 39 (from ~12900) |

f3 end state:
- CD is still rising, by +0.00016 per 5000 iterations.
- imbalance 0.03 %, jet 0.023721 kg/s
- bow shock on the axis at 16.075 mm

## Grid comparison
- CD: coarse s6 -0.12698, medium m3 -0.12439, fine f3 -0.12318.
- The apparent order is well below 1, so no GCI value is given. CD is still drifting on every mesh.
- Axis Mach minimum after the Mach disk: coarse 1.30, medium 0.87, fine 0.81.
- Axis-node Mach in the jet core has no local extrema on any mesh. The zig-zag in the line plots comes from
  plotting near-axis points.
- Analysis: `analyse_fine.sh`, output in `analysis/`.
