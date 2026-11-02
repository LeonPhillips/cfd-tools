# Case SA_p0j8.015_p0inf3.247 on meshMedium (srp_M2.5)

Same case, settings and confirmed assumptions as `../SA_p0j8.015_p0inf3.247_meshCoarse/` (see its `README.md` and
`README_case_assumptions.md`), on `../../mesh/srp_M2.5_meshMedium.su2`: meshCoarse plus 0.025 mm cells at
x 0.5-20 mm, r 0-9.5 mm. 468,874 cells.

## `0_interp/` (2026-10-28)
- The meshCoarse s6_cap3 solution (iteration 188239, CD -0.12698) interpolated onto meshMedium with SU2_CFD's own
  restart interpolation (user asked for SU2's method). `interp.cfg` = `s6_cap3.cfg` with the medium mesh,
  `SOLUTION_FILENAME` = the s6 restart, 1 iteration at CFL 1e-8 (no adapt); `interp.sh`, 24 ranks on m3f104.
  SU2 saw 233,282 restart points vs 464,934 mesh points and did its "recursive nearest neighbor interpolation".
- Result `restart.dat` (start point for a medium run), `flow.vtu`, `surface.vtu`.
- Check against s6 (coarse field resampled linearly onto the medium points, `diff_field.py`, `bad_points.py`):
  CD identical; mean |dMach| 0.003 in the refined box, 0.0001 elsewhere; 99th percentile 0.04.
  104 of 464,934 points differ by > 0.5 in Mach. Most sit on the bow shock and shear layers (sub-cell shift of a
  discontinuity), but there are real artefacts: 6 points on the axis at x 16.71-16.73 mm, r < 0.13 mm (just ahead
  of the bow shock at 16.06 mm) with Mach ~0 instead of 2.5, and 18 points inside the shock layer holding exact
  freestream values. They should wash out in the first iterations; start the medium run at a low CFL.
- Bow shock on the axis: 16.03-16.06 mm on both meshes (axis mesh points). `inspect_field.py` reports 16.75 mm for
  the interpolated field only because of the 6 bad points.
- Plot: `compare_coarse-s6_vs_medium-interp.png` (`plot_field.py`).

## Phase 2 on meshMedium (`2_jet_staged/`, 2026-10-28)
Fully second order, CFL floor 0.1, log-only watcher, 24 ranks on m3f104. Scripts `run_medium.sh` (m1, m2), `run_m3.sh`;
results `analyse_medium.sh` -> `analysis/`.

| Stage | CFL | Restarts from | Iterations | End state |
|---|---|---|---|---|
| `m1_cap0.5` | 0.1 -> cap 0.5 | `0_interp` | 5000 (6.6 min) | CD -0.12560, rms[Rho] -5.20, 0 CFL drops |
| `m2_cap1` | 0.5 -> cap 1 | m1 | 26673 of 45000 (stopped by the user with SIGTERM) | CD -0.12481, rms[Rho] -6.08, 0 CFL drops |
| `m3_cap3` | 1 -> cap 3 | m2 | 20000 (26 min, 0.079 s/it) | CD -0.12439, rms[Rho] -6.58, imbalance 0.03 %, 113 CFL drops |

- m3: CFL held 3 until ~6000 iterations, then cycled 3 -> 0.1 about every 500 iterations (as coarse s6), without harm:
  residuals kept falling, mass flows flat (jet 0.023720 kg/s, outlet 1.88964). CD still rising slowly and steadily
  (-0.12481 -> -0.12439, ~0.0001 per 5000 iterations); CD varies 0.15 % over the last 10000 (criterion < 1 % met).
- vs coarse s6 (CD -0.12698, itself still rising): medium CD 2.0 % higher (less negative); bow shock on the axis
  16.03-16.08 mm on both; the bow shock off the axis (r 5-10 mm) is thinner/sharper on medium.
- The interpolation artefacts are gone (axis ahead of the bow shock clean).
- Jet core on the axis (x 4.5-8 mm): the medium axis Mach shows a fast zig-zag (amplitude up to ~0.3) that the coarse
  field does not; the shock at x ~4.5 mm drops to Mach ~0.3 (coarse ~1.0). Not yet diagnosed.
- Plots: `analysis/m{1,2,3}_convergence.png`, `analysis/coarse-s6_vs_medium-m3.png`.
