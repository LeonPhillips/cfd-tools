# Case SA_p0j8.015_p0inf3.247 on geometry/mesh v4 (srp_M2.5)

Case SA_p0j8.015_p0inf3.247 (its confirmed assumptions are in `README_case_assumptions.md`; the first runs on mesh v2 were in `../SA_p0j8.015_p0inf3.247/`, deleted on 2026-10-28), run fresh on:

- **Geometry v4** `../../geometry/su2_shortDiff_wFillet_v4.STEP` (user, 2026-10-27): the duct ends at a vertical outlet
  at x = -308.1 mm (r 0-55 mm), the base corner fillet is 0.5 mm (was 0.1 mm), the slip wall rises r 45 -> 55 mm over
  250 mm. Marker map: `../../geometry/edge_map_v4.png`.
- **Mesh v4** `../../mesh/srp_M2.5_coarse_v4.su2`: 235.3k cells, first cell 1.8e-7 m, 17 wall layers (18/19 do not mesh
  on this geometry), 1.5 mm cells in the whole wake.
- **Outlet:** `MARKER_SUPERSONIC_OUTLET` (user, 2026-10-27) instead of `MARKER_OUTLET` at p_inf. Phase 1 is checked
  for Mach >= 1 on the whole outlet before Phase 2.
- **Phase 1** (`1_freestream/`): as v2 Phase 1 (CFL 1 -> 50 x1.2, 2-order drop of rms[Rho] from its peak once fully
  second order, cap 20000).
- **Phase 2, first try** (`2_jet/`, deleted 2026-10-28; CFL cycled and it did not converge): jet stepped on to 8.015 bar from the Phase 1 restart; first order for 500 iterations, second
  order by 1000; CFL starts at 1, x1.01 per iteration, cap 50, floor 0.5 (user), LINEAR_SOLVER_ITER 20; up to 100000
  iterations. Stopped by `watch_stop.sh` when CD varies < 1 % of |CD| over 10000 iterations
  and mass balance < 0.5 % (user), or on OSCILLATE / CORRUPT / DIVERGE.

Run: `run_v4.sh` (Phase 1 and `check_phase1.py`; its Phase 2 part now exits, since `2_jet/` was deleted). Log: `run_v4.log`.

## Phase 2 as run (`2_jet_staged/`)

Each stage restarts from the one before it. CFL floor 0.1, CFL start 1.0.

| Stage | CFL cap | Restarts from | Result |
|---|---|---|---|
| `s1_cap0.5` | 0.5 | `1_freestream` | 5000 iterations, clean, CD 1.38 -> 0.88 |
| `s2_cap1` | 1 | `s1_cap0.5` | 5000 iterations, CFL started cycling, CD 0.27 |
| `s3_cap2`, `s4_cap5` | 2, 5 | s2, s3 | deleted 2026-10-28 (s3 cycled, s4 unstable) |
| `s5_cap1_long` | 1 | `s2_cap1` | 200000 iterations, stable after ~20000; imbalance < 0.1 %, CD -0.184 still creeping |
| `s6_cap3` | 3 | `s5_cap1_long` | stopped on CONVERGED at 188240 (CD band 1.00 % over 10000, imbalance 0.04 %), CD -0.127, still drifting slowly |

Scripts: `run_staged.sh` (s1, s2), `run_s5.sh`, `run_s6.sh`. The watcher (`watch_stop.sh`) and its helpers
(`check_run.py`, `plot_history.py`) and the field tools (`inspect_field.py`, `plot_field.py`) were moved here from the old
folder on 2026-10-28. `run_su2.sh` now only has the `v4_phase1` stage.
