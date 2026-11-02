# k-omega SST, corrected jet inlet: p0j 7.4945 bar / p0inf 3.247 bar (srp_M2.5)

Runs: `../SST_p0j7.4945_p0inf3.247_mesh{Coarse,Medium,Fine,FineJet,FineJet2,VeryFine}/`. Log: `chain.log`.

## Confirmed assumptions (user, 2026-11-02)
- **The jet inlet boundary condition used until now (p0j = 8.015 bar) was wrong.** The user gave the jet inlet as
  static p = 7.0588 bar, V = 99.11 m/s. The inlet is the plenum face at x = -8.124 mm, r 0-2.8787 mm
  (area 2.6034e-5 m2).
- **Option A:** T0j = 288.15 K as before. Isentropic, gamma 1.4, R 287.058, cp 1004.703:
  - T = 288.15 - 99.11^2 / (2 * 1004.703) = 283.26 K
  - M = 99.11 / sqrt(1.4 * 287.058 * 283.26) = 0.2937
  - **p0j = 705880 * (1 + 0.2 M^2)^3.5 = 749,446 Pa = 7.4945 bar**
  - rho = 8.681 kg/m3, mass flow rho V A = 0.0224 kg/s (full 360 deg; the 8.015 bar runs carried 0.02372)
  - Check: the plenum/throat area ratio ~ (2.88/2.0)^2 = 2.07 gives an isentropic plenum Mach of ~0.294.
- **Boundary type unchanged:** `INLET_TYPE= TOTAL_CONDITIONS`, `MARKER_INLET= ( inlet_jet, 288.15, 749446.0, 1, 0, 0 )`.
  The choked throat sets the mass flow, so the inlet static p and V come out close to, not exactly at, the given values.
- **Everything else as the SST study** (`../SST_p0j8.015_p0inf3.247_tools/README.md`): freestream, gas, walls,
  numerics, SST V2003m with no compressibility correction, freestream TI 1 % / mu_t/mu 1, jet TI 2 % / mu_t/mu 10, meshes.
- **Start:** each mesh restarts from its own SST 8.015 bar final solution (coarse t4_cap3, others t3_cap3). Same mesh,
  so no interpolation; no Phase 1.
- **Stages per mesh:** t1 cap 0.5 (5000, CFL from 0.1), t2 cap 1 (25000, from 0.5), t3 cap 3 (20000, from 1), floor 0.1,
  x1.01 ramp. The coarse configs are the meshMedium ones with the coarse mesh file.
- **Nothing stops on convergence** (`CONV_RESIDUAL_MINVAL= -20`); the watchers only log.
- Order: coarse, medium, fine, fineJet, fineJet2, veryFine. Then the MASS_FLOW runs
  (`../SST_rhoj8.681_Uj99.11_p0inf3.247_tools/`).

## Tooling
- `setup_p0j7.4945.sh` built the configs from the SST 8.015 bar ones; only the inlet line, header and t1 source change.
- `test_inlets.sh`: short checks of this case and of the MASS_FLOW case before launch (log `test_inlets.log`).
- `run_p0j7.4945_chain.sh`: the chain, restartable. Per mesh, `analysis/conv_lines.*` (t2 vs t3) and
  `analysis/inlet_lines.*` (8.015 bar final vs 7.4945 bar t3). At the end, `../grid_compare/sst_p0j7.4945_grid_lines.*`
  and `../grid_compare/box_SST_p0j7.4945_0.025mm.npz`.
- `watch_fj.sh`, `box_snap.py`, `box_change.py`, `plot_conv_roi.py`, `stage_summary.py`: copies of the SST study tools.
- The run-tree copy of this folder was copied from the cfd-tools clone; the commit is in `COMMIT`.
