# k-omega SST, MASS_FLOW jet inlet: rho 8.681 kg/m3, V 99.11 m/s / p0inf 3.247 bar (srp_M2.5)

Runs: `../SST_rhoj8.681_Uj99.11_p0inf3.247_mesh{Coarse,FineJet}/`. Log: `chain.log`. Runs after the
`SST_p0j7.4945_p0inf3.247` chain (same tmux session).

## Confirmed assumptions (user, 2026-11-02)
- Jet inlet from the user: static p = 7.0588 bar, V = 99.11 m/s, with T0 = 288.15 K (option A):
  T = 283.26 K, **rho = 705880 / (287.058 * 283.26) = 8.681 kg/m3**, mass flow ~0.0224 kg/s (full 360 deg).
  Derivation in `../SST_p0j7.4945_p0inf3.247_tools/README.md`.
- **Boundary:** `INLET_TYPE= MASS_FLOW`, `MARKER_INLET= ( inlet_jet, 8.681, 99.11, 1, 0, 0 )`. SU2 reads density, then
  velocity magnitude, then direction (checked in `CEulerSolver::BC_Inlet`, v8.5.0). Pressure is extrapolated from the
  interior, so the inlet static p and T float.
- **Meshes: meshCoarse and meshFineJet only** (user).
- **meshCoarse from scratch, with the long staging of the original SST study:**
  - `1_freestream`: nozzle plugged, a fixed 20000 iterations, CFL 0.1 -> cap 0.5, floor 0.1 (as relaunched in the SST study).
  - `2_jet_staged`: t1 cap 0.3 (first order for 500 iterations, second order by 1000), t2 cap 0.5, t3 cap 1, t4 cap 3;
    50000 iterations each, floor 0.1, x1.01 ramp.
- **meshFineJet:** coarse t4 interpolated onto it (`0_interp`, SU2 restart interpolation), then t1 cap 0.5 (5000),
  t2 cap 1 (25000), t3 cap 3 (20000).
- Everything else as the SST study. **Nothing stops on convergence**; the watchers only log.

## Tooling
- `setup_rhoj8.681.sh` built the configs from the SST 8.015 bar coarse and fineJet ones.
- Tests: `../SST_p0j7.4945_p0inf3.247_tools/test_inlets.sh` (tests B and C).
- `run_rhoj8.681_chain.sh`: the chain, restartable. Per mesh, `analysis/conv_lines.*` (coarse t3 vs t4, fineJet t2 vs
  t3) and `analysis/inlet_lines.*` (total conditions p0j 7.4945 bar t3 vs MASS_FLOW, same mesh).
- The run-tree copy of this folder was copied from the cfd-tools clone; the commit is in `COMMIT`.
