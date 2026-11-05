# k-omega SST, jet total pressure p0j 6.6 bar / p0inf 3.247 bar (srp_M2.5)

Runs: `../SST_p0j6.6_p0inf3.247_mesh{Coarse,Medium,Fine,VeryFine}/`. Log: `chain.log`.

## Confirmed assumptions (user, 2026-11-05)
- **Jet inlet:** total conditions, **p0j = 6.6 bar** (660000 Pa), T0j = 288.15 K:
  `INLET_TYPE= TOTAL_CONDITIONS`, `MARKER_INLET= ( inlet_jet, 288.15, 660000.0, 1, 0, 0 )`.
- **Expected inlet values**, scaled from the 7.0 bar run by 6.6 / 7.0 = 0.943, since the choked throat fixes the
  plenum Mach:
  - static p ≈ 6.22 bar (7.0 bar run: ≈ 6.60)
  - V ≈ 98 m/s and T ≈ 283.4 K, both unchanged
  - rho ≈ 7.65 kg/m3 (7.0 bar run: ≈ 8.11)
  - jet mass flow ≈ 0.01952 kg/s (7.0 bar run: 0.020700)
- **Everything else as the 7.0 bar study** (`../SST_p0j7.0_p0inf3.247_tools/README.md`), which is the 7.4945 bar study
  with only p0j changed: freestream, gas, walls, numerics, SST V2003m with no compressibility correction, turbulence
  inflow, meshes.
- **Meshes:** coarse, medium, fine and veryFine. The configs are the 7.0 bar ones (veryFine: meshFine settings with
  `mesh/srp_M2.5_meshVeryFine.su2`).
- **Start:**
  - Coarse restarts from the **7.0 bar coarse final solution** (`t3_cap3`, same mesh).
  - Every other mesh is interpolated from the previous mesh's new 6.6 bar t3 (`0_interp`, 1 iteration): medium <- coarse,
    fine <- medium, veryFine <- fine.
- **Stages per mesh:** t1 cap 0.5 (5000, CFL from 0.1), t2 cap 1 (25000, from 0.5), t3 cap 3 (20000, from 1), floor 0.1,
  x1.01 ramp. Nothing stops on convergence (`CONV_RESIDUAL_MINVAL= -20`).
- **Convergence:** the live monitor (`monitor/watch_conv.sh`, every 1000 iterations, three figures per stage).
  - The Mach-disk probes are re-placed before each stage from its starting solution (`monitor/probes_md.py`).
- **Analysis:**
  - **After each mesh:** `analysis/inlet_lines.*`, Mach and density on the monitor lines, comparing the 7.0 bar t3 with
    the 6.6 bar t3 on the same mesh (all four meshes).
  - **At the end:** `../grid_compare/sst_p0j6.6_grid_{lines,features}.{png,csv}`.
  - Convergence plots: `monitor/plot_chain.py` per mesh (the default).
- **No separate test run** (user, 2026-11-05): the tooling and meshes are the ones that ran the 7.0 bar chain.

## Tooling
- `setup_p0j6.6.sh` builds the configs and copies this folder, `monitor/` and `grid_compare/feature_positions.py` into the
  run tree (`COMMIT` holds the commit).
- `run_p0j6.6_chain.sh`: the chain. It is restartable, and `MESHES="..."` limits the meshes.
