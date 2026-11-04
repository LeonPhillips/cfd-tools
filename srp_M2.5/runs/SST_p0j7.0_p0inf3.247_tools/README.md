# k-omega SST, jet total pressure p0j 7.0 bar / p0inf 3.247 bar (srp_M2.5)

Runs: `../SST_p0j7.0_p0inf3.247_mesh{Coarse,Medium,Fine,VeryFine}/`. Log: `chain.log`.

## Confirmed assumptions (user, 2026-11-04)
- **Jet inlet:** total conditions, **p0j = 7.0 bar** (700000 Pa), T0j = 288.15 K:
  `INLET_TYPE= TOTAL_CONDITIONS`, `MARKER_INLET= ( inlet_jet, 288.15, 700000.0, 1, 0, 0 )`.
- **Expected inlet values**, scaled from the 7.4945 bar run by 7.0 / 7.4945 = 0.934, since the choked throat fixes the
  plenum Mach:
  - static p ≈ 6.60 bar (7.4945 bar run: 7.0665)
  - V ≈ 98.0 m/s and T ≈ 283.4 K, both unchanged
  - rho ≈ 8.11 kg/m3
  - jet mass flow ≈ 0.0207 kg/s (7.4945 bar run: 0.022165)
- **Everything else as the 7.4945 bar study** (`../SST_p0j7.4945_p0inf3.247_tools/README.md`): freestream, gas, walls,
  numerics, SST V2003m with no compressibility correction, turbulence inflow, meshes.
- **Meshes:** coarse, medium, fine and veryFine (`mesh/srp_M2.5_meshVeryFine.su2`); no fineJet meshes.
  - veryFine has no 7.4945 bar configs, so it uses the meshFine ones with the veryFine mesh. The SST 8.015 bar veryFine
    configs differ from its fine ones only in the mesh file.
- **Start:**
  - Coarse restarts from the **7.4945 bar coarse final solution** (`t3_cap3`, same mesh).
  - Every other mesh is interpolated from the previous mesh's new 7.0 bar t3 (`0_interp`, 1 iteration): medium <- coarse,
    fine <- medium, veryFine <- fine.
- **Stages per mesh:** t1 cap 0.5 (5000, CFL from 0.1), t2 cap 1 (25000, from 0.5), t3 cap 3 (20000, from 1), floor 0.1,
  x1.01 ramp. Nothing stops on convergence (`CONV_RESIDUAL_MINVAL= -20`).
- **Convergence:** the live monitor (`monitor/watch_conv.sh`, every 1000 iterations, three figures per stage).
  - The Mach-disk probes are re-placed before each stage from its starting solution (`monitor/probes_md.py`).
- **Analysis:**
  - **After each mesh:** `analysis/inlet_lines.*`, Mach and density on the monitor lines, comparing the 7.4945 bar t3 with
    the 7.0 bar t3. veryFine is compared with the SST 8.015 bar veryFine t3 instead. There is no t2-vs-t3 comparison
    (user).
  - **At the end:** `../grid_compare/sst_p0j7.0_grid_{lines,features}.{png,csv}`. These are the monitor's `conv_lines`
    and `conv_features` plots with one curve or point per mesh (each mesh's final t3 snapshot) instead of per iteration.

## Tooling
- `setup_p0j7.0.sh` builds the configs and copies this folder, `monitor/` and `grid_compare/feature_positions.py` into the
  run tree (`COMMIT` holds the commit).
- `test_p0j7.0.sh`: short checks before launch (coarse t1 for 2000 iterations, the veryFine interpolation, the grid
  plots). Output goes to `../SST_p0j7.0_p0inf3.247_test/`.
- `run_p0j7.0_chain.sh`: the chain. It is restartable, and `MESHES="..."` limits the meshes.
