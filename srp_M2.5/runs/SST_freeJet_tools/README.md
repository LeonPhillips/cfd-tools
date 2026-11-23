# SST free jet: the srp_M2.5 nozzle and model without the counterflow (Gadi)

Set up on 2026-11-23. The same model, nozzle and jet inlet as the srp_M2.5 runs. The jet exhausts into still air
at a pressure just below its exit pressure, so it is slightly underexpanded. The case is run on three meshes, from
coarse to fine. Each mesh starts from the previous mesh's final solution, interpolated by SU2.

## Confirmed assumptions (user, 2026-11-23: "use the remainder of the defaults")

- **Jet inlet:** the same as the srp_M2.5 runs. `inlet_jet` (the plenum face) is a total-conditions inlet with
  p0j 6.6 bar and T0j 288.15 K, flow in +x, TI 2 % and mu_t/mu 10.
- **Ambient:** still air at **pa = 101,325 Pa and Ta = 288.15 K**. Air is an ideal gas with Sutherland viscosity.
  - The nozzle area ratio is (2.389/2.000)² = 1.427. The 1-D design exit Mach is 1.79, pe/p0 = 0.177, and the design
    NPR is 5.65.
  - At p0j 6.6 bar the 1-D pe = 1.169 bar, so **pe/pa ≈ 1.15 and NPR = 6.51**.
  - The CFD exit plane is not uniform (the axis Mach at the exit is 1.60), so `exit_plane.py` measures the area- and
    mass-averaged pe from the 6.6 bar meshFine t3 solution. The result is in `srp_M2.5_freeJet/exit_plane.txt`.
- **Geometry:** the v4 model and nozzle walls are unchanged (no-slip, adiabatic). The tunnel boundaries are replaced
  by an open domain: x from −308.1 mm (the v4 rear boundary) to 240 mm (~50 exit diameters), r from 0 to 120 mm
  (~25 exit diameters). The geometry and mesh script is `srp_M2.5/mesh/mesh_freejet.py`.
- **Outer boundaries (option A):**
  - `ambient_side` (r = 120 mm) and `ambient_rear` (x = −308.1 mm) are total-condition inlets at p0 = pa,
    T0 = Ta, with flow normal into the domain. They supply the entrained air.
  - `outlet` (x = 240 mm) is a static-pressure outlet at pa.
- **Meshes:** each level is a constant sqrt 2 finer than the one before. The fine level is the proposed mesh.

  | Mesh | Near field (x 0.5–20, r 0–4.5 mm) | Plume (x 20–100 mm, r 6.5 → 12 mm) |
  |---|---|---|
  | Coarse | 0.035 mm | 0.07 mm |
  | Medium | 0.025 mm | 0.05 mm |
  | Fine | 0.0177 mm | 0.035 mm |

  Every level shares the rest of the mesh:
  - Wall layers: first cell 1.8e-7 m, ratio 1.2, 17 layers. Cells are 0.02 mm at the top of the layers and 0.01 mm
    on the fillets.
  - Nozzle and plenum: 0.05 mm.
  - Far jet (x 100–240 mm, r 0–30 mm): 0.5 mm.
  - Around the body and jet: 0.6 mm. Behind the body: 1.5 mm. Elsewhere: up to 10 mm.
- **Solver:** SST (V2003m), axisymmetric. SLAU2 fluxes with MUSCL and the Venkatakrishnan-Wang limiter (coefficient
  0.05). First-order k-ω convection, implicit time stepping.
- **Start:**
  - Coarse starts cold, from uniform ambient air with a token Mach 0.01 in +x. It runs first order for 300
    iterations and ramps to second order by iteration 800.
  - Medium starts from the coarse t4 solution, and fine from the medium t4 solution. Each is interpolated in
    `0_interp` (1 iteration at CFL 1e-8), the method used in the srp_M2.5 grid study.
- **Stages on every mesh (50,000 iterations):**

  | Stage | Iterations | CFL |
  |---|---|---|
  | t1 | 1000 | 0.1 → cap 0.5 |
  | t2 | 1000 | 0.5 → cap 1 |
  | t3 | 8000 | 1 → cap 3 |
  | t4 | 40,000 | 1 → cap 3 |

  The CFL halves on trouble, with a floor of 0.1. Nothing stops on convergence: `CONV_RESIDUAL_MINVAL= -20`.
- **Live monitor:** every 1000 iterations, Mach on five lines:
  - r = 0.01, 1.2, 2.4 (the lip) and 4.8 mm, each from x = 1 to 100 mm
  - the shear layer, from (2.0, 2.45) to (40, 5.0) mm

  It also tracks two features: the x of the first four axis-Mach peaks (the end of each shock cell) and the
  supersonic core length. The figures are `convergence_chain.png`, `conv_lines_chain.png` and
  `conv_features_chain.png`, each also copied as `*_latest.png`. No pass threshold has been set.

## Gadi resources

| Job | Queue | Most it can charge |
|---|---|---|
| `mesh_fj.pbs` | normal, 4 cores, 1:30 | 12 SU |
| Test: coarse t1 + t2, 2000 iterations, `STOP_AFTER=t2_cap1` | normalsr, 104 cores, 0:15 | 52 SU |
| Chain: coarse, then medium, then fine, each `fj_sr.pbs` with `-W depend=afterok` | normalsr, 104 cores | 208 SU per hour of walltime requested |

The chain costs more than 100 SU, so it is submitted only after the user's OK.

## Files

| File | What it does |
|---|---|
| `../../mesh/mesh_freejet.py` | Geometry (v4 walls + open domain), markers by position, size fields, mesh. `GEOM_ONLY=1` checks the markers only |
| `mesh_fj.pbs` | `exit_plane.py`, then the three meshes and their views |
| `exit_plane.py` | pvbatch: area- and mass-averaged exit-plane pressure and Mach, and pe/pa |
| `plot_fj_mesh.py` | Mesh views and the boundary-marker map |
| `freejet_template.cfg`, `setup_freejet.sh` | Build the run folders under `cases/srp_M2.5_freeJet/runs/` |
| `run_fj_chain.sh`, `fj_sr.pbs` | The driver for one mesh (0_interp, then the stages) and its PBS job (`STOP_AFTER` for the test) |
| `watch_fj_gadi.sh`, `snap_fj.py`, `conv_fj.py`, `plot_fj.py` | The live monitor and the chain figures |
