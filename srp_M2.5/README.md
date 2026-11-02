# srp_M2.5: supersonic retropropulsion, M = 2.5, axisymmetric RANS (SU2)

## Geometry
- `geometry/su2_wFillets_v3.stp` (SolidWorks, one planar surface, metres; source on the laptop:
  `PhD/Modelling/CFD/RL_Design/geoms/su2/`). v2 was deleted on 2026-10-26.
- v3 fillets: 0.1 mm at body base/shoulder corners (edges 3, 5), the lip outer corner (8) and
  the lip inner corner (10); a small fillet where the nozzle contour meets the plenum wall (12).
  The body face / lip-top junction (edges 6/7) is a sharp concave corner. Fillets are only
  there to help RANS: the real part is sharp apart from machining radii.
- Body Ø18 mm x 58.6 mm; nozzle throat r ~ 2.0 mm, exit r ~ 2.27 mm; plenum r = 2.88 mm;
  domain 1.27 m long, r = 45 mm at the model (tunnel wall).

## Confirmed assumptions (user, 2026-10-26)
- Units metres. Freestream enters at edge 16 (x = +74 mm) and flows in -x. The jet exits in
  +x. Axis of symmetry is y = 0.
- Markers: axis = 1, 15; inlet_freestream = 16 (pressure inlet); inlet_jet = 14 (pressure
  inlet); wall_noslip = 2-13 (adiabatic); wall_slip = 17, 18, 19, 21; outlet = 20
  (wave-transmissive).
- Air, ideal gas (gamma 1.4, R 287.058), Sutherland viscosity, for freestream and jet.
- T0inf = T0j = 288.15 K. p0inf and p0j are set per case. The mesh was sized with
  p0inf = 6 bar, giving T_inf = 128.07 K, p_inf = 35.117 kPa, U = 567.2 m/s, Re/m = 6.11e7.
- Wall-resolved RANS (no wall functions). SA and SST are both to be run for comparison.
- The inlet distance and the bow-shock reflection off the r = 45 mm wall are intended
  (they match the experiment).

## Mesh: mesh/srp_M2.5_coarse.su2 (2026-10-26)
- *Mesh v1 files deleted on 2026-10-28; this description is kept for the record.*
- Script: `mesh/mesh_coarse.py` (Gmsh 4.15.2 in su2-venv). Run it in a job or on a held
  compute node, never on the login node.
- 234,777 cells (226,654 quads, 8,123 triangles), 233,681 nodes (job 61480195; Gmsh varies by a few cells run to run). Mostly quads (Gmsh
  algorithm 8 with simple recombination; blossom recombination crashes with wall layers).
- Wall layers on 2-13: first cell 3.875e-7 m (y+ ~ 1 from a flat-plate estimate at x = 18 mm),
  ratio 1.2, 16 layers (0.034 mm). 18 layers made Gmsh fail ("2 intersections"); 17 is the
  maximum that meshes. The layers end cleanly at the base/axis and plenum/jet-inlet
  junctions (PointsList).
- Sizes: fillets 0.01 mm, next to walls 0.02 mm, nozzle 0.05 mm, **plume box x 0..30 mm, r 0..9 mm:
  0.05 mm** (user, 2026-10-26), growing at 0.1 to a **0.6 mm max** in the jet/bow-shock box
  (x -10..75 mm, previously 0.4 mm), body side 0.6 mm, wake to x = -300 mm 1.5 mm, far 10 mm. Growth
  away from walls and boxes is limited to 0.2.
- Nozzle y+: the first cell is sized from the freestream. A throat at p0j = 6 bar would need
  ~2.4e-7 m, and the requirement scales roughly as 1/p0j^0.86, so check y+ in the nozzle
  once p0j is chosen.
- Quality (SU2 1-iteration load, `meshcheck.cfg`): orthogonality 72.1-90 deg; CV face aspect
  ratio <= 130; CV sub-volume ratio <= 48.7; all elements correctly oriented. Gmsh: min cell
  angle 7.5 deg (2 cells, at the sharp concave 6/7 corner); max neighbour area ratio 37.7.
- `dryrun.cfg` / `meshcheck.cfg` are only for checking the mesh; they are not case setups.
- Plots: `mesh/srp_M2.5_coarse_views.png`, `geometry/edge_map_v3.png`.
- ParaView (written each run): `mesh/srp_M2.5_coarse.vtk` (cells; quality_minSICN, cell_size, is_quad),
  `mesh/srp_M2.5_coarse_boundaries.vtk` (edges; marker_id 0=axis 1=inlet_freestream 2=inlet_jet
  3=wall_noslip 4=wall_slip 5=outlet), `mesh/srp_M2.5_coarse_su2quality.vtu` (SU2 Orthogonality,
  Aspect_Ratio, Volume_Ratio from meshcheck.cfg; the other fields are from 1 iteration and are meaningless).
  Coordinates are in metres.

## Mesh v2: mesh/srp_M2.5_coarse_v2.su2 (2026-10-27)
- *Mesh v2 and v3 files (and the v4 layer-count test meshes) deleted on 2026-10-28.*
- Same script and sizes as v1 (`mesh/mesh_coarse.py`; the v1 script is kept as `mesh_coarse_v1.py`),
  but with a first cell of **1.8e-7 m** on all walls. That gives a nozzle-throat y+ of about 1 at
  p0j = 8.015 bar and a body y+ of about 0.3 at p0inf = 3.247 bar (user approved). 20 layers failed in
  Gmsh, so 19 layers (0.028 mm) were used. Job script: `mesh/mesh_job_v2.slurm`.
- 245,605 cells (237,706 quads, 7,899 triangles), 244,583 nodes.
- SU2 quality: orthogonality 76.4-90 deg, aspect ratio <= 293, volume ratio <= 76.2, all elements
  correctly oriented. Worst cell angle is 3.7 deg at the sharp concave body-face/lip corner (x 0.5 mm, r 2.9 mm).

## Mesh v4: mesh/srp_M2.5_coarse_v4.su2 (2026-10-27)
- Geometry `geometry/su2_shortDiff_wFillet_v4.STEP` (vertical outlet at x = -308.1 mm, 0.5 mm base fillet, slip
  wall r 45 -> 55 mm). Script `mesh/mesh_coarse.py` (now the v4 version), job `mesh/mesh_job_v4.sh`.
- 235.3k cells, first cell 1.8e-7 m, 17 wall layers, 1.5 mm cells through the whole wake. SU2 check: `mesh/meshcheck_v4.log`.

## meshMedium: mesh/srp_M2.5_meshMedium.su2 (2026-10-28)
- The v4 (meshCoarse) mesh plus a **0.025 mm** box at **x 0.5..20 mm, r 0..9.5 mm** (user: the jet, Mach disk and bow
  shock region; half the 0.05 mm plume-box size there), growing at 0.1/mm back to 0.05 mm. Everything else unchanged.
- Built with `mesh/mesh_job_medium.sh` (`MED_BOX="0.5,20,0,9.5,0.025"` env var in `mesh_coarse.py`; unset gives meshCoarse).
- 468,874 cells (455,029 quads, 13,845 triangles), 464,934 nodes; meshCoarse 235,331. In the box 308,396 vs 81,972 cells.
- Quality as meshCoarse: SU2 aspect ratio <= 309, volume ratio <= 137; Gmsh min angle 2.39 deg (same concave corner).
- Comparison plot: `mesh/coarse_vs_medium.png` (`mesh/compare_coarse_medium.py`).

## meshFine: mesh/srp_M2.5_meshFine.su2 (2026-10-28)
- meshCoarse plus a **0.0177 mm** box at x 0.5..20 mm, r 0..9.5 mm (user: option B, medium refined by sqrt 2; levels
  0.05 / 0.025 / 0.0177 mm). 0.025/sqrt(2) = 0.017678 mm failed in Gmsh ("2 intersections in the 1D mesh" at the wall
  layers); 0.0177 mm meshed. Built with `mesh/mesh_job_fine.sh` (tries 0.0177, 0.0176, 0.0178, 0.0175 in turn), 4 min.
- An option D version (0.0125 mm, 1,382,137 cells) was built first and deleted at the user's request on 2026-10-28.
- 772,018 cells (753,121 quads, 18,897 triangles), 765,737 nodes. In the box 606,170 (medium 308,396, coarse 81,972).
- Quality: SU2 aspect ratio <= 309, volume ratio <= 32 (coarse/medium ~137); Gmsh min angle 6.49 deg (the 2.4 deg cell
  at the concave body-face/lip corner is gone). Mesh files 269 MB.
- Three-mesh comparison: `mesh/coarse_medium_fine.png` (`mesh/compare_meshes.py`).

## Runs
- `runs/SA_p0j8.015_p0inf3.247_meshCoarse/` (renamed from `_meshv4` by the user on 2026-10-28): SA case on mesh v4 = meshCoarse. Phase 1 `1_freestream/`, Phase 2 `2_jet_staged/`
  (see its README). The only run folder kept after the cleanup on 2026-10-28.
- Deleted on 2026-10-28: `runs/SA_p0j8.015_p0inf3.247/` (mesh v2: Phase 1, Phase 2 and continuations 1-5, test) and
  `runs/SA_p0j8.015_p0inf3.247_meshv3/` (interpolation and Phase 2 on mesh v3, and its diagnostic run). Their results
  are summarised in the laptop note `context/m3-environment.md`.
